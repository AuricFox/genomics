from flask import Flask, request, redirect, render_template, url_for, send_file, flash

import os
import sys

sys.path.append('./src/')
import utils
import bioinformatics as bio

app = Flask(__name__, static_folder='static')
app.secret_key = 'my_super_secret_totaly_unbreakable_key'

# ====================================================================
# Main Pages
# ====================================================================

@app.route("/")
@app.route("/home")
def home():
    return render_template('home.html', active='analysis-page')

@app.route("/analysis")
def analysis():
    return render_template('analysis.html', active='analysis-page')

@app.route("/alignment")
def alignment():
    return render_template('alignment.html', active='alignment-page')

@app.route("/variance")
def variance():
    return render_template('variance.html', active='variance-page')

@app.route("/assembly")
def assembly():
    return render_template('assembly.html', active='assembly-page')

@app.route("/phylogeny")
def phylogeny():
    return render_template('phylogeny.html', active='phylogeny-page')

# Custom page not found
@app.errorhandler(404)
def page_not_found(error):
    return render_template('404.html'), 404
'''
# Custom page not found
@app.errorhandler(500)
def server_error(error):
    return render_template('404.html'), 404
'''
# ====================================================================
# Sequence Processing Functions
# ====================================================================
# Handles sequence counting of codons, amino acids, and/or k-mers
@app.route("/sequence_analysis", methods=["POST"])
def sequence_analysis():
    try:
        # User submitted form data
        codon = request.form.get('codon', type=str)
        amino = request.form.get('amino', type=str)
        kmer = request.form.get('kmer', type=str)

        if codon is None and amino is None and kmer is None:
            flash('A Sequence Type Must Be Selected', 'error')
            return redirect(request.referrer)

        k = request.form.get('n-mer', type=int)
        file_type = request.form.get('file_type', type=str)

        file = request.files["file"]                    # Get user's submitted file
        file_path = utils.create_file(file)                   # Get temp file path

        if file_path is None:                           # Incorrect file was submitted
            flash(f'{file.filename} is an Invalid File or FileType', 'error')
            return redirect(request.referrer)

        data = bio.get_sequence_data(                   # Path of zip file
            file=file.filename, 
            codon=codon, 
            amino=amino, 
            kmer=kmer, 
            k=k,
            file_type=file_type
        )
        
        # Return zip file if found else return error message
        if 'zip_file' in data:
            response = send_file(data['zip_file'], as_attachment=True) 
        else:
            flash(f'Failed to create Zip File!', 'error')
            response = redirect(request.referrer)

    except Exception as e:
        flash(f'An Error occured: {str(e)}', 'error')
        response = redirect(request.referrer)

    finally:
        try:
            os.remove(file_path)
        except Exception as e:
            print(f'Error removing file: {str(e)}')
    
    return response

# ====================================================================
# Sequence Alignment Functions
# ====================================================================
# Handles sequence alignment of bases
@app.route("/sequence_alignment", methods=["POST"])
def sequence_alingment():

    try:
        match_point = request.form.get('match_point', type=int)
        match_penalty = request.form.get('match_penalty', type=int)
        gap_penalty = request.form.get('gap_penalty', type=int)
        ignore_gaps = request.form.get('ignore_gaps', type=bool)

        if match_point < 1:
            flash('Match Point Must Be Greater Than Zero!', 'error')
            return redirect(request.referrer)

        if match_penalty > -1:
            flash('Match Penalty Must Be Less Than Zero!', 'error')
            return redirect(request.referrer)

        if gap_penalty > -1:
            flash('Gap Penalty Must Be Less Than Zero!', 'error')
            return redirect(request.referrer)
        
        file_type = request.form.get('file_type', type=str)

        file1 = request.files['file1']
        file2 = request.files['file2']
        file_path1 = utils.create_file(file1)
        file_path2 = utils.create_file(file2)

        if file_path1 is None:
            flash(f'{file1.filename} is an Invalid File or FileType', 'error')
            return redirect(request.referrer)
        
        if file_path2 is None:
            flash(f'{file2.filename} is an Invalid File or FileType', 'error')
            return redirect(request.referrer)
        
        # Creating alignment data
        data = bio.get_alignment_data(
            file1=file_path1,
            file2=file_path2,
            gap_pen=gap_penalty,
            match_point=match_point,
            match_pen=match_penalty,
            ignore=ignore_gaps,
            file_type=file_type
        )

        # Return zip file if found else return error message
        if 'zip_file' in data:
            response = send_file(data['zip_file'], as_attachment=True) 
        else:
            flash(f'Failed to create Zip File!', 'error')
            response = redirect(request.referrer)
    
    except Exception as e:
        flash(f'An Error occured: {str(e)}', 'error')
        response = redirect(request.referrer)

    finally:
        try:
            utils.remove_files([file_path1, file_path2])
        except:
            print(f'Error removing file(s): {str(e)}')

    return response

# ====================================================================
# Run Main
# ====================================================================
if __name__ == "__main__":
    app.run()
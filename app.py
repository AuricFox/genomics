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
# Sequence Processing Function(s)
# ====================================================================
# Handles sequence counting of codons, amino acids, and/or k-mers
@app.route("/sequence_analysis", methods=["POST", "GET"])
def sequence_analysis():
    try:
        # User submitted form data
        codon = True if request.form.get('codon', type=str) is not None else False
        amino = True if request.form.get('amino', type=str) is not None else False
        kmer = True if request.form.get('kmer', type=str) is not None else False

        if not codon and not amino and not kmer:
            flash('A Sequence Type Must Be Selected!', 'error')
            return redirect(request.referrer)

        k = request.form.get('n-mer', type=int)
        file_type = request.form.get('file_type', type=str)

        # File that contains the sequence data being analyzed
        file = request.files["file"]
        file_path = utils.create_file(file)

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
            flash(f"{data['error']}", 'error')
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
# Sequence Alignment Function(s)
# ====================================================================
@app.route("/sequence_alignment", methods=["POST", "GET"])
def sequence_alingment():

    try:
        # User inputs that regulate the alignment
        match_point = request.form.get('match_point', type=int)
        match_penalty = request.form.get('match_penalty', type=int)
        gap_penalty = request.form.get('gap_penalty', type=int)
        ignore_gaps = True if request.form.get('ignore_gaps', type=str) is not None else False

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

        # Files that contain the sequences being aligned
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
            flash(f"{data['error']}", 'error')
            response = redirect(request.referrer)
    
    except Exception as e:
        flash(f'An Error Occured: {str(e)}', 'error')
        response = redirect(request.referrer)

    finally:
        try:
            utils.remove_files([file_path1, file_path2])
        except Exception as e:
            print(f'Error Removing File(s): {str(e)}')

    return response

# ====================================================================
# Sequence Variance Function(s)
# ====================================================================
@app.route("/sequence_variance", methods=["POST", "GET"])
def sequence_variance():

    try:
        # User inputs used for determining which files are returned
        plot = True if request.form.get('plot', type=str) is not None else False
        smooth = True if request.form.get('smooth', type=str) is not None else False
        vRegion = True if request.form.get('vRegion', type=str) is not None else False
        spacing = request.form.get('spacing', type=int)
        numV = request.form.get('numV', type=int)

        if spacing < 1 or spacing > 1000:
            flash('Spacing Must Be Between 1 and 1000!', 'error')
            return redirect(request.referrer)

        if numV < 1 or numV > 1000:
            flash('Number of V-Regions Must Be Between 1 and 1000!', 'error')
            return redirect(request.referrer)
        
        file_type = request.form.get('file_type', type=str)

        # File containing the series of sequences used for calculating variance
        file = request.files['file']
        file_path = utils.create_file(file)

        if file_path is None:
            flash(f'{file.filename} is an Invalid File or FileType', 'error')
            return redirect(request.referrer)
        
        # Creating variance data
        data = bio.get_variance_data(
            file=file_path,
            plot=plot,
            smooth=smooth,
            vRegion=vRegion,
            spacing=spacing,
            numV=numV,
            file_type=file_type
        )

        # Return zip file if found else return error message
        if 'zip_file' in data:
            response = send_file(data['zip_file'], as_attachment=True) 
        else:
            flash(f"{data['error']}", 'error')
            response = redirect(request.referrer)
    
    except Exception as e:
        flash(f'An Error Occured: {str(e)}', 'error')
        response = redirect(request.referrer)

    finally:
        try:
            utils.remove_files([file_path])
        except Exception as e:
            print(f'Error Removing File(s): {str(e)}')

    return response

# ====================================================================
# Sequence Phylogeny Function(s)
# ====================================================================
@app.route("/sequence_phylogeny", methods=["POST", "GET"])
def sequence_phylogeny():

    try:
        # File that contains the series of sequences used for creating phylogeny tree
        file = request.files['file']
        file_path = utils.create_file(file)

        if file_path is None:
            flash(f'{file.filename} is an Invalid File or FileType', 'error')
            return redirect(request.referrer)
        
        # Creating phylogeny data
        data = bio.get_phylogeny_data(file=file_path)

        # Return zip file if found else return error message
        if 'zip_file' in data:
            response = send_file(data['zip_file'], as_attachment=True) 
        else:
            flash(f"{data['error']}", 'error')
            response = redirect(request.referrer)
    
    except Exception as e:
        flash(f'An Error Occured: {str(e)}', 'error')
        response = redirect(request.referrer)

    finally:
        try:
            utils.remove_files([file_path])
        except Exception as e:
            print(f'Error Removing File(s): {str(e)}')

    return response

# ====================================================================
# Sequence De Bruijn Function(s)
# ====================================================================
@app.route("/sequence_assembly", methods=["POST", "GET"])
def sequence_assembly():
    files = []

    try:
        # User inputs for determining how substrings are processed in the graph
        k = request.form.get('k_size', type=int)
        cut = request.form.get('p_size', type=int)

        if k < 3 or k > 100:
            flash(f"The k-mer size must be between 3 and 100: K-mer size is {k}!")
            return redirect(request.referrer)
        if cut < 1 or cut > 99:
            flash(f"The size of the prefix/suffix must be between 1 and 99: size is {cut}!")
            return redirect(request.referrer)


        # File that contains the main sequence fragments
        seq_file = request.files['seq_file']
        seq_path = utils.create_file(seq_file)
        files.append(seq_path)

        if seq_path is None:
            flash(f'{seq_file.filename} is an Invalid File or FileType', 'error')
            return redirect(request.referrer)

        # File that contains the reference genome for checking accuracy
        ref_file = request.files['ref_file']
        ref_path = None

        # User wants an alignment check
        if ref_file is not None:
            ref_path = utils.create_file(ref_file)
            files.append(ref_path)

            # Make sure reference file is legit
            if ref_path is None:
                flash(f'{ref_file.filename} is an Invalid File or FileType', 'error')
                return redirect(request.referrer)
        
        # Creating assembled data
        data = bio.get_assembled_data(seq_file=seq_path, ref_file=ref_path, k=k, cut=cut)

        # Return zip file if found else return error message
        if 'zip_file' in data:
            response = send_file(data['zip_file'], as_attachment=True) 
        else:
            flash(f"{data['error']}", 'error')
            response = redirect(request.referrer)
    
    except Exception as e:
        flash(f'An Error Occured: {str(e)}', 'error')
        response = redirect(request.referrer)

    finally:
        try:
            utils.remove_files(files=files)
        except Exception as e:
            print(f'Error Removing File(s): {str(e)}')

    return response

# ====================================================================
# Run Main
# ====================================================================
if __name__ == "__main__":
    app.run()
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

# Accessing Home page
@app.route("/")
@app.route("/home")
def home():
    return render_template('home.html')

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

# --------------------------------------------------------------------
# Displays the results of the sequence analysis
@app.route("/sequence_results")
def sequence_results():
    return render_template('sequence_results.html')

# ====================================================================
# Sequence Alignment Functions
# ====================================================================
# Handles sequence alignment of bases
@app.route("/sequence_alignment", methods=["POST"])
def sequence_alingment():
    return redirect('/home')

# --------------------------------------------------------------------
# Displays the results of the sequence alignment
@app.route("/alignment_results")
def aligment_results():
    return render_template('aligment_results.html')

# ====================================================================
# Run Main
# ====================================================================
if __name__ == "__main__":
    app.run()
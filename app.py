from flask import Flask, request, redirect, render_template, url_for

import os
import sys
import re
import mimetypes

sys.path.append('./src/')
import bioinformatics as bio

app = Flask(__name__, static_folder='static')

# ====================================================================
# Main Pages
# ====================================================================

# Accessing Home page
@app.route("/")
@app.route("/home")
def home():
    return render_template('home.html')
'''
# Custom page not found
@app.errorhandler(404)
def page_not_found(error):
    return render_template('404.html'), 404

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
    # User submitted form data
    codon = request.form.get('codon', type=str)
    amino = request.form.get('amino', type=str)
    kmer = request.form.get('kmer', type=str)
    k = request.form.get('n-mer', type=int)
    file_type = request.form.get('file_type', type=str)
    
    file = request.files["file"]                    # Get user's submitted file
    file_path = create_file(file)                   # Get temp file path
    
    data = bio.get_sequence_data(                   # Path of zip file
        file=file.filename, 
        codon=codon, 
        amino=amino, 
        kmer=kmer, 
        k=k,
        file_type=file_type
    )

    print(data)
    os.remove(file_path)                            # File is no longer needed
    return redirect('/home')

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
# File Function(s)
# ====================================================================
def create_file(file):
    '''
    Takes in a file object, sanitizes, validates, and saves it to the temp directory

    Parameter(s):
        file: the user input file being saved

    Output(s):
        file_path (str): the path to the saved file
    '''

    # Replace special characters with underscores
    sanitized_name = re.sub(r'[\\/*?:"<>|]', '_', file.filename)
    # Remove leading and trailing whitespace
    file.filename = sanitized_name.strip()

    allowed_mime_types = ['image/jpeg', 'image/png', 'application/pdf', 'text/plain', None]
    allowed_extensions = ['.jpg', '.jpeg', '.png', '.pdf', '.fna', '.fastq', '.txt']

    # Get the file's MIME type and extension
    file_mime_type, _ = mimetypes.guess_type(file.filename)
    file_extension = os.path.splitext(file.filename)[1].lower()

    # Check if the file's MIME type or extension is allowed
    if file_mime_type not in allowed_mime_types or file_extension not in allowed_extensions:
        print(f'{file.filename} MIME type or extension is not supported! '
              f'MIME type: {file_mime_type}, Extension: {file_extension}')
        return None
    print('ksbvibeiornvornv')
    path = os.path.join(os.path.dirname(__file__), "src/temp")  # Path where file will be saved

    os.makedirs(path, exist_ok=True)                            # Create path if it doesn't exist
        
    file_path = os.path.join(path, file.filename)               # Creating saved file path
    file.save(file_path)                                        # Saving input file

    return file_path

# --------------------------------------------------------------------
def remove_file(filename:str):
    '''
    Takes in a file object and removes the file from the temp directory

    Parameter(s):
        file (str): the input file being removed

    Output(s): None
    '''

    path = os.path.join(os.path.dirname(__file__), "src/temp")  # Path where file is saved
    file_path = os.path.join(path, filename)                    # Creating saved file path

    try:
        os.remove(file_path)                                    # File is no longer needed
    except OSError as e:
        print(f'Error while removing file {filename}: {e}')

# ====================================================================
# Run Main
# ====================================================================
if __name__ == "__main__":
    app.run()
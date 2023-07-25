from flask import Flask, request, redirect, render_template, url_for, jsonify

import os
import sys
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
# Bioinformatics Processing Functions
# ====================================================================

# Accessing counting_codons Page
@app.route("/sequence_analysis", methods=["POST"])
def counting_codons():
    # User submitted form data
    codon = request.form.get('codon')
    amino = request.form.get('amino')
    kmer = request.form.get('kmer')
    k = request.form.get('n-mer', type=int)
    
    file = request.files["file"]                                # Get user's submitted file
    path = os.path.join(os.path.dirname(__file__), "src/temp")  # Path where file will be saved

    if not os.path.exists(path):                                # Checks if path exists
        os.makedirs(path)                                       # Create path if it doesn't exist
        
    file_path = os.path.join(path, file.filename)               # Creating saved file path
    file.save(file_path)                                        # Saving input file
    data = bio.get_sequence_data(file.filename, codon, amino, kmer, k)   # Get the totals

    print(data)
    os.remove(file_path)                                        # File is no longer needed
    return redirect('/home')

# Accessing codon_results Page
@app.route("/codon_results")
def codon_results():
    return render_template('codon_results.html',)

# ====================================================================
# Run Main
# ====================================================================
if __name__ == "__main__":
    app.run()
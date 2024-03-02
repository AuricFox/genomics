'''
Serves as the entry point for the web server
'''

from app import init_app
app = init_app()

if __name__ == "__main__":
    '''
    Run the Flask application
    Set host='0.0.0.0' to make the server externally visible
    Set port=5000 to use port 5000 (default for Flask)
    '''
    app.run(host='0.0.0.0')
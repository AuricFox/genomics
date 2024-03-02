'''
Flask App configuration

To set the configuration type, navigate to __init__.py file in the 
app directory and change the configuration in the init_app function

Using a production configuration:
app.config.from_object('config.ProdConfig')

Using a development configuration:
app.config.from_object('config.DevConfig')

Configuration Variables:

FLASK_DEBUG: provides logging for debugging purposes
SECRET_KEY: strings used to encrypt sensitive data
SERVER_NAME: app's domian name

More Info:
https://flask.palletsprojects.com/en/3.0.x/config/
'''

from os import environ, path
from dotenv import load_dotenv

BASEDIR = path.abspath(path.dirname(__file__))
load_dotenv(path.join(BASEDIR, '.env'))

class Config:
    '''
    Base config
    '''
    SECRET_KEY = environ.get('SECRET_KEY')
    SESSION_COOKIE_NAME = environ.get('SESSION_COOKIE_NAME')
    STATIC_FOLDER = 'static'
    TEMPLATES_FOLDER = 'templates'

    # Database
    #SQLALCHEMY_DATABASE_URI = environ.get('SQLALCHEMY_DATABASE_URI')
    #SQLALCHEMY_TRACK_MODIFICATIONS = False

    # AWS Secrets
    #AWS_SECRET_KEY = environ.get('AWS_SECRET_KEY')
    #AWS_KEY_ID = environ.get('AWS_KEY_ID')


class ProdConfig(Config):
    '''
    Production config
    '''
    FLASK_ENV = "production"
    FLASK_DEBUG = False
    DATABASE_URI = environ.get('PROD_DATABASE_URI')


class DevConfig(Config):
    '''
    Development config
    '''
    FLASK_ENV = "development"
    FLASK_DEBUG = True
    DATABASE_URI = environ.get('DEV_DATABASE_URI')
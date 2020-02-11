import logging
import os
from logging import NullHandler



logging.getLogger(__name__).addHandler(NullHandler())

def get_data_filename(filename):
    data_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')
    return os.path.join(data_dir, filename)




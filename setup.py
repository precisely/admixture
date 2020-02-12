"""
GENOMIX

Genetics utility code for Insitome products
"""
import re
from setuptools import setup, find_packages


def get_version():
    with open('version.py') as version_file:
        return re.search(r"""__version__\s+=\s+(['"])(?P<version>.+?)\1""",
                         version_file.read()).group('version')


install_requires = [
    "click==7.0",
    "PyVCF==0.6.8",
    "nose==1.3.7"
]


# FIXME: Is this necessary? Since this won't be distributed as a separate app,
# why not just combine all dependencies in install_requires? Also, `python
# setup.py test` is apparently deprecated.
test_requirements = [
    'mock',
    'pytest'
]

__NAME__ = "ancestry"
__author__ = "Gareth Highnam"
__license__ = "Copyright Precise.ly"
__copyright__ = "2020"

config = {
    'name': __NAME__,
    'license': __license__,
    'description': 'Ancestry and admixture code for Precisely products',
    'long_description': __doc__,
    'author': __author__,
    'url': 'https://github.com/precisely/admixture',
    'version': get_version(),
    'packages': find_packages(exclude=('tests',)),
    'include_package_data': True,
    'install_requires': install_requires,
    'zip_safe': False,
    'platforms': 'any',
    'entry_points': """
        [console_scripts]
        ancestry=ancestry.cli:cli
    """,
    'classifiers': [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Software Development :: Libraries :: Python Modules'
        ],
    'tests_require': install_requires + test_requirements,
}

setup(**config)

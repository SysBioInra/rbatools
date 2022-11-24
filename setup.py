# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.txt'), encoding='utf-8') as f:
    long_description = f.read()

with open(path.join(here, 'rbatools', '_version.py'), 'r') as f:
    version = f.readline().split("'")[1]

with open(path.join(here, 'rbatools', '_authors.py'), 'r') as f:
    authors = f.readline().split("'")[1]

setup(name='rbatools',

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version=version,
    description='Programming interface to resource allocation modelling with the Resource Balance Analysis (RBA) method.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    # The project's main homepage.
    url='https://rba.inrae.fr',
    # Author details
    author=authors,
    author_email='anne.goelzer@inrae.fr, wolfram.liebermeister@inrae.fr',

    # Choose your license
    license='GPLv3',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',

        'Programming Language :: Python :: 3',
    ],
    # What does your project relate to?
    keywords=[
        'metabolism',
        'resource allocation modelling'
        'Resource Balance Analysis',
        'molecular biology',
        'cell biology',
        'biochemistry',
        'systems biology',
        'computational biology',
        'mathematical modeling',
        'numerical simulation',
    ],

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(),
    #packages=find_packages(where="src"),
    install_requires=[
        'rbapy',
        'swiglpk',
        'sbtab>=1.0.6',
        'jxmlease',
        'urllib3',
    ],
    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    entry_points={
        'console_scripts': [
            'run-growth-rate-optimization=rbatools.cli.run_growth_rate_optimization:main',
            'generate-sbtab-of-model-for-html=rbatools.cli.generate_sbtab_of_model_for_html:main',
        ],
    },

)

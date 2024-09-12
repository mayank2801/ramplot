# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 17:42:20 2024

@author: mayank
"""

# setup.py

from setuptools import setup, find_packages

setup(
    name="ramplot",
    version="1.0.3",
    packages=find_packages(),
    package_data={
    'ramplot': ['RamBoundry/*'],  # Include all files in the `data` directory
    },
    include_package_data=True,  # Ensure package data is included
    entry_points={
        'console_scripts': [
            'ramplot=ramplot.cli:main',  # 'my-cli' is the command-line tool name
        ],

    },
    install_requires=[ 'numpy','pandas','biopython','matplotlib','mdanalysis','pytest-shutil','setuptools'
        # List your package dependencies here
    ],
    classifiers=[
    # Development Status
    'Development Status :: 5 - Production/Stable',  # Or '5 - Production/Stable' if your package is stable

    # Audience
    'Intended Audience :: Science/Research',
    'Intended Audience :: Education',
    'Intended Audience :: Healthcare Industry',

    # License (choose the appropriate license for your package)
    'License :: OSI Approved :: MIT License',  # Or any other license you choose

    # Programming Language
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'Programming Language :: Python :: 3.10',

    # Operating System
    'Operating System :: OS Independent',

    # Topics
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Scientific/Engineering :: Chemistry',
    'Topic :: Scientific/Engineering :: Visualization',
    'Topic :: Software Development :: Libraries :: Python Modules',

    # Natural Language
    'Natural Language :: English',
],
    keywords=['Ramachandran plot', 'Protein structure analysis', 'Structural biology tools', 'Bioinformatics visualization', 'Protein conformation', 'Amino acid dihedrals', 'Torsion angle plots', 'Protein folding analysis', 'Structural bioinformatics', 'Protein modeling', 'Ramachandran plot software', 'Visualize protein structure', 'Protein backbone', 'angles', 'Phi and Psi angles', 'Protein structure validation', 'Molecular dynamics analysis', 'Protein secondary structure', 'Protein dihedral angles', 'Structural biology software', 'Protein structure research', 'Bioinformatics researchers', 'Structural biologists', 'Computational biology tools', 'Protein data visualization', 'Molecular modeling software'],
    author="Mayank Kumar & R.S. Rathore",
    author_email="mayank2801@gmail.com, rsrathore@cusb.ac.in",
    description="A brief description of your package.",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://ramplot.in",  # Replace with your GitHub repo
    python_requires='>=3.6',
)

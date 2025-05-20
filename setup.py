from setuptools import setup, find_packages

setup(
    name='guacamole',
    version='0.1.0',
    author='Laurenz Holcik',
    author_email='laurenz.holcik@univie.ac.at',
    description='GC-aware species abundance estimation from metagenomic data.',
    long_description='This project implements GuaCAMOLE, a tool for species abundance estimation from metagenomic data.',
    url='https://github.com/cibiv/guacamole',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'numpy',
        'matplotlib',
        'seaborn',
        'biopython',
        'scipy',
        'qpsolvers',
    ],
    entry_points={
        'console_scripts': [
            'create-reference-dist=guacamole.create_reference_dist:main',
            'guacamole=guacamole.guacamole:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3.10',
        'License :: OSI Approved :: MIT License',
        'Operating System :: Linux',
    ],
    python_requires='>=3.10',
    )

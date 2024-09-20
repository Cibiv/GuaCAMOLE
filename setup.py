from setuptools import setup, find_packages

setup(
    name='guacamole',
    version='0.1.0',
    author='Your Name',  # Replace with your name
    author_email='your.email@example.com',  # Replace with your email
    description='GC-aware species abundance estimation from metagenomic data.',
    long_description='This project implements GuaCAMOLE, a tool for species abundance estimation from metagenomic data.',
    long_description_content_type='text/markdown',
    url='https://github.com/yourusername/guacamole',  # Replace with your GitHub repository
    packages=find_packages(),  # Automatically find packages in your project
    install_requires=[
        'pandas',
        'numpy==1.26.4',
        'matplotlib',
        'seaborn',
        'biopython',
        'scipy',
        'qpsolvers',
    ],
    entry_points={
        'console_scripts': [
            'create-reference-dist=guacamole.create_reference_dist:main',  # Entry point for create_reference_dist.py
            'guacamole=guacamole.guacamole:main',  # Entry point for guacamole.py
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3.10',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.10',
    )

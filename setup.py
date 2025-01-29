from setuptools import setup

setup(
    name='PoreAnalyser', #'PoreFinding',
    version='0.1.1',    
    description='A Python package for analysing (ion channel) pore profiles',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://poreanalyser.bioch.ox.ac.uk/',
    author='David Seiferth',
    author_email='david.seiferth@oriel.ox.ac.uk',
    python_requires='>=3.6,<3.13',
    #license='BSD 2-clause',
    packages=['PoreAnalyser'],
    install_requires=[
        'numpy>=1.0',  # 1.22.0
        'MDAnalysis>=2.0,<3.0', #2.0.0
        'matplotlib>=0.1', # 3.5.1 for base , for inline 0.1.3
        'pandas>=1.3', # 1.3.5
        'streamlit>=1.0', # 1.22.0
        'stmol==0.0.9 ', # 1.22.0
        'py3Dmol',      # 2.0.3
        'ipyspeck==0.6.1', # 0.6.1
        'ipywidgets==7.6.3', # 7.6.3
        'scipy', #.optimize #==1.7.3 # 1.7.3
        'altair<5',
        'nglview', #3.0.3
        'ipython_genutils==0.2.0',                    
                      ],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        #'License :: OSI Approved :: BSD License',  
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
    ],
)

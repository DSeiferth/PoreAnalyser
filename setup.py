from setuptools import setup

setup(
    name='PoreFinding',
    version='0.0.1',    
    description='A example Python package',
    url='https://huggingface.co/spaces/DSeiferth/PoreFinding_pdb',
    author='David Seiferth',
    author_email='david.seiferth@oriel.ox.ac.uk',
    #license='BSD 2-clause',
    packages=['PoreFinding'],
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
    ],
)
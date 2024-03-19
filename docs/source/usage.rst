Usage
=====

.. _installation:

Installation
------------

To use PoreFinding, first install it using pip:

.. code-block:: console

   (.venv) $ pip install PoreFinder
   
The PoreFinder package has been tested with python versions 3.9, 3.10 and 3.11.
Using a virtual environment is recommended:

.. code-block:: console

   $ conda create -n PoreFinder python=3.9 -y
   $ conda activate PoreFinder
   (PoreFinder) $ pip install PoreFinder

Creating Pore profiles
----------------

To analyse ion channel pores, you can initialise ``PoreFinding.PoreAnalysis()`` class
with an array of pdb structures that you want to analyse.


For example:

>>> import PoreFinder as pf
>>> p = 'PoreFinder/pdb_models/'
>>> pdb_array = [p+'8fe1.pdb']
>>> c = pf.PoreAnalysis(pdb_array, num_circle=20,)
>>> c.hole_analysis(plot_lines=True, legend_outside=False, title='', f_size=15, )
>>> c.hole_df 
>>> c.pathway_visualisation(index_model=0, f_end='_circle.pdb')
>>> c.ellipsoid_analysis(index_model=0)
>>> c.pathway_visualisation(0, f_end='_ellipsoid.pdb')

Running Streamlit app locally
------------------------------
After having installed PoreFinding locally, you can run the streamlit app yourself:

.. code-block:: console

   (.venv) $ streamlit run app.py

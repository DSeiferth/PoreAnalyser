Usage
=====

.. _installation:

Installation
------------

To use PoreAnalyser, first install it using pip:

.. code-block:: console

   (.venv) $ pip install PoreAnalyser
   
The PoreAnalyser package has been tested with python versions 3.9, 3.10 and 3.11.
Using a virtual environment is recommended:

.. code-block:: console

   $ conda create -n PoreAnalyser python=3.9 -y
   $ conda activate PoreAnalyser
   (PoreAnalyser) $ pip install PoreAnalyser

Creating Pore profiles
----------------

To analyse ion channel pores, you can initialise ``PoreAnalyser.PoreAnalysis()`` class
with an array of pdb structures that you want to analyse.


For example:

>>> import PoreAnalyser as pa
>>> p = 'PoreAnalyser/pdb_models/'
>>> pdb_array = [p+'8fe1.pdb']
>>> c = pa.PoreAnalysis(pdb_array, num_circle=20,)
>>> c.hole_analysis(plot_lines=True, legend_outside=False, title='', f_size=15, )
>>> c.hole_df 
>>> c.pathway_visualisation(index_model=0, f_end='_circle.pdb')
>>> c.ellipsoid_analysis(index_model=0)
>>> c.pathway_visualisation(0, f_end='_ellipsoid.pdb')
>>> g_hole_bulk, g_pa_bulk, g_hole, g_pa = c.conductance_estimation(index_model=0)
>>> print('Conductance of the pore is:', 'g_hole', g_hole_bulk, 'g_pa', g_pa_bulk, 'using the bulk conductivity')
>>> print('Conductance of the pore is:', 'g_hole', g_hole, 'g_pa', g_pa, 'using the pore conductivity model')

Example with trajectory:

>>> pdb_models = [fname+'.tpr', fname+'.xtc']
>>> pore_analysis = pa.PoreAnalysis(pdb_array=pdb_models, trajectory=True, traj_frames=10)
>>> pore_analysis.hole_analysis()
>>> pore_analysis.plt_trajectory_average(HOLE_profile=True)
>>> for i in range(10): pore_analysis.ellipsoid_analysis(index_model=i)
>>> pore_analysis.plt_trajectory_average(HOLE_profile=False)

Running Streamlit app locally
------------------------------
After having installed PoreAnalyser locally, you can run the streamlit app yourself:

.. code-block:: console

   (.venv) $ streamlit run app.py

Welcome to PoreFinding's documentation!
===================================

**PoreFinding** is a Python library for analysing (ion channel) 
pore profiles. 
Over the last two decades, advances in structural biology along with recent artificial intelligence–driven structure prediction algorithms, such as AlphaFold, have revealed a plethora of 3-D ion channel and nanopore structures in different conformational states. However, in nearly every case, these structures still require functional annotation. Different tools, such as HOLE and CHAP, allow the analysis of the physical dimensions of the pore running through an ion channel. Here, we present a package that allows users to calculate the pore profile of any input structure. Based on the well-established HOLE programme, we add a new feature to capture pore asymmetry by using an ellipsoidal probe particle.


It uses the `HOLE <https://github.com/osmart/hole2/>`_
programme wrapped by `MDAnalysis <https://github.com/MDAnalysis/mdanalysis>`_.

Check out the :doc:`usage` section for further information, including
how to :ref:`installation` the project.

.. note::

   This project is under active development.

Contents
--------

.. toctree::

   usage
   api
   
.. automodule:: package_name.module
   :members:

Interactive webservice (no installation needed)
================================================
You can try out the package without installing anything using a
`interactive web-service <https://poreanalyser.bioch.ox.ac.uk/>`_
hosted on a webserver of the Department of Biochemistry (University of Oxford). 

Output to download
-------------------

Conent and Usage of (zipped) HOLE output files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
example: pdb_name = '7tu9_aligned_z.pdb'
- pdb files
  - pdb_name: uploaded pdb file (aligned to z-axis)
  - pdb_name + '_circle.pdb': point cloud for pore surface
- vmd files
  - pdb_name + ".vmd": vmd surface for uploaded pdb file
  - "visualise_pathway_hole.tcl": vmd script for plotting the pore surface; the script can be used in the following way: "vmd -e visualise_pathway_hole.tcl -args  7tu9_aligned_z.pdb 7tu9_aligned_z.vmd"
- other files
  - "hole_pathway_profile.csv": A DataFrame containing the results of the hole analysis, with the following columns:
    - 'Label z [A]': the z-coordinate of each point along the pore axis.
    - 'Label Radius [A]': the radius of the pore at each point.
    - 'Label' corresponds to the labels provided in the `labels` parameter in the hole_analysis.analysis() function.
  - "README.md"
  - "hole.out": HOLE output
  - "hole_pathway_profile."+fig_format: pathway profile figure in desired format

Conent and Usage of (zipped) output files for pathfinding with an ellipsoidal probe particle
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
example: pdb_name = 7tu9_aligned_z
- vmd files
  - pdb_name+'.pdb_pathway_ellipse.vmd'
  - "visualise_pathway_hole.tcl"
  - the script can be used in the following way: "vmd -e visualise_pathway_hole.tcl -args  7tu9_aligned_z.pdb 7tu9_aligned_z.pdb_pathway_ellipse.vmd.vmd"
- pdb files 
  - pdb_name + '.pdb_ellipsoid.pdb'
  - pdb_name+'.pdb'
- other files
  - "README.md"
  - pdb_name + '.pdb_pathway_ellipse.txt': A DataFrame containing the results of the pathfinding analysis, with the following columns:
    - "x", "y", "z": position of ellipsoid center
    - "a", "b": larger and smaller radius of the ellipsoid
    - "theta": orientation of ellipsoid
  - figure plotting the spherical HOLE radius and the larger radius of the ellipsoid

PoreFinding
==============================

[//]: # (Badges)
[![pytest](https://github.com/simonlichtinger/PyMEMENTO/actions/workflows/run_pytest.yml/badge.svg)](https://github.com/simonlichtinger/PyMEMENTO/actions/workflows/run_pytest.yml)
[![Documentation Status](https://readthedocs.org/projects/porefinding/badge/?version=latest)](https://porefinding.readthedocs.io/en/latest/?badge=latest)
![Tests](https://github.com/github/docs/actions/workflows/python-package.yml/badge.svg)

---
title: PoreFinding_pdb
emoji: 📊
colorFrom: red
colorTo: green
sdk: streamlit
sdk_version: 1.19.0
app_file: app.py
pinned: false
license: mit
---

[Try out this protoype on HugginFace without installing anything](https://huggingface.co/spaces/DSeiferth/PoreFinding_pdb)

# What does this package add?
- Pore finding with new features: 
- Capture pore asymmetry.
  - Asymmetry of crystal/cryoEM structures due to heterogeneous subunit composition.
  - from crystal structure broken in simulations.
- Making a tool accessible without installation or download of any package.


# Path finding with ellipsoidal probe particle

1. Align principal axis to z-axis
2. Load HOLE output file with positions and radii of probes.
3. Loop through all spherical probe particles: 
    a) Ellipsoid initialized with spherical probe particle parameters from HOLE output. 
    b) First Nelder-Mead 4-dim optimization to insert ellipsoid with smaller bounds for parameters [x, y, r1, θ ]. 
    c) Second optimization with larger boundaries for parameters to further increase ellipsoid. The loop takes around 60s to complete...

# Conent and Usage of (zipped) HOLE output files
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

# Conent and Usage of (zipped) output files for pathfinding with an ellipsoidal probe particle
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

## Notes

[Hugging Face Spaces](https://huggingface.co/docs/hub/spaces) work as `git` repositories. To keep everything on GitHub but publish on Hugging Face, add the Hugging Face Space repository as a remote repository:

```bash
git remote add hf https://huggingface.co/spaces/DSeiferth/PoreFinding_pdb
```
## Acknowledgements
* Rocco Meli for pointing out streamlit and hugginface
* SBCB community for discussion

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

# Usage of (zipped) output files
example: uploaded_file.name = '7tu9.pdb'
- uploaded_file.name: uploaded pdb file
- uploaded_file.name+".pdb.vmd": vmd surface for uploaded pdb file
- "visualise_pathway_hole.tcl": vmd script for plotting the pore surface; the script can be used in the following way: "vmd -e visualise_pathway_hole.tcl -args  7tu9.pdb 7tu9.pdb.vmd"
- uploaded_file.name + '_circle.pdb': point cloud for pore surface
- "hole_pathway_profile.csv": A DataFrame containing the results of the hole analysis, with the following columns:
        - 'Label z [A]': the z-coordinate of each point along the pore axis.
        - 'Label Radius [A]': the radius of the pore at each point.
        'Label' corresponds to the labels provided in the `labels` parameter in the hole_analysis.analysis() function.
- "README.md"
- "hole.out": HOLE output

## Notes

[Hugging Face Spaces](https://huggingface.co/docs/hub/spaces) work as `git` repositories. To keep everything on GitHub but publish on Hugging Face, add the Hugging Face Space repository as a remote repository:

```bash
git remote add hf https://huggingface.co/spaces/DSeiferth/PoreFinding_pdb
```
## Acknowledgements
* Rocco Meli 

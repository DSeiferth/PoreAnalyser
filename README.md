---
title: PoreAnalyser
colorFrom: red
colorTo: green
sdk: streamlit
sdk_version: 1.19.0
app_file: app.py
pinned: false
license: LGPL 2.1
---

PoreAnalyser
==============================

[//]: # (Badges)
[![Documentation Status](https://readthedocs.org/projects/porefinding/badge/?version=latest)](https://porefinding.readthedocs.io/en/latest/?badge=latest) ![License: MIT](https://img.shields.io/badge/License-LGPL_2.1-blue) ![Unittests](https://github.com/DSeiferth/PoreAnalyser/actions/workflows/python-package.yml/badge.svg) ![package](https://github.com/DSeiferth/PoreAnalyser/actions/workflows/python-publish.yml/badge.svg) ![Docker](https://github.com/DSeiferth/PoreAnalyser/actions/workflows/docker-publish.yml/badge.svg)

[Try out this software without installing anything](https://poreanalyser.bioch.ox.ac.uk/)

[Read the accompanying preprint](https://doi.org/10.1101/2024.04.18.589791)

Recent advances in structural biology have led to a growing number of ion channel structures featuring heteromeric subunit assembly, exemplified by synaptic Glycine receptors ([GlyRs](https://www.nature.com/articles/s41467-023-37106-7)) and α4β2 nicotinic receptors. These structures exhibit inherent pore asymmetry, which has raised questions about the role of asymmetry in ion channel function.  Furthermore, molecular dynamics simulations performed on symmetrical homomeric channels often lead to thermal distortion that means conformations of the resulting ensemble are also asymmetrical. We introduce an algorithm that employs ellipsoidal probe particles, enabling a more comprehensive characterization of pore asymmetries. A constriction is more asymmetric for a larger difference between the smaller and larger radius of the ellipsoidal probe particle. 

#### Existing tools for pore pathfinding
- [HOLE](https://www.holeprogram.org/) uses Monte Carlo simulated annealing procedure to find the best route for a sphere with variable radius to squeeze through the channel.
- The Channel Annotation Package [CHAP](https://github.com/channotation/chap) combines  calculations of the pore radius, the hydrophobicity of a pore and water density in the pore to predict hydrophobic gates in ion channels.
- Other tools, such as MOLEonline and CAVER, do not use a probe based algorithm for path finding. Cavities are identified using Voronoi diagrams and molecular surfaces.

#### What does this package add?
- Adding new features to pore-path-finding tools to capture pore asymmetry.
- Capture pore asymmetry.
  - Asymmetry of crystal/cryoEM structures due to heterogeneous subunit composition.
  - From crystal structure broken in simulations.
- Making software tools accessible to the community via an interactive web-service. No installation needed when using the web-page. For python users, we publish an easy-to-install python package. 


# Path finding with ellipsoidal probe particle

1. Align principal axis to z-axis
2. HOLE analysis with spherical probe particle.
3. Load HOLE output file with positions and radii of probes.
4. Loop through all spherical probe particles: 
    a) Ellipsoid initialized with spherical probe particle parameters from HOLE output. 
    b) First Nelder-Mead 4-dim optimization to insert ellipsoid with smaller bounds for parameters [x, y, r1, θ ]. 
    c) Second optimization with larger boundaries for parameters to further increase ellipsoid. The loop takes around 60s to complete...
5. Plot pathway and render pore surface. 

### Installation
PoreAnalyser may be installed as the latest release from PyPI ( pip install PoreAnalyser ) or in the development version from this github repository. 
Detailed [installation instructions](https://porefinding.readthedocs.io/en/latest/usage.html#installation) can be found in the documentation.
After having installed PoreAnalyser locally, you can run the streamlit app yourself: streamlit run app.py

### Links to documentation
You can either upload your proteins of interest to the [webservice](https://poreanalyser.bioch.ox.ac.uk/) hosted on a webserver of the Department of Biochemistry (University of Oxford)
or you can [install](https://porefinding.readthedocs.io/en/latest/usage.html#installation) the PoreAnalyser python package on your machine. 
If you decide to use the webservice, you can download all output files and visualisation scripts to produce high quality figures. 
More information about the [output files](https://porefinding.readthedocs.io/en/latest/webservice.html) can be found in the documentation. 

To render 3d representations of the pore surface, you can use a variety of software ranging from py3Dmol, VMD to UCSF Chimera.
See [Visualisation tools](https://porefinding.readthedocs.io/en/latest/visualisation.html).

### Publication
[Seiferth and Biggin, 2024](https://doi.org/10.1016/j.bpj.2024.07.010)

### Preprint
![qrcode](https://connect.biorxiv.org/qr/qr_img.php?id=2024.04.18.589791)

Find out more about the influence of pore shape on conductance and permeation in our [preprint](https://doi.org/10.1101/2024.04.18.589791)

## Acknowledgements
* Rocco Meli for pointing out streamlit (and hugginface)
* SBCB community for discussion

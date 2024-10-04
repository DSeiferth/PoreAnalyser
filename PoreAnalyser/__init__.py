"""
Library for Pathfinding with an ellipsoidal probe particle
"""
# Import version info
#from .version_info import VERSION_INT, VERSION  # noqa

# Import main classes
from .poreanalyser import PoreAnalysis
#import sys
#from .hole_analysis import hole_analysis 
#sys.path.append('ProbeParticleEllipsoid/')
#from ellipsoid_optimisation import ellipsoid_optimisation as ellipsoid_analysis 
#from visualization import write_pdb_with_pore_surface, plt_ellipsoid_pathway, pathway_visu, st_write_ellipsoid, write_pdb_with_ellipsoid_surface, example_xy_plane, compare_volume, render_visu
#from download_files import download_output, download_Ellipsoid_output
#from ellipsoid_optimisation import ellipsoid_pathway

__version__ = "0.0.7"
__author__ = 'David Seiferth'

from .hole_analysis import hole_analysis

import numpy as np

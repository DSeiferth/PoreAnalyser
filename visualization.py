import MDAnalysis
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import py3Dmol

def write_pdb_with_pore_surface(path='', name='', end_radius=15, num_circle = 24):
    conf = path + name[:-4] + '.sph'
    top = conf
    sph = MDAnalysis.Universe(top, conf, topology_format='pdb', format='pdb') # tpr_resid_from_one=True
    radii = sph.atoms.occupancies
    resids = sph.atoms.resids
    sel = sph.select_atoms('resid '+str(min(resids[radii < end_radius ])) +':'+ str(max(resids[radii < end_radius ])) )
    n_residues = len(sel)
    n_atoms = num_circle * len(sel)
    ### create resindex list ###
    resindices = np.repeat(range(num_circle), len(sel))
    assert len(resindices) == n_atoms
    ### all water molecules belong to 1 segment ###
    segindices = [0] * n_residues
    ### create universe ###
    sol2 = MDAnalysis.Universe.empty(n_atoms,
                                n_residues=n_residues,
                                atom_resindex=resindices,
                                residue_segindex=segindices,
                                trajectory=True) # necessary for adding coordinates
    sol2.add_TopologyAttr('name', ['point']*n_residues*num_circle)
    sol2.add_TopologyAttr('resname', ['Pathway']*n_residues)
    sol2.add_TopologyAttr('resid', list(range(1, n_residues+1)))
    coordinates = []
    for count, atom in enumerate(sel):
        probe1 = sph.select_atoms('resname SPH and resid '+str(atom.resid))
        r = radii[np.where(resids==atom.resid)[0][0] ]
        for i in range(0,num_circle):
            p = probe1.positions[0] + r * np.array([np.cos(2*np.pi*i/num_circle), np.sin(2*np.pi*i/num_circle),0])
            coordinates.append(p)
    coord_array = np.array(coordinates)
    #assert coord_array.shape == (n_atoms, 1)
    sol2.atoms.positions = coord_array
    sel = sol2.select_atoms('name *')
    sel.write(path + name + '_circle.pdb')

def plt_ellipsoid_pathway(df_res, f_size=22, title='', end_radius=15):
    z = df_res['z']
    a = df_res['a']
    b = df_res['b']
    #theta = df_res['theta']
    #print(min(z), max(z))
    fig, ax = plt.subplots()
    ax.set_ylim([0, end_radius + 5])
    plt.plot(z, a, label='Larger radius of ellipsoid')
    plt.plot(z, b, label='Smaller radius of ellipsoid\n=radius of spherical probe')
    plt.ylabel(r"Pathway radii ($\AA$)", fontsize=f_size)
    plt.xlabel(r"z-coordinate ($\AA$)", fontsize=f_size)
    title_str = 'Pathfinding with ellipsoidal probe'
    if title!='':
        title_str = title_str + '\n' + title
    plt.title(title_str, fontsize=f_size)

    ax.tick_params(axis='both', which='major', labelsize=f_size)
    ax.legend(prop={'size': 12}) # loc='upper center'
    fig.tight_layout()
    #fig.savefig('GlyR/asymmetric_channel_ellipse_pathway.png', bbox_inches='tight')  
    #plt.show()
    return fig

def pathway_visu(path, name):
    with open(path+name) as ifile:
            system = "".join([x for x in ifile])
    xyzview = py3Dmol.view(height=800, width=800,) 
    xyzview.addModelsAsFrames(system)
    xyzview.setStyle({'model': -1}, {"cartoon": {'color': 'spectrum'}})
    with open(path+name + '_circle.pdb') as ifile:
        sph2 = "".join([x for x in ifile])
    xyzview.addModelsAsFrames(sph2)
    xyzview.addSurface(py3Dmol.SES,{'opacity':0.9,'color':'lightblue'}, {'model': -1})
    xyzview.zoomTo()
    return xyzview
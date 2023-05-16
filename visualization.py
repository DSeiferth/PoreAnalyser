import MDAnalysis
import numpy as np

def write_pdb_with_pore_surface(path='', name='', end_radius=15, num_circle = 24):
    conf = path + name + '.sph'
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

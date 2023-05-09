import numpy as np
import MDAnalysis
import nglview as nv
import matplotlib.pyplot as plt
import pandas as pd
from MDAnalysis.analysis import hole2


hole_exe = 'hole2//hole'
sph_proc = 'hole2/sph_process'
f_size = 18

import warnings; warnings.simplefilter('ignore')

from MDAnalysis.analysis import align
def align_and_write(universe, reference, names, out_path,
                   TMD_lower=0, TMD_higher=0):
    '''
    names=[names[i] ,names[0]]: 2nd name refers to reference
    '''
    try:
    #bla=0
    #if bla==0:
        #### all c-alphas around COM +/- 10 A ###

        CA_2align = universe.select_atoms('name CA ') #and chainID A and segid seg_0_PROA
        CA_ref = reference.select_atoms('name CA ') # and chainID A
        print('calphas model 1, reference (model 0)', len(CA_2align), len(CA_ref))
        align_IDs = np.unique(CA_2align.resids)
        ref_IDs = np.unique(CA_ref.resids)
        mask1 = np.isin(align_IDs, ref_IDs)
        mask2 = np.isin(ref_IDs, align_IDs)
        chains = '(chainID A )' # for glycine
        print(np.unique(align_IDs[mask1] == ref_IDs[mask2]),
              'number of common c-alphas', len(ref_IDs[mask2]),  len(align_IDs[mask1]))
        select = 'name CA and ('
        for i in ref_IDs[mask2]:
            select = select + 'resid '+str(i) + ' or '
        select = select + 'resid '+str(ref_IDs[mask2][0])+' ) ' #'and ' + chains
        #print(select)
        sel1 = universe.select_atoms(select)
        sel2 = reference.select_atoms(select)
        print('selection in reference', len(sel2), 'selection in structure to be aligned',
              len(sel1))

        rmsds = align.alignto(universe.atoms,  # mobile
                          reference.atoms,  # reference
                          select=select, # selection to operate on
                          #match_atoms=True # whether to match atoms
                         )
        print('returns (old_rmsd, new_rmsd)',rmsds)
        #mer = MDAnalysis.Merge(universe.atoms, reference.atoms).select_atoms('protein') #name *
        #mer.write(out_path+'merged_'+name+'.pdb')


        sel1 = 'protein and prop  z<'+str(TMD_higher)
        sel2 = ' and prop z>'+str(TMD_lower)

        aligned = MDAnalysis.Merge(universe.atoms).select_atoms('protein')
        if TMD_higher!=0 and TMD_lower!=0 :
            TM = aligned.select_atoms(sel1+sel2)
            TM.write(out_path+  names[0]+'.pdb')
        else:
            aligned.atoms.write(out_path+  names[0]+'.pdb')

        ref = MDAnalysis.Merge(reference.atoms).select_atoms('protein')
        if TMD_higher!=0 and TMD_lower!=0 :
            TM = ref.select_atoms(sel1+sel2)
            TM.write(out_path+  names[1]+'.pdb')
        else:
            ref.atoms.write(out_path+names[1]+'.pdb')
        print('aligned', names[0], 'with', names[1], 'as reference')
    except:
        print('ERROR: Alignment did not work')
        print('usually the following error:')
        print('electionError: Reference and trajectory atom selections do not contain the same number of atoms: ')
        print('Printing the orignal models as pdb')
        universe.atoms.write(out_path+  names[0]+'.pdb')
        reference.atoms.write(out_path+names[1]+'.pdb')
        return names, [0,0,0]

    return names, rmsds[1]

def visualise_aligned(names, path):
    #conf =  path + 'aligned_' + names[0] + '.pdb'
    conf =  path +  names[0] + '.pdb'
    top = conf
    u = MDAnalysis.Universe(top, conf, topology_format='pdb', format='pdb')

    #conf =  path + 'reference_' + names[1] + '.pdb'
    conf =  path +  names[1] + '.pdb'
    top = conf
    ref = MDAnalysis.Universe(top, conf, topology_format='pdb', format='pdb')

    mer = MDAnalysis.Merge(u.atoms, ref.atoms).select_atoms('protein') #name *
    v = nv.show_mdanalysis(mer.atoms)
    return v

def hole_analysis(name, path, typ='pdb', end_radius=20, sel='protein'):
    tmpdir = path #+ 'tmpdir/'
    conf = path + name
    top = conf
    sys = MDAnalysis.Universe(top, conf, topology_format=typ, format=typ)
    ha2 = hole2.HoleAnalysis(
                                sys, select=sel,
                                cpoint='center_of_geometry',
                                executable=hole_exe,
                                tmpdir=tmpdir,
				                sph_process=sph_proc,
				                #sphpdb_file=path+name+'.sph',
                                #sphpdb=path+name+'.sph',
                                end_radius=end_radius,
				                #keep_files=True
				)
    ha2.run(random_seed=31415)
    #print(ha2.results.profiles)
    #print(ha2.results.sphpdbs)
    #print(ha2.results.outfiles)
    try:
        #http://minium.com.au/UserGuide/stable/examples/analysis/polymers_and_membranes/hole.html
        ha2.create_vmd_surface(filename=path+name+'.vmd') #source hole1.vmd
        print('vmd surface created for', name)
    except:
        print('ERROR: No vmd surface created for', name)

    radii2, edges2 = ha2.bin_radii(bins=100, range=None)
    means2, edges2 = ha2.histogram_radii(bins=100, range=None,
                                          aggregator=np.mean)
    midpoints2 = 0.5*(edges2[1:]+edges2[:-1])
    return midpoints2, means2

def visualisation(name , path='/biggin/b198/orie4254/Documents/CHAP/', out=0,
                  num_circle=24, TMD_higher=0,  TMD_lower=0, end_radius=15):
    conf = path + name + '.sph'
    top = conf
    sph = MDAnalysis.Universe(top, conf, topology_format='pdb', format='pdb') # tpr_resid_from_one=True
    radii = sph.atoms.occupancies
    resids = sph.atoms.resids
    print([min(resids[radii < end_radius ]), max(resids[radii < end_radius ])])

    sel = sph.select_atoms('resid '+str(min(resids[radii < end_radius ])) +':'+ str(max(resids[radii < end_radius ])) )

    conf =  path + name + '.pdb' # '.pdb1'
    top = conf
    u = MDAnalysis.Universe(top, conf, topology_format='pdb', format='pdb')

    ### create new universe with center ###
    n_residues = len(sel)
    n_atoms = 1 * len(sel)
    # create resindex list
    resindices = np.repeat(range(len(sel)), 1)
    if out: print(len(resindices) , n_atoms )
    assert len(resindices) == n_atoms
    if out: print("resindices:", resindices[:10])
    # all water molecules belong to 1 segment
    segindices = [0] * n_residues
    if out: print("segindices:", segindices[:10])
    # create universe
    sol = MDAnalysis.Universe.empty(n_atoms,
                             n_residues=n_residues,
                             atom_resindex=resindices,
                             residue_segindex=segindices,
                             trajectory=True) # necessary for adding coordinates
    sol.add_TopologyAttr('name', ['center']*n_residues)
    if out: print(sol.atoms.names)
    sol.add_TopologyAttr('resname', ['Pathway']*n_residues)
    if out: print(sol.atoms.resnames)
    sol.add_TopologyAttr('resid', list(range(1, n_residues+1)))
    if out: print(sol.atoms.resids)
    coordinates = []
    for count, atom in enumerate(sel):
        probe1 = sph.select_atoms('resname SPH and resid '+str(atom.resid))
        if out: print(probe1.positions[0], 'radius', radii[np.where(resids==atom.resid)[0][0] ])
        coordinates.append(probe1.positions[0])
    if out: print(coordinates[:10])
    coord_array = np.array(coordinates)
    if out: print(coord_array.shape, n_atoms)
    sol.atoms.positions = coord_array

    mer = MDAnalysis.Merge(u.atoms, sol.atoms).select_atoms('name *')

    ### create new universe with circle ###
    # num_circle = 24
    n_residues = len(sel)
    n_atoms = num_circle * len(sel)
    # create resindex list
    resindices = np.repeat(range(num_circle), len(sel))
    if out: print(len(resindices) , n_atoms )
    assert len(resindices) == n_atoms
    if out: print("resindices:", resindices[:10])
    # all water molecules belong to 1 segment
    segindices = [0] * n_residues
    if out: print("segindices:", segindices[:10])
    # create universe
    sol2 = MDAnalysis.Universe.empty(n_atoms,
                             n_residues=n_residues,
                             atom_resindex=resindices,
                             residue_segindex=segindices,
                             trajectory=True) # necessary for adding coordinates
    sol2.add_TopologyAttr('name', ['point']*n_residues*num_circle)
    if out: print(sol.atoms.names)
    sol2.add_TopologyAttr('resname', ['Pathway']*n_residues)
    if out: print(sol.atoms.resnames)
    sol2.add_TopologyAttr('resid', list(range(1, n_residues+1)))
    if out: print(sol.atoms.resids)
    coordinates = []
    for count, atom in enumerate(sel):
        probe1 = sph.select_atoms('resname SPH and resid '+str(atom.resid))
        r = radii[np.where(resids==atom.resid)[0][0] ]
        if out: print(probe1.positions[0], 'radius', r)
        for i in range(0,num_circle):
            p = probe1.positions[0] + r * np.array([np.cos(2*np.pi*i/num_circle), np.sin(2*np.pi*i/num_circle),0])
            coordinates.append(p)
    if out: print(coordinates[:10])
    coord_array = np.array(coordinates)
    if out: print(coord_array.shape, n_atoms)
    #assert coord_array.shape == (n_atoms, 1)
    sol2.atoms.positions = coord_array

    mer2 = MDAnalysis.Merge(mer.atoms, sol2.atoms).select_atoms('name *')
    if TMD_higher!=0 and TMD_lower!=0:
        sel1 = 'prop  z<'+str(TMD_higher)
        sel2 = ' and prop z>'+str(TMD_lower)
        TMD = mer2.select_atoms(sel1+sel2)
        v = nv.show_mdanalysis(TMD)
    else:
        v = nv.show_mdanalysis(mer2.atoms)
    return v

def analysis(names,labels, path='/biggin/b198/orie4254/Documents/CHAP/', end_radius=15,
            TMD_lower=0, TMD_higher=0, save='', title='', typ='gro', sel='protein', legend_outside=False):
    ### Downloading pdb models ###
    #for name in names:
    #    ! wget 'https://files.rcsb.org/download/'$name'.pdb1'
    #    ! cp $name'.pdb1' $outpath
    #    ! cp $name'.pdb1' $path

    ### align model ###
    conf =  path + names[0] +'.'+ typ
    print('conf', conf)
    #! ls $conf
    top = conf
    ref = MDAnalysis.Universe(top, conf, topology_format=typ, format=typ)

    if len(names) == 1 :
        print('ref', len(ref.atoms),ref)
        if TMD_lower!=0 and TMD_higher!=0:
            sel1 = 'protein and prop  z<'+str(TMD_higher)
            sel2 = ' and prop z>'+str(TMD_lower)
            TM = ref.select_atoms(sel1+sel2)
            TM.write(path+names[0]+'.pdb')
        else:
            ref.atoms.write(path+names[0]+'.pdb')
    else:
        for i in range(1, len(names)):
            conf =  path + names[i] + '.'+typ
            top = conf
            u = MDAnalysis.Universe(top, conf, topology_format=typ, format=typ)
            print('align', names[i], 'with', names[0], 'as reference')
            print('TMD_lower', TMD_lower, 'TMD_higher', TMD_higher)
            align_and_write(universe=u, reference=ref, names=[names[i] ,names[0]],
                            out_path=path, TMD_lower=TMD_lower, TMD_higher=TMD_higher)
            #aligned = visualise_aligned(names, path)

    #! cp $outpath$name'.pdb' $path
    ### hole analysis ###
    aligned_path = path
    fig, ax = plt.subplots()
    plt.title(r'HOLE ' + title, fontsize=f_size)
    colors = ['black', 'blue', 'orange', 'purple','green','red','gray', 'brown',
              'cyan', 'violet', 'olive', 'peru', 'slategray',
    ]
    for i in range(len(names)):
        midpoints, means = hole_analysis(name=names[i]+'.pdb', path=aligned_path,
                                                  typ='pdb', end_radius=end_radius, sel=sel)
        rmin = min(means)
        ax.plot(midpoints, means, color=colors[i],
                label=labels[i] + r' R$_{min}$ = '+str(round(rmin,2)) + r'$\AA$')
        if i==0:
                    d = {labels[i]+' z [A]': midpoints,
                         labels[i]+' Radius [A]': means,
                        }
                    df = pd.DataFrame(data=d)
        else:
            df[labels[i]+' z [A]'] = midpoints
            df[labels[i]+' Radius [A]'] = means
    plt.ylabel(r"HOLE radius $R$ ($\AA$)", fontsize=f_size)
    plt.xlabel(r"Pore coordinate $\zeta$ ($\AA$)", fontsize=f_size)

    y1 = 0.0
    y2 = 20
    ax.set_ylim([y1,y2])
    xlim = ax.get_xlim()
    ax.plot(xlim, [1.15,1.15], '--',color='red') # label=r'r < 1.15 $\AA$'
    ax.plot(xlim, [2.3,2.3], '--',color='green') # label=r'r < 2.30 $\AA$'

    ax.tick_params(axis='both', which='major', labelsize=f_size)
    if legend_outside:
        plt.legend(prop={'size': 10},  loc='upper center', bbox_to_anchor=(1.4,1.0), frameon=False)
    else:
        ax.legend(prop={'size': 10}, loc='upper left') # loc='upper center'
        fig.tight_layout()
    fig.savefig(path + save[:-1] +'HOLE_pathwayprofile.png', bbox_inches='tight')
    #df.to_csv(outpath + '/Pathway_HOLE_comparison_pdb.csv', sep=',', index=False, header=True)
    plt.show()
    ### visualise pathway ###
    pathways = []
    for count, name in enumerate(names):
        try:
            print(name, '### visualise pathway ###')
            if count<2:
                #! echo $path/{name}.pdb
                #!cp $path/{name}.pdb .
                print('Not copying ...')
            else:
                print('last one..., hard coded')
                path = path + '/PROD/'
                #!cp $path/prod.pdb .
                name = 'prod'
            ha2 =  hole2.hole(
                    pdbfile=name +'.pdb',
                    #cpoint='center_of_geometry',
                    executable=hole_exe,
                    #tmpdir=path,
                    #sph_process=sph_proc,
                    sphpdb_file=name+'.sph',
                    end_radius=end_radius,
                    keep_files=True
            )
            #pathway = visualisation(name = names[0], path=path, out=0, end_radius=end_radius )
            #pathways.append(pathway)
            #! mv {name}.sph $path
            #! rm {name}.pdb
        except:
            print('ERROR with', name, 'no SPH file generated')
    return pathways

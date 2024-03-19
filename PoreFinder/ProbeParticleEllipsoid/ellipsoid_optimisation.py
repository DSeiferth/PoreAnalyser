import numpy as np
import MDAnalysis
import matplotlib.pyplot as plt
import pandas as pd

from ellipse_lib import atom, ellipse, distance_ellipse, dist_ellipse_vdwSphere, assign_radius

from scipy.optimize import minimize
from scipy.optimize import NonlinearConstraint

import multiprocessing 
import time
import os


def penalty_overlap_4dim(x, args):
    '''
    Parameters:
    - x (list): A 4-dimensional vector representing the parameters of the ellipse.
        - x[0]: Radius to expand.
        - x[1]: Angle of rotation.
        - x[2]: x-coordinate of the center.
        - x[3]: y-coordinate of the center.
    - args (tuple): A tuple of arguments.
        - args[0]: Constant radius used in the creation of the ellipse.
        - args[1]: List of Sphere objects representing the surroundings of the ellipse.
        - args[2] (optional): Boolean flag indicating whether to stop the loop upon detecting an overlap. Default is True.
            - if True: return a high penalty if any overlap is detected.
            - if False: return the absolute value of the minimum distance between the ellipse and the spheres if we have an overlap.
    Returns:
    float: The penalty score based on the conditions specified.

    Notes:
    This function creates an Ellipse object using the input parameters and evaluates the overlap penalty with a set of vdW spheres.
    The penalty is calculated based on the minimum distance between the ellipse and the spheres. If the stop_loop flag is set to True,
    the function returns a high (positive) penalty if any overlap is detected. 
    '''
    e1 = ellipse(a=x[0], b=args[0], theta=x[1], cx=x[2], cy=x[3])
    a_vec = args[1]
    if len(args)!=3:
        stop_loop = True 
    else:
        stop_loop = args[2]
    #print('ellipse: a', e1.a, 'b', e1.b, 'theta', e1.theta, 'center', e1.cx, e1.cy)
    ### spheres ###        
    dist = np.zeros(len(a_vec))
    for i, atom in enumerate(a_vec):
        dist[i] = dist_ellipse_vdwSphere(ellipse=e1, sphere=atom)
        if stop_loop: ### make overlap costly and avoid calculating minimum ###
            if dist[i] < 0: ### overlap ###
                return 1e9 #abs(dist[i]) 
        else: ### calculate minimum distance ###
            if i == 0 :
                min_dist = dist[i]
            elif dist[i]< min_dist:
                min_dist = dist[i]
    #min_dist = min(dist)
    #print('min_dist', min_dist)
    if stop_loop:
        return - x[0]   ### negative score == larger radius ###
    else:
        if min_dist<0: ### overlap ###
            #print('min_dist<0', abs(min_dist))
            score = abs(min_dist) ### positive penalty ###
        else:
            score = - x[0]   ### negative score == larger radius ###
        #print(score)
        return score


def neighbor_vec(universe, probe, probe1, n_xy_fac, out=0, call=0, pathway_sel='protein'):
    """
    Identify neighboring atoms within a specified spatial range around a probe in a molecular system.

    Parameters:
    - universe (MDAnalysis.universe): The molecular dynamics universe representing the system.
    - probe (Sphere): A Sphere object representing the central probe.
    - probe1 (MDAnalysis AtomGroup): AtomGroup representing the probe atoms from initial HOLE run (resname SPH).
    - n_xy_fac (float): The factor used to determine the spatial range in the xy-plane around the probe.
    - out (int, optional): An integer flag indicating whether to print debugging information. Default is 0 (no printing).
    - call (int, optional): An integer indicating the recursive call level. Default is 0.
    - pathway_sel (str, optional): The selection string to identify the pathway in the molecular system. Default is 'protein'.

    Returns:
    tuple: A tuple containing:
        - list: A list of Atom objects representing neighboring atoms within the specified range.
        - list: A list of strings representing labels for the neighboring atoms.
        - float: The spatial range (n_xy) used for the neighbor search.

    Notes:
    This function identifies neighboring atoms around a central probe within the specified spatial range.
    The spatial range is determined in the xy-plane based on the probe's radius and the given factor (n_xy_fac).
    The selection is performed within the specified pathway (default is 'protein').

    If the number of neighboring atoms is too large (greater than 150) and the recursive call level is less than 4,
    the function reduces the spatial range (n_xy_fac) and makes a recursive call to refine the selection.
    If the number of neighboring atoms is too small (less than 30) and the recursive call level is less than 4,
    the function increases the spatial range and makes a recursive call to expand the selection.
    """ 
    protein = universe.select_atoms(pathway_sel)
    if out: print('number of atoms to find pathway', len(protein))
    start_neighbors = time.time()
    n_xy = n_xy_fac*probe.r 
    if out: print('n_xy', n_xy, 'probe.r', probe.r)
    n_z = 2
    string_x = 'prop x >'+str(probe.x-n_xy) + ' and prop x<'+str(probe.x+n_xy)
    string_y = 'prop y >'+str(probe.y-n_xy) + ' and prop y<'+str(probe.y+n_xy)
    string_z = 'prop z >'+str(probe.z-n_z) + ' and prop z<'+str(probe.z+n_z)
    select = string_x + ' and ' + string_y + ' and ' + string_z 
    if out: print(select)
    neighbors = universe.select_atoms('group id1 and '+ select, id1=protein)
    end = time.time()
    if out: print("TIME(start_neighbors)=",end - start_neighbors)
    if out:  print('number of neighbors' ,len(neighbors))
    a_vec = []
    neighbour_labels = []
    neighbour_string = ''
    neighbour_resnum_resid = []
    z_slice = probe1.positions[0][2]
    for n in neighbors:
        #print(n.type,  assign_radius(n.type))
        R0 = assign_radius(n.type)
        z0 = n.position[2]
        R_projected = np.sqrt(R0*R0 - (z0-z_slice)*(z0-z_slice))
        if R_projected > 0:
            a_vec.append(atom(n.position[0],n.position[1], r=R_projected ))
            neighbour_labels.append(n.resname + ' ' + str(n.resid)) # + ' ' + n.type
            neighbour_resnum_resid.append(n.resname + ' ' + str(n.resid))
            neighbour_string = neighbour_string + n.resname + ' ' + str(n.resid) + ' ' + n.type + '\n'
    if out: print('number of neighbors with R>0' ,len(a_vec), 'n_xy_fac', n_xy_fac)
    if len(a_vec) > 150 and call < 4:
        if out: print('DECREASE n_xy_fac')
        call += 1
        return neighbor_vec(universe, probe, probe1, 0.75*n_xy_fac,
                                              call=call, out=out, pathway_sel=pathway_sel)
        print('number of neighbors with R>0' ,len(a_vec), 'n_xy_fac', 0.75*n_xy_fac)
    elif len(a_vec)<30 and call<4:
        if out: 
            print('INCREASE n_xy_fac')
        call += 1 
        return neighbor_vec(universe, probe, probe1, 1.25*n_xy_fac,
                                              call=call, out=out, pathway_sel=pathway_sel)
        if out: print('number of neighbors with R>0' ,len(a_vec), 'n_xy_fac', 1.25*n_xy_fac)
    elif len(a_vec)<30 and call>4:
        n_xy_fac = 1.5
        n_xy = n_xy_fac*probe.r 
        print('len(a_vec)', len(a_vec), 'call', call, 'set back to n_xy = n_xy_fac*probe.r ', n_xy , n_xy_fac,probe.r )
        return a_vec, neighbour_labels, n_xy
    return a_vec, neighbour_labels, n_xy


def insert_ellipse(index, dataframe, universe, 
                   out=0, plt_path='', rmax=50, 
                   show=0, label=0, n_xy_fac=1.6,
                   timing=0,
                   f_size=22,
                   pathway_sel='protein',
                   opt_method='nelder-mead',
                  ):
    """
    Inserts an ellipse into a plot with specified parameters.

    Parameters:
    index : int
    The index to locate the ellipse within the dataframe.
    dataframe : pandas.DataFrame
    A dataframe that contains the x, y, z, r and resid columns.
    universe : MDAnalysis.Universe
    A molecular dynamics universe that contains the coordinates of atoms.
    out : int, default 0
    A flag to control output print statements.
    plt_path : str, default ''
    A path to save the plot to.
    rmax : float, default 50
    The maximum radius for the ellipsoid (to be deleted, set rmax to n_xy)

    show : int, default 0
    A flag to control whether the plot is displayed.
    label : int, default 0
    A flag to control whether labels are displayed on the plot.
    n_xy_fac: float, default 1.6
    when calculating the neighbor vector, atoms within hole_radius*n_xy_fac are considered
    
    timing: int, default 0
    A flag to control whether timings for sub tasts shoudl be printed out
    Returns:
    None
    """

    start = time.time()
    rID = dataframe['resid'].loc[index]
    probe1 = universe.select_atoms('resid '+str(rID) + ' and resname SPH')
    if out: print('probe1', len(probe1))
    if out: 
        print('resid of probe', rID, probe1.positions, probe1.resnames)
        print('dataframe.loc[index]',dataframe.loc[index]) # iloc != loc
    if dataframe['z'].loc[index] - probe1.positions[0][2] > 0.01:
        print('ERROR in index')
        return -1, -1

    probe = atom(dataframe['x'].loc[index],dataframe['y'].loc[index],
                 z=dataframe['z'].loc[index],
                 r=dataframe['r'].loc[index])

    z_slice = probe1.positions[0][2]
    #Radius = min( 10*dataframe['r'].loc[index], R_N)
    #initial_radius = dataframe['r'].loc[index]

    ### neighbor vector ###
    a_vec, neighbour_labels, n_xy = neighbor_vec(universe, probe, probe1, n_xy_fac, out=out, pathway_sel=pathway_sel)
    assert len(a_vec)>2, "in function 'insert_ellipse': len(a_vec)<2="+str(len(a_vec))

    ### does HOLE have overlap already ????
    r = probe.r 
    while penalty_overlap_4dim([r, 0, probe.x, probe.y], [r, a_vec]) > 0:
        print('PROBE overlap', r, '->', 0.9*r, 'z', probe.z )
        r = 0.95*r
    probe.r = r

    fig, ax = plt.subplots()
    plt.gca().set_aspect('equal', adjustable='box')
    
    ### probe ###
    p0 = ellipse(a=probe.r, b=probe.r, theta=0, cx=probe.x, cy=probe.y)
    x0, y0 = p0.draw()
    ax.plot(x0, y0, color='blue')
    ax.plot(probe.x,probe.y, '-x', color='blue')
    ### vdw particles ###
    start_neighbors_loop = time.time()
    prev = neighbour_labels[0]
    for count, a in enumerate(a_vec):
        vdw = ellipse(a=a.r, b=a.r, theta=0, cx=a.x, cy=a.y)
        x0, y0 = vdw.draw()
        ax.plot(x0, y0, color='black')
        if label:
            if count == 0 or neighbour_labels[count]!=prev:
                center2vdw = np.array([-probe.x+a.x, -probe.y+a.y])
                center2vdw = a.r/np.linalg.norm(center2vdw) * center2vdw
                ax.text(a.x + center2vdw[0],a.y+center2vdw[1],  
                        neighbour_labels[count], fontsize=10, 
                        color='blue'
                        #verticalalignment='top'
                       )
            prev = neighbour_labels[count]
        else:
            ax.plot(a.x,a.y, '-x', color='black')
    end = time.time()
    if timing: print("TIME(neighbors_loop)=",end - start_neighbors_loop)

    ### Test if closest atom has overlap ###
    #ind, dist = find_nearest_atom(probe.x,probe.y, a_vec)
    #if out: print('atom ', ind, 'has minimal distance to center of probe', dist, 
    #      'probe radius=', probe.r)

    ### Nelder mead optimisation ###
    def optimisation_ellipsoid(ax, p0, a_vec, dx, rad_fac, pt, opt, bnds, timing, opt_method, col='green' ):
        start_opt1 = time.time()
        result = minimize(penalty_overlap_4dim, pt, 
                        args = [rad_fac*p0.b, a_vec], ### set smaller radius to 90% to allow for more flexibility
                        method=opt_method,
                        bounds=bnds,
                        #constraints=constr,
                        options = opt
                        #tol=1e-10
                        )
        end = time.time()
        if timing: print("TIME(start_opt1)=",end - start_opt1)

        sol = result['x']
        if out: print(result) #result
        p1 = ellipse(a=sol[0], b=p0.b, theta=sol[1],
                     cx=sol[2], cy=sol[3], 
                     #cx= p0.cx, cy= p0.cy
                     )
        x1, y1 = p1.draw()
        ax.plot(p1.cx, p1.cy, '-x', color=col)
        if result['fun'] > 0:
            print('ERROR: overlap', result['fun'])
            ax.plot(x1, y1, '--', color=col)
            return p0 #-1, -1
        else:
            ax.plot(x1, y1, color=col)
            if out: print('4d optimisation small variation of center (r, theta, x, y)', sol)
        return p1

    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.NonlinearConstraint.html#scipy.optimize.NonlinearConstraint
    ### 4d with small variation ###
    dx = 0.1*p0.a ### boundary for center of ellipse ###
    #rad_fac = 0.8 # set radius to 95% of HOLE output to allow for more flexibility 

    opt = {}
    if 'nelder' in opt_method:
        #'xatol': 0.00001, # Absolute error in xopt between iterations that is acceptable for convergence.
        #'fatol': 0.00001, # Absolute error in func(xopt) between iterations that is acceptable for convergence.
        rad_fac = [0.9, 0.99]
        dr = 0.15 # for second optimisation
        opt['maxiter'] = 800 # Maximum allowed number of iterations and function evaluations. Will default to N*200, where N is the number of variables
        opt['adaptive'] = True # Adapt algorithm parameters to dimensionality of problem. Useful for high-dimensional minimization 
        opt_init = 'initial_simplex'
    elif 'powell' in opt_method:
        rad_fac = [0.8, 0.9]
        dr = 1
        opt_init = 'direc'
        opt['maxiter'] = 4*1000 # Will default to N*1000, where N is the number of variables
    else:
        opt_init = ''
    pt = [rad_fac[0]*p0.a, 0,  p0.cx, p0.cy] # starting point =  HOLE output
    # (p0.a, rmax)
    bnds = ((0, n_xy), (-np.pi, np.pi), (p0.cx-dx, p0.cx+dx), (p0.cy-dx, p0.cy+dx))

    opt_HOLE = opt.copy()
    opt_HOLE[opt_init] = [pt, [p0.a+0.1, np.pi/4, p0.cx, p0.cy], 
                                [p0.a+0.11, np.pi/2, p0.cx, p0.cy],
                                [p0.a+0.12, np.pi/4, p0.cx+dx/2, p0.cy+dx/2],
                                [p0.a+0.13, np.pi/2, p0.cx-dx/2, p0.cy-dx/2]]
    
    p1_HOLE = optimisation_ellipsoid(ax=ax, p0=p0, a_vec=a_vec, dx=dx, rad_fac=rad_fac[0], 
                                pt=pt, opt=opt_HOLE, bnds=bnds, timing=timing, opt_method=opt_method,
                                 col='green' )
    

    COG = [0,0]
    len_a_vec = len(a_vec)
    for a in a_vec:
        COG[0] += a.x
        COG[1] += a.y
    COG[0] = COG[0]/len_a_vec
    COG[1] = COG[1]/len_a_vec
    print('COG', COG, p0.cx, p0.cy )
    p0_COG = ellipse(a=probe.r, b=probe.r, theta=0, cx=COG[0], cy=COG[1])
    pt_COG = [rad_fac[0]*p0.a, 0,  COG[0], COG[1] ]

    opt_COG = opt.copy()
    opt_COG[opt_init] = [pt_COG, [p0_COG.a+0.13, -np.pi/4, p0_COG.cx, p0_COG.cy], 
                                [p0_COG.a+0.12, np.pi/2, p0_COG.cx, p0_COG.cy],
                                [p0_COG.a+0.11, -np.pi/4, p0_COG.cx+dx/2, p0_COG.cy+dx/2],
                                [p0_COG.a+0.10, np.pi/2, p0_COG.cx-dx/2, p0_COG.cy-dx/2]]

    p1_COG = optimisation_ellipsoid(ax, p0_COG, a_vec, dx, rad_fac[0], pt_COG, opt_COG, bnds, timing, opt_method,
                                 col='yellow' )
    print('p1_HOLE',p1_HOLE.a)
    print('p1_COG', p1_COG.a)
    # comparing COG and HOLE starting point increases ratio from 2.0 to 2.3 for CNT 10x50 ratio 1/2
    # rad_fac 0.9 and 0.95
    if p1_COG.a > p1_HOLE.a:
        p1 = p1_COG
    else:
        p1 = p1_HOLE

    ### 4d optimisation works better without initial simplex ###
    vec_long_axis = np.array([p1.cx-p0.cx, p1.cy-p0.cy]) # TO DO: theta
    ### claculate vector along long axis ###


    dx = 0.5*p1.a # 3 ### boundary for center of ellipse ###
    bnds =  ((p1.a, n_xy), (-np.pi, np.pi), (p0.cx-dx, p0.cx+dx), (p0.cy-dx, p0.cy+dx))
    #rad_fac = 0.9 # set radius to 97% of previous optimisation to allow for more flexibility 
    pt = [rad_fac[1]*p1.a, p1.theta, p1.cx, p1.cy ] #p0.cx, p0.cy -1,-1
    axis1 = [p1.a*np.cos(p1.theta), p1.a*np.sin(p1.theta)]
    axis2 = [-p1.b*np.sin(p1.theta), p1.b*np.cos(p1.theta)]
    opt2 = opt.copy()
    opt2[opt_init] = [pt, 
                             [rad_fac[1]*p1.a, p1.theta, p1.cx+dr/p1.a*axis1[0], p1.cy+dr/p1.a*axis1[1]],
                             [rad_fac[1]*p1.a, p1.theta, p1.cx-dr/p1.a*axis1[0], p1.cy-dr/p1.a*axis1[1]],
                             [rad_fac[1]*p1.a, p1.theta, p1.cx+dr/p1.b*axis2[0], p1.cy+dr/p1.b*axis2[1]],
                             [rad_fac[1]*p1.a, p1.theta, p1.cx-dr/p1.b*axis2[0], p1.cy-dr/p1.b*axis2[1]],
                            ]
    p2 = optimisation_ellipsoid(ax=ax, p0=p1, a_vec=a_vec, dx=dx, rad_fac=rad_fac[1], 
                                pt=pt, opt=opt2, bnds=bnds, timing=timing , opt_method=opt_method,
                                col='brown')

    ax.set_xlabel(r"x ($\AA$)", fontsize=f_size)
    ax.set_ylabel(r"y ($\AA$)", fontsize=f_size)
    ax.tick_params(axis='both', which='major', labelsize=f_size)
    plt.title(r'xy-plane at z = '+str(z_slice)+r' $\AA$', fontsize=f_size)

    ### check whether center moved ###
    dist_prev = np.linalg.norm( np.array([p2.cx,p2.cy]) -  np.array([p1.cx,p1.cy]) )
    if dist_prev > max(p0.b, 7) or p2.a/p1.a > 1.5:
        print('ATTENTION', 'dist_prev > p0.b', dist_prev , p0.b, 'at z=',z_slice, '\n')
        p2 = p1
        #ax.plot(x2, y2, '--',color='brown')
    if len(a_vec) < 30:
        p2 = p1
        print('ATTENTION', 'number of neighbors low', len(a_vec), 'at z=',z_slice, '\n')
        #ax.plot(x2, y2, '--',color='brown')
        if p1.a > 3*p0.a:
            print('ERROR: 1st opimisation p1.a', p1.a, '>3*p0.a, p0.a=',p0.a,  'at z=',z_slice, '\n' )
            return -1, -1
    #else:
        #ax.plot(x2, y2, color='brown')
    
    if label: 
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        plt.subplots_adjust(right=0.35)
        # neighbour_resnum_resid = np.unique(neighbour_resnum_resid)
        string = ''
        # for info in neighbour_resnum_resid:
        #    string = string + info +'\n'
        #ax.text(xlim[1],ylim[0],  neighbour_string, fontsize=8) #verticalalignment='top'
        ax.text(1.015*xlim[1],ylim[0],  string, fontsize=12)
    
    fig.tight_layout()
    if plt_path!='': 
        fig.savefig(plt_path+'/z=' +str(z_slice) +'A.png', bbox_inches='tight')  
    if show: plt.show()
    else: plt.close()
    end = time.time()
    if timing: print("TIME(insert_ellipse)=",end - start, '\n')
    return p2, z_slice


def insert_ellipse_async(index, dataframe, universe, out=0, plt_path='', rmax=50, 
                         show=0, label=0, n_xy_fac=1.6,num_processes=None, timeout=20, f_size=22,
                         pathway_sel='protein', opt_method='nelder-mead'):
    result_queue = multiprocessing.Queue()  # Queue to store the result
    processed_indices = multiprocessing.Manager().list()  # Shared list to track processed indices
    process_times = multiprocessing.Manager().list()  # Shared list to store process times
    
    if num_processes is None:
        num_processes = multiprocessing.cpu_count()
    print('Number of processes used', num_processes)

    def insert_ellipse_worker(index):
        if index in processed_indices:
            return  # Skip already processed indices
        start_time = time.time()  # Start time of the process
        result = insert_ellipse(index, dataframe, universe, out, plt_path, rmax,
                                show, label, n_xy_fac, f_size=f_size, pathway_sel=pathway_sel, opt_method=opt_method)
        elapsed_time = time.time() - start_time  # Elapsed time of the process
        result_queue.put(result)  # Store the result in the queue
        processed_indices.append(index)  # Track the processed index
        process_times.append(elapsed_time)  # Store the process time

    processes = []
    for idx in index:
        if idx in processed_indices:
            continue  # Skip already processed indices
        process = multiprocessing.Process(target=insert_ellipse_worker, args=(idx,))
        process.start()
        process._start_time = time.time()  # Store the process start time
        processes.append(process)
        if len(processes) >= num_processes:
            for p in processes:
                p.join(timeout)  # Wait for the processes to finish or timeout
                if p.is_alive():  # Check if the process is still alive after the timeout
                    p.terminate()  # Terminate the process
                    print('Terminated idx', idx)
            processes = []

    # Wait for the remaining processes to finish
    for process in processes:
        process.join(timeout)
        if process.is_alive():  # Check if the process is still alive after the timeout
            process.terminate()  # Terminate the process
            print('Terminated')

    results = []
    while not result_queue.empty():
        result = result_queue.get()  # Retrieve the result from the queue
        results.append(result)
    times_list = list(process_times)
    print('statistics for times for each process:')
    print('min', min(times_list), 'max', max(times_list))
    print('mean', np.mean(times_list), 'median', np.median(times_list),
          'std', np.std(times_list))
    return results


def ellipsoid_pathway(p, pdb_name, sph_name, 
                      slice_dz=1,
                      parallel=False,
                      end_radius=15,                    
                      num_processes=None, 
                      timeout=20,
                      start_index=50,
                      f_size=22,
                      out=0,
                      n_xy_fac=1.6,
                      pathway_sel='protein',
                      opt_method='nelder-mead',

                     ):
    """Generate ellipsoids to represent the pore path of a biomolecule.

    Given the path to a directory (`p`), a PDB file name (`pdb_name`), and a
    SPH file name (`sph_name`), this function generates ellipsoids to represent
    the pore path of a biomolecule.

    Parameters
    ----------
    p : str
        The path to the directory containing the PDB and SPH files.
    pdb_name : str
        The name of the PDB file (without the extension).
    sph_name : str
        The name of the SPH file (without the extension).

    Returns
    -------
    None
        This function does not return anything, but generates a text file and
        saves plots of the ellipsoids.

    Raises
    ------
    IOError
        If the PDB or SPH file cannot be read.

    Notes
    -----
    This function uses the MDAnalysis library to read and analyze the PDB and
    SPH files, and the pandas and matplotlib libraries to plot the data and
    save the output.

    The output of this function is a text file that contains the parameters of
    the ellipsoids that represent the pore path of the biomolecule. The format
    of the text file is as follows:

        #x, y, z, a, b, theta
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0
        ...

    The `x`, `y`, and `z` values represent the center of the ellipsoid, while
    the `a`, `b`, and `theta` values represent the semi-axes and orientation of
    the ellipsoid.

    The plots of the ellipsoids are saved in the directory
    `p+pdb_name+'_pathway_slices/'`, where `p` is the path to the directory
    containing the PDB and SPH files, and `pdb_name` is the name of the PDB
    file.
    """
    start = time.time()
    print(pdb_name, sph_name)
    try:
        conf = p + sph_name + '.sph'
        top = conf
        sph = MDAnalysis.Universe(top, conf, topology_format='pdb', format='pdb') # tpr_resid_from_one=True
    except:
        print(' ERROR could not find SPH file', p + sph_name + '.sph')
        conf = 'hole000.sph'
        top = conf
        sph = MDAnalysis.Universe(top, conf, topology_format='pdb', format='pdb') # tpr_resid_from_one=True
    

    conf =  p + pdb_name 
    top = conf
    u = MDAnalysis.Universe(top, conf, topology_format='pdb', format='pdb')

    mer = MDAnalysis.Merge(u.atoms, sph.atoms).select_atoms('name *')
    print('len(mer)', len(mer.atoms), 'len(sph)', len(sph.atoms), 'len(u)', len(u.atoms))

    radii = sph.atoms.occupancies
    resids = sph.atoms.resids
    x_coordinates = []
    y_coordinates = []
    z_coordinates = []
    for i, pos in enumerate(sph.atoms.positions):
        x_coordinates.append(pos[0])
        y_coordinates.append(pos[1])
        z_coordinates.append(pos[2])

    x_coordinates = np.array(x_coordinates)
    y_coordinates = np.array(y_coordinates)
    z_coordinates = np.array(z_coordinates)

    ind_sort = np.argsort(z_coordinates)
    print('len(ind_sort), len(x_coordinates)', len(ind_sort), len(x_coordinates))
    x_coordinates = x_coordinates[ind_sort]
    y_coordinates = y_coordinates[ind_sort]
    z_coordinates = z_coordinates[ind_sort]
    print('min(z)', z_coordinates[0], 'max(z)', z_coordinates[-1])
    radii = radii[ind_sort]
    resids = resids[ind_sort]

    d = {'x': x_coordinates, 'y': y_coordinates,'z': z_coordinates, 'r': radii, 'resid': resids}
    df = pd.DataFrame(data=d)
    df2 = df[(df['r']<end_radius) ]
    
    print('df2.columns',df2.columns)
    print('df2.index',df2.index)
    print('len(df)',len(df), 'len(df2)', len(df2))

    fig, ax = plt.subplots()
    ax.set_ylim([0, end_radius])
    plt.plot(df2['z'], df2['r'])
    plt.ylabel(r"HOLE quantities ($\AA$)", fontsize=f_size)
    plt.xlabel(r"Pore coordinate $\zeta$ ($\AA$)", fontsize=f_size)

    ax.tick_params(axis='both', which='major', labelsize=f_size)
    fig.tight_layout()
    plt.show()

    if not os.path.exists(p+pdb_name+'_pathway_slices/'):
        os.makedirs(p+pdb_name+'_pathway_slices/')
    f = open(p+pdb_name+'_pathway_ellipse.txt','w')
    f.write('#x, y, z, a, b, theta\n')
    failed = []
    vec = np.array(df2.index)
    print('int(start_index)',int(start_index),':-int(start_index):',
          -int(start_index),'slice_dz',slice_dz)
    vec = vec[int(start_index):-int(start_index):slice_dz]
    print('z-coordinated processes: min',df2['z'].loc[vec[0]], 
          'max', df2['z'].loc[vec[-1]])
    print('vec',vec)
    if parallel:
        n = multiprocessing.cpu_count()
        print('multiprocessing.cpu_count()', n)

        start_time = time.time()
        results = insert_ellipse_async(index=vec, dataframe=df2, universe=mer, 
                                       out=out, plt_path=p+pdb_name+'_pathway_slices/',
                                       num_processes=num_processes, timeout=timeout,
                                       n_xy_fac=n_xy_fac,
                                       pathway_sel=pathway_sel, opt_method=opt_method
                                      )   
        # print(results)
        
        print("--- %s seconds --- for pool.starmap_async" % (time.time() - start_time))
        start_time = time.time()
        for count, result in enumerate(results):
            try:
                e, z = result
                position_center = str(e.cx) + ', ' +  str(e.cy) + ', ' + str(z) + ', '
                param = str(e.a) + ', ' +  str(e.b) + ', ' + str(e.theta) + '\n'
                f.write(position_center + param)
            except:
                failed.append(count)
        print("--- %s seconds --- for processing results" % (time.time() - start_time))
    else:
        start_time = time.time()
        for i in vec:
            # print(df2.iloc[i])
            # try:
                print('try', i ,df2['resid'].loc[i])
                e, z = insert_ellipse(index=i, dataframe=df2, universe=mer, 
                                      out=0, show=1,
                                    plt_path=p+pdb_name+'_pathway_slices/',
                                    pathway_sel=pathway_sel,
                                    opt_method=opt_method,
                                    #R_N=25
                                    )
                position_center = str(e.cx) + ', ' +  str(e.cy) + ', ' + str(z) + ', '
                param = str(e.a) + ', ' +  str(e.b) + ', ' + str(e.theta) + '\n'
                f.write(position_center + param)
            # except:
            #    failed.append(i)
        print("--- %s seconds --- for insert_ellipse with one worker (not parallel)" % (time.time() - start_time))
    print('failed slices', failed)
    f.close()
    end = time.time()
    print("TIME(ellipsoid_pathway)=",end - start)

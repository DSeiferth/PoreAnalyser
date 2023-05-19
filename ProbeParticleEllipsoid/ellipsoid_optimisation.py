import numpy as np
import MDAnalysis
import nglview as nv
import matplotlib.pyplot as plt
import pandas as pd

from ellipse_lib import atom, ellipse, distance_ellipse, dist_ellipse_vdwSphere, assign_radius

from scipy.optimize import minimize
from scipy.optimize import NonlinearConstraint

import multiprocessing 
import time
import os

def penalty_overlap_4dim(x, args, stop_loop=True):
    '''
    x: 4-dim vec
    0) radius to expand 
    1) angle
    2)-3) position of centre
    args: 
    0) constant radius
    1) vec with surroundings vdw-spheres
    '''
    e1 = ellipse(a=x[0], b=args[0], theta=x[1], cx=x[2], cy=x[3])
    a_vec = args[1]
    
    #print('ellipse: a', e1.a, 'b', e1.b, 'theta', e1.theta, 'center', e1.cx, e1.cy)
    ### spheres ###        
    dist = np.zeros(len(a_vec))
    for i, atom in enumerate(a_vec):
        dist[i] = dist_ellipse_vdwSphere(ellipse=e1, sphere=atom)
        if stop_loop: ### make overlap costly and avoid calculating minimum ###
            if dist[i]<0 : ### overlap ###
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
    
def neighbor_vec(universe, probe, probe1, n_xy_fac, out=0, call=0):
    protein = universe.select_atoms('protein')
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
        if R_projected>0:
            a_vec.append(atom(n.position[0],n.position[1], r=R_projected ))
            neighbour_labels.append(n.resname + ' ' + str(n.resid)) # + ' ' + n.type
            neighbour_resnum_resid.append(n.resname + ' ' + str(n.resid))
            neighbour_string = neighbour_string + n.resname + ' ' + str(n.resid) + ' ' + n.type + '\n'
    if out: print('number of neighbors with R>0' ,len(a_vec), 'n_xy_fac', n_xy_fac)
    if len(a_vec)>150 and call<4:
        if out: print('DECREASE n_xy_fac')
        call += 1
        return neighbor_vec(universe, probe, probe1, 0.75*n_xy_fac,
                                              call=call, out=out)
        print('number of neighbors with R>0' ,len(a_vec), 'n_xy_fac', 0.75*n_xy_fac)
    elif len(a_vec)<30 and call<4:
        if out: print('INCREASE n_xy_fac')
        call += 1 
        return neighbor_vec(universe, probe, probe1, 1.25*n_xy_fac,
                                              call=call, out=out)
        if out: print('number of neighbors with R>0' ,len(a_vec), 'n_xy_fac', 1.25*n_xy_fac)
        #return -1, -1
    
    return a_vec, neighbour_labels

def insert_ellipse(index, dataframe, universe, 
                   out=0, plt_path='', rmax=50, 
                   show=0, label=0, n_xy_fac=1.6,
                   timing = 0,
                   f_size = 22
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
    The maximum radius for the ellipsoid

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
    initial_radius = dataframe['r'].loc[index]

    ### neighbor vector ###
    a_vec, neighbour_labels = neighbor_vec(universe, probe, probe1, n_xy_fac, out=out)
    assert len(a_vec)>2, "in function 'insert_ellipse': len(a_vec)<2="+str(len(a_vec))

    fig, ax = plt.subplots()
    plt.gca().set_aspect('equal', adjustable='box')
    
    ### probe ###
    p0 = ellipse(a=probe.r, b=probe.r, theta=0, cx=probe.x, cy=probe.y)
    x0, y0 = p0.draw()
    ax.plot(x0, y0, color='blue')
    ax.plot(probe.x,probe.y, '-x', color='blue')
    ### vdw particles ###
    start_neighbors_loop = time.time()
    #prev = neighbour_labels[0]
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
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.NonlinearConstraint.html#scipy.optimize.NonlinearConstraint
    ### 4d with small variation ###
    dx = 0.1*p0.a ### boundary for center of ellipse ###
    pt = [p0.a, 0,  p0.cx, p0.cy]
    opt = {'maxiter': 800, # Maximum allowed number of iterations and function evaluations. Will default to N*200, where N is the number of variables
           #'xatol': 0.00001, # Absolute error in xopt between iterations that is acceptable for convergence.
           #'fatol': 0.00001, # Absolute error in func(xopt) between iterations that is acceptable for convergence.
           'adaptive': True, #Adapt algorithm parameters to dimensionality of problem. Useful for high-dimensional minimization 
           'initial_simplex':[pt, [p0.a+0.1, np.pi/4, p0.cx, p0.cy], 
                              [p0.a+0.11, np.pi/2, p0.cx, p0.cy],
                              [p0.a+0.12, np.pi/4, p0.cx+dx/2, p0.cy+dx/2],
                              [p0.a+0.13, np.pi/2, p0.cx-dx/2, p0.cy-dx/2]]
          }
    # (p0.a, rmax)
    bnds = ((0, rmax), (-np.pi, np.pi), (p0.cx-dx, p0.cx+dx), (p0.cy-dx, p0.cy+dx))
    start_opt1 = time.time()
    result = minimize(penalty_overlap_4dim, pt, 
                      args = [p0.b, a_vec],
                      method='nelder-mead',
                      bounds=bnds,
                      #constraints=constr,
                      options = opt
                      #tol=1e-10
                     )
    end = time.time()
    if timing: print("TIME(start_opt1)=",end - start_opt1)

    sol = result['x']
    if result['fun'] > 0:
        print('ERROR: overlap', result['fun'])
        return -1, -1
    if out: print(result)
    if out: print('4d optimisation small variation of center (r, theta, x, y)', sol) #result
    p1 = ellipse(a=sol[0], b=p0.b, theta=sol[1], cx= p0.cx, cy= p0.cy)
    x1, y1 = p1.draw()
    ax.plot(x1, y1, color='green')
    ax.plot(p1.cx, p1.cy, '-x', color='green')

    ### 4d optimisation works better without initial simplex ###
    dx = 0.5*p1.a # 3 ### boundary for center of ellipse ###
    bnds = ((p1.a, rmax), (-np.pi, np.pi), (p0.cx-dx, p0.cx+dx), (p0.cy-dx, p0.cy+dx))
    pt = [p1.a, p1.theta, p1.cx, p1.cy ] #p0.cx, p0.cy -1,-1
    opt = {'maxiter': 800, # Maximum allowed number of iterations and function evaluations. Will default to N*200, where N is the number of variables
           'adaptive': True, #Adapt algorithm parameters to dimensionality of problem. Useful for high-dimensional minimization 
           #'initial_simplex':[pt, [p1.a+0.1,  p1.theta+np.pi/4, p1.cx, p1.cy], 
           #                   [p1.a+0.15,  p1.theta+np.pi/2, p1.cx, p1.cy],
           #                   [p1.a+0.2,  p1.theta+np.pi/6, p1.cx, p1.cy],
           #                   [p1.a+0.25,  p1.theta+np.pi/3, p1.cx, p1.cy]]
          }
    start_opt2 = time.time()
    result = minimize(penalty_overlap_4dim, pt, 
                      args = [p1.b, a_vec],
                      method='nelder-mead',
                      bounds=bnds,
                      #constraints=constr,
                      options = opt
                      #tol=1e-10
                     )
    end = time.time()
    if timing:  print("TIME(start_opt2)=",end - start_opt2)
    sol = result['x']
    if out: print('4d optimisation large variation of center (r, theta, x, y)', sol)

    p2 = ellipse(a=sol[0], b=p0.b, theta=sol[1], cx=sol[2], cy=sol[3])
    x2, y2 = p2.draw()
    ax.plot(x2, y2, color='brown')
    ax.plot(p2.cx, p2.cy, '-x', color='brown')
    ax.set_xlabel(r"x ($\AA$)", fontsize=f_size)
    ax.set_ylabel(r"y ($\AA$)", fontsize=f_size)
    ax.tick_params(axis='both', which='major', labelsize=f_size)
    plt.title(r'xy-plane at z = '+str(z_slice)+r' $\AA$', fontsize=f_size)
    
    if label: 
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        plt.subplots_adjust(right=0.35)
        neighbour_resnum_resid = np.unique(neighbour_resnum_resid)
        string = ''
        for info in neighbour_resnum_resid:
            string = string + info +'\n'
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
                         show=0, label=0, n_xy_fac=1.6,num_processes=None, timeout=20, f_size=22):
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
                                show, label, n_xy_fac, f_size=f_size)
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
                      end_radius = 15,
                      
                      num_processes=None, 
                      timeout=20,
                      start_index = 50,
                      f_size = 22,
                      out = 0,
                      n_xy_fac = 1.6

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
    conf = p + sph_name + '.sph'
    top = conf
    sph = MDAnalysis.Universe(top, conf, topology_format='pdb', format='pdb') # tpr_resid_from_one=True

    conf =  p + pdb_name + '.pdb'
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
    #ax.legend(prop={'size': 11}) # loc='upper center'
    fig.tight_layout()
    plt.show()

    if not os.path.exists(p+pdb_name+'_pathway_slices_parallel2/'):
        os.makedirs(p+pdb_name+'_pathway_slices_parallel2/')
    f = open(p+pdb_name+'_pathway_ellipse_parallel2.txt','w')
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
                                       out=out, plt_path=p+pdb_name+'_pathway_slices_parallel2/',
                                       num_processes=num_processes, timeout=timeout,
                                       n_xy_fac=n_xy_fac
                                      )   
        #print(results)
        
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
            #print(df2.iloc[i])
            try:
                print('try', i ,df2['resid'].loc[i])
                e, z = insert_ellipse(index=i, dataframe=df2, universe=mer, 
                                      out=0, show=1,
                                    plt_path=p+pdb_name+'_pathway_slices_parallel2/',
                                    R_N=25)
                position_center = str(e.cx) + ', ' +  str(e.cy) + ', ' + str(z) + ', '
                param = str(e.a) + ', ' +  str(e.b) + ', ' + str(e.theta) + '\n'
                f.write(position_center + param)
            except:
                failed.append(i)
        print("--- %s seconds --- for insert_ellipse with one worker (not parallel)" % (time.time() - start_time))
    print('failed slices', failed)
    f.close()
    end = time.time()
    print("TIME(ellipsoid_pathway)=",end - start)

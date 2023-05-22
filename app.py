import streamlit as st
import hole_analysis as hole_analysis
import os
import io
import MDAnalysis
from stmol import showmol
import py3Dmol
import numpy as np
from visualization import write_pdb_with_pore_surface, plt_ellipsoid_pathway, pathway_visu, st_write_ellipsoid, write_pdb_with_ellipsoid_surface
from download_files import download_output, download_Ellipsoid_output

try:
    import multiprocessing 
    print("Number of processors: ", multiprocessing.cpu_count())
    #st.write("Number of processors: ", multiprocessing.cpu_count())
    parallel = True
except:
    parallel = False
    st.write("Could not import multiprocessing library, => multiprocessing disabled")
import sys
sys.path.append('ProbeParticleEllipsoid/')
from ellipsoid_optimisation import ellipsoid_pathway
import pandas as pd
import matplotlib.pyplot as plt


st.title("Extending the capabilities of the HOLE Package for Annotation of Ion Channels")
#st.latex(r''' a+a r^1+a r^2+a r^3 ''')

str1 = 'Over the last two decades, advances in structural biology along with recent artificial intelligence–driven structure prediction '
str2 = 'algorithms, such as AlphaFold, have revealed a plethora of 3-D ion channel and nanopore structures in different conformational states. '
str3 = 'However, in nearly every case, these structures still require functional annotation.'
str4 = 'Different tools, such as HOLE and CHAP, allow the analysis of the physical dimensions of the pore running through an ion channel.'  
str5 = 'Here, we present an interactive web-page that allows users to calculate the pore profile of any input structure. '
str6 = 'Based on the well-established HOLE programme, we add a new feature to capture pore asymmetry by using an ellipsoidal probe particle.'
st.write(str1+str2+str3+str4+str5+str6)

st.subheader("Pathway Finding Settings")
string1 = 'Radius in Å, which is considered to be the end of the pore. '
string2 = 'This keyword can be used to specify the radius above '
string3 = 'which the program regards a result as indicating that the end of the pore has been reached. '
string4 = 'This may need to be increased for large channels, or reduced for small channels. '
end_radius = st.text_input(label=r'end_radius in $\AA$', value='15', max_chars=3,
              help=string1+string2+string3+string4)
end_radius = int(end_radius)
st.write('The current end_radius is', end_radius, r'$\AA$')

st.subheader("Plotting options")
fig_format = st.text_input(label='Format to download pathway figure', value='png', max_chars=4,
              help='default: png, other options: jpeg, tiff, eps, pdf, ...')
string1 = r'The dashed red line indicates where the pore radius is to tight for a water molecule (r < 1.15 $\AA$). '
string2 = r'The dashed green line indicates where there is room for a single water (r < 2.30 $\AA$).'
st.write(string1+string2)
plot_lines = st.text_input(label='Plot red and green lines (default: True)', value='True', max_chars=5,
              help=string1+string2)
if plot_lines == 'True':
    plot_lines = True
else:
    plot_lines = False
title = st.text_input(label='Write a title for your plot', value='', help='Title string for plot')
f_size = st.text_input(label='Font size for figure', value='22', help='default=22')
f_size = int(f_size)

st.subheader("Upload pdb file(s)")
uploaded_files = st.file_uploader("Choose a file", label_visibility="visible",  accept_multiple_files=True )

st.subheader("Pathfinding with a spherical probe particle (HOLE)")
string1 = "HOLE is a program that allows the analysis and visualisation of the pore dimensions of the holes "
string2 = "through molecular structures of ion channels (Smart et al., 1996)."
st.write(string1+string2)
string1 = "Here, we use the MDAnalysis interface for HOLE to analyse an ion channel pore or transporter pathway. "
string2 = "The original HOLE documentation can be found here: https://www.holeprogram.org"
st.write(string1+string2)

st.write('First, we align the principal axis to the z-axis.')

labels = []
names = []
names_aligned = []
if uploaded_files:
    for uploaded_file in uploaded_files:
        #st.write("Filename: ", uploaded_file.name)
        labels.append(uploaded_file.name[:-4])
        names.append(uploaded_file.name)
        names_aligned.append(uploaded_file.name[:4]+'_aligned_z.pdb')
        with open(uploaded_file.name,"wb") as f:
            f.write(uploaded_file.getbuffer())
    #st.write('Uploaded', names)
    fig , df = hole_analysis.analysis(names, labels=labels, path='', end_radius=end_radius, save='Uploaded', title=title,
                                            legend_outside=True, plot_lines=plot_lines, f_size=f_size)
    st.pyplot(fig)
    path_save = ''
    st.write("Pathway visualisation for ", names[0])
    write_pdb_with_pore_surface(path=path_save, name=names[0], end_radius=end_radius, num_circle = 24)
    xyzview = pathway_visu(path=path_save, name=names[0])
    showmol(xyzview, height=800, width=800)

    st.subheader("Download HOLE output files")
    df.to_csv('hole_pathway_profile.csv',sep=',')
    fn ="hole_pathway_profile."+fig_format
    fig.savefig(fn, format=fig_format, bbox_inches='tight')
    ### Download ###
    download_output(names_aligned[0][:-4], fn, df, fig, fig_format, path_save, names )

    ### Elipsoid ###
    st_write_ellipsoid()
    ellipsoid_pathway(p=path_save, 
                        pdb_name = names_aligned[0], 
                        sph_name = names_aligned[0][:-4], 
                        slice_dz=4, parallel=True, 
                        num_processes=None, timeout=6, 
                        start_index = 1, end_radius=end_radius-1,
                        out = 0,
                        n_xy_fac = 1.6
                    )
    res = np.loadtxt(path_save + names_aligned[0]+ '_pathway_ellipse.txt', 
                    comments='#', delimiter=',')
    df_res = pd.DataFrame(data=res, columns=['x', 'y', 'z', 'a', 'b', 'theta'])
    df_res.sort_values('z', inplace=True)
    fig = plt_ellipsoid_pathway(df_res, f_size=f_size, title=title, end_radius=end_radius)
    st.pyplot(fig)
    ### visualization fo ellipsoidal surface ###
    write_pdb_with_ellipsoid_surface(p='', pdbname=names_aligned[0], 
                                     fname=names_aligned[0]+'_pathway_ellipse.txt', num_circle = 24)
    xyzview = pathway_visu(path='', name=names_aligned[0], f_end='_ellipsoid.pdb')
    showmol(xyzview, height=800, width=800)

    ### Download Ellipsoid output###
    st.subheader("Download files for pathway with ellipsoidal probe particle")
    fn ="ellipsoid_pathway_profile."+fig_format
    fig.savefig(fn, format=fig_format, bbox_inches='tight')
    download_Ellipsoid_output(names_aligned[0][:-4], fn, path_save,  )

    #st.write('ERROR with', names)
else:
    st.markdown("Example application with 7tu9")
    st.write("Example Filename: ", "pdb_models/7tu9.pdb")
    path_save = 'pdb_models/'
    titles = [r'$\alpha$$\beta$ Heteromeric Glycine Receptor']

    labels = [
        'GlyR-Stry', 
        #'GlyR-Gly',
        #'GlyR-Gly-Ivm'
            ]
    names = [
        '7tu9.pdb', #'GlyR-Stry', 
        #'7tvi.pdb', #'GlyR-Gly',
        #'8fe1.pdb', # 'GlyR-Gly-Ivm'
    ]
    fig ,csv = hole_analysis.analysis(names, labels=labels, path=path_save, 
                             end_radius=end_radius, 
                       #TMD_lower=59, TMD_higher=97,
                        save='', title=titles[0], 
                       legend_outside=True,
                       plot_lines=plot_lines,
                       f_size=f_size
                       )

    st.pyplot(fig)
    write_pdb_with_pore_surface(path=path_save, name=names[0], end_radius=end_radius, num_circle = 24)
    # https://github.com/napoles-uach/stmol
    # https://william-dawson.github.io/using-py3dmol.html
    xyzview = pathway_visu(path=path_save, name=names[0])
    showmol(xyzview, height=800, width=800)

    ### Ellipsoidal probe particle ###
    st_write_ellipsoid()
    res = np.loadtxt('pdb_models/7tu9_aligned_z.pdb_pathway_ellipse.txt', 
                 comments='#', delimiter=',')
    df_res = pd.DataFrame(data=res, columns=['x', 'y', 'z', 'a', 'b', 'theta'])
    df_res.sort_values('z', inplace=True)
    fig = plt_ellipsoid_pathway(df_res, f_size=f_size, title=title, end_radius=end_radius)
    st.pyplot(fig)
    write_pdb_with_ellipsoid_surface(p='pdb_models/', pdbname='7tu9_aligned_z.pdb',
                                     fname='7tu9_aligned_z.pdb_pathway_ellipse.txt', num_circle = 24)
    xyzview = pathway_visu(path='pdb_models/', name='7tu9_aligned_z.pdb', f_end='_ellipsoid.pdb')
    showmol(xyzview, height=800, width=800)


#st.write(os.listdir())

st.write("Smart, O.S., Neduvelil, J.G., Wang, X., Wallace, B.A., Sansom, M.S.P., 1996. HOLE: A program for the analysis of the pore dimensions of ion channel structural models. Journal of Molecular Graphics 14, 354–360. https://doi.org/10.1016/S0263-7855(97)00009-X")
st.write("Gowers, R., Linke, M., Barnoud, J., Reddy, T., Melo, M., Seyler, S., Domański, J., Dotson, D., Buchoux, S., Kenney, I., Beckstein, O., 2016. MDAnalysis: A Python Package for the Rapid Analysis of Molecular Dynamics Simulations. Presented at the Python in Science Conference, Austin, Texas, pp. 98–105. https://doi.org/10.25080/Majora-629e541a-00e")


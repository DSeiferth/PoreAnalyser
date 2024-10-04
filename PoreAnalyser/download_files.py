from zipfile import ZipFile
import streamlit as st
import pandas as pd
import io
import glob

@st.cache_data
def convert_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv().encode('utf-8')

def download_output(pdb_name, fn, df, fig, fig_format,
                    path_save, names ):
    # To zip entire directory use command “shutil.make_archive(“name”,”zip”, root_dir)
    # To select the files to zip use command “ZipFile.write(filename)”
    # shutil.make_archive('', 'zip', dir_name)
    with ZipFile("poreFinding_HOLE_output.zip", "w") as newzip:
            newzip.write(pdb_name+'.vmd')
            newzip.write("visualise_pathway_hole.tcl")
            newzip.write("pymol_pore_visu.py")
            newzip.write("chimera_pore.py")
            newzip.write("chimeraX_pore.py")
            newzip.write(path_save + names[0] + '_circle.pdb')
            newzip.write(pdb_name+'.pdb')
            newzip.write("README.md")
            newzip.write("hole.out")
            newzip.write('hole_pathway_profile.csv')
            newzip.write(fn)
    with open("poreFinding_HOLE_output.zip", "rb") as file:
            st.download_button(
                label="Download ZIP",
                data=file,
                file_name="poreFinding_HOLE_output.zip",
                help="usage: vmd -e visualise_pathway_hole.tcl -args  "+pdb_name+" "+pdb_name+".vmd"
                #mime='text/csv',
            )

    csv = convert_df(df)
    st.download_button(
            label="Download pathway profile as CSV",
            data=csv,
            file_name='hole_pathway_profile.csv',
            mime='text/csv',
        )
    fn ="hole_pathway_profile."+fig_format
    img = io.BytesIO() # Create an in-memory buffer
    fig.savefig(img, format=fig_format, bbox_inches='tight')
    st.download_button(
            label="Download figure",
            data=img,
            file_name=fn,
            mime="image/"+fig_format
        )
    ### download vmd file ###
    with open(pdb_name+'.vmd', "rb") as file:
        st.download_button(
              label="Download vmd pathway visualisation",
                data=file,
                file_name=pdb_name+'.vmd',
                help="usage: vmd -e visualise_pathway_hole.tcl -args  "+pdb_name+" "+pdb_name+".vmd"
        )
    with open('visualise_pathway_hole.tcl', "rb") as file:
            st.download_button(
                label="Download vmd visualisation TCL script (you need the corresponding pdb and vmd file)",
                data=file,
                file_name="visualise_pathway_hole.tcl",
                help="usage: vmd -e visualise_pathway_hole.tcl -args  fname.pdb fname.pdb.vmd"
                #mime='text/csv',
            )

def download_Ellipsoid_output(pdb_name, fn, path_save,  ):
    # To zip entire directory use command “shutil.make_archive(“name”,”zip”, root_dir)
    # To select the files to zip use command “ZipFile.write(filename)”
    # shutil.make_archive('', 'zip', dir_name)
    with ZipFile("poreFinding_Ellipsoid_output.zip", "w") as newzip:
            ### vmd ###
            newzip.write(pdb_name+'.pdb_pathway_ellipse.vmd')
            newzip.write("visualise_pathway_hole.tcl")
            newzip.write("pymol_pore_visu.py")
            newzip.write("chimera_pore.py")
            newzip.write("chimeraX_pore.py")
            ### pdb files ### 
            newzip.write(path_save + pdb_name + '.pdb_ellipsoid.pdb')
            newzip.write(pdb_name+'.pdb')

            newzip.write("README.md")
            newzip.write(pdb_name + '.pdb_pathway_ellipse.txt')
            newzip.write(fn)
            for file in glob.glob(pdb_name+'.pdb_pathway_slices/*png'):
                  newzip.write(file)
    with open("poreFinding_Ellipsoid_output.zip", "rb") as file:
            st.download_button(
                label="Download ZIP",
                data=file,
                file_name="poreFinding_Ellipsoid_output.zip",
                help="usage: vmd -e visualise_pathway_hole.tcl -args  "+pdb_name+".pdb "+pdb_name+".pdb_pathway_ellipse.vmd"
                #mime='text/csv',
            )

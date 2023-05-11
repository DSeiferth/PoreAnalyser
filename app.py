import streamlit as st
import hole_analysis as hole_analysis
#import os

st.title("Pore Analysis with HOLE")
#st.header("this is the markdown")
#st.latex(r''' a+a r^1+a r^2+a r^3 ''')

#x = st.slider('Select a value')
#st.write(x, 'squared is', x * x)

st.write("HOLE is a program that allows the analysis and visualisation of the pore dimensions of the holes through molecular structures of ion channels Smart et al., 1996.")

st.subheader("Upload pdb file(s)")
uploaded_files = st.file_uploader("Choose a file", label_visibility="visible",  accept_multiple_files=True )

labels = []
names = []
if uploaded_files:
    for uploaded_file in uploaded_files:
        st.write("Filename: ", uploaded_file.name)
        labels.append(uploaded_file.name)
        names.append(uploaded_file.name)
        with open(uploaded_file.name,"wb") as f:
            f.write(uploaded_file.getbuffer())
    st.write('Uploaded', names)
    try:
        fig = hole_analysis.analysis(names, labels=labels, path='', end_radius=15, save='Uploaded', title='',legend_outside=True)
        st.pyplot(fig)
    except:
        st.write('ERROR with', names)
else:
    st.markdown("Example application with 7tu9")
    st.write("Example Filename: ", "pdb_models/7tu9.pdb")
    path_save = 'pdb_models/'
    dirs = ['']
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
    fig = hole_analysis.analysis(names, labels=labels, path=path_save+dirs[0], 
                             end_radius=15, 
                       #TMD_lower=59, TMD_higher=97,
                        save='', title=titles[0], 
                       legend_outside=True
                       )

    st.pyplot(fig)


st.write("Smart, O.S., Neduvelil, J.G., Wang, X., Wallace, B.A., Sansom, M.S.P., 1996. HOLE: A program for the analysis of the pore dimensions of ion channel structural models. Journal of Molecular Graphics 14, 354–360. https://doi.org/10.1016/S0263-7855(97)00009-X")
import streamlit as st
import hole_analysis as hole_analysis
import os

x = st.slider('Select a value')
st.write(x, 'squared is', x * x)


path_save = 'pdb_models/'
dirs = ['']

titles = ['Alpha1/BetaB Heteromeric Glycine Receptor']

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

uploaded_files = st.file_uploader("Choose a file", label_visibility="visible",  accept_multiple_files=True )

labels = []
names = []
if uploaded_files:
    for uploaded_file in uploaded_files:
        st.write("Filename: ", uploaded_file.name)
        labels.append(uploaded_file.name)
        names.append(uploaded_file.name)
        with open(os.path.join("tempDir",uploaded_file.name),"wb") as f:
            f.write(uploaded_file.getbuffer())
    fig = hole_analysis.analysis(names, labels=labels, path='tempDir/', end_radius=15, save='Uploaded', title='',legend_outside=True)
    st.pyplot(fig)

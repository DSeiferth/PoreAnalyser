import streamlit as st
import hole_analysis as hole_analysis
#import os
import io
import plotly.graph_objects as go

@st.cache_data
def convert_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv().encode('utf-8')

st.title("Pore Analysis with HOLE")
#st.latex(r''' a+a r^1+a r^2+a r^3 ''')

#x = st.slider('Select a value')
#st.write(x, 'squared is', x * x)

st.write("HOLE is a program that allows the analysis and visualisation of the pore dimensions of the holes through molecular structures of ion channels Smart et al., 1996.")

st.header("Pathway Finding Settings")
string1 = 'Radius in Å, which is considered to be the end of the pore. '
string2 = 'This keyword can be used to specify the radius above '
string3 = 'which the program regards a result as indicating that the end of the pore has been reached. '
string4 = 'This may need to be increased for large channels, or reduced for small channels. '
end_radius = st.text_input(label=r'end_radius in $\AA$', value='15', max_chars=3,
              help=string1+string2+string3+string4)
end_radius = int(end_radius)
st.write('The current end_radius is', end_radius, r'$\AA$')
fig_format = st.text_input(label='Format to download pathway figure', value='png', max_chars=4,
              help='default png')

st.subheader("Upload pdb file(s)")
uploaded_files = st.file_uploader("Choose a file", label_visibility="visible",  accept_multiple_files=True )

labels = []
names = []
if uploaded_files:
    for uploaded_file in uploaded_files:
        #st.write("Filename: ", uploaded_file.name)
        labels.append(uploaded_file.name)
        names.append(uploaded_file.name)
        with open(uploaded_file.name,"wb") as f:
            f.write(uploaded_file.getbuffer())
    #st.write('Uploaded', names)
    bla=1
    if bla:
    #try:
        fig , df = hole_analysis.analysis(names, labels=labels, path='', end_radius=end_radius, save='Uploaded', title='',legend_outside=True)
        st.pyplot(fig)
        #st.write("pathway ", df)
        csv = convert_df(df)
        st.download_button(
            label="Download pathway profile as CSV",
            data=csv,
            file_name='hole_pathway_profile.csv',
            mime='text/csv',
        )
        fn ="hole_pathway_profile."+fig_format
        img = io.BytesIO() # Create an in-memory buffer
        fig.savefig(img, format=fig_format)
        #byte_im = buffer.getvalue()
        # Save the figure as a pdf to the buffer
        #fig.write_image(file=buffer, format=fig_format) # 'Figure' object has no attribute 'write_image'
        # Download the pdf from the buffer
        st.download_button(
            label="Download figure",
            data=img,
            file_name=fn,
            mime="image/"+fig_format
        )
    #except:
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
    fig ,csv = hole_analysis.analysis(names, labels=labels, path=path_save+dirs[0], 
                             end_radius=end_radius, 
                       #TMD_lower=59, TMD_higher=97,
                        save='', title=titles[0], 
                       legend_outside=True
                       )

    st.pyplot(fig)


st.write("Smart, O.S., Neduvelil, J.G., Wang, X., Wallace, B.A., Sansom, M.S.P., 1996. HOLE: A program for the analysis of the pore dimensions of ion channel structural models. Journal of Molecular Graphics 14, 354–360. https://doi.org/10.1016/S0263-7855(97)00009-X")
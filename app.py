import streamlit as st
import hole_analysis as hole_analysis
#import os
import io

@st.cache_data
def convert_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv().encode('utf-8')

st.title("Pore Analysis with HOLE")
#st.latex(r''' a+a r^1+a r^2+a r^3 ''')

string1 = "HOLE is a program that allows the analysis and visualisation of the pore dimensions of the holes"
string2 = "through molecular structures of ion channels Smart et al., 1996."
st.write(string1+string2)

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
              help='default png')
string1 = r'The dashed red line indicates where the pore radius is to tight for a water molecule (r < 1.15 $\AA$). '
string2 = r'The dashed green line indicates where there is room for a single water (r < 2.30 $\AA$).'
st.write(string1+string2)
plot_lines = st.text_input(label='Plot red and green lines (default: True)', value='True', max_chars=5,
              help=string1+string2)
plot_lines = bool(plot_lines)
title = st.text_input(label='Write a title for your plot', value='', help='Title string for plot')
f_size = st.text_input(label='Font size for figure', value='22', help='default=22')
f_size = int(f_size)

st.header("Upload pdb file(s)")
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
    try:
        fig , df = hole_analysis.analysis(names, labels=labels, path='', end_radius=end_radius, save='Uploaded', title=title,
                                          legend_outside=True, plot_lines=plot_lines, f_size=f_size)
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
        fig.savefig(img, format=fig_format, bbox_inches='tight')
        st.download_button(
            label="Download figure",
            data=img,
            file_name=fn,
            mime="image/"+fig_format
        )
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
    fig ,csv = hole_analysis.analysis(names, labels=labels, path=path_save+dirs[0], 
                             end_radius=end_radius, 
                       #TMD_lower=59, TMD_higher=97,
                        save='', title=titles[0], 
                       legend_outside=True,
                       plot_line=plot_line,
                       f_size=f_size
                       )

    st.pyplot(fig)


st.write("Smart, O.S., Neduvelil, J.G., Wang, X., Wallace, B.A., Sansom, M.S.P., 1996. HOLE: A program for the analysis of the pore dimensions of ion channel structural models. Journal of Molecular Graphics 14, 354–360. https://doi.org/10.1016/S0263-7855(97)00009-X")
st.write("Gowers, R., Linke, M., Barnoud, J., Reddy, T., Melo, M., Seyler, S., Domański, J., Dotson, D., Buchoux, S., Kenney, I., Beckstein, O., 2016. MDAnalysis: A Python Package for the Rapid Analysis of Molecular Dynamics Simulations. Presented at the Python in Science Conference, Austin, Texas, pp. 98–105. https://doi.org/10.25080/Majora-629e541a-00e")
import streamlit as st
import pydynamo

st.title("Radial Distribution Tool")

output_file = st.file_uploader("Choose a Output file")

if output_file is not None:
    of = pydynamo.OutputFile(output_file)
    of.tree.findall('.//RadialDistribution/Species')
    print(of)
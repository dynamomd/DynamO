import streamlit as st
import pydynamo

st.title("Radial Distribution Tool")

output_file = st.file_uploader("Choose a Output file")

if output_file is not None:
    of = pydynamo.OutputFile(output_file)
    grs = of.tree.findall('.//RadialDistributionMoments/Species')
    all_species = set(s.attrib["Name1"] for s in grs).union(set(s.attrib["Name2"] for s in grs))
    print(all_species)

    for species in grs:
        for moment in species.findall('./Moment'):
            for line in moment.text.split("\n")[:3]:
                print([(species.attrib["Name1"], species.attrib["Name2"], moment.attrib['Order'], *line.split())])
    print("Done!")

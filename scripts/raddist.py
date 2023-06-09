import streamlit as st
import pydynamo
import math
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

st.title("Radial Distribution Tool")

output_file = st.file_uploader("Choose a Output file")

def getGrData(of : pydynamo.OutputFile):
    grs = of.tree.findall('.//RadialDistributionMoments/Species')
    N=int(of.tree.find('.//ParticleCount').attrib['val'])
    #Loop over all species
    for species in grs:
        for moment in species.findall('./Moment'):
            for line in moment.text.strip().split("\n"):
                r, v = map(float, line.split())
                yield (N, species.attrib["Name1"], species.attrib["Name2"], int(moment.attrib['Order']), of.numdensity(), r, v)

if output_file is None:
    st.error("Please load a data file to process")
else:
    of = pydynamo.OutputFile(output_file)

    
    #Grab the data into a big pandas dataframe
    gr_df = pd.DataFrame(getGrData(of), columns=["N", "Species 1", "Species 2", "Order","N/V", "R", "Value"])

    all_species = set(gr_df["Species 1"])
    all_densities = set(gr_df["N/V"])
    print("Found the following species ", all_species)
    print("Found the following densities ", all_densities)
    species1 = st.selectbox("Species 1", all_species)
    species2 = st.selectbox("Species 2", all_species)
    density = st.selectbox("Density", all_densities)
    st.text("V*="+str(math.sqrt(2)/density))
    gr_df = gr_df.set_index(["N", "Species 1", "Species 2", "Order","N/V", "R"])

    #Lets grab Series for each average/order taken
    avg_pow1 = gr_df.xs(0, level="Order")["Value"]
    avg_pow2 = gr_df.xs(1, level="Order")["Value"]
    avg_pow3 = gr_df.xs(2, level="Order")["Value"]

    #Now make a new dataframe with the A values, first is just a simple rename
    df = gr_df.xs(0, level="Order").rename(columns={"Value":"A1"}, inplace=False)
    df["A2"] = (avg_pow2 - avg_pow1**2) / 2 #This is the formula for the standard deviation
    df["A3"] = (avg_pow3 + 2 * (avg_pow1**3) - 3 * avg_pow1 * (avg_pow2))

    #Plot the values
    N=1372
    dfplot = df.loc[(N, species1, species2, density)]
    fig = go.Figure(
        data=[
            go.Scatter(x=dfplot.index, y=dfplot["A1"] / N, name="A1"),
            go.Scatter(x=dfplot.index, y=dfplot["A2"] / N, name="A2"),
            go.Scatter(x=dfplot.index, y=dfplot["A3"] / N, name="A3"),
        ],
        layout=go.Layout(
            title=go.layout.Title(text="Moments of the cumulative distribution function"),
        )
    )
    fig.update_layout(xaxis_title=r"$r / sigma$")
    st.plotly_chart(fig, use_container_width=True)

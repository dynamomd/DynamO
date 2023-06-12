import streamlit as st
import pydynamo
import math
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

import jax.numpy as jnp
import jax

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

def try_jit(f):
    try:
        ftmp = jax.jit(ftmp)
        ftmp(1.0, 1.0)
        return ftmp
    except:
        return f

class A_EOS:
    def __init__(self, Aexpr):
        self.a = try_jit(Aexpr)
        #dA/dV = -p
        self.p = try_jit(jax.grad(lambda V, kT : -A_FCC_HS_Young_Alder(V, kT), argnums=0))
        #dA/dT = -S
        self.s = try_jit(jax.grad(lambda V, kT : -A_FCC_HS_Young_Alder(V, kT), argnums=1))
        

def A_FCC_HS_Young_Alder(V, kT):
    V = V * math.sqrt(2)
    return kT * (-3 * jnp.log((V-1)/V)+5.124*jnp.log(V)-20.78*V + 9.52*V**2 -1.98*V**3+15.05)
EOS_FCC_HS_Young_Alder = A_EOS(A_FCC_HS_Young_Alder)

def A_Liq_HS_Young_Alder(V, kT):
    V = V * math.sqrt(2)
    y = math.pi * math.sqrt(2)/6/V
    return  kT * ((4 * y - 3 * y**2) / (1-y)**2 + jnp.log(y) + math.log(6 / math.pi / math.e))
EOS_Liq_HS_Young_Alder = A_EOS(A_Liq_HS_Young_Alder)


def A_HCP_HS_Young_Alder(V, kT):
    V = V * math.sqrt(2)
    if 1.0 <= V < 1.1:
        return A_FCC_HS_Young_Alder(V, kT) + 0.05 * (1.1 - V - jnp.log(1.1 / V)) + 0.005 * math.log(1.5 / 1.1)
    elif 1.1 <= V:
        return A_FCC_HS_Young_Alder(V, kT) + 0.005 * math.log(1.5 / V)
    else:
        raise RuntimeError("Out of range")
    
EOS_HCP_HS_Young_Alder = A_EOS(A_HCP_HS_Young_Alder)


for EOS in [EOS_FCC_HS_Young_Alder, EOS_HCP_HS_Young_Alder, EOS_Liq_HS_Young_Alder]:
    print(EOS.a(1.1, 1.0))
    print(EOS.p(1.1, 1.0))
    print(EOS.s(1.1, 1.0))


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

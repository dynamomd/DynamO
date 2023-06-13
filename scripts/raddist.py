import streamlit as st
import pydynamo
import math
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import scipy.special
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
    return kT * (-3 * jnp.log((V - 1) / V)+5.124 * jnp.log(V) - 20.78 * V + 9.52 * V**2 - 1.98 * V**3 + 15.05)
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
    elif 1.1 <= V < 1.5:
        return A_FCC_HS_Young_Alder(V, kT) + 0.005 * math.log(1.5 / V)
    else:
        raise RuntimeError("Out of range")
    
EOS_HCP_HS_Young_Alder = A_EOS(A_HCP_HS_Young_Alder)

for EOS in [EOS_FCC_HS_Young_Alder, EOS_HCP_HS_Young_Alder, EOS_Liq_HS_Young_Alder]:
    rho=0.9718
    V=1.0/rho
    kT=1.0
    print(EOS.a(V, kT))
    print(EOS.p(V, kT))
    print(EOS.s(V, kT))

if output_file is None:
    st.error("Please load a data file to process")
else:
    of = pydynamo.OutputFile(output_file)

    #Grab the data into a big pandas dataframe
    gr_df = pd.DataFrame(getGrData(of), columns=["N", "Species 1", "Species 2", "Order","N/V", "R", "Value"])

    all_species = set(gr_df["Species 1"])
    all_densities = set(gr_df["N/V"])
    all_N = set(gr_df["N"])
    max_moment = max(gr_df["Order"])+1

    print("Found the following species ", all_species)
    print("Found the following densities ", all_densities)
    print("Found the max moment ", max_moment)

    gr_df = gr_df.set_index(["N", "Species 1", "Species 2", "Order","N/V", "R"])

    # The Helmholtz free energy is related to the partition function, which is related to the configuration integral
    # -F/(k_BT) = ln Z = ln ⟨exp(-βU)⟩
    # This has the form of a cumulant expansion generating function K(t)=ln⟨exp(tX)⟩, thus we can perform the cumulant expansion

    # -F/(k_BT) = Σ_{n=1}^∞ κ_n (-β)^n / (n!)
    #The cumulants κ_n can be written in terms of the central moments #https://en.wikipedia.org/wiki/Cumulant#The_first_several_cumulants_as_functions_of_the_moments
    # κ_1 = ⟨U⟩
    # κ_2 = ⟨(U-⟨U⟩)^2⟩
    # κ_3 = ⟨(U-⟨U⟩)^3⟩
    # κ_4 = ⟨(U-⟨U⟩)^4⟩ - 3(⟨(U-⟨U⟩)^2⟩)^2 

    # In the simulation we have collected the moments of, N(r), the number of pairs below a radius r. 
    # These moments are collected about an origin/offset N₀(r), i.e. ⟨(N(r)-N₀(r))^n⟩
    # The origin is a numerical trick to reduce the size of the moments to reduce precision issues, its an
    # approximation of the mean, by taking the initial value of N(r) at t=0.

    # We eventually want to convert to the central moments ⟨(N(r)-⟨N(r)⟩)^n⟩, To
    # do this we need to determine the first moment ⟨N(r)⟩. There is a general
    # expression for transforming moments:
    #
    # ⟨(x-b)^n⟩ = Σ_{i=0}^n comb(n, i) ⟨(x-a)^i⟩ (a-b)^{n-i}
    #
    # comb is combination function in scipy.special (binomial coefficient).
    # But the answer for this first step is simple, x=N(r), b=0, a=N₀(r), n=1 gives the straightforward identity
    #
    # ⟨N(r)⟩ = ⟨N(r)-N₀(r)⟩ + N₀(r)

    # For all other steps n>2, b=⟨N(r)⟩, and a = N₀(r).
    # ⟨(N(r)-⟨N(r)⟩)^n⟩ = Σ_{i=0}^n comb(n, i) ⟨(N(r)-N₀(r))^i⟩ (N₀(r)-⟨N(r)⟩)^{n-i}

    #First, grab the moments ⟨(N(r)-N₀(r))^n⟩
    N_N0_n = [gr_df.xs(i, level="Order")["Value"] for i in range(max_moment)]
    #And the offset N₀
    N0 = gr_df.xs(-1, level="Order")["Value"]

    subs = ["₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉"]
    pwrs = ["⁰","¹", "²", "³", "⁴", "⁵", "⁶", "⁷","⁸", "⁹"]

    # Now make a new dataframe with all the data
    #
    # Start with the initial offset/origin to build the dataframe
    df = gr_df.xs(-1, level="Order").rename(columns={"Value":"N₀"}, inplace=False)
    df["⟨(N-N₀)"+pwrs[0]+"⟩"] = 1
    # Put all the simulation moments in too
    for i in range(1, max_moment+1):
        df["⟨(N-N₀)"+pwrs[i]+"⟩"] = gr_df.xs(i-1, level="Order")["Value"]
    # Calculate the first moment
    df["⟨N⟩"] = df["⟨(N-N₀)¹⟩"] + df["N₀"]

    ##Now make the other central moments
    for n in range(2, max_moment+1):
        df["⟨(N-⟨N⟩)"+pwrs[n]+"⟩"] = sum(scipy.special.comb(n, i) * df["⟨(N-N₀)"+pwrs[i]+"⟩"] * (df["N₀"]-df["⟨N⟩"])**(n-i) for i in range(n+1))

    if max_moment >= 1:
        df["κ₁"] = df["⟨N⟩"]
    if max_moment >= 2:
        df["κ₂"] = df["⟨(N-⟨N⟩)²⟩"]
    if max_moment >= 3:
        df["κ₃"] = df["⟨(N-⟨N⟩)³⟩"]
    if max_moment >= 4:
        df["κ₄"] = df["⟨(N-⟨N⟩)⁴⟩"] - 3 * (df["⟨(N-⟨N⟩)²⟩"] ** 2)
    if max_moment >= 5:
        df["κ₅"] = df["⟨(N-⟨N⟩)⁵⟩"] - 10 * df["⟨(N-⟨N⟩)³⟩"] * df["⟨(N-⟨N⟩)²⟩"]

    st.dataframe(df)

    ##Plot the values
    species1 = st.selectbox("Species 1", all_species)
    species2 = st.selectbox("Species 2", all_species)
    density = st.selectbox("Density", all_densities)
    st.text("V*="+str(math.sqrt(2)/density))
    N = st.selectbox("N", all_N)

    dfplot = df.loc[(N, species1, species2, density)]
    fig = go.Figure(
        data=[
            go.Scatter(x=dfplot.index, y=dfplot["κ₄"] / N, name="κ₄"),
            go.Scatter(x=dfplot.index, y=dfplot["κ₃"] / N, name="κ₃"),
            go.Scatter(x=dfplot.index, y=dfplot["κ₂"] / N, name="κ₂"),
            go.Scatter(x=dfplot.index, y=dfplot["κ₁"] / N, name="κ₁"),
        ],
        layout=go.Layout(
            title=go.layout.Title(text="Cumulants"),
        )
    )
    fig.update_layout(xaxis_title=r"<i>r / σ</i>", yaxis_title="<i>κ<sub>n</sub></i>")
    st.plotly_chart(fig, use_container_width=True)

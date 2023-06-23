import streamlit as st
import pydynamo
import math
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import scipy.special
import jax.numpy as jnp
import numpy as np
import jax
import glob
import datastat
import time 
from pandas.api.types import (
    is_categorical_dtype,
    is_datetime64_any_dtype,
    is_numeric_dtype,
    is_object_dtype,
)

import shapely
from shapely.geometry import LineString

from jax_spline import InterpolatedUnivariateSpline

subs = ["₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉"]
pwrs = ["⁰","¹", "²", "³", "⁴", "⁵", "⁶", "⁷","⁸", "⁹"]

def filter_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds a UI on top of a dataframe to let viewers filter columns

    Args:
        df (pd.DataFrame): Original dataframe

    Returns:
        pd.DataFrame: Filtered dataframe
    """
    modify = st.checkbox("Add filters")

    if not modify:
        return df

    df = df.reset_index(inplace=False)

    # Try to convert datetimes into a standard format (datetime, no timezone)
    for col in df.columns:
        if is_object_dtype(df[col]):
            try:
                df[col] = pd.to_datetime(df[col])
            except Exception:
                pass

        if is_datetime64_any_dtype(df[col]):
            df[col] = df[col].dt.tz_localize(None)

    modification_container = st.container()

    with modification_container:
        to_filter_columns = st.multiselect("Filter dataframe on", df.columns)
        for column in to_filter_columns:
            left, right = st.columns((1, 20))
            # Treat columns with < 10 unique values as categorical
            if is_categorical_dtype(df[column]) or df[column].nunique() < 10:
                user_cat_input = right.multiselect(
                    f"Values for {column}",
                    df[column].unique(),
                    default=list(df[column].unique()),
                )
                df = df[df[column].isin(user_cat_input)]
            elif is_numeric_dtype(df[column]):
                _min = float(df[column].min())
                _max = float(df[column].max())
                step = (_max - _min) / 100
                user_num_input = right.slider(
                    f"Values for {column}",
                    min_value=_min,
                    max_value=_max,
                    value=(_min, _max),
                    step=step,
                )
                df = df[df[column].between(*user_num_input)]
            elif is_datetime64_any_dtype(df[column]):
                user_date_input = right.date_input(
                    f"Values for {column}",
                    value=(
                        df[column].min(),
                        df[column].max(),
                    ),
                )
                if len(user_date_input) == 2:
                    user_date_input = tuple(map(pd.to_datetime, user_date_input))
                    start_date, end_date = user_date_input
                    df = df.loc[df[column].between(start_date, end_date)]
            else:
                user_text_input = right.text_input(
                    f"Substring or regex in {column}",
                )
                if user_text_input:
                    df = df[df[column].astype(str).str.contains(user_text_input)]

    return df


def getGrData(of : pydynamo.OutputFile):
    grs = of.tree.findall('.//RadialDistributionMoments/Species')
    N=int(of.tree.find('.//ParticleCount').attrib['val'])
    #Loop over all species
    for species in grs:
        for moment in species.findall('./Moment'):
            for line in moment.text.strip().split("\n"):
                r, v = map(float, line.split())
                yield (N, species.attrib["Name1"], species.attrib["Name2"], int(moment.attrib['Order']), of.numdensity(), r, v)

@st.cache_data
def get_df():
    return pd.concat(pd.DataFrame(getGrData(pydynamo.OutputFile(file)), columns=["N", "Species 1", "Species 2", "Order","N/V", "R", "Value"]) for file in glob.glob("/home/mjki2mb2/dynamo-repo/scripts/HS_TPT/*/1.data.xml.bz2", recursive=True))


def try_jit(f):
#    return f
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
        self.p = try_jit(jax.grad(lambda V, kT : -Aexpr(V, kT), argnums=0))
        #dA/dT = -S
        self.s = try_jit(jax.grad(lambda V, kT : -Aexpr(V, kT), argnums=1))
        
def A_FCC_HS_Young_Alder(V, kT):
    V = V * math.sqrt(2)
    return kT * (-3 * jnp.log((V - 1) / V)+5.124 * jnp.log(V) - 20.78 * V + 9.52 * V**2 - 1.98 * V**3 + 15.05)
EOS_FCC_HS_Young_Alder = A_EOS(A_FCC_HS_Young_Alder)

def A_Fluid_HS_Young_Alder(V, kT):
    V = V * math.sqrt(2)
    y = math.pi * math.sqrt(2)/6/V
    return  kT * ((4 * y - 3 * y**2) / (1-y)**2 + jnp.log(y) + math.log(6 / math.pi / math.e))
EOS_Fluid_HS_Young_Alder = A_EOS(A_Fluid_HS_Young_Alder)

def A_HCP_HS_Young_Alder(V, kT):
    Vstar = V * math.sqrt(2)
    if 1.0 <= Vstar < 1.1:
        return A_FCC_HS_Young_Alder(V, kT) + 0.05 * (1.1 - Vstar - jnp.log(1.1 / Vstar)) + 0.005 * math.log(1.5 / 1.1)
    elif 1.1 <= Vstar < 1.5:
        return A_FCC_HS_Young_Alder(V, kT) + 0.005 * jnp.log(1.5 / Vstar)
    else:
        raise RuntimeError("Out of range rho="+str(math.sqrt(2)/V))
    
EOS_HCP_HS_Young_Alder = A_EOS(A_HCP_HS_Young_Alder)

#for EOS in [EOS_FCC_HS_Young_Alder, EOS_HCP_HS_Young_Alder, EOS_Fluid_HS_Young_Alder]:
#    rho=0.9718
#    V=1.0/rho
#    kT=1.0
#    print(EOS.a(V, kT))
#    print(EOS.p(V, kT))
#    print(EOS.s(V, kT))
#
class A_TPT_EOS:
    def __init__(self, df, N, A_ref, s):
        import scipy.interpolate as si

        self.A_ref = A_ref
        self.terms = []
        self.dtermsdV = []
        self.ddtermsdV = []
        order = 1
        while 'A'+subs[order] in df:
            #Fit a spline in terms of V
            x=np.concatenate(([0], df.index.to_numpy()))
            y=np.concatenate(([0], df['A'+subs[order]].to_numpy() / N))
            #print(x,y)
            spline = si.UnivariateSpline(x=x, y=y, s=s)
            self.terms.append(spline)
            self.dtermsdV.append(spline.derivative())
            self.ddtermsdV.append(spline.derivative(2))
            order +=1

    def a(self, V, kT, nterms=None):
        beta = 1/kT
        rho = 1/V
        return self.A_ref.a(V, kT) + sum(term(rho) * beta**(idx+1) for idx, term in enumerate(self.terms[:nterms])) 

    def p(self, V, kT, nterms=None):
        beta = 1/kT
        #-rho**2 is the jacobian from V to rho
        rho = 1/V
        return self.A_ref.p(V, kT) + (rho**2) * sum(term(rho) * beta**(idx+1) for idx, term in enumerate(self.dtermsdV[:nterms])) 

    def s(self, V, kT, nterms=None):
        rho = 1/V
        return self.A_ref.s(V, kT) - sum(term(rho) * (-(idx+1)) * kT**(-(idx+2)) for idx, term in enumerate(self.terms[:nterms]))
    
    def mu(self, V, kT, nterms=None):
        return self.a(V, kT, nterms) + self.p(V, kT, nterms) * V

st.title("Radial Distribution Tool")
#output_file = st.file_uploader("Choose a Output file")
#if output_file is None:
#    st.error("Please load a data file to process")
#else:

tab_df, tab_Aterms, tab_phase = st.tabs(["Dataframe", "A terms", "Phase"])
import collections

if True:
    print("Loading data files")
    df = get_df()
    all_species = sorted(set(df["Species 1"]))
    all_densities = sorted(set(df["N/V"]))
    all_N = sorted(set(df["N"]))
    all_R = sorted(set(df["R"]))
    max_moment = max(df["Order"])+1
    df.set_index(["N", "Species 1", "Species 2", "Order","N/V", "R"], inplace=True)
    gr_df = df
    #print("Found the following species ", all_species)
    #print("Found the following densities ", all_densities)
    #print("Found the max moment ", max_moment)

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

    print("Processing moments")
    #First, grab the moments ⟨(N(r)-N₀(r))^n⟩
    N_N0_n = [gr_df.xs(i, level="Order")["Value"] for i in range(max_moment)]
    #And the offset N₀
    N0 = gr_df.xs(-1, level="Order")["Value"]


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

    well_depth = -1
    for n in range(1, max_moment+1):
        #The -1 arises from the negative sign on the beta in the cumulant expansion being moved in here to agree with Young_Alder's notation
        df["A"+subs[n]] = (-1) ** (n+1) * (well_depth)**n * df["κ"+subs[n]] / math.factorial(n)

    print("Creating dataframe tab")
    with tab_df:
        st.dataframe(filter_dataframe(df))

    print("Creating Aterms tab")
    with tab_Aterms:
        ##Plot the values
        N = st.selectbox("N", all_N)
        species1 = st.selectbox("Species 1", all_species)
        species2 = st.selectbox("Species 2", all_species)

    format_func = lambda x : "ρ="+str(datastat.roundSF(x, 12))+",V*="+str(datastat.roundSF(math.sqrt(2)/x, 12))
    with tab_Aterms:
        densities = st.multiselect("Density", all_densities, format_func=format_func)

    data = []
    for rho in densities:
        dfplot = df.loc[(N, species1, species2, rho)]
        data = data + [
            go.Scatter(x=dfplot.index, y=dfplot["A₄"] / N, name="A₄/N,"+format_func(rho)),
            go.Scatter(x=dfplot.index, y=dfplot["A₃"] / N, name="A₃/N,"+format_func(rho)),
            go.Scatter(x=dfplot.index, y=dfplot["A₂"] / N, name="A₂/N,"+format_func(rho)),
            go.Scatter(x=dfplot.index, y=dfplot["A₁"] / N, name="A₁/N,"+format_func(rho)),
        ]        
    fig = go.Figure(
        data=data,
        layout=go.Layout(
            title=go.layout.Title(text="Cumulants"),
        )
    )
    fig.update_layout(xaxis_title=r"<i>r / σ</i>", yaxis_title="<i>A<sub>n</sub>/N</i>")
    with tab_Aterms:
        st.plotly_chart(fig, use_container_width=True)


    print("Creating Phase tab")
    with tab_phase:
        Rval = st.selectbox("Lambda", all_R, format_func=lambda x: str(datastat.roundSF(x, 12)), index=149)

        #This is a slice of all the densities for a particular Rval. 
        print("Slicing and sorting the data")
        dfplot = df.loc[(N, species1, species2, slice(None), Rval)].sort_index()
        st.dataframe(dfplot)

        s = st.slider("Spline smoothing, 0=none", min_value=0.0, max_value=0.1, step=0.001, value=0.05)
        print("Building the equation of state")
        EOS_Fluid = A_TPT_EOS(dfplot, N, EOS_Fluid_HS_Young_Alder, s)
        EOS_FCC = A_TPT_EOS(dfplot, N, EOS_FCC_HS_Young_Alder, s)
        EOS_HCP = A_TPT_EOS(dfplot, N, EOS_HCP_HS_Young_Alder, s)

        print("Plotting the spline fit to An")
        n = st.slider("Theory order n", min_value=1, max_value=max_moment, step=1, value=2)
        xs=np.linspace(0.001, 1.4, 100)
        fig = go.Figure(
            data=[
                go.Scatter(x=dfplot.index, y=dfplot["A"+subs[n]]/N, name="A"+subs[n], mode="markers"),
                go.Scatter(x=xs, y=EOS_Fluid.terms[n-1](xs), name="fit A"+subs[n], mode="lines"),
                ],
            layout=go.Layout(
            title=go.layout.Title(text="A terms"),
            )
        )
        fig.update_layout(xaxis_title=r"<i> N σ³ / V</i>", yaxis_title="<i>A<sub>n</sub>/N</i>")
        st.plotly_chart(fig, use_container_width=True)
#
        print("Plotting the pressure ", end="")
        start = time.process_time()
        kT = st.slider("kT", min_value=0.1, max_value=5.0, step=0.1, value=1.0)
        xs=np.linspace(0.001, 1.4, 200)
        xs_solid=np.linspace(1.0, 1.4, 200)
        fig = go.Figure(
            data=[
                go.Scatter(x=xs, y=[EOS_Fluid.p(1/rho, kT, nterms=n) for rho in xs], name="Fluid", mode="lines"),
                go.Scatter(x=xs_solid, y=[EOS_FCC.p(1/rho, kT, nterms=n) for rho in xs_solid], name="FCC", mode="lines"),
                go.Scatter(x=xs_solid, y=[EOS_HCP.p(1/rho, kT, nterms=n) for rho in xs_solid], name="HCP", mode="lines"),
                ],
            layout=go.Layout(
            title=go.layout.Title(text="Pressure"),
            )
        )
        fig.update_layout(xaxis_title=r"<i> N σ³ / V</i>", yaxis_title="<i>p</i>")
        st.plotly_chart(fig, use_container_width=True)
        print(f"{time.process_time()-start:.2f} seconds")

        print("Calculating the mu-p loop ", end="")
        start = time.process_time()
        FluidData = np.column_stack(([EOS_Fluid.p(1 / rho, kT, nterms=n) for rho in xs],       [EOS_Fluid.mu(1/rho, kT, nterms=n) for rho in xs]))
        FCCData = np.column_stack(([EOS_FCC.p(1 / rho, kT, nterms=n) for rho in xs_solid], [EOS_FCC.mu(1/rho, kT, nterms=n) for rho in xs_solid]))
        HCPData = np.column_stack(([EOS_HCP.p(1 / rho, kT, nterms=n) for rho in xs_solid], [EOS_HCP.mu(1/rho, kT, nterms=n) for rho in xs_solid]))
        print(f"{time.process_time()-start:.2f} seconds")

        print("Solving for intersections ", end="")
        start = time.process_time()
        curves = ((FluidData, "Fluid"), (FCCData, "FCC"), (HCPData, "HCP"))
        for ci in range(len(curves)):
            Data1, name1 = curves[ci]
            for cj in range(ci, len(curves)):
                Data2, name2 = curves[cj]
                if name1 == name2:
                    for i in range(1, Data1.shape[0]):
                        for j in range(i+2, Data1.shape[0]):
                            l1 = LineString([Data1[i], Data1[i-1]])
                            l2 = LineString([Data1[j], Data1[j-1]])
                            if l1.intersects(l2):
                                print(name1,"-", name2, l1.intersection(l2))
                else:
                    for i in range(1, Data1.shape[0]):
                        for j in range(1, Data2.shape[0]):
                            l1 = LineString([Data1[i], Data1[i-1]])
                            l2 = LineString([Data2[j], Data2[j-1]])
                            if l1.intersects(l2):
                                print(name1,"-", name2, l1.intersection(l2))
        
        print(f"{time.process_time()-start:.2f} seconds")

        print("Plotting mu-p", end="")
        fig = go.Figure(
            data=[go.Scatter(x=Data[:,0], y=Data[:,1], name=name, mode="lines") for Data, name in ((FluidData, "Fluid"), (FCCData, "FCC"), (HCPData, "HCP"))],
            layout=go.Layout(
            title=go.layout.Title(text="Loop diagram"),
            )
        )
        fig.update_layout(xaxis_title=r"<i> p </i>", yaxis_title="<i>μ-μ_fluid</i>")
        st.plotly_chart(fig, use_container_width=True)
        print(f"{time.process_time()-start:.2f} seconds")

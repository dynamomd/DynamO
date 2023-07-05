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
import uncertainties
import uncertainties.unumpy
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
        self.dpdV = try_jit(jax.grad(self.p, argnums=0))
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

from functools import partial

class A_TPT_EOS:
    def __init__(self, df, A_ref, s):
        import scipy.interpolate as si

        self.A_ref = A_ref
        self.terms = []
        self.dtermsdV = []
        self.ddtermsdV = []
        order = 1
        while 'A'+subs[order] in df:
            #Fit a spline in terms of V
            x=np.concatenate(([0.0], df.index.to_numpy()))
            y=np.concatenate(([0.0], uncertainties.unumpy.nominal_values(df['A'+subs[order]])))
            yerr=np.concatenate(([0.0], uncertainties.unumpy.std_devs(df['A'+subs[order]])))
            #Some value have zero uncertainty which doesn't work with fitting, so here we give them the minimum uncertainty
            yerr_min = np.ma.masked_equal(yerr, 0.0, copy=False).min()
            yerr[yerr == 0] = yerr_min
            #Weights are 1/std_dev, so we follow scipy.interpolate.UnivariateSpline guidance for s which is equal to data count.
            spline = si.UnivariateSpline(x=x, y=y, w=1/yerr, s= s * 4 * y.shape[0])
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

    def dpdV(self, V, kT, nterms=None):
        beta = 1/kT
        rho = 1/V
        return self.A_ref.p(V, kT) - 2 * (rho**3) * sum(term(rho) * beta**(idx+1) for idx, term in enumerate(self.dtermsdV[:nterms]))  - (rho**4) * sum(term(rho) * beta**(idx+1) for idx, term in enumerate(self.ddtermsdV[:nterms])) 

    def s(self, V, kT, nterms=None):
        rho = 1/V
        return self.A_ref.s(V, kT) - sum(term(rho) * (-(idx+1)) * kT**(-(idx+2)) for idx, term in enumerate(self.terms[:nterms]))
    
    def mu(self, V, kT, nterms=None):
        return self.a(V, kT, nterms) + self.p(V, kT, nterms) * V
    
    def u(self, V, kT, nterms=None):
        return self.a(V, kT, nterms) + kT * self.s(V, kT, nterms)

def fixup_values(values):
   #Get a unique list of values without small numerical differences (10s.f.) in sorted order
   return list(sorted(set(float(f'{float(f"{v:.10g}"):g}') for v in values)))


@st.cache_resource
def get_df():
    import pickle
    max_moment=5
    print("  Running pickle load")

    df = pickle.load(open('/home/mjki2mb2/dynamo-repo/scripts/HS_TPT.pkl', 'rb'))
    df.set_index(["N", "ndensity", "InitState"], inplace=True)
    all_N = fixup_values(df.index.get_level_values("N"))

    #Selecting N
    print("  Dropping unneeded columns")
    df.drop(["Lambda", "PhiT", "Rso", "kT"], axis=1, inplace=True)

    all_densities = fixup_values(df.index.get_level_values("ndensity"))
    all_InitState = sorted(set(df.index.get_level_values("InitState")))

    #We need to unpack the moments into a new dataframe
    def moment_to_df(index, row):
        #row by row, grab the moments as a numpy array
        moments = row['RadialDistribution'].ufloat()
        #Create the array of R values that index the rows of the moments
        Rvals = np.linspace(bin_width, bin_width*moments.shape[1], moments.shape[1], endpoint=True)
        #Stack this together with the moment values
        data = np.hstack((Rvals.reshape((-1,1)), np.transpose(moments)))
        #We then build a dataframe with this information
        df_m = pd.DataFrame(data, columns=["R"]+["⟨(N-⟨N⟩)"+pwrs[i]+"⟩" for i in range(1, moments.shape[0]+1)])
        #We add the index now individually, as smarter tricks with numpy will promote the float of ndensity to string to match InitState, or something equally stupid.
        for k, v in zip(df.index.names, index):
            df_m[k] = v
        df_m.drop(df_m[df_m['R'] <= 1.0].index, inplace=True)
        df_m.drop(df_m[df_m['R'] > 2.0].index, inplace=True)
        df_m.set_index(df.index.names+["R"], inplace=True)
        return df_m
    
    print("  Unpacking moments")
    df_moments = pd.concat([moment_to_df(index, row) for index, row in df.iterrows()])

    all_Rs = fixup_values(df_moments.index.get_level_values("R"))

    # The Helmholtz free energy is related to the partition function, which is related to the configuration integral
    # -F/(k_BT) = ln Z = ln ⟨exp(-βU)⟩
    # This has the form of a cumulant expansion generating function K(t)=ln⟨exp(tX)⟩, thus we can perform the cumulant expansion

    # -F/(k_BT) = Σ_{n=1}^∞ κ_n (-β)^n / (n!)
    #The cumulants κ_n can be written in terms of the central moments #https://en.wikipedia.org/wiki/Cumulant#The_first_several_cumulants_as_functions_of_the_moments
    # κ_1 = ⟨U⟩
    # κ_2 = ⟨(U-⟨U⟩)^2⟩
    # κ_3 = ⟨(U-⟨U⟩)^3⟩
    # κ_4 = ⟨(U-⟨U⟩)^4⟩ - 3(⟨(U-⟨U⟩)^2⟩)^2 


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

    print("  Generating cumulants")
    if "⟨(N-⟨N⟩)¹⟩" in df_moments:
        df_moments["κ₁"] = df_moments["⟨(N-⟨N⟩)¹⟩"]
    if "⟨(N-⟨N⟩)²⟩" in df_moments:
        df_moments["κ₂"] = df_moments["⟨(N-⟨N⟩)²⟩"]
    if "⟨(N-⟨N⟩)³⟩" in df_moments:
        df_moments["κ₃"] = df_moments["⟨(N-⟨N⟩)³⟩"]
    if "⟨(N-⟨N⟩)⁴⟩" in df_moments:
        df_moments["κ₄"] = df_moments["⟨(N-⟨N⟩)⁴⟩"] - 3 * (df_moments["⟨(N-⟨N⟩)²⟩"] ** 2)
    if "⟨(N-⟨N⟩)⁵⟩" in df_moments:
        df_moments["κ₅"] = df_moments["⟨(N-⟨N⟩)⁵⟩"] - 10 * df_moments["⟨(N-⟨N⟩)³⟩"] * df_moments["⟨(N-⟨N⟩)²⟩"]

    print("  Generating Helmholtz terms")
    well_depth = -1
    for n in range(1, max_moment+1):
        if "κ"+subs[n] not in df_moments:
            break
        #The -1 arises from the negative sign on the beta in the cumulant expansion being moved in here to agree with Young_Alder's notation
        df_moments["A"+subs[n]] = (-1) ** (n+1) * (well_depth)**n * df_moments["κ"+subs[n]] / math.factorial(n) / df_moments.index.get_level_values("N")

    return df, all_N, all_densities, all_InitState, all_Rs, df_moments

st.title("Radial Distribution Tool")
tab_df, tab_Aterms, tab_phase = st.tabs(["Dataframe", "A terms", "Phase"])
import collections

if True:
    max_moment = 5
    bin_width = 0.01
    with tab_df:
        print("Loading data files")
        df, all_N, all_densities, all_InitState, all_Rs, df_moments = get_df()

        print("Creating dataframe tab", flush=True)
        start = time.process_time()
        st.dataframe(filter_dataframe(df_moments))
        print(f"{time.process_time()-start:.2f} seconds")

    print("Creating Aterms tab")
    start = time.process_time()
    format_func = lambda x : "ρ="+str(datastat.roundSF(x, 12))+",V*="+str(datastat.roundSF(math.sqrt(2)/x, 12))
    with tab_Aterms:
        densities = st.multiselect("Density", all_densities, format_func=format_func)
    data = []
    for rho in densities:
        dfplot_FCC = df_moments.loc[(rho, 'FCC')]
        dfplot_HCP = df_moments.loc[(rho, 'HCP')]
        data = data + [go.Scatter(x=dfplot_FCC.index, y=uncertainties.unumpy.nominal_values(dfplot_FCC["A"+subs[i]]), error_y=dict(type="data", array=uncertainties.unumpy.std_devs(dfplot_FCC["A"+subs[i]]), visible=True), name="A"+subs[i]+"/N,FCC,"+format_func(rho)) for i in range(1,5)]
    
    fig = go.Figure(
        data=data,
        layout=go.Layout(
            title=go.layout.Title(text="Cumulants"),
        )
    )
    fig.update_layout(xaxis_title=r"<i>r / σ</i>", yaxis_title="<i>A<sub>n</sub>/N</i>")

    with tab_Aterms:
        st.plotly_chart(fig, use_container_width=True)
    print(f"{time.process_time()-start:.2f} seconds")

    print("Creating Phase tab")
    with tab_phase:
        Rval = st.selectbox("Lambda", all_Rs, format_func=lambda x: str(datastat.roundSF(x, 12)), index=50)

        #This is a slice of all the densities for a particular Rval. 
        print("Slicing and sorting the data")
        dfplot = df_moments.reset_index()
        dfplot = dfplot[dfplot['R'].between(Rval-bin_width*0.5, Rval+bin_width*0.5)]
        
        st.dataframe(dfplot)
        
        rho_minsolid = 0.95
        rho_maxfluid = 1.0

        s = st.slider("Spline smoothing", min_value=0.0, max_value=1.0, step=0.01, value=0.5, format="%g")
        print("Building the equations of state")
        start = time.process_time()

        phases = {}
        for name, InitState, solid, BaseEOS in [("Fluid", "SC", False, EOS_Fluid_HS_Young_Alder), ("FCC", "FCC", True, EOS_FCC_HS_Young_Alder), ("HCP", "HCP", True, EOS_HCP_HS_Young_Alder)]:
            dfphase = dfplot[dfplot["InitState"]==InitState]
            if solid:
                filter_points = dfphase['ndensity'] > rho_minsolid
                xs = np.linspace(rho_minsolid, 1.4, 201)
            else:
                filter_points = dfphase['ndensity'] < rho_maxfluid
                xs = np.linspace(1e-12, rho_maxfluid, 201)
            datapoints = dfphase[filter_points].set_index("ndensity").sort_index()
            phases[name] = dict(InitState=InitState, solid=solid, datapoints=datapoints, EOS=A_TPT_EOS(datapoints, BaseEOS, s), xs=xs)
        print(f"{time.process_time()-start:.2f} seconds")

        print("Plotting the spline fit to An")
        start = time.process_time()
        n = st.slider("Theory order n", min_value=1, max_value=max_moment, step=1, value=2)
        fig = go.Figure(
            data=[go.Scatter(x=d['datapoints'].index, y=uncertainties.unumpy.nominal_values(d['datapoints']["A"+subs[n]]), error_y=dict(type="data", array=uncertainties.unumpy.std_devs(d['datapoints']["A"+subs[n]]), visible=True), name="A"+subs[n]+","+name, mode="markers") for name, d in phases.items()]\
                +[go.Scatter(x=d['xs'], y=d['EOS'].terms[n-1](d['xs']), name=name+" A"+subs[n], mode="lines") for name, d in phases.items()],
            layout=go.Layout(title=go.layout.Title(text="A terms"),),
        )
        fig.update_layout(xaxis_title=r"<i> N σ³ / V</i>", yaxis_title="<i>A<sub>n</sub>/N</i>")
        st.plotly_chart(fig, use_container_width=True)
        print(f"{time.process_time()-start:.2f} seconds")

        kT = st.slider("kT", min_value=0.1, max_value=5.0, step=0.01, value=1.1)

        def find_tielines(kT, timeit=False):
            #Calculate all the segments for each EOS
            def generator():
                for name, d in phases.items():
                    last = None
                    for idx,rho in enumerate(d['xs']):
                        next = (rho, (float(d['EOS'].p(1 / rho, kT, nterms=n)), float(d['EOS'].mu(1/rho, kT, nterms=n)), float(d['EOS'].dpdV(1/rho, kT, nterms=n))))
                        if last != None:
                            if last[1][0] < next[1][0]:
                                yield (name, last[0], next[0], *(last[1]), *(next[1]), idx)
                            else:
                                yield (name, next[0], last[0], *(next[1]), *(last[1]), idx)
                        last = next
            if timeit:
                print("Calculating the segments ", flush=True)
                start = time.process_time()
            #Create the dataframe so this appears easier to read
            df = pd.DataFrame(generator(),
                              columns=["phase", "ρ1", "ρ2", "p1", "μ1", "dp/dv1", "p2", "μ2", "dp/dv2", "phaseindex"])
            if timeit:
                print(f"{time.process_time()-start:.2f} seconds")

            #Here we also sort by pressure and chemical potential, so that when
            #we iterate, we pull points of increasing pressure.
            df.sort_values(by=["p1","μ1"], inplace=True)
            #Filter out any segments including unstable points. 
            df_scan = df[(df['dp/dv1'] > 0) & (df['dp/dv2'] > 0)]

            if timeit:
                print("Scanning segments for transitions", flush=True)
                start = time.process_time()
            #We scan up in segments, with "active" segments being ones that had
            #upper pressure limits above the previous segment's lower limit.
            active_segments = [] 
            active_tielines = []
            tielines = []
            for index, seg1 in df_scan.iterrows():
                #print("Checking index", index, list(seg1))
                #Check if any transitions have been determined to be stable or unstable
                temp_active_stable_transitions = active_tielines
                active_tielines = []
                for transition in temp_active_stable_transitions:
                    p,mu,density1, density2, phase1, phase2 = transition
                    if p < seg1['p1']:
                        #We've moved passed this transition in pressure, so it must be stable
                        tielines.append(transition + [True])
                    else:
                        mu_seg1 = seg1['μ1'] + (seg1['μ2'] - seg1['μ1']) * (p - seg1['p1']) / (seg1['p2'] - seg1['p1'])
                        if mu_seg1 > mu:
                            #Its still stable, so keep it
                            active_tielines.append(transition)
                        else:
                            #Its not stable, so discard
                            tielines.append(transition + [False])

                #Now we check if this new segment brings in a transition. This
                #will happen exactly once for each crossing.
                l1 = LineString([(seg1["p1"], seg1["μ1"]), (seg1["p2"], seg1["μ2"])])
                for aseg_idx in active_segments:
                    seg2 = df.loc[aseg_idx]
                    if seg1['phase'] == seg2['phase'] and abs(seg1['phaseindex'] - seg2['phaseindex']) == 1:
                        #We skip adjacent line segments of the same phase, as they always intersect at their common points.
                        continue
                    l2 = LineString([(seg2["p1"], seg2["μ1"]), (seg2["p2"], seg2["μ2"])])
                    if l1.intersects(l2):
                        #We have an intersection/tie line. Calculate its properties
                        intersection = l1.intersection(l2)
                        p = intersection.x
                        mu = intersection.y
                        frac1 = (p - seg1['p1']) / (seg1['p2'] - seg1['p1'])
                        frac2 = (p - seg2['p1']) / (seg2['p2'] - seg2['p1'])
                        density1 = frac1 * (seg1['ρ2'] - seg1['ρ1']) + seg1['ρ1']
                        density2 = frac2 * (seg2['ρ2'] - seg2['ρ1']) + seg2['ρ1']
                        #Check that this tieline is on sensible parts of the EOS
                        if phases[seg1['phase']]['EOS'].dpdV(1/density1, kT, nterms=n) < 0:
                            continue
                        if phases[seg2['phase']]['EOS'].dpdV(1/density2, kT, nterms=n) < 0:
                            continue
                        if phases[seg1['phase']]['EOS'].p(1/density1, kT, nterms=n) < 0:
                            continue
                        if phases[seg2['phase']]['EOS'].p(1/density2, kT, nterms=n) < 0:
                            continue
                        #We need to check that mu is below all other tie lines (skip the intersecting one)
                        stable=True
                        for seg3_idx in active_segments:
                            if seg3_idx == aseg_idx:
                                continue
                            seg3 = df.loc[seg3_idx]
                            if seg3['p1'] <= p <= seg3['p2']:
                                #OK, check if the transition is above it
                                mu_seg3 = seg3['μ1'] + (seg3['μ2'] - seg3['μ1']) * (p - seg3['p1']) / (seg3['p2'] - seg3['p1'])
                                if mu_seg3 < mu:
                                    stable = False
                                    break
                        if stable:
                            if density1 < density2:
                                active_tielines.append([p,mu, density1, density2, seg1['phase'], seg2['phase']])
                            else:
                                active_tielines.append([p,mu, density2, density1, seg2['phase'], seg1['phase']])
                
                #First, only keep segments that are overlapping with the current one
                active_segments = list(filter(lambda idx: df.loc[idx]['p2'] >= seg1['p1'], active_segments))
                active_segments.append(index)

            if timeit:
                print(f"{time.process_time()-start:.2f} seconds")
            #Resort the curve values for plotting as curves
            df.sort_values(by=["phase", "phaseindex"], inplace=True)
            return pd.DataFrame(tielines, columns=['p', 'μ', 'ρ1', 'ρ2', 'name1', 'name2', 'stable']), df

        print("Solving for the tielines and data ", flush=True)
        start = time.process_time()
        tielines, curve_points = find_tielines(kT, timeit=True)
        print(f"{time.process_time()-start:.2f} seconds")
        #print("Discovered tielines\n", tielines)

        print("Plotting the pressures")
        fig = go.Figure(
            data=[go.Scatter(x=curve_points[curve_points['phase'] == name]["ρ1"], y=curve_points[curve_points['phase'] == name]["p1"], name=name, mode="lines") for name, d in phases.items()],
            layout=go.Layout(
                title=go.layout.Title(text="Pressure"),
            )
        )
        fig.update_layout(xaxis_title=r"<i> N σ³ / V</i>", yaxis_title="<i>p</i>")
        st.plotly_chart(fig, use_container_width=True)

        stable_tielines = tielines[tielines['stable'] == True]
        unstable_tielines = tielines[tielines['stable'] == False]
        fig = go.Figure(
            data=[go.Scatter(x=curve_points[curve_points['phase'] == name]["p1"], y=curve_points[curve_points['phase'] == name]["μ1"], name=name, mode="lines") for name, d in phases.items()]
                +[go.Scatter(x=stable_tielines['p'], 
                             y=stable_tielines['μ'],
                             text=stable_tielines['name1']+"→"+stable_tielines['name2'],
                             mode="markers+text",
                             name="Stable Transitions")]
                +[go.Scatter(x=unstable_tielines['p'], 
                             y=unstable_tielines['μ'],
                             text=unstable_tielines['name1']+"→"+unstable_tielines['name2'],
                             mode="markers+text",
                             name="Unstable Transitions")],
            layout=go.Layout(
            title=go.layout.Title(text="Loop diagram"),
            )
        )
        fig.update_layout(xaxis_title=r"<i> p </i>", yaxis_title="<i>μ</i>")
        st.plotly_chart(fig, use_container_width=True)

        #print("Solving for the phase diagram", flush=True)
        #tielines = pd.concat([find_tielines(kT)[0].assign(kT=kT) for kT in np.linspace(0.25, 2.5, 50)])
        #groupd_ties = tielines.groupby(['name1', 'name2'])
        #fig = go.Figure(
        #    data=[go.Scatter(x=list(group['ρ1'])+list(group['ρ2']),y=list(group['kT'])+list(group['kT']), mode='markers', name=key[0]+'⇌'+key[1]) for key, group in groupd_ties],
        #    layout=go.Layout(
        #    title=go.layout.Title(text="Phase diagram"),
        #    )
        #)
        #fig.update_layout(xaxis_title=r"<i> ρ</i>", yaxis_title="<i>k<sub>B</sub>T</i>")
        #st.plotly_chart(fig, use_container_width=True)

        print("DONE!")
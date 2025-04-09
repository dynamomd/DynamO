import math
import sys

import numpy
import uncertainties

def roundSF(x, n):
    """Return 'x' rounded to 'n' significant digits."""
    y=abs(x)
    if y <= sys.float_info.min: return 0.0
    return round( x, int( n-math.ceil(math.log10(y)) ) )

def print_to_14sf(f):
    """Utility function to print a variable to 14 significant figures"""
    if isinstance(f, str):
        return f
    return '{:.{p}g}'.format(f, p=14)

def conv_to_14sf(f):
    if isinstance(f, float):
        return float(print_to_14sf(f))
    else:
        return f

def simpson_impl(x, f):
    """Simpson rule for irregularly spaced data. Implementation is to
    allow uncertainties as the points for integration.

        Parameters
        ----------
        x : list or np.array of floats
                Sampling points for the function values
        f : list or np.array of floats
                Function values at the sampling points

        Returns
        -------
        float : approximation for the integral

    """
    N = len(x) - 1
    h = numpy.diff(x)

    result = 0.0
    for i in range(1, N, 2):
        hph = h[i] + h[i - 1]
        result += f[i] * ( h[i]**3 + h[i - 1]**3
                           + 3. * h[i] * h[i - 1] * hph )\
                     / ( 6 * h[i] * h[i - 1] )
        result += f[i - 1] * ( 2. * h[i - 1]**3 - h[i]**3
                              + 3. * h[i] * h[i - 1]**2)\
                     / ( 6 * h[i - 1] * hph)
        result += f[i + 1] * ( 2. * h[i]**3 - h[i - 1]**3
                              + 3. * h[i - 1] * h[i]**2)\
                     / ( 6 * h[i] * hph )

    if (N + 1) % 2 == 0:
        result += f[N] * ( 2 * h[N - 1]**2
                          + 3. * h[N - 2] * h[N - 1])\
                     / ( 6 * ( h[N - 2] + h[N - 1] ) )
        result += f[N - 1] * ( h[N - 1]**2
                           + 3*h[N - 1]* h[N - 2] )\
                     / ( 6 * h[N - 2] )
        result -= f[N - 2] * h[N - 1]**3\
                     / ( 6 * h[N - 2] * ( h[N - 2] + h[N - 1] ) )
    return result

def trapz_impl(x, f):
    """Just a normal trapezium rule, which supports uncertainties.ufloat
    as points for integration."""
    retval = 0
    if len(x) < 2:
        raise RuntimeError("Cannot integrate less than 2 points")
    N = len(x) - 1
    for i in range(N):
        retval += (x[i+1] - x[i]) * (f[i] + f[i+1]) / 2
    return retval

def simpson_integrate(x,f):
    """A simpson integration rule which propogates error and adds an
    estimate of the integration truncation error."""
    trapval = trapz_impl(x, f)
    simsval = simpson_impl(x, f)
    int_error_estimate = abs(simsval - trapval).nominal_value
    overall_error = (simsval.std_dev**2 + int_error_estimate**2)**0.5
    return uncertainties.ufloat(simsval.nominal_value, overall_error), int_error_estimate, simsval.std_dev

def linear_interp(xdata, ydata, xval = 0):
    """A linear interpolation function that utilises the uncertainties of
    the ydata to weight each point. This also propogates the
    uncertainty in the fit parameters too.
    """
    
    from scipy.optimize import curve_fit
    f1 = lambda x,a,b : b * x + a
    f2 = lambda x,a,b,c : c * x * x + b * x + a
    from uncertainties import ufloat
    linfit_popt,  linfit_pcov  = curve_fit(f = f1, xdata = xdata, ydata = [y.nominal_value for y in ydata], sigma = [y.std_dev for y in ydata], absolute_sigma=True)
    quadfit_popt, quadfit_pcov = curve_fit(f = f2, xdata = xdata, ydata = [y.nominal_value for y in ydata], sigma = [y.std_dev for y in ydata], absolute_sigma=True)

    lincoeffs = [ufloat(linfit_popt[i], linfit_pcov[i,i]) for i in range(2)]
    quadcoeffs = [ufloat(quadfit_popt[i], quadfit_pcov[i,i]) for i in range(3)]

    linval = f1(xval, *lincoeffs)
    quadval = f2(xval, *quadcoeffs)
    diff = linval.nominal_value - quadval.nominal_value
    return ufloat(linval.nominal_value, math.sqrt(diff**2 + linval.std_dev**2)), lincoeffs
    
def split_unc(df):
    """Method to split uncertainties values inside pandas dataframes into
    their uncertainty and average (for easier interfacing with external
    tools, like CSV)"""
    import pandas as pd
    import uncertainties
    odf = pd.DataFrame()
    for (colname, coldata) in df.iteritems():
        anyunc=False
        def convert_to_val(val):
            if isinstance(val, uncertainties.UFloat):
                anyunc = True
                return val.nominal_value
            else:
                return 0

        def convert_to_unc(val):
            if isinstance(val, uncertainties.UFloat):
                return val.std_dev
            else:
                return 0

        newdata = coldata.map(convert_to_val)
        if anyunc:
            odf[colname] = newdata
            odf[colname+' unc'] = coldata.map(convert_to_unc)
        else:
            odf[colname] = coldata
    return odf

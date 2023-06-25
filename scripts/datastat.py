import uncertainties, uncertainties.unumpy, pandas, math, numpy, sys

def roundSF(x, n):
    """Return 'x' rounded to 'n' significant digits."""
    y=abs(x)
    if y <= sys.float_info.min: return 0.0
    return round( x, int( n-math.ceil(math.log10(y)) ) )

class WeightedFloat():
    '''This class implements weighted arithmetic means along with an
    estimate of the standard error for the mean
    https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Bootstrapping_validation
    that is unbiased
    https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Reliability_weights.

    It eventually returns via its ufloat() method, a uncertainty float
    (ufloat), containing the average and error estimate for the
    average.
    '''

    def __init__(self, value = 0, weight = 0):
        """Initialise the weighted value."""
        self._ww_vv_sum = weight * weight * value * value
        self._ww_v_sum = weight * weight * value
        self._ww_sum = weight * weight
        self._w_vv_sum = weight * value * value
        self._w_v_sum = weight * value
        self._w_sum = weight
        self._count = 1

    def __add__(self, v):
        if not isinstance(v, WeightedFloat):
            raise RuntimeError("Cannot add non-WeightedFloat to WeightedFloat")

        retval = v
        retval._ww_vv_sum += self._ww_vv_sum
        retval._ww_v_sum += self._ww_v_sum
        retval._ww_sum += self._ww_sum
        retval._w_vv_sum += self._w_vv_sum
        retval._w_v_sum += self._w_v_sum
        retval._w_sum += self._w_sum
        retval._count += self._count

        return retval

    def avg(self):
        if self._w_sum == 0:
            return 0
        return self._w_v_sum / self._w_sum

    def std_error(self):
        avg = self.avg()
        if self._w_sum == 0:
            #We actually have no data, no error with this estimate!
            return float('nan')
        bias = (1 - self._ww_sum / (self._w_sum * self._w_sum))
        if bias == 0:
            #We have one sample, so infinite error?
            return float('nan')
        unbiased_mean_stderror_sq = (self._ww_vv_sum - 2 * self._ww_v_sum * avg + self._ww_sum * avg * avg) / (self._w_sum * self._w_sum * bias)
        if unbiased_mean_stderror_sq < 0:
            # Sometimes, round-off error for identical values gives almost zero, but negative, error
            unbiased_mean_stderror_sq = 0
        return math.sqrt(unbiased_mean_stderror_sq)

    def std_dev(self):
        return self.std_error() * math.sqrt(self._w_sum)
    
    def ufloat(self):
        if self._w_sum == 0:
            return uncertainties.ufloat(0, float('nan'))
        else:
            return uncertainties.ufloat(self.avg(), self.std_error())

    def __str__(self):        
        return str(self.ufloat())

    def __repr__(self):
        return repr(self.ufloat())

import numpy
class WeightedArray():
    '''This class implements weighted arithmetic means along with an
    estimate of the standard error for the mean
    https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Bootstrapping_validation
    that is unbiased
    https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Reliability_weights.

    It eventually returns via its ufloat() method, a uncertainty float
    (ufloat), containing the average and error estimate for the
    average.
    '''

    def __init__(self, value = numpy.array([0]), weight = numpy.array([0])):
        """Initialise the weighted value."""
        ww = weight * weight
        vv = numpy.multiply(value, value)
        self._ww_vv_sum = ww * vv
        self._ww_v_sum = ww * value
        self._ww_sum = ww
        self._w_vv_sum = weight * vv
        self._w_v_sum = weight * value
        self._w_sum = weight
        self._count = 1

    def __add__(self, v):
        if not isinstance(v, WeightedArray):
            raise RuntimeError("Cannot add non-WeightedArray to WeightedArray")

        retval = v
        if self._w_sum == 0:
            #Shortcut if this has no value, then just take the value passed in (and its shape!)
            return retval

        retval._ww_vv_sum += self._ww_vv_sum
        retval._ww_v_sum += self._ww_v_sum
        retval._ww_sum += self._ww_sum
        retval._w_vv_sum += self._w_vv_sum
        retval._w_v_sum += self._w_v_sum
        retval._w_sum += self._w_sum
        retval._count += self._count

        return retval

    def avg(self):
        if self._w_sum == 0:
            return 0
        return self._w_v_sum / self._w_sum

    def std_error(self):
        avg = self.avg()
        if self._w_sum == 0:
            #We actually have no data, no error with this estimate!
            return float('nan')
        bias = (1 - self._ww_sum / (self._w_sum * self._w_sum))
        if bias == 0:
            #We have one sample, so infinite error?
            return float('nan')
        unbiased_mean_stderror_sq = (self._ww_vv_sum - 2 * numpy.multiply(self._ww_v_sum, avg) + self._ww_sum * numpy.multiply(avg, avg)) / (self._w_sum * self._w_sum * bias)
    
        # Sometimes, round-off error for identical values gives almost zero, but negative, error
        return numpy.sqrt(numpy.maximum(unbiased_mean_stderror_sq, numpy.zeros_like(unbiased_mean_stderror_sq)))

    def std_dev(self):
        return self.std_error() * math.sqrt(self._w_sum)
    
    def ufloat(self):
        return uncertainties.unumpy.uarray(self.avg(), self.std_error())

    def __str__(self):
        return str(self.ufloat())

    def __repr__(self):
        return repr(self.ufloat())


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

if __name__ == "__main__":
    from scipy.optimize import curve_fit
    x = [2,4,8]
    yavg = [2.2, 2.4, 2.8]
    ysig = [1, 2, 3]
    
    results = []
    for i in range(20000):
        y = [yavg[j] + numpy.random.normal(0, ysig[j]) for j in range(len(yavg))]
        popt, pcov = curve_fit(lambda x,a,b: a + b * x, x, y)
        results.append(popt[0])

    from statistics import mean, stdev
    print(mean(results),"+/-",stdev(results))
    print(linear_interp(x, [uncertainties.ufloat(y,dy) for y,dy in zip(yavg, ysig)]))
    print("Tests")

    a = WeightedArray(numpy.array([1, 2]), 1)
    b = WeightedArray(numpy.array([2, 4]), 0.1)
    c = a+b
    c1 = WeightedFloat(1, 1) + WeightedFloat(2, 0.1)
    c2 = WeightedFloat(2, 1) + WeightedFloat(4, 0.1)
    print(c,c1,c2)

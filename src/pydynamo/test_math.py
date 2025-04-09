from pydynamo.math import *


def test_math():
    from statistics import mean, stdev

    from scipy.optimize import curve_fit
    x = [2,4,8]
    yavg = [2.2, 2.4, 2.8]
    ysig = [1, 2, 3]
    
    results = []
    for i in range(200):
        y = [yavg[j] + numpy.random.normal(0, ysig[j]) for j in range(len(yavg))]
        popt, pcov = curve_fit(lambda x,a,b: a + b * x, x, y)
        results.append(popt[0])

    result, coeffs = linear_interp(x, [uncertainties.ufloat(y,dy) for y,dy in zip(yavg, ysig)])
    # Currently I don't check anything here. I'm not even sure what this code was meant to do?
    # TODO - Write meaningful tests for the math functions?
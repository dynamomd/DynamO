import math

import numpy
import uncertainties
import uncertainties.unumpy


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

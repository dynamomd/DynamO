import math

import numpy
import uncertainties
import uncertainties.unumpy


class WeightedFloat():
    '''This class implements weighted arithmetic means along with an
    estimate of the standard error for the mean.
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

    def stats(self):
        if self._w_sum == 0:
            # We actually have no weighted data!
            return float('nan'), float('nan'), float('nan')

        avg = self._w_v_sum / self._w_sum

        var_per_time = (1 / self._count) * (self._w_sum * avg * avg - 2 * avg * self._w_v_sum + self._w_vv_sum)

        avg_var = var_per_time / self._w_sum
        if avg_var < 0:
            # Sometimes, round-off error for identical values gives almost zero, but negative, error
            avg_var = 0
        
        return avg, math.sqrt(avg_var), var_per_time

    def avg(self):
        if self._w_sum == 0:
            return float('nan')
        return self._w_v_sum / self._w_sum

    def std_dev(self):
        '''Return the standard deviation of the mean.'''
        avg, std_dev, _ = self.stats()
        return std_dev
    
    def std_error(self):
        '''Return the standard error of the mean.'''
        avg, std_dev, _ = self.stats()
        return std_dev / math.sqrt(self._w_sum)

    def ufloat(self):
        '''Return the average and error as a ufloat.'''
        avg, std_dev, _ = self.stats()
        return uncertainties.ufloat(avg, std_dev)

    def __str__(self):        
        return str(self.ufloat())

    def __repr__(self):
        return repr(self.ufloat())

class WeightedArray():
    '''This class implements weighted arithmetic means along with an
    estimate of the standard error for the mean.
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

    def stats(self):
        if self._w_sum == 0:
            v = numpy.empty_like(self._w_v_sum)
            v.fill(float('nan'))
            return v, v, v
        
        avg = self._w_v_sum / self._w_sum

        var_per_time = (1 / self._count) * (self._w_sum * numpy.multiply(avg, avg) - 2 * numpy.multiply(avg, self._w_v_sum) + self._w_vv_sum)

        # Sometimes, round-off error for identical values gives almost zero, but negative, error, so we clamp this away
        avg_var = numpy.maximum(var_per_time / self._w_sum, numpy.zeros_like(var_per_time))
        
        return avg, numpy.sqrt(avg_var), var_per_time

    def avg(self):
        '''Return the average of the weighted values.'''
        avg, _, _ = self.stats()
        return avg
    
    def std_dev(self):
        '''Return the standard deviation of the mean.'''
        avg, std_dev, _ = self.stats()
        return std_dev
    
    def ufloat(self):
        '''Return the average and error as a ufloat.'''
        avg, std_dev, _ = self.stats()
        return uncertainties.unumpy.uarray(avg, std_dev)

    def __str__(self):
        return str(self.ufloat())

    def __repr__(self):
        return repr(self.ufloat())

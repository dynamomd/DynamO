import math
from collections import defaultdict

import numpy
import uncertainties
import uncertainties.unumpy


class KeyedArray():
    """A key-value store of values that has element-wise addition and multiplication.
       Any missing values are assumed to be zero. This is needed for observations that do not appear in some simulations.
    """
    def __init__(self, values = None, type = float):
        self.type = type
        self.store = defaultdict(type)
        if values is not None:
            if isinstance(values, KeyedArray):
                for k, v in values.items():
                    self.store[k] = v
            elif isinstance(values, dict):
                for k, v in values.items():
                    self.store[k] = v
            else:
                raise RuntimeError("Cannot create KeyedArray from non-dict or non-KeyedArray")
    
    def __getitem__(self, key):
        """Get the value for a key, or zero if it does not exist."""
        return self.store[key]
    
    def __setitem__(self, key, value):
        """Set the value for a key."""
        self.store[key] = value
        return self.store[key]
    
    def __delitem__(self, key):
        """Delete the value for a key."""
        del self.store[key]
        return self.store[key]
    
    def __contains__(self, key):
        """Check if a key exists."""
        return key in self.store
    
    def __iter__(self):
        """Iterate over the keys."""
        return iter(self.store)
    
    def __len__(self):
        """Get the number of keys."""
        return len(self.store)
    
    def items(self):
        """Get the items as a list of tuples."""
        return self.store.items()
    
    def keys(self):
        """Get the keys as a list."""
        return self.store.keys()

    def values(self):
        """Get the values as a list."""
        return self.store.values()
    
    def __add__(self, rhs):
        if not isinstance(rhs, KeyedArray):
            raise RuntimeError("Cannot add non-KeyedArray to KeyedArray")
        
        retval = KeyedArray()
        retval.store = self.store.copy()

        for k, v in rhs.items():
            retval[k] = retval[k] + v

        return retval

    def __neg__(self):
        retval = KeyedArray()
        for k, v in self.items():
            retval[k] = -v
        return retval

    def __sub__(self, rhs):
        return self + (-rhs)

    def __mul__(self, rhs):
        import copy
        retval = KeyedArray()

        if isinstance(rhs, KeyedArray):
            # Only the intersection of the two KeyedArrays is multiplied, as everything else has at least one zero term.
            for k in set(self.keys()).intersection(rhs.keys()):
                retval[k] = copy.copy(self[k]) * rhs[k]
        else:
            for k, v in self.items():
                retval[k] = copy.copy(v) * rhs

        return retval
    
    def __rmul__(self, lhs):
        # Assume commutative multiplication
        return self.__mul__(lhs)
    

    def __truediv__(self, rhs):
        import copy
        retval = KeyedArray()

        if isinstance(rhs, KeyedArray):
            raise RuntimeError("Cannot divide KeyedArray by KeyedArray, what to do about missing values which are zeros?")
        else:
            for k, v in self.items():
                retval[k] = copy.copy(v) / rhs

        return retval
    
    def __repr__(self):
        return str(self.store)
    
    def __str__(self):
        return str(self.store)

def element_wise_multiply(a, b):
    if isinstance(a, numpy.ndarray) or isinstance(b, numpy.ndarray):
        return numpy.multiply(a, b)
    else:
        return a * b

def zeros_like(a):
    if isinstance(a, numpy.ndarray):
        return numpy.zeros_like(a)
    elif isinstance(a, KeyedArray):
        return a * 0.0
    else:
        return 0

def maximum(a, b):
    if isinstance(a, KeyedArray) and isinstance(b, KeyedArray):
        return KeyedArray({k: maximum(a[k], b[k]) for k in set(a.keys()).union(b.keys())})
    else:
        return numpy.maximum(a,b)
    
def sqrt(a):
    if isinstance(a, KeyedArray):
        return KeyedArray({k: sqrt(a[k]) for k in a.keys()})
    else:
        return numpy.sqrt(a)

class WeightedType():

    '''This class implements weighted arithmetic means along with an
    estimate of the standard error for the mean.
    '''
    def __init__(self, value = 0, weight = 0):
        """Initialise the weighted value."""
        ww = weight * weight
        vv = element_wise_multiply(value, value)
        self._ww_vv_sum = ww * vv
        self._ww_v_sum = ww * value
        self._ww_sum = ww
        self._w_vv_sum = weight * vv
        self._w_v_sum = weight * value
        self._w_sum = weight
        self._count = 1

    def __add__(self, v):
        if not isinstance(v, WeightedType):
            if not isinstance(v._w_v_sum, type(self._w_v_sum)):
                raise RuntimeError("Cannot add non-WeightedType to WeightedType")
            raise RuntimeError("Cannot add non-WeightedType to WeightedType")

        import copy
        retval = copy.copy(v)
        if self._w_sum == 0:
            #Shortcut if this has no weight, then just take the value passed in (and its shape!)
            return retval
    
        # We merge the two running sums of the weighted values 
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
            if isinstance(self._w_v_sum, numpy.ndarray):
                v = numpy.empty_like(self._w_v_sum)
                v.fill(float('nan'))
                return v, v, v
            else:
                return float('nan'), float('nan'), float('nan')

        avg = self._w_v_sum / self._w_sum

        var_per_time = (1 / self._count) * (self._w_sum * element_wise_multiply(avg, avg) - 2 * element_wise_multiply(avg, self._w_v_sum) + self._w_vv_sum)

        # Sometimes, round-off error for identical values gives almost zero, but negative, error, so we clamp this away
        avg_var = maximum(var_per_time / self._w_sum, zeros_like(var_per_time))
        
        return avg, sqrt(avg_var), var_per_time

    def avg(self):
        avg, _, _ = self.stats()
        return avg

    def std_dev(self):
        '''Return the standard deviation of the mean.'''
        _, std_dev, _ = self.stats()
        return std_dev
    
    def std_error(self):
        '''Return the standard error of the mean.'''
        _, std_dev, _ = self.stats()
        return std_dev / math.sqrt(self._w_sum)

    def ufloat(self):
        '''Return the average and error as a ufloat.'''
        avg, std_dev, _ = self.stats()
        if isinstance(avg, numpy.ndarray):
            # If the average is an array, we need to convert it to a ufloat array
            return uncertainties.unumpy.uarray(avg, std_dev)
        elif isinstance(avg, KeyedArray):
            return KeyedArray({k: uncertainties.ufloat(avg[k], std_dev[k]) for k in avg.keys()})
        else:
            return uncertainties.ufloat(avg, std_dev)

    def __str__(self):        
        return str(self.ufloat())

    def __repr__(self):
        return repr(self.ufloat())

from collections import defaultdict

import pytest

from pydynamo.weighted_types import *

def test_weighted_array_equals_weighted_float():
    # Check that array operations are equivalent to float operations
    a = WeightedArray(numpy.array([1, 2]), 1)
    b = WeightedArray(numpy.array([2, 4]), 0.1)
    c = a+b
    c1 = WeightedFloat(1, 1) + WeightedFloat(2, 0.1)
    c2 = WeightedFloat(2, 1) + WeightedFloat(4, 0.1)
    
    assert c.avg()[0] == pytest.approx(c1.avg())
    assert c.std_error()[0] == pytest.approx(c1.std_error())
    assert c.std_dev()[0] == pytest.approx(c1.std_dev())
    assert c.avg()[1] == pytest.approx(c2.avg())
    assert c.std_error()[1] == pytest.approx(c2.std_error())
    assert c.std_dev()[1] == pytest.approx(c2.std_dev())

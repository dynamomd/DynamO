import numpy
import pytest

from pydynamo.weighted_types import *


def test_weighted_estimators():
    # Create a process we know/understand
    true_avg = 5 # Center of the distribution
    true_stddev = 1 # Standard deviation of the distribution for a single sample
    time_step = 30 # How much time each sample represents
    N_pop = 10000 # Number of samples in the population
    values = numpy.random.normal(true_avg, true_stddev, N_pop)
    
    final_uncertainty = true_stddev / math.sqrt(N_pop)
    final_stddev_per_time = true_stddev * time_step
    limit = 5 # 5 sigma limit == 99.99994% sure its wrong

    # Compute estimates of the process properties
    sum = WeightedType(0, 0)
    for v in values:
        sum += WeightedType(v, time_step)
    avg, avg_std_dev, std_dev_time = sum.stats()

    #print("\n", true_avg, final_uncertainty, true_stddev)
    #print("\n!!!", avg, avg_std_dev, std_dev_time)

    # Check that the average is within 5 sigma of the true average
    assert abs((avg - true_avg) / (final_uncertainty)) < limit

    # Check that the estimated standard deviation is within 10% of the true standard deviation
    assert avg_std_dev == pytest.approx(final_uncertainty, rel=0.1)
    assert std_dev_time == pytest.approx(final_stddev_per_time, rel=0.1)

def test_weighted_array_equals_weighted_float():
    # Check that array operations are equivalent to float operations
    a = WeightedType(numpy.array([1, 2]), 1)
    b = WeightedType(numpy.array([2, 4]), 0.1)
    c = a+b
    c1 = WeightedType(1, 1) + WeightedType(2, 0.1)
    c2 = WeightedType(2, 1) + WeightedType(4, 0.1)

    assert c.avg()[0] == pytest.approx(c1.avg())
    assert c.avg()[1] == pytest.approx(c2.avg())
    for i in range(3):
        assert c.stats()[i][0] == pytest.approx(c1.stats()[i])
        assert c.stats()[i][1] == pytest.approx(c2.stats()[i])


def test_keyed_array():
    a = KeyedArray()
    a["a"] = 1
    a["b"] = 2
    a["c"] = 3
    b = KeyedArray()
    b["a"] = 2
    b["b"] = 4

    c = a + b
    assert c["a"] == pytest.approx(3)
    assert c["b"] == pytest.approx(6)
    assert c["c"] == pytest.approx(3)
    c = a - b
    assert c["a"] == pytest.approx(-1)
    assert c["b"] == pytest.approx(-2)
    assert c["c"] == pytest.approx(3)
    c = a * b
    assert c["a"] == pytest.approx(2)
    assert c["b"] == pytest.approx(8)
    assert "c" not in c
    c = a * 1.5
    assert c["a"] == pytest.approx(1.5)
    assert c["b"] == pytest.approx(3)
    assert c["c"] == pytest.approx(4.5)

    # Check division of two keyed arrays fails
    with pytest.raises(Exception):
       c = a / b
    
    c = a / 2
    assert c["a"] == pytest.approx(0.5)
    assert c["b"] == pytest.approx(1)
    assert c["c"] == pytest.approx(1.5)

    c = sqrt(a)
    assert c["a"] == pytest.approx(1)
    assert c["b"] == pytest.approx(math.sqrt(2))

    # Weighted types
    a = WeightedType(a, 1)
    b = WeightedType(b, 0.1)    
    c = a+b
    c1 = WeightedType(1, 1) + WeightedType(2, 0.1)
    c2 = WeightedType(2, 1) + WeightedType(4, 0.1)
    assert c.avg()["a"] == pytest.approx(c1.avg())
    assert c.avg()["b"] == pytest.approx(c2.avg())
    for i in range(3):
        assert c.stats()[i]["a"] == pytest.approx(c1.stats()[i])
        assert c.stats()[i]["b"] == pytest.approx(c2.stats()[i])

def test_keyed_keyed_array():
    a = KeyedArray(type=KeyedArray)
    a["a"]["1"] = 1
    a["a"]["2"] = 3
    a["b"]["2"] = 3
    a["c"]["3"] = 4

    # Check items that are there
    assert a["a"]["1"] == 1
    assert a["a"]["2"] == 3
    assert a["b"]["2"] == 3
    assert a["c"]["3"] == 4

    # Check items that are not there
    assert a["a"]["3"] == 0

    b = KeyedArray(type=KeyedArray)
    b["a"] = KeyedArray(float, {"1":2})

    c = a + b
    assert c["a"]["1"] == pytest.approx(3)
    assert c["b"]["2"] == pytest.approx(3)
    assert c["c"]["3"] == pytest.approx(4)

    c = a - b
    assert c["a"]["1"] == pytest.approx(-1)
    assert c["b"]["2"] == pytest.approx(3)
    assert c["c"]["3"] == pytest.approx(4)

    c = a * b
    assert c["a"]["1"] == pytest.approx(2)
    assert "b" not in c
    assert "c" not in c
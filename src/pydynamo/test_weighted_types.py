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
    sum = WeightedFloat(0, 0)
    for v in values:
        sum += WeightedFloat(v, time_step)
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
    a = WeightedArray(numpy.array([1, 2]), 1)
    b = WeightedArray(numpy.array([2, 4]), 0.1)
    c = a+b
    c1 = WeightedFloat(1, 1) + WeightedFloat(2, 0.1)
    c2 = WeightedFloat(2, 1) + WeightedFloat(4, 0.1)

    assert c.avg()[0] == pytest.approx(c1.avg())
    assert c.avg()[1] == pytest.approx(c2.avg())
    for i in range(3):
        assert c.stats()[i][0] == pytest.approx(c1.stats()[i])
        assert c.stats()[i][1] == pytest.approx(c2.stats()[i])

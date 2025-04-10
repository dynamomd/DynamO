import math

import pytest

import pydynamo


def test_pressure():
    """Here we test basic functionality of the SingleAttribute class and the Pressure property."""
    from subprocess import PIPE, run

    # Create a hard sphere config
    createlog = run(["dynamod", "-m", "0"], check=True, stdout=PIPE, stderr=PIPE)
    runlog = run(["dynarun", "config.out.xml.bz2", "-c", "10000"], check=True, stdout=PIPE, stderr=PIPE)

    outputfile = pydynamo.OutputFile("output.xml.bz2")
    outputplugin = pydynamo.OutputFile.output_props["p"]
    prop = outputplugin.init()
    result = outputplugin.result(state=[], outputfile=outputfile, configfilename="config.out.xml.bz2", counter=1, manager=None, output_dir=None)

    assert result.avg() == pytest.approx(outputfile.p())
    #assert math.isnan(result.std_error())

    runlog = run(["dynarun", "config.out.xml.bz2", "-c", "10000"], check=True, stdout=PIPE, stderr=PIPE)
    outputfile2 = pydynamo.OutputFile("output.xml.bz2")
    result += outputplugin.result(state=[], outputfile=outputfile2, configfilename="config.out.xml.bz2", counter=2, manager=None, output_dir=None)

    time_weighted_avg = outputfile.p() * outputfile.t() + outputfile2.p() * outputfile2.t()
    time_weighted_avg /= outputfile.t() + outputfile2.t()
    
    assert result.avg() == pytest.approx(time_weighted_avg)
    #print("\n",repr(result))
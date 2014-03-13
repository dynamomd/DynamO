#!/usr/bin/env python
# -*- coding: utf-8 -*-
####################################################
# Sample script that shows how to save result data #
####################################################
import urllib
import urllib2
from datetime import datetime
import commands
import time
import os
import sys
import math

# You need to enter the real URL and have the server running
CODESPEED_URL = 'http://localhost:8000/'
Ntests=3
testPath="build-dir"
commitID = commands.getoutput('git rev-parse HEAD')
current_date = datetime.today()
# Mandatory fields

def avg(vallist):
    return sum(vallist) / len(vallist)

def stddev(vallist):
    avgval=avg(vallist)
    return math.sqrt(sum(val * val - avgval * avgval for val in vallist) / len(vallist))

def RunAndAddData(benchmarkname, executable):
    print "Running",benchmarkname,
    sys.stdout.flush()
    timings=[]
    for i in range(3):
        print i,
        sys.stdout.flush()
        start=float(time.time())
        if os.system(executable+" > run.log 2>&1") != 0:
            continue
        timings.append(float(time.time()) - start)
    if not timings:
        return
    print "Done!"
    data = {
        'commitid': commitID,
        'branch': 'master',
        'project': 'DynamO',
        'executable': 'default',
        'benchmark': benchmarkname,
        'environment': "i7 Laptop",
        'result_value': avg(timings),
        'min': min(timings),
        'max': max(timings),
        'std_dev': stddev(timings),
        'result_date': current_date,
        }

    params = urllib.urlencode(data)
    response = "None"
    print "Saving result for executable %s, revision %s, benchmark %s" % (
        data['executable'], data['commitid'], data['benchmark'])
    try:
        f = urllib2.urlopen(CODESPEED_URL + 'result/add/', params)
    except urllib2.HTTPError as e:
        print str(e)
        print e.read()
        return
    response = f.read()
    f.close()
    print "Server (%s) response: %s\n" % (CODESPEED_URL, response)

def findJamFiles(path):
    '''Searches the provided path for jamfile files and returns them
    in a dictionary, with the folder as a key.'''
    jamfiles={}
    for root, dirs, files in os.walk(path):
        if "jamfile" in files:
            jamfiles[root.split(os.path.sep)[-1]] = root
    return jamfiles

def findUnitTests(jamfile):
    '''Given a path to a jamfile, parse the jamfile for all unit-test
    tags, their name, and their sources.'''
    f=open(jamfile, 'r')
    line=""
    unittests=[]
    while True:
        newline = f.readline()
        if not newline:
            return unittests
        if newline[0] == '#':
            continue
        line = line.strip()+" "+newline.strip()
        if ";" not in line:
            continue
        keyword, args = line.split(None, 1)
        if keyword == "unit-test":
            name, sources = args.split(";")[0].split(":")[0:2]
            unittests.append({"name":name, "sources":sources.split()})
        line=""
    return unittests

def touch(fname):
    with file(fname, 'a'):
        os.utime(fname, None)

class TestBuildFail(Exception):
    pass

def builtTestAndExtractPath(projectpath, unittest):
    print "Building",unittest['name']
    import subprocess
    touch(os.path.join(projectpath, unittest["sources"][0]))
    command=['bjam','-j8',projectpath+'//'+unittest["name"],'testing.launcher="echo TESTEXELOCATION"']
    try:
        output=subprocess.check_output(' '.join(command), shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError, e:
        raise TestBuildFail("Build of test failed")
        
    for line in output.split('\n'):
        if line[:15] == "TESTEXELOCATION":
            return line.split()[-1]
    raise Exception("Failed to locate exe")

if __name__ == "__main__":
    for project, projectpath in findJamFiles(".").iteritems():
        print "##### PROJECT ",project, projectpath
        for unittest in findUnitTests(os.path.join(projectpath, "jamfile")):
            try:
                RunAndAddData(project+':'+unittest['name'], builtTestAndExtractPath(projectpath, unittest))
            except TestBuildFail, e:
                print " Build failed!"

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
import math

# You need to enter the real URL and have the server running
CODESPEED_URL = 'http://localhost:8000/'
Ntests=3
testPath="build-dir"
commitID = commands.getoutput('git rev-parse HEAD')
current_date = datetime.today()
# Mandatory fields

def add(executable, benchmark, timing, minVal=None, maxVal=None, dev=None):
    data = {
        'commitid': commitID,
        'branch': 'default',
        'project': 'DynamO',
        'executable': executable,
        'benchmark': benchmark,
        'environment': "i7 Laptop",
        'result_value': timing,
        'result_date': current_date
        }

    if minVal is not None:
        data.update({'min': minVal})
    if maxVal is not None:
        data.update({'max': minVal})
    if dev is not None:
        data.update({'std_dev': dev})
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

if __name__ == "__main__":
    for line in open('src/dynamo/jamfile'):
        data = line.split()
        if len(data) == 0:
            continue
        if data[0] == 'unit-test':
            testname=data[1]
            for root, dirs, files in os.walk("./build-dir"):
                if testname+'.passed' in files:
                    print "Running unit-test "+testname
                    timings=[]
                    for i in range(Ntests):
                        start=float(time.time())
                        os.system(os.path.join(root, testname)+" > /dev/null 2>&1")
                        timings.append(float(time.time()) - start)
                    avg=sum(timings)/len(timings)
                    minval=min(timings)
                    maxval=max(timings)
                    dev=math.sqrt(sum((t - avg)*(t - avg) for t in timings) / len(timings))
                    print "Avg =",avg,"Min =",minval,"Max =",maxval,"Dev =",dev
                    add('dynamo', testname, sum(timings)/len(timings), minVal=min(timings), maxVal=max(timings), dev=dev)

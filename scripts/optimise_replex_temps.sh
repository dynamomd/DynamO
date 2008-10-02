#!/bin/bash
../../replex_opt | gawk '{print "bzcat replex."$1".xml.bz2 | xml ed -u '"'/DYNAMOconfig/Dynamics/SystemEvents/System/@Temperature'"' -v \""$2"\" | bzip2 > replex."$1".xml.bz2"}' | /bin/bash

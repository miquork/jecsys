#!/usr/bin/env python

import sys
import json

f = open(sys.argv[1])
map = json.load(f)
f.close()

for run, lumis in map.iteritems():
   for lumi in lumis:
       print run, " ".join([str(num) for num in lumi])

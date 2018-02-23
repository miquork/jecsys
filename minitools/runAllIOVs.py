#! /usr/bin/python
import os

IOV_list= ['B','C','D','E','F','BCDEF']
for iov in IOV_list:
    os.system("mkdir pdf/"+iov)
    os.system("root -b -q 'mk_reprocess_RunEpoch.C(\""+iov+"\")'")

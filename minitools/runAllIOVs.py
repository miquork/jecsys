#! /usr/bin/python
import os

IOV_list= ['A','B','C','ABC']
for iov in IOV_list:
    os.system("mkdir pdf/"+iov)
    os.system("root -b -q 'mk_reprocess_RunEpoch.C(\""+iov+"\")'")

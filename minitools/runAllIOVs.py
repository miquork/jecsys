#! /usr/bin/python
import os

# Run BCDEF first in case using global fit 'dofsrcombo' option
#IOV_list= ['BCDEF','B','C','D','E','F']
IOV_list= ['2018ABCD','2018A','2018B','2018C','2018D']
for iov in IOV_list:
    os.system("mkdir pdf/"+iov)
    os.system("root -b -q 'mk_reprocess_RunEpoch.C(\""+iov+"\")'")

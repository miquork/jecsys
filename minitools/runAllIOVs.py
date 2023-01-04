#! /usr/bin/python
import os

# Run BCDEF first in case using global fit 'dofsrcombo' option
#IOV_list= ['BCDEF','B','C','D','E','F']
#IOV_list= ['2018ABCD','2018A','2018B','2018C','2018D']
#IOV_list= ['2016GH','2016BCDEF','2016BCD','2016EF']
#IOV_list= ['2016BCDEF','2016GH','2016BCD','2016EF']
IOV_list= ['2016BCDEF','2016GH','2017BCDEF','2017H','2018ABCD','Run2Test']
for iov in IOV_list:
    os.system("mkdir pdf/"+iov)
    os.system("root -b -q 'mk_reprocess_RunEpoch.C(\""+iov+"\")'")

os.system("root -b -q 'recombine.C'")
os.system("root -b -q 'mk_reprocess_RunEpoch.C(\"Run2Test\")'")

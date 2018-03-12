#! /usr/bin/python
import os

IOV_list= ['BCDEFGH','BCD','EF','GH','G','H']
for iov in IOV_list:
    os.system("mkdir pdf/"+iov)
    os.system("root -b -q 'mk_reprocess_RunEpoch.C(\""+iov+"\")'")

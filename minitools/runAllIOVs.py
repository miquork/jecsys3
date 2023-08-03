#! /usr/bin/python
import os

IOV_list= ['Run22C','Run22D','Run22CD','Run22E','Run22F','Run22G','Run22FG',
           'Run23B','Run23C1','Run23C2','Run23C2','Run23C3','Run23BC123',
           'Run23C4']
for iov in IOV_list:
    os.system("mkdir pdf/"+iov)
    os.system("root -b -q 'mk_reprocess_RunEpoch.C(\""+iov+"\")'")

#! /usr/bin/python
import os

#IOV_list= ['2016ULAPV','2016UL','2018UL','Run2']
#IOV_list= ['Run22C','Run22D','Run22CD','Run22E','Run22F','Run22G','Run22FG',
#           'Run23B','Run23C1','Run23C2','Run23C2','Run23C3','Run23BC123',
#           'Run23C4']
#IOV_list= ['Run22CD','Run22E','Run22FG','Run23BC123','Run23C4D','Run3']
#IOV_list= ['Run22CD','Run22E','Run22FG','Run23C123',
#           'Run23C4','Run23D']
#IOV_list= ['Run23C123','Run23C4','Run23D']
#          'Run23C4D']
#IOV_list= ['Run24BCD','Run24E','Run24BCDE','Run24CP','Run24CR','Run24CS']
#IOV_list= ['Run24F','Run24E','Run24BCD','Run24E']
IOV_list= ['Run24G','Run24F','Run24E','Run24BCD']
for iov in IOV_list:
#os.system("mkdir pdf/"+iov)
    os.system("root -b -q 'mk_reprocess_RunEpoch.C(\""+iov+"\")'")
#os.system("root -b -q 'recombine.C'")
#os.system("root -b -q 'mk_reprocess_RunEpoch.C(\"Run3\")'")

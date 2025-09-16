#! /usr/bin/python3
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
#IOV_list= ['Run24G','Run24F','Run24E','Run24BCD']
#IOV_list= ['Run24I','Run24H','Run24G','Run24F','Run24E','Run24BCD']
#IOV_list= ['2024B_nib1','2024C_nib1','2024D_nib1','2024Ev1_nib1','2024Ev2_nib1',
#           '2024F_nib1','2024F_nib2','2024F_nib3','2024G_nib1','2024G_nib2',
#           '2024H_nib1','2024I_nib1']
#`IOV_list= ['2024C_nib1','2024D_nib1','2024Ev1_nib1','2024Ev2_nib1',
#IOV_list= ['2024C_nib1','2024D_nib1','2024E_nib1',
#           '2024F_nib1','2024F_nib2','2024F_nib3','2024G_nib1','2024G_nib2',
#           '2024H_nib1','2024I_nib1',
#           '2024_nib']
#           '2025C']
IOV_list= ['2025C','2025D','2025E','2025CDE']
#IOV_list= ['2025C']

for iov in IOV_list:
#os.system("mkdir pdf/"+iov)
    os.system("root -b -q 'mk_reprocess_RunEpoch.C(\""+iov+"\")'")
#os.system("root -b -q 'recombine.C'")
#os.system("root -b -q 'mk_reprocess_RunEpoch.C(\"Run3\")'")

#! /usr/bin/python3
import os

IOV_list= ['2024_nib',
           '2024C_nib1','2024D_nib1','2024E_nib1','2024CDE_nib',
           '2024F_nib1','2024F_nib2','2024F_nib3','2024G_nib1','2024G_nib2',
           '2024H_nib1','2024I_nib1','2024FGHI_nib','2024_nib',
           '2025C','2025D','2025E','2025F','2025G','2025CDEFG',
           '2026B','2026C','2026D']

# remove old L2Res files, but request permission (-i) for it first
os.system("rm -i rootfiles/JERSF.root");
os.system("rm -i rootfiles/JERSF_closure.root");

# startup run to force recompile JERSF to include latest Config.C changes
os.system("root -b -q 'JERSF.C++g(false,\""+IOV_list[0]+"\",\"Summer24MG_NOJERSF\")'")

# derive JER SF in loop without recompiling to go faster
for iov in IOV_list:
    os.system("root -l -b -q 'JERSF.C+g(false,\""+iov+"\",\"Summer24MG_NOJERSF\")'")

# test closure
#for iov in IOV_list:
#    os.system("root -l -b -q 'JERSF.C+g(true,\""+iov+"\",\"Summer24MG_NOJERSF\")'")

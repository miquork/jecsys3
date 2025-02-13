// Purpose: common configurations for L2Res.C, L3Res.C/globalFit, JERSF.C
#ifndef __CONFIG_C__
#define __CONFIG_C__

#include <map>
#include <string>
std::map<std::string, std::string> mlum;

mlum["2024B"] = "0.13 fb^{-1}";
mlum["2024C"] = "7.5 fb^{-1}";
mlum["2024D"] = "8.3 fb^{-1}";
mlum["2024BC"] = "7.7 fb^{-1}";
mlum["2024BCD"] = "15.3 fb^{-1}";
mlum["2024CP"] = "7.5 fb^{-1}";
mlum["2024CR"] = "7.5 fb^{-1}";
mlum["2024CS"] = "7.5 fb^{-1}";
mlum["2024E"] = "11.3 fb^{-1}";
mlum["2024F"] = "27.8 fb^{-1}";
mlum["2024FG"] = "65.5 fb^{-1}";
mlum["2024G"] = "37.8 fb^{-1}";
mlum["2024H"] = "5.4 fb^{-1}";
mlum["2024I"] = "11.6 fb^{-1}";

mlum["2024B_nib1"] = "0.13 fb^{-1}";
mlum["2024C_nib1"] = "7.5 fb^{-1}";
mlum["2024D_nib1"] = "8.3 fb^{-1}";
//mlum["2024E"] = "11.3 fb^{-1}";
mlum["2024Ev1_nib1"] = "X.X fb^{-1}";
mlum["2024Ev2_nib1"] = "X.X fb^{-1}";
//mlum["2024F"] = "27.8 fb^{-1}";
mlum["2024F_nib1"] = "0.88 fb^{-1}";
mlum["2024F_nib2"] = "X.X fb^{-1}";
mlum["2024F_nib3"] = "X.X fb^{-1}";
//mlum["2024G"] = "37.8 fb^{-1}";
mlum["2024G_nib1"] = "X.X fb^{-1}";
mlum["2024G_nib2"] = "X.X fb^{-1}";
mlum["2024H_nib1"] = "5.4 fb^{-1}";
mlum["2024I_nib1"] = "11.6 fb^{-1}";

// rootfiles/brilcalc/sumfibs.py

// NIB-level:
mlum["2022B_nib1"] = "0.097 fb^{-1}";
mlum["2022C_nib1"] = "5.0 fb^{-1}";
mlum["2022D_nib1"] = "2.97 fb^{-1}";
mlum["2022E_nib1"] = "5.8 fb^{-1}";
mlum["2022F_nib1"] = "17.8 fb^{-1}";
mlum["2022G_nib1"] = "3.08 fb^{-1}";
mlum["2023B_nib1"] = "0.64 fb^{-1}";
mlum["2023C_nib1"] = "7.2 fb^{-1}";
mlum["2023Cv4_nib1"] = "10.4 fb^{-1}";
mlum["2023Cv4_nib2"] = "0.407 fb^{-1}";
mlum["2023D_nib1"] = "9.7 fb^{-1}";
mlum["2024B_nib1"] = "0.130 fb^{-1}";
mlum["2024C_nib1"] = "7.2 fb^{-1}";
mlum["2024D_nib1"] = "8.0 fb^{-1}";
mlum["2024Ev1_nib1"] = "6.3 fb^{-1}";
mlum["2024Ev2_nib1"] = "5.0 fb^{-1}";
mlum["2024F_nib1"] = "0.88 fb^{-1}";
mlum["2024F_nib2"] = "14.4 fb^{-1}";
mlum["2024F_nib3"] = "12.5 fb^{-1}";
mlum["2024G_nib1"] = "16.9 fb^{-1}";
mlum["2024G_nib2"] = "20.9 fb^{-1}";
mlum["2024H_nib1"] = "5.4 fb^{-1}";
mlum["2024I_nib1"] = "11.5 fb^{-1}";

// ERA-level:
mlum["2022B"] = "0.097 fb^{-1}";
mlum["2022C"] = "5.0 fb^{-1}";
mlum["2022D"] = "2.97 fb^{-1}";
mlum["2022E"] = "5.8 fb^{-1}";
mlum["2022F"] = "17.8 fb^{-1}";
mlum["2022G"] = "3.08 fb^{-1}";
mlum["2023B"] = "0.64 fb^{-1}";
mlum["2023C"] = "18.1 fb^{-1}";
mlum["2023D"] = "9.7 fb^{-1}";
mlum["2024B"] = "0.130 fb^{-1}";
mlum["2024C"] = "7.2 fb^{-1}";
mlum["2024D"] = "8.0 fb^{-1}";
mlum["2024E"] = "11.3 fb^{-1}";
mlum["2024F"] = "27.8 fb^{-1}";
mlum["2024G"] = "37.8 fb^{-1}";
mlum["2024H"] = "5.4 fb^{-1}";
mlum["2024I"] = "11.5 fb^{-1}";

// YEAR-level:
mlum["2022"] = "34.8 fb^{-1}";
mlum["2023"] = "28.4 fb^{-1}";
mlum["2024"] = "109 fb^{-1}";

#endif

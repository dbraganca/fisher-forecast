(* ::Package:: *)

SetDirectory[NotebookDirectory[]]


(*Calculate interpolators for power spectrum*)
jmatpath="../GitHub/python-integer-power-project/2. Jmat_loopvals/";
LoopTablesPath="tabs/";

p13jfunc[k_]:=Import[jmatpath<>"P13_bias_Jmat/P13_Jfunc_"<>ToString[Round[k,0.00001]]<>"_.h5","jmat"];
p22jfunc[k_]:=Import[jmatpath<>"P22_bias_Jmat/P22_Jfunc_"<>ToString[Round[k,0.00001]]<>"_.h5","jmat"];

(*ktabTot includes both CMASS and LOWZ*)
ktabTot={
0.0158`,0.02139`,0.0258`,0.03139`,0.0358`,0.04139`,0.0458`,0.05139`,
0.0558`,0.06139`,0.0658`,0.07139`,0.0758`,0.08139`,0.0858`,0.09139`,
0.0958`,0.10139`,0.1058`,0.11139`,0.1158`,0.12139`,0.1258`,0.13139`,
0.1358`,0.14139`,0.1458`,0.15139`,0.1558`,0.16139`,0.1658`,0.17139`,
0.1758`,0.18139`,0.1858`,0.19139`,0.1958`,0.2058`,0.2158`};

p13jfuncTab=Transpose[{ktabTot,p13jfunc/@ktabTot}];
P13interp=Interpolation[p13jfuncTab];

p22jfuncTab=Transpose[{ktabTot,p22jfunc/@ktabTot}];
P22interp=Interpolation[p22jfuncTab];


Export[LoopTablesPath<>"P22interp.m",P22interp]
Export[LoopTablesPath<>"P13interp.m",P13interp]


(*make ktab*)
ktabCMASS = Import["C:\\Users\\diogo\\Dropbox\\FFTLog\\Diogo_work\\GitHub\\python-integer-power-project\\3. Ctabs\\CMASS_ks.csv","CSV"]//Flatten;
ktabLOWZ = Import["C:\\Users\\diogo\\Dropbox\\FFTLog\\Diogo_work\\GitHub\\python-integer-power-project\\3. Ctabs\\LOWZ_ks.csv","CSV"]//Flatten;


ktabLOWZ


Export["covsks/ktabCMASS.m",ktabCMASS]
Export["covsks/ktabLOWZ.m",ktabLOWZ]




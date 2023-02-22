(* ::Package:: *)

SetDirectory[NotebookDirectory[]]


(* ::Subsection:: *)
(*Make P13 and P22 interpolators*)


(*Calculate interpolators for power spectrum*)
jmatpath="../GitHub/python-integer-power-project/2. Jmat_loopvals/";
LoopTablesPath="tabs/";
(*make ktab*)
ktabCMASS = Import["C:\\Users\\diogo\\Dropbox\\FFTLog\\Diogo_work\\GitHub\\python-integer-power-project\\3. Ctabs\\CMASS_ks.csv","CSV"]//Flatten;
ktabLOWZ = Import["C:\\Users\\diogo\\Dropbox\\FFTLog\\Diogo_work\\GitHub\\python-integer-power-project\\3. Ctabs\\LOWZ_ks.csv","CSV"]//Flatten;
Export["covsks/ktabCMASS.m",ktabCMASS]
Export["covsks/ktabLOWZ.m",ktabLOWZ]

p13jfunc[k_]:=Import[jmatpath<>"P13_bias_Jmat/P13_Jfunc_"<>ToString[Round[k,0.00001]]<>"_.h5","jmat"];
p22jfunc[k_]:=Import[jmatpath<>"P22_bias_Jmat/P22_Jfunc_"<>ToString[Round[k,0.00001]]<>"_.h5","jmat"];

(*ktabTot includes both CMASS and LOWZ*)
ktabTot=Sort[Join[ktabLOWZ,ktabCMASS]];


ktabTot={0.0158`,0.02139`,0.0258`,0.03139`,0.0358`,0.04139`,0.0458`,0.05139`,0.0558`,
		 0.06139`,0.0658`,0.07139`,0.0758`,0.08139`,0.0858`,0.09139`,0.0958`,0.10139`,
		 0.1058`,0.11139`,0.1158`,0.12139`,0.1258`,0.13139`,0.1358`,0.14139`,0.1458`,
		 0.15139`,0.1558`,0.16139`,0.1658`,0.17139`,0.1758`,0.18139`,0.1858`,0.19139`,
		 0.1958`,0.20139`,0.2058`,0.21139`,0.2158`,0.22139`};


p13jfuncTab=Transpose[{ktabTot,p13jfunc/@ktabTot}];
P13interp=Interpolation[p13jfuncTab];

p22jfuncTab=Transpose[{ktabTot,p22jfunc/@ktabTot}];
P22interp=Interpolation[p22jfuncTab];


Export[LoopTablesPath<>"P22interp.m",P22interp]
Export[LoopTablesPath<>"P13interp.m",P13interp]


(* ::Subsubsection:: *)
(*k list for MegaMapper*)


ktabraw = Table[0.001`+i 0.005,{i,0,300}];


(* ::Subsection:: *)
(*Make mu list from bispectrum diagrams coefficients (this is essential to speed up angular integration)*)


coefpath="tabs/";
BLoop\[Mu]List={1,u^2,u^4,u^6,u^8,u^10,u^12,u \[Mu]1,u^3 \[Mu]1,u^5 \[Mu]1,u^7 \[Mu]1,u^9 \[Mu]1,u^11 \[Mu]1,\[Mu]1^2,u^2 \[Mu]1^2,u^4 \[Mu]1^2,u^6 \[Mu]1^2,u^8 \[Mu]1^2,u^10 \[Mu]1^2,u \[Mu]1^3,u^3 \[Mu]1^3,u^5 \[Mu]1^3,u^7 \[Mu]1^3,u^9 \[Mu]1^3,\[Mu]1^4,u^2 \[Mu]1^4,u^4 \[Mu]1^4,u^6 \[Mu]1^4,u^8 \[Mu]1^4,u \[Mu]1^5,u^3 \[Mu]1^5,u^5 \[Mu]1^5,u^7 \[Mu]1^5,\[Mu]1^6,u^2 \[Mu]1^6,u^4 \[Mu]1^6,u^6 \[Mu]1^6,u \[Mu]1^7,u^3 \[Mu]1^7,u^5 \[Mu]1^7,\[Mu]1^8,u^2 \[Mu]1^8,u^4 \[Mu]1^8,u \[Mu]1^9,u^3 \[Mu]1^9,\[Mu]1^10,u^2 \[Mu]1^10,u \[Mu]1^11,\[Mu]1^12};
(*check it matches Yaniv*)
BLoop\[Mu]List=={1,u^2,u^4,u^6,u^8,u^10,u^12,u \[Mu]1,u^3 \[Mu]1,u^5 \[Mu]1,u^7 \[Mu]1,u^9 \[Mu]1,u^11 \[Mu]1,\[Mu]1^2,u^2 \[Mu]1^2,u^4 \[Mu]1^2,u^6 \[Mu]1^2,u^8 \[Mu]1^2,u^10 \[Mu]1^2,u \[Mu]1^3,u^3 \[Mu]1^3,u^5 \[Mu]1^3,u^7 \[Mu]1^3,u^9 \[Mu]1^3,\[Mu]1^4,u^2 \[Mu]1^4,u^4 \[Mu]1^4,u^6 \[Mu]1^4,u^8 \[Mu]1^4,u \[Mu]1^5,u^3 \[Mu]1^5,u^5 \[Mu]1^5,u^7 \[Mu]1^5,\[Mu]1^6,u^2 \[Mu]1^6,u^4 \[Mu]1^6,u^6 \[Mu]1^6,u \[Mu]1^7,u^3 \[Mu]1^7,u^5 \[Mu]1^7,\[Mu]1^8,u^2 \[Mu]1^8,u^4 \[Mu]1^8,u \[Mu]1^9,u^3 \[Mu]1^9,\[Mu]1^10,u^2 \[Mu]1^10,u \[Mu]1^11,\[Mu]1^12}

(*creates Function that obtains the coefficient of each term in biaslist in each ctab coefficient "coef" *)
biaslisterfn[vars_,biaslist_]:=Module[{tab0,biasexps,biaslister},
biasexps=Flatten[CoefficientRules[biaslist,vars][[All,All,1]],1];
tab0=Table[_,Length[vars]];
biaslister[coef_]:=(biasexps/.CoefficientRules[coef,vars])/.{tab0:>0};
biaslister
];


BLoop\[Mu]lister=biaslisterfn[{\[Mu]1,u},BLoop\[Mu]List];


(*Cosine of angle between k1 and k2*)
ys[k1_,k2_,k3_]=(-k1^2-k2^2+k3^2)/(2 k1 k2);
(*Cosine of angle between k1 and k3*)
y13s[k1_,k2_,k3_]=(-k1^2-k3^2+k2^2)/(2 k1 k3);
(*Use u variable instead of Sqrt[1-\[Mu]1^2]Sin[\[Phi]] to speed up angular integration later*)
\[Mu]2ey= (y \[Mu]1+Sqrt[1-y^2] u);
\[Mu]3ey=-Sqrt[1-y13^2]u+y13 \[Mu]1;


(* ::Subsubsection:: *)
(*Tree level*)


Btreeimp=Import[coefpath<>"BtreePerms.m"];


BTree=Btreeimp/.{\[Mu]2->\[Mu]2ey,\[Mu]3->\[Mu]3ey}//Expand;


BTreeCoefs=BLoop\[Mu]lister[BTree];


(*Check decomposition works*)
BTreeCoefs.BLoop\[Mu]List-BTree//Expand


Export[coefpath<>"Btreecoefslistmu.m",BTreeCoefs]


(* ::Subsubsection:: *)
(*Diagrams*)


b222coeflist=Import[coefpath<>"B222simpcoefslist.m"];
b222coeflistmu=BLoop\[Mu]lister[b222coeflist/.{\[Mu]2->\[Mu]2ey,\[Mu]3->\[Mu]3ey}];
(*Check decomposition works*)
Total[Abs[b222coeflistmu.BLoop\[Mu]List-b222coeflist/.{\[Mu]2->\[Mu]2ey,\[Mu]3->\[Mu]3ey}//Expand]]


Export[coefpath<>"b222coeflistmu.m",b222coeflistmu]


(* ::Subsubsection:: *)
(**)


b3211coeflist=Import[coefpath<>"B3211simpcoefslist.m"];
b3211coeflistmu=BLoop\[Mu]lister[b3211coeflist/.{\[Mu]2->\[Mu]2ey,\[Mu]3->\[Mu]3ey}];
(*Check decomposition works*)
Total[Abs[b3211coeflistmu.BLoop\[Mu]List-b3211coeflist/.{\[Mu]2->\[Mu]2ey,\[Mu]3->\[Mu]3ey}//Expand]]


Export[coefpath<>"b3211coeflistmu.m",b3211coeflistmu]


b3212coeflist=Import[coefpath<>"B3212simpcoefslist.m"];
b3212coeflistmu=BLoop\[Mu]lister[b3212coeflist/.{\[Mu]2->\[Mu]2ey,\[Mu]3->\[Mu]3ey}];
(*Check decomposition works*)
Total[Abs[b3212coeflistmu.BLoop\[Mu]List-b3212coeflist/.{\[Mu]2->\[Mu]2ey,\[Mu]3->\[Mu]3ey}//Expand]]


Export[coefpath<>"b3212coeflistmu.m",b3212coeflistmu]


b411coeflist=Import[coefpath<>"B411simpcoefslist.m"];
b411coeflistmu=BLoop\[Mu]lister[b411coeflist/.{\[Mu]2->\[Mu]2ey,\[Mu]3->\[Mu]3ey}];


(*Check decomposition works*)
Total[Abs[b411coeflistmu.BLoop\[Mu]List-b411coeflist/.{\[Mu]2->\[Mu]2ey,\[Mu]3->\[Mu]3ey}//Expand]]


b411coeflistmu//Variables


Export[coefpath<>"b411coeflistmu.m",b411coeflistmu]

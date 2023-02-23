(* ::Package:: *)

SetDirectory[NotebookDirectory[]]


(* ::Subsection:: *)
(*Make P13 and P22 interpolators*)


(*Calculate interpolators for power spectrum*)
jmatpath="../GitHub/python-integer-power-project/2. Jmat_loopvals/";
LoopTablesPath="tabs/";


(*make ktab*)
ktabCMASS = Import["covsks/ktabCMASS.m"];
ktabLOWZ = Import["covsks/ktabLOWZ.m"];
ktabraw = Import["covsks/ktabBase.m"];




p13jfunc[k_]:=Import[jmatpath<>"P13/"<>ToString[Round[k,0.00001]]<>"_.h5","jmat"];
p22jfunc[k_]:=Import[jmatpath<>"P22/"<>ToString[Round[k,0.00001]]<>"_.h5","jmat"];

(*ktabTot includes both CMASS and LOWZ, as well as base*)
ktabTot=Sort[Join[ktabLOWZ,ktabCMASS,ktabraw]];


ktabTot={0.001`,0.006`,0.011`,0.0158`,0.016`,0.021`,0.02139`,0.0258`,0.026000000000000002`,0.031`,
0.03139`,0.0358`,0.036000000000000004`,0.041`,0.04139`,0.0458`,0.046`,0.051000000000000004`,0.05139`,
0.0558`,0.056`,0.061`,0.06139`,0.0658`,0.066`,0.07100000000000001`,0.07139`,0.0758`,0.076`,0.081`,0.08139`,
0.0858`,0.08600000000000001`,0.091`,0.09139`,0.0958`,0.096`,0.101`,0.10139`,0.1058`,0.106`,0.111`,0.11139`,
0.1158`,0.116`,0.121`,0.12139`,0.1258`,0.126`,0.131`,0.13139`,0.1358`,0.136`,0.14100000000000001`,0.14139`,
0.1458`,0.146`,0.151`,0.15139`,0.1558`,0.156`,0.161`,0.16139`,0.1658`,0.166`,0.171`,0.17139`,0.1758`,0.17600000000000002`,
0.181`,0.18139`,0.1858`,0.186`,0.191`,0.19139`,0.1958`,0.196`,0.201`,0.20139`,0.2058`,0.20600000000000002`,0.211`,
0.21139`,0.2158`,0.216`,0.221`,0.22139`,0.226`,0.231`,0.23600000000000002`,0.241`,0.246`,0.251`,0.256`,0.261`,0.266`,
0.271`,0.276`,0.281`,0.28600000000000003`,0.291`,0.296`,0.301`,0.306`,0.311`,0.316`,0.321`,0.326`,0.331`,0.336`,0.341`,
0.34600000000000003`,0.35100000000000003`,0.356`,0.361`,0.366`,0.371`,0.376`,0.381`,0.386`,0.391`,0.396`,0.401`,0.406`,
0.41100000000000003`,0.41600000000000004`,0.421`,0.426`,0.431`,0.436`,0.441`,0.446`,0.451`,0.456`,0.461`,0.466`,
0.47100000000000003`,0.47600000000000003`,0.481`,0.486`,0.491`,0.496`,0.501`,0.506`,0.511`,0.516`,0.521`,0.526`,
0.531`,0.536`,0.541`,0.546`,0.551`,0.556`,0.561`,0.5660000000000001`,0.5710000000000001`,0.5760000000000001`,0.581`,
0.586`,0.591`,0.596`,0.601`,0.606`,0.611`,0.616`,0.621`,0.626`,0.631`,0.636`,0.641`,0.646`,0.651`,0.656`,0.661`,
0.666`,0.671`,0.676`,0.681`,0.686`,0.6910000000000001`,0.6960000000000001`,0.7010000000000001`,0.706`,0.711`,0.716`,
0.721`,0.726`,0.731`,0.736`,0.741`,0.746`,0.751`,0.756`,0.761`,0.766`,0.771`,0.776`,0.781`,0.786`,0.791`,0.796`,0.801`,
0.806`,0.811`,0.8160000000000001`,0.8210000000000001`,0.8260000000000001`,0.8310000000000001`,0.836`,0.841`,0.846`,
0.851`,0.856`,0.861`,0.866`,0.871`,0.876`,0.881`,0.886`,0.891`,0.896`,0.901`,0.906`,0.911`,0.916`,0.921`,0.926`,0.931`,
0.936`,0.9410000000000001`,0.9460000000000001`,0.9510000000000001`,0.9560000000000001`,0.961`,0.966`,0.971`,0.976`,0.981`,
0.986`,0.991`,0.996`,1.001`,1.006`,1.011`,1.016`,1.021`,1.0259999999999998`,1.031`,1.0359999999999998`,1.041`,
1.0459999999999998`,1.051`,1.0559999999999998`,1.061`,1.0659999999999998`,1.071`,1.0759999999999998`,1.081`,
1.0859999999999999`,1.091`,1.0959999999999999`,1.101`,1.1059999999999999`,1.111`,1.1159999999999999`,1.121`,
1.126`,1.131`,1.136`,1.141`,1.146`,1.151`,1.156`,1.1609999999999998`,1.166`,1.1709999999999998`,1.176`,
1.1809999999999998`,1.186`,1.1909999999999998`,1.196`,1.2009999999999998`,1.206`,1.2109999999999999`,1.216`,
1.2209999999999999`,1.226`,1.2309999999999999`,1.236`,1.2409999999999999`,1.246`,1.251`,1.256`,1.261`,1.266`,
1.271`,1.276`,1.281`,1.2859999999999998`,1.291`,1.2959999999999998`,1.301`,1.3059999999999998`,1.311`,1.3159999999999998`,
1.321`,1.3259999999999998`,1.331`,1.3359999999999999`,1.341`,1.3459999999999999`,1.351`,1.3559999999999999`,1.361`,
1.3659999999999999`,1.371`,1.376`,1.381`,1.386`,1.391`,1.396`,1.401`,1.406`,1.4109999999999998`,1.416`,1.4209999999999998`,
1.426`,1.4309999999999998`,1.436`,1.4409999999999998`,1.446`,1.4509999999999998`,1.456`,1.4609999999999999`,1.466`,
1.4709999999999999`,1.476`,1.4809999999999999`,1.486`,1.4909999999999999`,1.496`,1.501`};


p13jfuncTab=Transpose[{ktabTot,p13jfunc/@ktabTot}];
P13interp=Interpolation[p13jfuncTab];

p22jfuncTab=Transpose[{ktabTot,p22jfunc/@ktabTot}];
P22interp=Interpolation[p22jfuncTab];


Export[LoopTablesPath<>"P22interp.m",P22interp]
Export[LoopTablesPath<>"P13interp.m",P13interp]


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

getMonomialsSingle[vars_,expr_]:=Module[{biaslist,coefrules},
coefrules=CoefficientRules[expr,vars];
biaslist=Sort[Times@@(vars^#)&/@coefrules[[All,1]]];
biaslist
]


BLoop\[Mu]lister=biaslisterfn[{\[Mu]1,u},BLoop\[Mu]List];


(*Cosine of angle between k1 and k2*)
ys[k1_,k2_,k3_]=(-k1^2-k2^2+k3^2)/(2 k1 k2);
(*Cosine of angle between k1 and k3*)
y13s[k1_,k2_,k3_]=(-k1^2-k3^2+k2^2)/(2 k1 k3);
(*Use u variable instead of Sqrt[1-\[Mu]1^2]Sin[\[Phi]] to speed up angular integration later*)
\[Mu]2ey= (y \[Mu]1+Sqrt[1-y^2] u);
\[Mu]3ey=-Sqrt[1-y13^2]u+y13 \[Mu]1;


(* ::Subsubsection::Closed:: *)
(*Tree level*)


Btreeimp=Import[coefpath<>"BtreePerms.m"];


BTree=Btreeimp/.{\[Mu]2->\[Mu]2ey,\[Mu]3->\[Mu]3ey}//Expand;


BTreeCoefs=BLoop\[Mu]lister[BTree];


(*Check decomposition works*)
BTreeCoefs.BLoop\[Mu]List-BTree//Expand


Export[coefpath<>"Btreecoefslistmu.m",BTreeCoefs]


(* ::Subsubsection::Closed:: *)
(*Diagrams*)


b222coeflist=Import[coefpath<>"B222simpcoefslist.m"];
b222coeflistmu=BLoop\[Mu]lister[b222coeflist/.{\[Mu]2->\[Mu]2ey,\[Mu]3->\[Mu]3ey}];
(*Check decomposition works*)
Total[Abs[b222coeflistmu.BLoop\[Mu]List-b222coeflist/.{\[Mu]2->\[Mu]2ey,\[Mu]3->\[Mu]3ey}//Expand]]

b3211coeflist=Import[coefpath<>"B3211simpcoefslist.m"];
b3211coeflistmu=BLoop\[Mu]lister[b3211coeflist/.{\[Mu]2->\[Mu]2ey,\[Mu]3->\[Mu]3ey}];
(*Check decomposition works*)
Total[Abs[b3211coeflistmu.BLoop\[Mu]List-b3211coeflist/.{\[Mu]2->\[Mu]2ey,\[Mu]3->\[Mu]3ey}//Expand]]

b3212coeflist=Import[coefpath<>"B3212simpcoefslist.m"];
b3212coeflistmu=BLoop\[Mu]lister[b3212coeflist/.{\[Mu]2->\[Mu]2ey,\[Mu]3->\[Mu]3ey}];
(*Check decomposition works*)
Total[Abs[b3212coeflistmu.BLoop\[Mu]List-b3212coeflist/.{\[Mu]2->\[Mu]2ey,\[Mu]3->\[Mu]3ey}//Expand]]

b411coeflist=Import[coefpath<>"B411simpcoefslist.m"];
b411coeflistmu=BLoop\[Mu]lister[b411coeflist/.{\[Mu]2->\[Mu]2ey,\[Mu]3->\[Mu]3ey}];
(*Check decomposition works*)
Total[Abs[b411coeflistmu.BLoop\[Mu]List-b411coeflist/.{\[Mu]2->\[Mu]2ey,\[Mu]3->\[Mu]3ey}//Expand]]


Export[coefpath<>"b222coeflistmu.m",b222coeflistmu]
Export[coefpath<>"b3211coeflistmu.m",b3211coeflistmu]
Export[coefpath<>"b3212coeflistmu.m",b3212coeflistmu]
Export[coefpath<>"b411coeflistmu.m",b411coeflistmu]


(* ::Subsection::Closed:: *)
(*Make mu list for bispectrum covariance*)


PTreeCov[\[Mu]_]=(b1+f1 \[Mu]^2)^2 pk+Be1/nb;
invPTree[\[Mu]_]=1/PTreeCov[\[Mu]];


DiogoExpansion[func_,degree_:6]:=Module[{nderivapprox,alist,funcapprox,solacoef,numDoF,eqList},
funcapprox[\[Mu]_]=Sum[a[i] \[Mu]^i,{i,0,degree,2}];
alist = Table[a[i],{i,0,degree,2}];
nderivapprox[\[Mu]_,n_]= Derivative[n][funcapprox][\[Mu]];

equalAt0=funcapprox[0]==func[0];
equalAt1=funcapprox[1]==func[1];
secDerivAt0= nderivapprox[0,2]==Derivative[2][func][0]//Simplify;
equalAtMid= funcapprox[3/4]==func[3/4];
fourthDerivAt0=nderivapprox[0,4]==Derivative[4][func][0]//Simplify;

eqList={equalAt0,equalAt1,secDerivAt0,equalAtMid,fourthDerivAt0};

numDoF=Floor[degree/2+1];
solacoef=Solve[eqList[[;;numDoF]],alist][[1]]//Simplify;

funcapprox[\[Mu]]/.solacoef
];


(*Check*)
result=DiogoExpansion[invPTree];//AbsoluteTiming
subs={pk->10^5,f1->0.7,nb->10^-4,b1->1.,Be1->1};
fourierPtree=Series[invPTree[\[Mu]]/.subs,{\[Mu],0,20}]//Normal;
plotRatio=Plot[Evaluate[{PTreeCov[\[Mu]]result,PTreeCov[\[Mu]] fourierPtree,1}/.subs],{\[Mu],-1,1},
		PlotRange->{0.99,1.01},
		ImageSize->Medium];
plotvalue=Plot[Evaluate[{result,fourierPtree,invPTree[\[Mu]]}/.subs],{\[Mu],-1,1},
		PlotRange->Automatic,
		ImageSize->Medium];
Row[{plotRatio,plotvalue}]


(* ::Text:: *)
(*We get a very good approximation at all \[Mu] already at degree 6 *)


invPTreeApprox[pk_,\[Mu]_]=DiogoExpansion[invPTree];


invCovBtot=Collect[invPTreeApprox[pk1,\[Mu]1]invPTreeApprox[pk2,\[Mu]2ey]invPTreeApprox[pk3,\[Mu]3ey],
				{\[Mu]1,u}];


invCovBtot//Variables


(*invcovnum=invCovBtot/.b1\[Rule]1.9/.f1\[Rule]0.77/.nb\[Rule]4 10^-4/.pk1\[Rule]1 10^4/.pk2\[Rule]1.1 10^4/.pk3\[Rule]1.2 10^4/.Be1\[Rule]1.01;*)
(*covcoefs=getMonomialsSingle[{\[Mu]1,u},invcovnum];*)
covcoefs={1,u^2,u^4,u^6,u^8,u^10,u^12,u \[Mu]1,u^3 \[Mu]1,u^5 \[Mu]1,u^7 \[Mu]1,u^9 \[Mu]1,u^11 \[Mu]1,\[Mu]1^2,u^2 \[Mu]1^2,u^4 \[Mu]1^2,u^6 \[Mu]1^2,u^8 \[Mu]1^2,u^10 \[Mu]1^2,
		u^12 \[Mu]1^2,u \[Mu]1^3,u^3 \[Mu]1^3,u^5 \[Mu]1^3,u^7 \[Mu]1^3,u^9 \[Mu]1^3,u^11 \[Mu]1^3,\[Mu]1^4,u^2 \[Mu]1^4,u^4 \[Mu]1^4,u^6 \[Mu]1^4,u^8 \[Mu]1^4,u^10 \[Mu]1^4,u^12 \[Mu]1^4,
		u \[Mu]1^5,u^3 \[Mu]1^5,u^5 \[Mu]1^5,u^7 \[Mu]1^5,u^9 \[Mu]1^5,u^11 \[Mu]1^5,\[Mu]1^6,u^2 \[Mu]1^6,u^4 \[Mu]1^6,u^6 \[Mu]1^6,u^8 \[Mu]1^6,u^10 \[Mu]1^6,u^12 \[Mu]1^6,u \[Mu]1^7,
		u^3 \[Mu]1^7,u^5 \[Mu]1^7,u^7 \[Mu]1^7,u^9 \[Mu]1^7,u^11 \[Mu]1^7,\[Mu]1^8,u^2 \[Mu]1^8,u^4 \[Mu]1^8,u^6 \[Mu]1^8,u^8 \[Mu]1^8,u^10 \[Mu]1^8,u \[Mu]1^9,u^3 \[Mu]1^9,u^5 \[Mu]1^9,
		u^7 \[Mu]1^9,u^9 \[Mu]1^9,\[Mu]1^10,u^2 \[Mu]1^10,u^4 \[Mu]1^10,u^6 \[Mu]1^10,u^8 \[Mu]1^10,u \[Mu]1^11,u^3 \[Mu]1^11,u^5 \[Mu]1^11,u^7 \[Mu]1^11,\[Mu]1^12,u^2 \[Mu]1^12,u^4 \[Mu]1^12,
		u^6 \[Mu]1^12,u \[Mu]1^13,u^3 \[Mu]1^13,u^5 \[Mu]1^13,\[Mu]1^14,u^2 \[Mu]1^14,u^4 \[Mu]1^14,u \[Mu]1^15,u^3 \[Mu]1^15,\[Mu]1^16,u^2 \[Mu]1^16,u \[Mu]1^17,\[Mu]1^18};


covlister=biaslisterfn[{\[Mu]1,u},covcoefs];
BCovList=covlister[invCovBtot];
(*BCovList//Variables*)
(*{b1,Be1,fb,nb,pk1,pk2,pk3,y,y13}*)
BCovList//Dimensions


Export[LoopTablesPath<>"BCovList.m",BCovList]


(*check*)
(BCovList.covcoefs-invCovBtot)/.{b1->1.9,f1->0.77,nb->4 10^-4,
								pk1->1 10^4,pk2->1.1 10^4,pk3->1.2 10^4,
								Be1->1.01,y->0.4,y13->0.5}//Expand

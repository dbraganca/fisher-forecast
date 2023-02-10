(* ::Package:: *)

BeginPackage["FisherFunctions`"]

changeb1::usage = "Change bias values based on b1 rescaling"
Evolve::usage = "Adjust biases with new time assuming 1/D dependence"

fi::usage = "Number format for jmat imports"; 
Tfi::usage = "String format for jmat imports"; 
EV::usage = "Evaluates function on a table"; 
EV2::usage = "Evaluates function on a table of 3-vectors"; 
BAllcoefs
stochp
nonstoch
IRls
PlinSubs
PlinsubsT
Prior4sky


Begin["`Private`"]

(* Parameters: background cosmology, shift sizes, k-scales *)
lnAsback=3.044; 
nsback = 0.965; 
hback=0.673; 
\[Omega]bback = 0.02237;
\[Omega]cback = 0.1203;
m\[Nu]back = 0.3;
\[Omega]mback = \[Omega]cback+\[Omega]bback+m\[Nu]back/93.14;
\[CapitalOmega]cback = \[Omega]cback/hback^2;
\[CapitalOmega]mback = \[Omega]mback/hback^2;
\[CapitalOmega]bback = \[Omega]bback/hback^2;
CosmoBackground = {lnAsback,nsback,hback,\[CapitalOmega]mback,\[Omega]bback,m\[Nu]back}; 

(*Shift sizes*)
\[CapitalDelta]lnAs = 0.1; 
\[CapitalDelta]h = 0.02; 
\[CapitalDelta]\[CapitalOmega]m=0.04; 
\[CapitalDelta]\[Omega]b = 0.02;
\[CapitalDelta]ns = 0.02; 
\[CapitalDelta]m\[Nu]=0.25;
CosmoShifts = {\[CapitalDelta]lnAs,\[CapitalDelta]ns nsback,\[CapitalDelta]h hback,\[CapitalDelta]\[CapitalOmega]m \[CapitalOmega]mback,\[CapitalDelta]\[Omega]b \[Omega]bback,\[CapitalDelta]m\[Nu] m\[Nu]back};

(* Nmax log k bins, form kmin to kmax *)
kBins[Nmax_,kmin_,kmax_]:=
	Module[{\[CapitalDelta],result},
		\[CapitalDelta]=(1/(Nmax-1))Log[kmax/kmin]//N;
		result=Table[kmin Exp[(i-1) \[CapitalDelta]],{i,1,Nmax}];
		result
];

(*Extend some power spectrum to high values*)
highval[plindat_]:=
Module[{nshigh,Ashigh,pextrahigh,kmax=1000},
	nshigh =(Log[plindat[[-1]]]-Log[plindat[[-2]]])[[2]]/( Log[plindat[[-1]]]-Log[plindat[[-2]]])[[1]];
	Ashigh = plindat[[-1]][[2]]/plindat[[-1]][[1]]^nshigh;

	pextrahigh = Table[{Log10[plindat[[-1]][[1]]]/0.01+0.01+10^((i-1) 0.01), 
						Ashigh ((Log10[plindat[[-1]][[1]]]/0.01+0.01+10^((i-1) 0.01))/plindat[[-1]][[1]])^nshigh},
						{i,1,IntegerPart[Log10[kmax]/0.01]+1}
					];
					
	Join[plindat,pextrahigh]
];



(*Cosmological parameters and fNL parameters*)
cosmopars ={lnAs,ns,h,\[CapitalOmega]m,\[Omega]b,m\[Nu]};
fnlcoefs={fNLloc,fNLeq,fNLorth};

(* All EFT parameters + biases*)
BAllcoefs = {b1,c2,b3,b4,c4,b6,b7,b8,b9,b10,
			 b11,b12,b13,b14,b15,Bc1,Bc2,Bc3,
			 Bc4,Be1,Be2,ce2,Bc5,Bc6,Bc7,Bc8,
			 Bc9,Bc10,Bc11,Bc12,Bc13,Bc14,Bd1,
			 Bd2,Bd3,Be3,Be4,Be5,Be6,Be7,Be8,
			 Be9,Be10,Be11,Be12,Tst};
(*Stochastic terms*)
stochp = {Bd1,Bd2,Bd3,Be1,Be2,ce2,Be3,Be4,Be5,
		  Be6,Be7,Be8,Be9,Be10,Be11,Be12};
(*Non-stochastic terms*)
nonstoch = Complement[BAllcoefs,stochp];


(*Bias transformation*)
(*Change bestfit to new linear bias (this is an approximation, where we assume all biases rescale like the linear one)*)
changeb1[biasfix_,newb1_] := Thread[biasfix[[All,1]]-> biasfix[[All,2]] newb1/biasfix[[1,2]]];
(* Adjust biases with new time assuming 1/D dependence *)
Evolve[Dgold_,Dgnew_,biasfix_]:=Join[Thread[nonstoch->(nonstoch/.biasfix) Dgold/Dgnew],
									 Thread[stochp->(stochp/.biasfix)]];


(* IR resummation functions (using E-Hu)*)

Ho=1/2998 ;
Tcmb=2.728;
alph[\[Omega]mo_,\[Omega]bo_]:=1-0.328 Log[431 \[Omega]mo] \[Omega]bo/\[Omega]mo+0.38 Log[22.3  \[Omega]mo ](\[Omega]bo/\[Omega]mo)^2;
ss[h_,\[Omega]mo_,\[Omega]bo_]:=44.5 Log[9.83/(\[Omega]mo)]/Sqrt[1+10(\[Omega]bo)^0.75] h;
gamma[k_,h_,\[Omega]mo_,\[Omega]bo_]:=\[Omega]mo/h  (alph[\[Omega]mo,\[Omega]bo]+(1-alph[\[Omega]mo,\[Omega]bo])/(1+(0.43 ss[h,\[Omega]mo,\[Omega]bo] k)^4));
qq[k_,h_,\[Omega]mo_,\[Omega]bo_]:= (Tcmb/2.7)^2 k/gamma[k,h,\[Omega]mo,\[Omega]bo];
L0[k_,h_,\[Omega]mo_,\[Omega]bo_]:=Log[2 E+1.8 qq[k,h,\[Omega]mo,\[Omega]bo]];
C0[k_,h_,\[Omega]mo_,\[Omega]bo_]:=14.2+731/(1+62.5 qq[k,h,\[Omega]mo,\[Omega]bo]);
Tnw[k_,h_,\[Omega]mo_,\[Omega]bo_]:=L0[k,h,\[Omega]mo,\[Omega]bo]/(L0[k,h,\[Omega]mo,\[Omega]bo]+C0[k,h,\[Omega]mo,\[Omega]bo]qq[k,h,\[Omega]mo,\[Omega]bo]^2)
Dlin[a_,h_,\[Omega]mo_]:=-(a/(3 Sqrt[\[Omega]mo/h^2 (a^3 (1-\[Omega]mo/h^2)+\[Omega]mo/h^2)])) (-8 \[Omega]mo/h^2 Hypergeometric2F1[-(1/2),5/6,11/6,-((a^3 (1-\[Omega]mo/h^2))/(\[Omega]mo/h^2))]+(2 a^3 (1-\[Omega]mo/h^2)+5 \[Omega]mo/h^2) Hypergeometric2F1[1/2,5/6,11/6,-((a^3 (1-\[Omega]mo/h^2))/(\[Omega]mo/h^2))]);
Dlinz[z_,h_,\[Omega]mo_]:=Dlin[1/(1+z),h,\[Omega]mo];
kP[h_]:=0.05 h^-1;
PLnw[k_,z_,h_,\[CapitalOmega]mo_,\[Omega]bo_,lnAs_,ns_]:=2 \[Pi]^2 4/25 (Exp[lnAs]*(10^-10) k )/(Ho^4 \[CapitalOmega]mo^2) Dlinz[z,h,\[CapitalOmega]mo h^2]^2 (k /kP[h])^(ns-1) Tnw[k,h,\[CapitalOmega]mo h^2,\[Omega]bo]^2;

pwnwdat[dat_,z_,h_,\[CapitalOmega]mo_,\[Omega]bo_,lnAs_,ns_]:=
Module[{plindat, iplin, klin, plinEHtab, iR,fk,pnw,pw,ipnw,ipw,\[CapitalSigma]2,\[Delta]\[CapitalSigma]2,\[Lambda]=0.25},
	plindat=dat;
	klin=plindat[[All,1]];
	plinEHtab = Table[{klin[[i]],PLnw[klin[[i]],z,h,\[CapitalOmega]mo,\[Omega]bo,lnAs,ns]},{i,1,Length[klin]}];
	iR=Interpolation[Transpose@{klin,plindat[[All,2]]/plinEHtab[[All,2]]}];
	fk[k_]:=NIntegrate[iR[10^x] PDF[NormalDistribution[Log10[k],\[Lambda]]][x],{x,Log10[k]-4 \[Lambda],Log10[k]+4 \[Lambda]},
			PrecisionGoal->5,Method->{"GlobalAdaptive","SymbolicProcessing"->0}];	
	
	pnw=Table[{klin[[i]],plinEHtab[[i,2]] fk[klin[[i]]]},{i,1,Length[klin]}];
	pw=Table[{klin[[i]],plindat[[i,2]]-pnw[[i,2]]},{i,1,Length[klin]}];
	ipnw=Interpolation[pnw];
	ipw=Interpolation[pw];
	\[CapitalSigma]2=1/(6 \[Pi]^2) NIntegrate[(1-SphericalBesselJ[0,q losc]+2 SphericalBesselJ[2,q losc])ipnw[q]/.losc->110,{q,klin[[1]],0.2},
			PrecisionGoal->5,Method->{"GlobalAdaptive","SymbolicProcessing"->0}];
	\[Delta]\[CapitalSigma]2=1/(2 \[Pi]^2) NIntegrate[SphericalBesselJ[2,q losc]ipnw[q]/.losc->110,{q,klin[[1]],0.2},
		    PrecisionGoal->5,Method-> {"GlobalAdaptive","SymbolicProcessing"->0}];
		    
	{pnw,pw,\[CapitalSigma]2,\[Delta]\[CapitalSigma]2}
];

pwnwavg[dat_,f_,z_,h_,\[CapitalOmega]mo_,\[Omega]bo_,lnAs_,ns_]:=
Module[{full,klin,plindat,pnw,pw,ipnw,ipw,\[CapitalSigma]2,\[Delta]\[CapitalSigma]2,\[CapitalSigma]tot2avg},
	{pnw,pw,\[CapitalSigma]2,\[Delta]\[CapitalSigma]2} = pwnwdat[dat,z,h,\[CapitalOmega]mo,\[Omega]bo,lnAs,ns];
	plindat=dat;
	klin=plindat[[All,1]];

	\[CapitalSigma]tot2avg = -(1/15) f^2 (2 \[Delta]\[CapitalSigma]2-5 \[CapitalSigma]2)+\[CapitalSigma]2+(2 f \[CapitalSigma]2)/3;

	Table[{klin[[i]],pnw[[i,2]]+ Exp[-klin[[i]]^2\[CapitalSigma]tot2avg]pw[[i,2]]},{i,1,Length[klin]}]
];

pwnwavgTree[dat_,f_,z_,h_,\[CapitalOmega]mo_,\[Omega]bo_,lnAs_,ns_]:=
Module[{full,klin,plindat,pnw,pw,ipnw,ipw,\[CapitalSigma]2,\[Delta]\[CapitalSigma]2,\[CapitalSigma]tot2avg},
	{pnw,pw,\[CapitalSigma]2,\[Delta]\[CapitalSigma]2}= pwnwdat[dat,z,h,\[CapitalOmega]mo,\[Omega]bo,lnAs,ns];
	plindat=dat;
	klin=plindat[[All,1]];
	\[CapitalSigma]tot2avg = -(1/15) f^2 (2 \[Delta]\[CapitalSigma]2-5 \[CapitalSigma]2)+\[CapitalSigma]2+(2 f \[CapitalSigma]2)/3;
	Table[{klin[[i]],pnw[[i,2]]+(1+klin[[i]]^2 \[CapitalSigma]tot2avg) Exp[-klin[[i]]^2\[CapitalSigma]tot2avg]pw[[i,2]]},{i,1,Length[klin]}]
];

IRrep[splitwnw_] := {\[CapitalSigma]2-> splitwnw[[3]],\[Delta]\[CapitalSigma]2-> splitwnw[[4]]};
\[CapitalSigma]tot2[\[Mu]_,f_]:=(1+f \[Mu]^2 (2+f))\[CapitalSigma]2+f^2 \[Mu]^2 (\[Mu]^2-1)\[Delta]\[CapitalSigma]2;


(*Fitting power spectrum*)
genPlinCoef[imax_,dat_]:=
Module[{kmaxFit = 0.6, kminFit = 0.001, knFit, k, 
		MyPsmoothtab,errors,Parameters,
		MyFitFunctions,nlm,Plin,Pinteger,coefnowig},
	Plin=Interpolation[dat];
	knFit=kBins[100,kminFit,kmaxFit];
	MyPsmoothtab = Table[{knFit[[i]],Plin[knFit[[i]]]},{i,1,Length[knFit]}];
	errors=Flatten[Table[{MyPsmoothtab[[i,2]] 10^-6},{i,1,Length[knFit]}]];

	kuv0=0.0001;
	kuv1=0.069;
	kuv2=0.0082;
	kuv3=0.0013;
	kuv4=0.0000135;
	kpeak1=-0.034;
	kpeak2=-0.001;
	kpeak3=-0.000076;
	kpeak4=-0.0000156;

	Parameters=Flatten[Table[{a1[i],a2[i],a3[i],a4[i]},{i,0,imax}]];
	MyFitFunctions=Delete[Flatten[
				Join[{1/(1+k^2/kuv0)},
					Table[{ 1/(1+(k^2-kpeak1)^2/kuv1^2)^(i), (k/(1/20))^2/(1+(k^2-kpeak2)^2/kuv2^2)^(i+1),1/(1+(k^2-kpeak3)^2/kuv3^2)^(i+2),1/(1+(k^2-kpeak4)^2/kuv4^2)^(i+1)},{i,0,imax}]]
					],
				2];
	nlm=LinearModelFit[MyPsmoothtab,MyFitFunctions, k, Weights->1/errors^2, IncludeConstantBasis->False]; 
	Pinteger[k_]=nlm[k];
	coefnowig=N[nlm["BestFitParameters"], 50];
	
	{Pinteger[q],coefnowig}
];


(*IRls[dat_,hfactor_,\[CapitalOmega]mfacor_,\[Omega]bfacor_,Asfactor_,nsfactor_,zpk_,f1_]:=
Module[{iplin,splitwnw,ls,lsnw,iplinnw,iplinw,sig,iplinresum,iplinresumTree},
	iplin =highval[dat]//Interpolation;
	splitwnw=pwnwdat[dat,zpk,hback*hfactor,\[CapitalOmega]mback*\[CapitalOmega]mfacor,\[Omega]bback*\[Omega]bfacor,lnAsback+Asfactor,nsback*nsfactor];
	ls=genPlinCoef[3,dat][[2]];
lsnw=genPlinCoef[3,splitwnw[[1]]][[2]];
iplinnw = highval[splitwnw[[1]]]//Interpolation;
iplinw = highval[splitwnw[[2]]]//Interpolation;
sig=\[CapitalSigma]tot2[\[Mu],f1]/.IRrep[splitwnw];
iplinresum = pwnwavg[dat,f1,zpk,hback*hfactor,\[CapitalOmega]mback*\[CapitalOmega]mfacor,\[Omega]bback*\[Omega]bfacor,lnAsback+Asfactor,nsback*nsfactor]//Interpolation; 
iplinresumTree = pwnwavgTree[dat,f1,zpk,hback*hfactor,\[CapitalOmega]mback*\[CapitalOmega]mfacor,\[Omega]bback*\[Omega]bfacor,lnAsback+Asfactor,nsback*nsfactor]//Interpolation; 
{ls,lsnw,iplin,iplinnw,iplinw,sig,iplinresum,iplinresumTree,f1}
]*)

(*Final loading function for single power spectrum*)
IRls[dat_,{h_,\[CapitalOmega]m_,\[Omega]b_,lnAs_,ns_},zpk_,f1_]:=
Module[{iplin,splitwnw,ls,lsnw,iplinnw,iplinw,sig,iplinresum,iplinresumTree},
	iplin =highval[dat]//Interpolation;
	splitwnw=pwnwdat[dat,zpk,h,\[CapitalOmega]m,\[Omega]b,lnAs,ns];
	ls=genPlinCoef[3,dat][[2]];
	lsnw=genPlinCoef[3,splitwnw[[1]]][[2]];
	iplinnw = highval[splitwnw[[1]]]//Interpolation;
	iplinw = highval[splitwnw[[2]]]//Interpolation;
	sig=\[CapitalSigma]tot2[\[Mu]1,f1]/.IRrep[splitwnw];
	iplinresum = pwnwavg[dat,f1,zpk,h,\[CapitalOmega]m,\[Omega]b,lnAs,ns]//Interpolation; 
	iplinresumTree = pwnwavgTree[dat,f1,zpk,h,\[CapitalOmega]m,\[Omega]b,lnAs,ns]//Interpolation; 
	
	{ls,lsnw,iplin,iplinnw,iplinw,sig,iplinresum,iplinresumTree,f1}
];


(*Load power spectra and transfer functions*)
CosmoPksFactors={
(*Filename \[Rule] {h,\[CapitalOmega]m,\[Omega]b,lnAs,ns}*)
"pk_fid.dat"-> {1,1,1,0,1},

"pk_lnAs_pl.dat"->{1,1,1,\[CapitalDelta]lnAs,1},
"pk_lnAs_min.dat"->{1,1,1,-\[CapitalDelta]lnAs,1},

"pk_ns_pl.dat"->{1,1,1,0,1+\[CapitalDelta]ns},
"pk_ns_min.dat"->{1,1,1,0,1-\[CapitalDelta]ns},

"pk_h_pl.dat"->{1+\[CapitalDelta]h,1,1,0,1},
"pk_h_min.dat"->{1-\[CapitalDelta]h,1,1,0,1},

"pk_Om_pl.dat"->{1,1+\[CapitalDelta]\[CapitalOmega]m,1,0,1},
"pk_Om_min.dat"->{1,1-\[CapitalDelta]\[CapitalOmega]m,1,0,1},

"pk_ob_pl.dat"->{1,1,1+\[CapitalDelta]\[Omega]b,0,1},
"pk_ob_min.dat"->{1,1,1-\[CapitalDelta]\[Omega]b,0,1},

"pk_nu_pl.dat"-> {1,1,1,0,1},
"pk_nu_min.dat"-> {1,1,1,0,1}
};

PlinSubs[surveypath_, zpk_, cosmoback_: CosmoBackground]:=
Module[{cosmonames,cosmopaths,hb,\[CapitalOmega]mb,\[Omega]bb,lnAsb,nsb,m\[Nu]b,factors,cosmoparams,importPS,
	    flist,fback,f\[CapitalOmega]mhigh,f\[CapitalOmega]mlow,f\[Omega]bhigh,f\[Omega]blow,fhhigh,fhlow,fm\[Nu]high,fm\[Nu]low},

(*output is list of {ls,lsnw,iplin,iplinnw,iplinw,sig} for each cosmology*)
(*No need to include neutrino mass because this is only for purposes of IR resummation*)
{lnAsb,nsb,hb,\[CapitalOmega]mb,\[Omega]bb,m\[Nu]b}=cosmoback;

(*Get f's*)
{fback,f\[CapitalOmega]mhigh,f\[CapitalOmega]mlow,f\[Omega]bhigh,f\[Omega]blow,fhhigh,fhlow,fm\[Nu]high,fm\[Nu]low}=Flatten[Import[surveypath<>"allfs.dat"]];
flist={
	fback,
	fback,fback,(*corresponding to lnAs*)
	fback,fback,(*corresponding to ns*)
	fhhigh,fhlow,
	f\[CapitalOmega]mhigh,f\[CapitalOmega]mlow,
	f\[Omega]bhigh,f\[Omega]blow,
	fm\[Nu]high,fm\[Nu]low
};

cosmonames = CosmoPksFactors[[All,1]];
factors = CosmoPksFactors[[All,2]];
cosmopaths=(surveypath<>#)&/@cosmonames;

cosmoparams={hb #[[1]],\[CapitalOmega]mb #[[2]],\[Omega]bb #[[3]],lnAsb+#[[4]],nsb #[[5]]}&/@factors;

importPS=Import[#]&/@cosmopaths;

Table[IRls[importPS[[i]],cosmoparams[[i]],zpk,flist[[i]]],{i,1,Length[cosmoparams]}]
];

GetT[surveypath_]:=
Module[{TkFiles = {"Tk_fid.dat",
				   "Tk_lnAs_pl.dat","Tk_lnAs_min.dat",
				   "Tk_ns_pl.dat","Tk_ns_min.dat",
				   "Tk_h_pl.dat","Tk_h_min.dat",
				   "Tk_Om_pl.dat","Tk_Om_min.dat",
				   "Tk_ob_pl.dat","Tk_ob_min.dat",
				   "Tk_nu_pl.dat","Tk_nu_min.dat"}},
	Table[Interpolation[Import[surveypath<>TkFiles[[i]]]],{i,1,Length[TkFiles]}]
	];
	
PlinsubsT[surveypath_, zpk_]:=Module[{Ts,Pks},
	Pks = PlinSubs[surveypath,zpk];
	Ts = GetT[surveypath];
	Table[Join[Pks[[i]],{Ts[[i]]}],{i,1,Length[Ts]}]
];


(*Utility functions*)

(*Formating jmat imports*)
fi[a_] := NumberForm[a, {5, 5}]; 
Tfi[a_] := ToString[fi[a]];

(*Function evaluation*) 
EV[f_, tab_] := Monitor[Table[f[tab[[s]]], {s, 1, Length[tab]}], s]; 
EV2[f_, tab_] := Monitor[Table[f[tab[[s, 1]], tab[[s, 2]], tab[[s, 3]]], {s, 1, Length[tab]}], s]; 

(*Functions for efficient matrix manipulation*)
Upperright[mat_]:= UpperTriangularize[mat];
Mirror[mat_]:=mat+Transpose[mat]-DiagonalMatrix[Diagonal[mat]];

(*Removes Row and Col n from matrix A*)
RemoveRowCol[A_,n_]:=Delete[Transpose[Delete[Transpose[A], n]],n];
(*Removes Rows and Cols in removelist from matrix A*)
RemoveRowCols[A_,removelist_]:=Module[{keeplist},
keeplist=Complement[Range[Length[A]],removelist];
A[[keeplist,keeplist]]
];

(*Generates triangles where each side is greater than kmin and less than kmax*)
makeTri[kmin_,kmax_,dk_]:=Module[{kpairs,kpairstriangle,triangles,ktab},
ktab = Table[k,{k,kmin,kmax,dk}];
triangles = Select[Subsets[ktab, {3}], 
    #[[1]] + #[[2]] > #[[3]] &];
 Sort[triangles]
];

(*Get keff from triangle and bin size*)
Getkeff[{k1_,k2_,k3_},\[CapitalDelta]k_]:=Module[{k1eff,k2eff,k3eff,VT},
VT = NIntegrate[q1 q2 q3 Boole[q1<=q2+q3]Boole[q2<=q1+q3]Boole[q3<=q1+q2],
	{q1,k1-\[CapitalDelta]k/2,k1+\[CapitalDelta]k/2},{q2,k2-\[CapitalDelta]k/2,k2+\[CapitalDelta]k/2},{q3,k3-\[CapitalDelta]k/2,k3+\[CapitalDelta]k/2}
	];
k1eff = 1/VT NIntegrate[q1 q2 q3 Min[{q1,q2,q3}] Boole[q1<=q2+q3]Boole[q2<=q1+q3]Boole[q3<=q1+q2],
		{q1,k1-\[CapitalDelta]k/2,k1+\[CapitalDelta]k/2},{q2,k2-\[CapitalDelta]k/2,k2+\[CapitalDelta]k/2},{q3,k3-\[CapitalDelta]k/2,k3+\[CapitalDelta]k/2}
		];
k2eff =1/VT NIntegrate[q1 q2 q3 Median[{q1,q2,q3}] Boole[q1<=q2+q3]Boole[q2<=q1+q3]Boole[q3<=q1+q2],
		{q1,k1-\[CapitalDelta]k/2,k1+\[CapitalDelta]k/2},{q2,k2-\[CapitalDelta]k/2,k2+\[CapitalDelta]k/2},{q3,k3-\[CapitalDelta]k/2,k3+\[CapitalDelta]k/2}
		];
k3eff = 1/VT NIntegrate[q1 q2 q3 Max[{q1,q2,q3}] Boole[q1<=q2+q3]Boole[q2<=q1+q3]Boole[q3<=q1+q2],
		{q1,k1-\[CapitalDelta]k/2,k1+\[CapitalDelta]k/2},{q2,k2-\[CapitalDelta]k/2,k2+\[CapitalDelta]k/2},{q3,k3-\[CapitalDelta]k/2,k3+\[CapitalDelta]k/2}
		];
		
{k1eff,k2eff,k3eff}
];

(*Calculate keff for all triangles generated by makeTri*)
makeTrikeff[kmin_,kmax_,dk_]:=Module[{tris},
tris=makeTri[kmin,kmax,dk];
Monitor[Table[Getkeff[tris[[i]],dk],{i,1,Length[tris]}],(i*100.0)/Length[tris]]
];

(*Get powerspectrum multipoles*)
mono[f_]:=(f//Expand)/.\[Mu]1^n_.-> (1+(-1)^n)/(2 (1+n));
quad[f_]:=(f-(f/.\[Mu]1->0)//Expand)/.\[Mu]1^n_.-> (5 (1+(-1)^n) n)/(2 (1+n) (3+n));

(* Make list of derivatives *)
(* The order of derivatives here is {lnAs, ns, h, \[CapitalOmega]m, \[Omega]b, \[Nu]m, log(b1), BAllcoefs[[2;;]]}*)
Makederiv[prederiv_]:=Module[{EFTderiv, cosmoderiv, derivs, fnlderiv},
(*prederiv is an array of the resummed power spectra for all shifted parameters,
prederiv\[LeftDoubleBracket]1\[RightDoubleBracket] is the fiducial power spectrum*)
EFTderiv=Table[D[prederiv[[1]],BAllcoefs[[i]]],{i,1,Length[BAllcoefs]}];
EFTderiv[[1]]=EFTderiv[[1]]b1;(*corresponding to log b1*)
cosmoderiv = (prederiv[[2;; ;;2]]-prederiv[[3;; ;;2]])/(2 CosmoShifts);
fnlderiv=Table[D[prederiv[[1]],fnlcoefs[[i]]],{i,1,Length[fnlcoefs]}];

derivs= Join[cosmoderiv,EFTderiv];
derivs
];


(*Fisher matrix manipulations*)

numCosmoFnlParams = Length[cosmoparams]+Length[fnlcoefs];

(*Combine Fishers from different skies*)
CombineMs[FisherMs_]:=Module[{MAu,MBu,double,l,Msu},
lenOtherParams = Length[FisherMs[[1]]]-numCosmoFnlParams;
numSkies = Length[FisherMs];
Msu = UpperTriangularize/@FisherMs;

double = Table[0,{cosmopar+numSkies*lenOtherParams},{cosmopar+numSkies*lenOtherParams}];

(*sum over cosmological parameters*)
double[[1;;numCosmoFnlParams,1;;numCosmoFnlParams]]=Total[Msu[[All,1;;numCosmoFnlParams,1;;numCosmoFnlParams]]];


(*fill remaining entries*)
For[i=1,i<numCosmoFnlParams+1,i++,
	For[j=numCosmoFnlParams+1,j<Length[double],j+=l,
		For[s=0,s<l,s++,
			double[[i,j+s]]=Msu[[s+1,i,(j-(numCosmoFnlParams+1))/l+(numCosmoFnlParams+1)]];
		]
	]
];
For[i=cosmopar+1,i<Length[double],i+=l,
	For[j=cosmopar+1,j<Length[double],j+=l,
		For[s=0,s<l,s++,
			double[[i+s,j+s]]=Msu[[s+1,(i-(cosmopar+1))/l+cosmopar+1,(j-(cosmopar+1))/l+cosmopar+1]];
		]
	]
];
Mirror[double]
];

(*Some Fisher matrix  manipulations *)

AllParams =Join[cosmopars,fnlcoefs,BAllcoefs];

(*Get index of specific parameters in AllParams*)
Getpos[params_]:=Table[Position[AllParams,params[[i]]],{i,1,Length[params]}]//Flatten;

(*Get Sub fisher corresponding to specific parameters*)
SubFisher[Fisher_,params_]:=Fisher[[Getpos[params],Getpos[params]]];


BAllcoefs4sky = Join[{BAllcoefs},{BAllcoefs},{BAllcoefs},{BAllcoefs}]//Transpose//Flatten;
AllParams4sky =Join[cosmopars,fnlcoefs,BAllcoefs4sky];

(*Get index of specific parameters in AllParams4sky*)
Getpos4sky[params_]:=Table[Position[AllParams4sky,params[[i]]],{i,1,Length[params]}]//Flatten;

(*Get Sub fisher corresponding to specific parameters*)
SubFisher4sky[Fisher_,params_]:=Fisher[[Getpos4sky[params],Getpos4sky[params]]];

(*Get std dev from fisher matrix*)
GetErrors[Fisher_]:=Quiet[(Diagonal[(Inverse[Fisher])])^(1/2)];

(*surveySpec = {zlist,nblist,Vlist,Nlist,Dglist,fslist}*)

(*zeff is given by the zs average weighted by N *)
Getzeffs[surveySpec_,cut_]:={Total[surveySpec[[1,;;cut]] surveySpec[[4,;;cut]]]/Total[surveySpec[[4,;;cut]]],
							Total[surveySpec[[1,cut+1;;]] surveySpec[[4,cut+1;;]]]/Total[surveySpec[[4,cut+1;;]]]};
							
							
fixsurv[surveySpec_,dkk_] := Table[{dk->dkk[[2]],nb->surveySpec[[2,i]],Vs->surveySpec[[3,i]]},{i,1,Length[surveySpec[[1]]]}];


Getnbeff[surveySpec_,cut_]:={Total[surveySpec[[2,;;cut]]surveySpec[[4,;;cut]]]/Total[surveySpec[[4,;;cut]]],
					        Total[surveySpec[[2,cut+1;;]]surveySpec[[4,cut+1;;]]]/Total[surveySpec[[4,cut+1;;]]]};
(* fundamental mode of the survey*)
getkf[surveySpec_]:=2\[Pi]/surveySpec[[3]]^(1/3);

(* effective fundamental mode of the survey*)
getkfeff[surveySpec_,cut_]:={Total[getkf[surveySpec][[;;cut]]surveySpec[[4,;;cut]]]/Total[surveySpec[[4,;;cut]]],
							Total[getkf[surveySpec][[cut+1;;]]surveySpec[[4,cut+1;;]]]/Total[surveySpec[[4,cut+1;;]]]};


(* Priors *)
BBN\[Nu]cosmoprior = DiagonalMatrix[Join[
	{0,0,0,0,1/ (\[Sigma]\[Omega]b^2),0},(*these are the cosmological parameters*)
	{0,0,0},(*these are the fnl parameters*)
	Table[0,{i,1,Length[BAllcoefs]}](*these are the EFT parameters*)
	]];
	
EFTpriorDiag={1/\[Sigma]b1^2,1/\[Sigma]c2^2,1/\[Sigma]b3^2,1/\[Sigma]b4^2,1/\[Sigma]c4^2,1/\[Sigma]b6^2,1/\[Sigma]b7^2,1/\[Sigma]b8^2,1/\[Sigma]b9^2,1/\[Sigma]b10^2,
			 1/\[Sigma]b11^2,1/\[Sigma]b12^2,1/\[Sigma]b13^2,1/\[Sigma]b14^2,1/\[Sigma]b15^2,1/\[Sigma]Bc1^2,1/\[Sigma]Bc2^2,1/\[Sigma]Bc3^2,1/\[Sigma]Bc4^2,
			 1/\[Sigma]Be1^2,1/\[Sigma]Be2^2,1/\[Sigma]ce2^2,1/\[Sigma]Bc5^2,1/\[Sigma]Bc6^2,1/\[Sigma]Bc7^2,1/\[Sigma]Bc8^2,1/\[Sigma]Bc9^2,1/\[Sigma]Bc10^2,
			 1/\[Sigma]Bc11^2,1/\[Sigma]Bc12^2,1/\[Sigma]Bc13^2,1/\[Sigma]Bc14^2,1/\[Sigma]Bd1^2,1/\[Sigma]Bd2^2,1/\[Sigma]Bd3^2,1/\[Sigma]Be3^2,1/\[Sigma]Be4^2,
			 1/\[Sigma]Be5^2,1/\[Sigma]Be6^2,1/\[Sigma]Be7^2,1/\[Sigma]Be8^2,1/\[Sigma]Be9^2,1/\[Sigma]Be10^2,1/\[Sigma]Be11^2,1/\[Sigma]Be12^2,1/\[Sigma]Tst^2};
			 
EFTprior=DiagonalMatrix[Join[
	{0,0,0,0,0,0},(*these are the cosmological parameters*)
	{0,0,0},(*these are the fnl parameters*)
	EFTpriorDiag(*these are the EFT parameters*)
]];

priorrep=Thread[{\[Sigma]\[Omega]b, 
				\[Sigma]b1,\[Sigma]c2,\[Sigma]b3,\[Sigma]b4,\[Sigma]c4,\[Sigma]b6,\[Sigma]b7,\[Sigma]b8,\[Sigma]b9,\[Sigma]b10,
				\[Sigma]b11,\[Sigma]b12,\[Sigma]b13,\[Sigma]b14,\[Sigma]b15,\[Sigma]Bc1,\[Sigma]Bc2,\[Sigma]Bc3,\[Sigma]Bc4,
				\[Sigma]Be1,\[Sigma]Be2,\[Sigma]ce2,\[Sigma]Bc5,\[Sigma]Bc6,\[Sigma]Bc7,\[Sigma]Bc8,\[Sigma]Bc9,\[Sigma]Bc10,
				\[Sigma]Bc11,\[Sigma]Bc12,\[Sigma]Bc13,\[Sigma]Bc14,\[Sigma]Bd1,\[Sigma]Bd2,\[Sigma]Bd3,\[Sigma]Be3,\[Sigma]Be4,
				\[Sigma]Be5,\[Sigma]Be6,\[Sigma]Be7,\[Sigma]Be8,\[Sigma]Be9,\[Sigma]Be10,\[Sigma]Be11,\[Sigma]Be12,\[Sigma]Tst}->
				{0.000036,
				Sqrt[0.8],2,2,2,2,2,2,2,2,2,
				2,2,2,2,2,4,4,4,2,2,
				4,2,2,2,2,2,2,2,2,2,
				2,2,2,2,2,2,2,2,2,2,
				2,2,2,2,2,2}];

priors =(BBN\[Nu]cosmoprior+EFTprior)/.priorrep;

(*Create a block Sparse Array whose non-zero blocks (the input matrices) are in the diagonal *)
blockArray[mat_]:=SparseArray[(Tuples[Range@#-{1,0,0}].{Rest@#,{1,0},{0,1}}&@Dimensions@mat)-> Flatten@mat];

\[Rho]sub={\[Rho]12 -> (1-0.1^2/2),\[Rho]13->(1-0.2^2/2)};

buildCorrPrior[priorlist_, cosmolen_]:=Module[{fullPriorpre,fullPrior,corrmat},
corrmat={{1, \[Rho]12, \[Rho]13,\[Rho]12 \[Rho]13},
		{\[Rho]12,1, \[Rho]12 \[Rho]13, \[Rho]13},
		{\[Rho]13, \[Rho]12 \[Rho]13, 1, \[Rho]12},
		{\[Rho]12 \[Rho]13, \[Rho]13, \[Rho]12, 1}}/.\[Rho]sub;
fullPriorpre=blockArray[Table[Diagonal[priorlist][[i]]Inverse[corrmat],{i,cosmolen+1, Length[priorlist]}]];
fullPrior=ArrayPad[fullPriorpre,{{cosmolen,0},{cosmolen,0}}]//Normal;
fullPrior[[5,5]]=1/\[Sigma]\[Omega]b^2;
fullPrior[[6,6]]=1;
fullPrior];

(*Dimensions of prior4sky should be 4*45+6 - 1 bias prior for each sky*)
Prior4sky = buildCorrPrior[priors, numCosmoFnlParams]/.priorrep;


(* Powerspectrum Fisher *)
(*Preliminaries*)

(*Counterterm definition*)
CounterP13=(b1+f1 \[Mu]^2)plin[k](-((Bc1 k^2)/km^2)+(f1 (Bc2-(Bc4 f1)/2) k^2 \[Mu]^2)/kr^2-(Bc3 f1^2 k^2 \[Mu]^4)/(2 kr^2));
StochasticP22=1/nb (Be1+(Be2 k^2)/km^2+(3 ce2 k^2 \[Mu]^2)/(2 km^2));
PCounter[k_,plin_,f1_]=2 CounterP13+StochasticP22;

(*Loop definitions and imports*)
DumpGet["tabs/P13interp.mx"]; 
DumpGet["tabs/P22interp.mx"];

P13coefs[f1_] = Import["tabs/P13simpcoefs.m"]; 
P13kers[k_] = Import["tabs/P13simppermks.m"];
P13n[k_]:=P13kers[k].P13interp[k];
P13ev[P13_,ls_,f_,Plink_]:=P13.ls.P13coefs[f] Plink;

P22coefs[f1_] = Import["tabs/P22simpcoefs.m"]/.b2->(c2 + c4)/Sqrt[2]/.b5->(c2 - c4)/Sqrt[2]; 
P22kers[k_] = Import["tabs/P22simppermks.m"];
P22n[k_]:=P22kers[k].P22interp[k];
P22ev[P22_,ls_,f_]:=P22.ls.ls.P22coefs[f];

(*fNL definitions*)
PNG[k_,{fnl_,\[Beta]s_},{plin_,f1_,T\[Alpha]_}]:=(4 (-1+b1) fnl (k/knl)^\[Beta]s \[Delta]c (b1+f1 \[Mu]^2) plin[k])/T\[Alpha][k];
PNGfull[k_,{plin_,f1_,T\[Alpha]_}]:=(PNG[k, {fNLloc,0},{plin,f1,T\[Alpha]}]
							 +PNG[k,{fNLeq,2},{plin,f1,T\[Alpha]}]
							 +PNG[k,{fNLorth,2},{plin,f1,T\[Alpha]}]);
fNLparamfix={\[Delta]c->1.68,p->8.52};


(*Derivatives Resummed*)

(*Tree level derivatives*)
PTreeRes[k_,{pnw_,pw_,sig_,f1_}]:=(b1+f1 \[Mu]1^2)^2 (pnw[k]+(1+k^2 sig)Exp[-k^2sig]pw[k]);

PTreeDerivs[k_,plinsubs_]:=Module[{prederiv},
	(*plinsubs\[LeftDoubleBracket]i\[RightDoubleBracket] is {ls,lsnw,iplin,iplinnw,iplinw,sig,iplinresum,iplinresumTree,f1,Tk}*)
	prederiv = Table[PTreeRes[k,plinsubs[[i]][[{4,5,6,9}]]],{i,1,Length[plinsubs]}];
	Makederiv[prederiv]];

(*Counter-term derivatives*)
PCounterRes[k_,{plin_,plinnw_,sig_,f1_}]:=(PCounter[k,plinnw,f1](1-Exp[-k^2sig])
										   +Exp[-k^2sig]PCounter[k,plin,f1]);
PCounterDerivs[k_,plinsubs_]:=Module[{prederiv},
	prederiv = Table[PCounterRes[k,plinsubs[[i]][[{3,4,6,9}]]],{i,1,Length[plinsubs]}];
	Makederiv[prederiv]
];

(*Loop derivatives*)
PLoopDerivs[k_,pks_]:=Module[{P13mat,P22mat,Plint,Ploopeval,Ploopresum,PTreeresum,prederiv,EFTderiv,cosmoderiv,derivs},
	P13mat=P13n[k]; 
	P22mat=P22n[k]; 
	Ploopeval[ls_,plin_,f_]:=P13ev[P13mat,ls,f,plin[k]]+P22ev[P22mat,ls,f];
	
	
	(*function to do IR resummation*)
	Ploopresum[{plin_,plinnw_,fitcoefs_,fitcoefsnw_,sig_,f1_}]:=
	Ploopeval[fitcoefsnw,plinnw,f1](1-Exp[-k^2sig])+Exp[-k^2sig]Ploopeval[fitcoefs,plin,f1];
		
	prederiv = Table[Ploopresum[pks[[i]][[{3,4,1,2,6,9}]]],{i,1,Length[pks]}];
	Makederiv[prederiv]];

(*fNL derivatives*)
PNGDerivs[k_,pks_]:=Module[{PTreeresum,prederiv,EFTderiv,cosmoderiv,derivs},
	prederiv = Table[PNGfull[k,pks[[i]][[{3,9,10}]]],{i,1,Length[pks]}];
	Makederiv[prederiv]];

(* Full P1loop derivatives*)	
P1Loopderivs[k_,pks_]:=Series[PLoopDerivs[k,pks]+PTreeDerivs[k,pks]
							 +PCounterDerivs[k,pks]+PNGDerivs[k,pks],
							 {\[Mu],0,10}]//Normal;


(*Derivatives non-resummed*)

(*WIP*)


End[]
EndPackage[]

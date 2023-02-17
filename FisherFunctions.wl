(* ::Package:: *)

(* ::Section:: *)
(*Public section*)


BeginPackage["FisherFunctions`"]

(*public variables to be used outside*)
{b1,c2,b2,Bb3,b3,Bb4,b4,b5,c4,Bb6,b6,Bb7,b7,Bb8,b8,Bb9,b9,Bb10,b10,
Bb11,b11,b12,b13,b14,b15,Bc1,Bc2,Bc3,
Bc4,Be1,Be2,ce2,Bc5,Bc6,Bc7,Bc8,
Bc9,Bc10,Bc11,Bc12,Bc13,Bc14,Bd1,
Bd2,Bd3,Be3,Be4,Be5,Be6,Be7,Be8,
Be9,Be10,Be11,Be12,Tst,
nb,knl,
hback};

fi::usage = "Number format for jmat imports"; 
Tfi::usage = "String format for jmat imports"; 
EV::usage = "Evaluates function on a table"; 
EV2::usage = "Evaluates function on a table of 3-vectors"; 

mono::usage = "integrates over the angles. Equivalent to taking the monopole.";
quad::usage = "Extracts quadrupole contribution.";
removelg2::usage = "Extracts only the monopole and quadrupole contributions from a given quantity";

BAllcoefs::usage = "lists all cosmological+EFT+fNL parameters used";
changeb1::usage = "Change bias values based on b1 rescaling";
Evolve::usage = "Adjust biases with new time assuming 1/D dependence";
shiftedbfixss::usage = "shift biases assuming same relation to b1";
normbias

fnlfix;
paramfix;



stochp;
nonstoch;
IRls;
PlinSubs;
PlinsubsT::usage = "loads pks for a given survey";
Prior4sky;

Getnbeff::usage = "get nb effectives for survey";
Getzeffs::usage = "get z effectives for survey";

P1LoopDerivs::usage = "Calculates an array of 1-loop power-spectrum derivatives with respect to all parameters"; 
BTreeDerivs::usage = "Calculates an array of tree-level bi-spectrum derivatives with respect to all parameters"; 
B1LoopDerivs::usage = "Calculates an array of 1-loop bi-spectrum derivatives with respect to all parameters"; 

GetFullzCovP::usage = "Get covariances for power-spectrum for different z's";
GetFullzCovPNS::usage = "Get covariances for power-spectrum with no shot noise for different z's";
GetFullzCovPn::usage = "Get covariances for power-spectrum with shot noise given by (n x Be1)for different z's";

GetFullzCovB::usage = "Get covariances for bi-spectrum for different z's";
GetFullzCovBNS::usage = "Get covariances for bi-spectrum with no shot noise for different z's";

GetFisherPnew::usage = "calculate power-spectrum full Fisher matrix using analytical covariance and derivatives";
GetFisherPnewNum::usage = "calculates Fisher matrix using numerical covariance and derivatives";

BTreemaster::usage = "Master angular integral tensor for tree level bispectrum";
BLoopmaster::usage = "Master angular integral tensor for 1-loop bispectrum";
BLoopmastermono::usage = "Master angular integral tensor for 1-loop bispectrum monopole";
GetFisherB::usage = "calculates bi-spectrum Fisher matrix using analytical covariance, derivatives, and precomputed master integral tensor";
GetFisherBNum::usage = "calculates bi-spectrum monopole Fisher matrix using numerical covariance and derivatives";

jmatPath;
loopTablesPath;
covsPath;

jmatimp;

TabExport::usage = "Run to generate and save k-dependent coefficients to be used later";
SubFisher;
SubFisher4sky;

B222kers;

PTreeShot;

Bcovlist;


(* ::Section:: *)
(*Private section*)


(* ::Text:: *)
(*TODO :*)
(* Careful: biaslist not matching with Yaniv*)
(* I cannot import .mx: Should use either .m or .wl (or plain.dat)*)


Begin["`Private`"]


(* ::Subsubsection:: *)
(*Preliminaries: paths and imports*)


jmatPath = "C:\\Users\\diogo\\Dropbox\\FFTLog\\Diogo_work\\GitHub\\python-integer-power-project\\2. Jmat_loopvals\\";
loopTablesPath = "tabs/";
saveLoopCoefsPath = "saveloops/";
covsPath = "covsks/";


(* ::Subsubsection::Closed:: *)
(*Background cosmology and shifts (including fNL)*)


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

(*Cosmological parameters and fNL parameters*)
cosmopars ={lnAs,ns,h,\[CapitalOmega]m,\[Omega]b,m\[Nu]};
(*fNL parameters*)
fnlcoefs={fNLloc,fNLeq,fNLorth};
fnlfix = {fNLloc->0,fNLeq->0,fNLorth->0};
paramfix = {\[Delta]c->1.68,p->8.52};

(* k-scales *)
kr=knl/2.8`;
km=knl;


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



(* ::Subsubsection::Closed:: *)
(*EFT-parameters and related functions*)


(* All EFT parameters + biases*)
BAllcoefs = {b1,c2,b3,b4,c4,b6,b7,b8,b9,b10,
			 b11,b12,b13,b14,b15,Bc1,Bc2,Bc3,
			 Bc4,Be1,Be2,ce2,Bc5,Bc6,Bc7,Bc8,
			 Bc9,Bc10,Bc11,Bc12,Bc13,Bc14,Bd1,
			 Bd2,Bd3,Be3,Be4,Be5,Be6,Be7,Be8,
			 Be9,Be10,Be11,Be12,Tst};
			 
			 
b25Toc24Sub= {b2->(c2+c4)/Sqrt[2],b5->(c2-c4)/Sqrt[2]};

biases = {b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15};
normbias =Thread[{Bb1,Bb2,Bb3,Bb4,Bb5,Bb6,Bb7,Bb8,Bb9,Bb10,Bb11,Bb12,Bb13,Bb14,Bb15}->biases];

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
									 
(*Change bestfit to new redshift*)
cfixresc[Dgold_,Dgnew_,biasfix_] := Join[Thread[stochp->(stochp/.biasfix)],Thread[biasfix[[All,1]]->biasfix[[All,2]] Dgold/Dgnew],fnlfix,paramfix]
shiftedbfix[spec_,Dg_,biasfix_]:=Table[cfixresc[Dg,spec[[5,i]],biasfix],{i,1,Length[spec[[5]]]}];
cfixrescss[Dgold_,Dgnew_,biasfix_] := Join[Thread[biasfix[[All,1]]->biasfix[[All,2]] Dgold/Dgnew],fnlfix,paramfix]
shiftedbfixss[spec_,Dg_,biasfix_]:=Table[cfixrescss[Dg,spec[[5,i]],biasfix],{i,1,Length[spec[[5]]]}];

(*Change bestfit to new redshift and no shotnoise*)
cfixrescNS[Dgold_,Dgnew_,biasfix_] := Join[{Be1->0},Thread[stochp->(stochp/.biasfix)],Thread[biasfix[[All,1]]->biasfix[[All,2]] Dgold/Dgnew],fnlfix,paramfix]//Flatten;
shiftedbfixNS[spec_,Dg_,biasfix_]:=Table[cfixrescNS[Dg,spec[[5,i]],biasfix],{i,1,Length[spec[[5]]]}];

(*Change bestfit to new redshift and no shotnoise*)
cfixrescn[Dgold_,Dgnew_,biasfix_] := Join[{Be1->n (Be1/.biasfix)},Thread[stochp->(stochp/.biasfix)],Thread[biasfix[[All,1]]->biasfix[[All,2]] Dgold/Dgnew],fnlfix,paramfix]//Flatten;
shiftedbfixn[spec_,Dg_,biasfix_]:=Table[cfixrescn[Dg,spec[[5,i]],biasfix],{i,1,Length[spec[[5]]]}];


(* ::Subsubsection::Closed:: *)
(*IR-resummation using Eisenstein-Hu*)


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
\[CapitalSigma]tot2[\[Mu]1_,f_]:=(1+f \[Mu]1^2 (2+f))\[CapitalSigma]2+f^2 \[Mu]1^2 (\[Mu]1^2-1)\[Delta]\[CapitalSigma]2;


(* ::Subsubsection::Closed:: *)
(*Plin decomposition*)


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


(* ::Subsubsection::Closed:: *)
(*Final functions to get fitting coefficients from IR-resummed power spectrum*)


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

PlinSubs[surveypath_String, zpk_, cosmoback_: CosmoBackground]:=
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


(* ::Subsubsection:: *)
(*Utility functions for evaluating tables and Fisher matrix manipulations*)


(*Utility functions*)

(*Formating jmat imports*)
fi[k_] := Round[k, N[10^-5]]; 
Tfi[k_] := ToString[fi[k]];

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
cosmoderiv[[6]]=m\[Nu]back cosmoderiv[[6]];(*corresponding to log m\[Nu]*)
fnlderiv=Table[D[prederiv[[1]],fnlcoefs[[i]]],{i,1,Length[fnlcoefs]}];

derivs= Join[cosmoderiv,EFTderiv];
derivs
];


(*define rules for very fast integration - especially important for bispectrum*)
usub={u->Sqrt[1-\[Mu]1^2]Sin[\[Phi]]};
(*angint[m_,n_]=1/(4\[Pi])Integrate[\[Mu]1^nu^m/.usub,{\[Mu]1,-1,1},{\[Phi],0,2\[Pi]},Assumptions\[Rule]{-1<\[Mu]1<1,0<\[Phi]<2\[Pi],m>0,n>0,m\[Element]Integers,n\[Element]Integers}]*)
angint[m_,n_]=((1+(-1)^m) (1+(-1)^n) Gamma[(1+m)/2] Gamma[(1+n)/2])/(8 Sqrt[\[Pi]] Gamma[1/2 (3+m+n)]);
rules\[Mu]1u={\[Mu]1^n_. u^m_.->angint[m,n],\[Mu]1^n_.->angint[0,n],u^m_.->angint[m,0]};


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


(* effective nb for the survey*)
Getnbeff[surveySpec_,cut_]:={Total[surveySpec[[2,;;cut]]surveySpec[[4,;;cut]]]/Total[surveySpec[[4,;;cut]]],
					        Total[surveySpec[[2,cut+1;;]]surveySpec[[4,cut+1;;]]]/Total[surveySpec[[4,cut+1;;]]]};
(* fundamental mode of the survey*)
getkf[surveySpec_]:=2\[Pi]/surveySpec[[3]]^(1/3);

(* effective fundamental mode of the survey*)
getkfeff[surveySpec_,cut_]:={Total[getkf[surveySpec][[;;cut]]surveySpec[[4,;;cut]]]/Total[surveySpec[[4,;;cut]]],
							Total[getkf[surveySpec][[cut+1;;]]surveySpec[[4,cut+1;;]]]/Total[surveySpec[[4,cut+1;;]]]};


(* ::Subsubsection::Closed:: *)
(*Priors for cosmological and EFT parameters for 1 and 4 skies*)


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


(* ::Subsection:: *)
(*Power-spectrum Derivatives up to 1-loop - Resummed and Non-resummed*)


(* Powerspectrum Fisher *)
(*Preliminaries*)

(*Counterterm definition*)
CounterP13=(b1+f1 \[Mu]1^2)plin[k](-((Bc1 k^2)/km^2)+(f1 (Bc2-(Bc4 f1)/2) k^2 \[Mu]1^2)/kr^2-(Bc3 f1^2 k^2 \[Mu]1^4)/(2 kr^2));
StochasticP22=1/nb (Be1+(Be2 k^2)/km^2+(3 ce2 k^2 \[Mu]1^2)/(2 km^2));
PCounter[k_,plin_,f1_]=2 CounterP13+StochasticP22;

(*Loop definitions and imports*)
P13interp = Import[loopTablesPath<>"P13interp.m"];
P22interp = Import[loopTablesPath<>"P22interp.m"];

P13coefs[f1_] = Import[loopTablesPath<>"P13simpcoefs.m"]/.\[Mu]->\[Mu]1; 
P13kers[k_] = Import[loopTablesPath<>"P13simppermks.m"];


P13n[k_]:=P13kers[k].P13interp[k];
P13ev[P13_,ls_,f_,Plink_]:=P13.ls.P13coefs[f] Plink;

P22coefs[f1_] = Import[loopTablesPath<>"P22simpcoefs.m"]/.b2->(c2 + c4)/Sqrt[2]/.b5->(c2 - c4)/Sqrt[2]/.\[Mu]->\[Mu]1; 
P22kers[k_] = Import[loopTablesPath<>"P22simppermks.m"];
P22n[k_]:=P22kers[k].P22interp[k];
P22ev[P22_,ls_,f_]:=P22.ls.ls.P22coefs[f];

(*fNL definitions*)
PNG[k_,{fnl_,\[Beta]s_},{plin_,f1_,T\[Alpha]_}]:=(4 (-1+b1) fnl (k/knl)^\[Beta]s \[Delta]c (b1+f1 \[Mu]1^2) plin[k])/T\[Alpha][k];
PNGfull[k_,{plin_,f1_,T\[Alpha]_}]:=(PNG[k, {fNLloc,0},{plin,f1,T\[Alpha]}]
							 +PNG[k,{fNLeq,2},{plin,f1,T\[Alpha]}]
							 +PNG[k,{fNLorth,2},{plin,f1,T\[Alpha]}]);
fNLparamfix={\[Delta]c->1.68,p->8.52};


(*Check dimensions compatibility to do dot product*)
Dimensions[P13kers[0.1]]=={Dimensions[P13coefs[0.7]][[1]],
						  Dimensions[P13interp[0.1]][[1]]}
Dimensions[P22kers[0.1]]=={Dimensions[P22coefs[0.7]][[1]],
						  Dimensions[P22interp[0.1]][[1]]}


(*fNL derivatives*)
PNGDerivs[k_,pks_]:=Module[{PTreeresum,prederiv,EFTderiv,cosmoderiv,derivs},
	prederiv = Table[PNGfull[k,pks[[i]][[{3,9,10}]]],{i,1,Length[pks]}];
	Makederiv[prederiv]];


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

(* Full P1loop derivatives*)	
P1LoopDerivs[k_,pks_]:=Series[PLoopDerivs[k,pks]+PTreeDerivs[k,pks]
							 +PCounterDerivs[k,pks]+PNGDerivs[k,pks],
							 {\[Mu]1,0,10}]//Normal;


(*Derivatives non-resummed*)

PTreeNR[k_,{plin_,f1_}]:=(b1+f1 \[Mu]1^2)^2 plin[k];

(*not resumed*)
PTreeDerivsNR[k_,plinsubs_]:=Module[{prederiv},
prederiv = Table[PTreeNR[k,plinsubs[[i]][[{3,9}]]],{i,1,Length[plinsubs]}];
Makederiv[prederiv]]

PLoopDerivsNR[k_,pks_]:=Module[{P13mat,P22mat,Plint,Ploopeval,PloopNR,prederiv,EFTderiv,cosmoderiv,derivs},
P13mat=P13n[k]; 
P22mat=P22n[k]; 
Plint[Plin_]=Plin[k];
Ploopeval[ls_,plin_,f_]:= P13ev[P13mat,ls,f,plin[k]]+ P22ev[P22mat,ls,f];
PloopNR[{ls_,plin_,f_}]:=Ploopeval[ls, plin, f];
prederiv = Table[PloopNR[pks[[i]][[{1,3,9}]]],{i,1,Length[pks]}];
Makederiv[prederiv]];


P1LoopDerivsNR[k_,pks_]:=Series[PLoopDerivsNR[k,pks]+PTreeDerivsNR[k,pks]+
								PCounterDerivsNR[k,pks]+PNGDerivs[k,pks],
								{\[Mu]1,0,10}]//Normal;


(* ::Subsubsection:: *)
(*Covariance*)


(*Covariance*)
PTreeShot[k_,plin_,f1_]:=(b1+f1 \[Mu]1^2)^2 plin[k]+Be1/nb;
ExpandInvPsq[k_,plin_,f1_,bias_]:=Normal[Series[1/(PTreeShot[k,plin,f1]/.bias)^2,{\[Mu]1,0,30}]];

InvCovP[k_,pks_,bias_]:=( (4\[Pi]^2)/(k^2 dk Vs))^-1 ExpandInvPsq[k,pks[[1,3]],pks[[1,9]],bias];
InvCovP2[k_,pks_,bias_,Survey_,f1_]:=( (4\[Pi]^2)/(k^2 dk Vs))^-1 ExpandInvPsq[k,pks[[1,3]],f1, bias]/.Survey;

CovPmulti[k_,pks_,biastab_,Surveytab_,fs_]:=Monitor[Sum[InvCovP2[k,pks,biastab[[i]],Surveytab[[i]],fs[[i]]],{i,1,Length[biastab]}],i]

GetFullzCovP[ks_,pks_,spec_,Dg_,biasfix_,dkk_]:=CovPmulti[ks,pks,shiftedbfix[spec,Dg,biasfix],fixsurv[spec,dkk],spec[[-1]]]//Expand;
GetFullzCovPNS[ks_,pks_,spec_,Dg_,biasfix_,dkk_]:=CovPmulti[ks,pks,shiftedbfixNS[spec,Dg,biasfix],fixsurv[spec,dkk],spec[[-1]]]//Expand;
GetFullzCovPn[ks_,pks_,spec_,Dg_,biasfix_,dkk_]:=CovPmulti[ks,pks,shiftedbfixn[spec,Dg,biasfix],fixsurv[spec,dkk],spec[[-1]]]//Expand;



(* ::Subsubsection:: *)
(*Fisher Matrices*)


(*Calculate derivatives and covariances for each k*)
GetDerivs[ktab_,pks_]:=P1LoopDerivs[#,pks]&/@ktab;
GetDerivsNR[ktab_,pks_]:=P1LoopDerivsNR[#,pks]&/@ktab;

GetInvCovs[ktab_,pks_,bias_]:=InvCovP[#,pks,bias]&/@ktab;

(*Fisher with Numerical covariance*)
GetFisherPNum[ktab_,biasfix_,pks_,Survey_,Covnum_]:=Module[{hartlapFactor,derivs,monoquad},
	derivs = Table[P1Loopderivs[ktab[[i]],pks],{i,1,Length[ktab]}]/.biasfix/.Survey;
	monoquad = Join[mono[derivs],quad[derivs]];
	hartlapFactor=(2048 - Length[Covnum] - 2)/(2048 - 1);
	hartlapFactor*(Transpose[monoquad].Inverse[Covnum].monoquad)];
	
(*Fisher with Numerical covariance - Non-Resummed*)
GetFisherPNumNR[ktab_,biasfix_,pks_,Survey_,Covnum_]:=Module[{hartlapFactor,derivs,monoquad},
	derivs = Table[P1LoopderivsNR[ktab[[i]],pks],{i,1,Length[ktab]}]/.biasfix/.Survey;
	monoquad = Join[mono[derivs],quad[derivs]];
	hartlapFactor=(2048 - Length[Covnum] - 2)/(2048 - 1);
	hartlapFactor*(monoquad//Transpose).Inverse[Covnum].monoquad];
	
(*Fisher with Analytical covariance monopole and quadrupole*)

(*only keep first two multiploes of expression a (function of \[Mu]1)*)
removelg2[a_]:=mono[a]+quad[a]LegendreP[2,\[Mu]1];

GetFisherP02[ktab_,biasfix_,pks_,Survey_]:=Module[{derivred,derivs,invcovs,Fpre,ktabnew},
	derivs = GetDerivs[ktab, pks]/.Dispatch[Join[biasfix,Survey]];
	derivred = removelg2[derivs];
	
	invcovs= GetInvCovs[ktab,pks,biasfix]/.Dispatch[Join[biasfix,Survey]];
	
	Fpre =Monitor[Sum[(KroneckerProduct[derivs[[i]],derivs[[i]]])invcovs[[i]],{i,1,Length[derivs]}],i];
	Mirror[mono[Upperright[Fpre]]]];
	
(*Fisher with Analytical covariance monopole and quadrupole no shotnoise*)
GetFisherP02nos[ktab_,biasfix_,biasfixNoShot_,pks_,Survey_]:=Module[{derivred,derivs,invcovs,Fpre,ktabnew},
	(*calculate derivatives*)
	derivs = GetDerivs[ktab, pks]/.Dispatch[Join[biasfix,Survey]];
	derivred = removelg2[derivs];
	(*calcualte inverse covariances*)
	invcovs= GetInvCovs[ktab,pks,biasfixNoShot]/.Dispatch[Join[biasfixNoShot,Survey]];
	(*calculate Fisher*)
	Fpre =Monitor[Sum[(KroneckerProduct[derivs[[i]],derivs[[i]]])invcovs[[i]],{i,1,Length[derivs]}],i];
	Mirror[mono[Upperright[Fpre]]]];

(*Fisher with Analytical covariance*)
GetFisherP[ktab_,biasfix_,pks_,Survey_]:=Module[{derivred,derivs,invcovs,Fpre,ktabnew},
	derivs = GetDerivs[ktab, pks]/.Dispatch[Join[biasfix,Survey]];
	invcovs= GetInvCovs[ktab,pks,biasfix]/.Dispatch[Join[biasfix,Survey]];
	
	Fpre =Monitor[Sum[(KroneckerProduct[derivs[[i]],derivs[[i]]])cov[[i]],{i,1,Length[derivs]}],i];
	Mirror[mono[Upperright[Fpre]]]];

(*Get Fisher List*)
GetFisherPnew[derivs_,covs_] := Mirror[mono[Upperright[Monitor[Sum[(KroneckerProduct[derivs[[i]],derivs[[i]]])covs[[i]],{i,1,Length[derivs]}],i]]]];
GetFisherPnewNum[derivs_,covs_] := (2048 - Length[covs] - 2)/(2048 - 1)*(derivs//Transpose).Inverse[covs].derivs
GetFisherlistP[derivs_, covs_]:=Monitor[Table[Mirror[mono[Upperright[(KroneckerProduct[derivs[[i]],derivs[[i]]])covs[[i]]]]],{i,1,Length[derivs]}],i];



(* ::Subsection:: *)
(*Bi-spectrum derivatives up to 1-loop*)


jmatimp[triang_,diagram_String]:=Import[jmatPath<>diagram<>"/"<>
								 Tfi[triang[[1]]]<>"_"<>Tfi[triang[[2]]]<>"_"<>Tfi[triang[[3]]]<>"_.h5","jmat"];

saveLoopCoefPath[triangle_,diagram_String]:=(saveLoopCoefsPath <> diagram <>"/"
										<>Tfi[triangle[[1]]]<>"_"<>Tfi[triangle[[2]]]<>"_"<>Tfi[triangle[[3]]]<>"_.m");
ExportKer[file_,func_,triang_]:=If[Not[FileExistsQ[file]],Export[file,func[triang]]];
										
(*import coefficients (k-dependent)*)
B222kers[{k1_,k2_,k3_}] = Transpose[Import[loopTablesPath<>"B222simppermks.m"]];
B3211kers[{k1_,k2_,k3_}] = TensorTranspose[Import[loopTablesPath<>"B3211simppermks.m"],{3,1,2}];
B3212kers[{k1_,k2_,k3_}] = TensorTranspose[Import[loopTablesPath<>"B3212simppermks.m"],{3,1,2}];
B411kers[{k1_,k2_,k3_}] = TensorTranspose[Import[loopTablesPath<>"B411simppermks.m"],{3,1,2}];

(*Evaluate coefficients at the necessary ks*)
TabExport[triangs_]:=Module[{B222precomp,file1,file2,file3,file4},
	Monitor[
		Do[
		ExportKer[saveLoopCoefPath[triangs[[i]],"B222"],B222kers,triangs[[i]]];
		ExportKer[saveLoopCoefPath[triangs[[i]],"B3211"],B3211kers,triangs[[i]]];
		ExportKer[saveLoopCoefPath[triangs[[i]],"B3212"],B3212kers,triangs[[i]]];
		ExportKer[saveLoopCoefPath[triangs[[i]],"B411"],B411kers,triangs[[i]]],
		{i,1,Length[triangs]}]
		,i]
		];													


(* ::Subsubsection:: *)
(*Preliminary functions and definitions*)


(*Angle definitions*)
ys[k1_,k2_,k3_]=(-k1^2-k2^2+k3^2)/(2 k1 k2); 
y13s[k1_,k2_,k3_]=(-k1^2+k2^2-k3^2)/(2 k1 k3);
(*Tree term definition*)
BTree[k1_,k2_,k3_,Plin_,f1_,y_,y13_]= Import[loopTablesPath<>"Treelist.m"]/.normbias/.b25Toc24Sub;

(*Counter term definition*)
counterlist2=Import[loopTablesPath<>"counterlistG.m"]/.b25Toc24Sub;

BCounter[k1_,k2_,k3_,Plin_,f1_,y_,y13_]=counterlist2;

(*Loop definitions and imports B[loop]coefs is the decomposition of the bias terms*)
B222coefs[f_,y_,y13_]=Import[loopTablesPath<>"B222simpcoefslist.m"]/.normbias/.b25Toc24Sub;
B3211coefs[f_,y_,y13_]=Import[loopTablesPath<>"B3211simpcoefslist.m"]/.normbias/.b25Toc24Sub;
B3212coefs[f_,y_,y13_]=Import[loopTablesPath<>"B3212simpcoefslist.m"]/.normbias/.b25Toc24Sub;
B411coefs[f_,y_,y13_]=Import[loopTablesPath<>"B411simpcoefslist.m"]/.{b12->Bb11,b11->Bb12}/.normbias/.b25Toc24Sub;

(*Combine jmat and k-dependent coef so that we get a quantity for each permutation, bias combination, fitting coefficients*)
B222n[k1_,k2_,k3_]:=Import[saveLoopCoefPath[{k1,k2,k3},"B222"]].jmatimp[{k1,k2,k3},"B222"];

B3211n[k1_,k2_,k3_]:=Module[{savedCoefs,jmat},
		savedCoefs = Import[saveLoopCoefPath[{k1,k2,k3},"B3211"]]; 
		jmat = TensorTranspose[jmatimp[{k1,k2,k3},"B3211"],{1,3,4,2}];
		Table[(savedCoefs[[i,All]]).(jmat[[All,i]]),{i,1,6}]
		];
		
B3212n[k1_,k2_,k3_]:=Module[{savedCoefs,jmat},
		savedCoefs =Import[saveLoopCoefPath[{k1,k2,k3},"B3212"]]; 
		jmat=TensorTranspose[jmatimp[{k1,k2,k3},"B3212"],{1,3,2}];
		Table[(savedCoefs[[i,All]]).(jmat[[All,i]]),{i,1,6}]
		];
		
B411n[k1_,k2_,k3_]:=Module[{saveker,jmat},
		saveker = Import[saveLoopCoefPath[{k1,k2,k3},"B411"]]; 
		jmat = TensorTranspose[jmatimp[{k1,k2,k3},"B411"],{1,3,2}];
		Table[(saveker[[i,All]]).(jmat[[All,i]]),{i,1,3}]
		];
		
B222ev[B222_,ls_, biascoefs_] := B222.ls.ls.ls.biascoefs;
B3211ev[B3211_,ls_,biascoefs_,PTab_] := Flatten[B3211,1].ls.ls.Flatten[biascoefs PTab,1];
B3212ev[B3212_,ls_,biascoefs_,PTab_] := Flatten[B3212,1].ls.Flatten[biascoefs PTab,1];
B411ev[B411_,ls_,biascoefs_,PTab_] := Flatten[B411,1].ls.Flatten[biascoefs PTab,1];

(*Subscript[f, NL] definitions and imports*)
Bloc[k1_,k2_,k3_]:=2(1/(k1^3 k2^3)+1/(k2^3 k3^3)+1/(k3^3 k1^3));
Bequil[k1_,k2_,k3_]:= 6(-1/(k1^3 k2^3)-1/(k1^3 k3^3)-1/(k2^3 k3^3)-2/(k1^2 k2^2 k3^2)+(1/(k1 k2^2 k3^3)+1/(k1 k3^2 k2^3)+1/(k2 k1^2 k3^3)+1/(k2 k3^2 k1^3)+1/(k3 k1^2 k2^3)+1/(k3 k2^2 k1^3)));
Borth[k1_,k2_,k3_]:=6 (((-k1+k2+k3)(k1-k2+k3)(k1+k2-k3))/(k1^3 k2^3 k3^3) (1+p)-p 1/(k1^4 k2^4 k3^4) (2/3 (k1 k2+k2 k3+k3 k1)-1/3 (k1^2+k2^2+k3^2))^3);
BNG[{k1_,k2_,k3_},{fnl_,\[Beta]s_,FNG_},{plin_,f1_,T\[Alpha]_}]=fnl Import[loopTablesPath<>"BNGlist.m"];
(*BNGFull[{k1_,k2_,k3_},{plin_,f1_,T\[Alpha]_}]:=BNG[{k1,k2,k3},{fNLloc,0,Bloc},{plin,f1,T\[Alpha]}]+BNG[{k1,k2,k3},{fNLeq,2,Bequil},{plin,f1,T\[Alpha]}]+BNG[{k1,k2,k3},{fNLorth,2,Borth},{plin,f1,T\[Alpha]}];*)
BNGFull[{k1_,k2_,k3_},{plin_,f1_,T\[Alpha]_}]=Import[loopTablesPath<>"BNGsimp.m"]/.y13->(-k1^2+k2^2-k3^2)/(2 k1 k3)/.y->(-k1^2-k2^2+k3^2)/(2 k1 k2)/.paramfix;


(* ::Subsubsection::Closed:: *)
(*fNL terms derivatives*)


BNGDerivs[{k1_,k2_,k3_},pks_]:=Module[{prederiv,BNGev},
	BNGev[{plin_,f1_,T\[Alpha]_}]=BNGFull[{k1,k2,k3},{plin,f1,T\[Alpha]}]/.paramfix;
	prederiv = Table[BNGev[pks[[i]][[{3,9,10}]]],{i,1,Length[pks]}];
	Makederiv[prederiv]];


(* ::Subsubsection::Closed:: *)
(*Derivatives: re-summed*)


(* Tree level *)

BTreeDerivs[{k1_,k2_,k3_},pks_]:=Module[{y,y13,Btreeev,prederiv},
	y=ys[k1,k2,k3];
	y13=y13s[k1,k2,k3];
	Btreeev[{Plin_,f_}] = BTree[k1,k2,k3,Plin,f,y,y13];
	prederiv = Table[Btreeev[pks[[i]][[{8,9}]]],{i,1,Length[pks]}];
	Makederiv[prederiv]];
	
(* counterterm derivatives*)
BCounterDerivs[{k1_,k2_,k3_},pks_]:=Module[{y,y13,Bcounterev,prederiv},
	y=ys[k1,k2,k3];
	y13=y13s[k1,k2,k3];
	Bcounterev[{Plin_,f_}]=BCounter[k1,k2,k3,Plin,f,y,y13];
	prederiv = Table[Bcounterev[pks[[i]][[{7,9}]]],{i,1,Length[pks]}];
	Makederiv[prederiv]
	];

(*Loop derivatives*)
BLoopDerivs[{k1_,k2_,k3_},pks_]:=Module[
	{B222,B3211,B3212,B411,
	B222coef,B3211coef,B3212coef,B411coef,
	Bloopeval,y,y13,PTab3211,PTab3212,PTab411,prederiv},
	y=ys[k1,k2,k3];
	y13=y13s[k1,k2,k3];

	B222=B222n[k1,k2,k3];
	B3211=B3211n[k1,k2,k3];
	B3212=B3212n[k1,k2,k3];
	B411=B411n[k1,k2,k3];
	
	B222coef[f_]=B222coefs[f,y,y13];
	B3211coef[f_]=B3211coefs[f,y,y13];
	B3212coef[f_]=B3212coefs[f,y,y13];
	B411coef[f_]=B411coefs[f,y,y13];
	PTab3211[P11_] = {P11[k1],P11[k1],P11[k2],P11[k2],P11[k3],P11[k3]};
	PTab3212[P11_] = {P11[k1] P11[k2],P11[k1] P11[k3],P11[k1] P11[k2],
					  P11[k2] P11[k3],P11[k1] P11[k3],P11[k2] P11[k3]};
	PTab411[P11_] = {P11[k1] P11[k2],P11[k2] P11[k3],P11[k1] P11[k3]};
	Bloopeval[{ls_,plin_,f_}]:= (B222ev[B222, ls, B222coef[f]]
							   + B3211ev[B3211, ls, B3211coef[f], PTab3211[plin]]
							   + B3212ev[B3212, ls, B3212coef[f], PTab3212[plin]]
							   + B411ev[B411, ls, B411coef[f], PTab411[plin]]);


	(*The order of derivatives here is {lnAs, ns, h, \[CapitalOmega]m, \[Omega]b, \[Nu]m, log(b1), BAllcoefs[[2;;]]}*)
	prederiv = Table[Bloopeval[pks[[i]][[{1,7,9}]]],{i,1,Length[pks]}];
	Makederiv[prederiv]];
	
(*Full B1loop derivatives*)
B1LoopDerivs[triangle_,pks_]:= (BLoopDerivs[triangle,pks]+BTreeDerivs[triangle,pks]
							  + BCounterDerivs[triangle,pks]+BNGDerivs[triangle,pks]);
BLoopDerivs[triangle_,pks_]:=BLoopDerivs[triangle,pks]+BCounterDerivs[triangle,pks];
BTreeDerivs[triangle_,pks_]:=BTreeDerivs[triangle,pks]+BNGDerivs[triangle,pks];


(* ::Subsubsection::Closed:: *)
(*Derivatives: not re-summed*)


(*Tree derivatives*)
BTreeDerivsNR[{k1_,k2_,k3_},pks_]:=Module[{y,y13,Btreeevpre,Btreeev,prederiv},
	y=ys[k1,k2,k3];
	y13=y13s[k1,k2,k3];
	Btreeev[{Plin_,f_}] = BTree[k1,k2,k3,Plin,f,y,y13];
	prederiv = Table[Btreeev[pks[[i]][[{3,9}]]],{i,1,Length[pks]}];
	Makederiv[prederiv]]

(* counterterm derivatives*)
BCounterDerivsNR[{k1_,k2_,k3_},pks_]:=Module[{y,y13,Bcounterev,prederiv},
	y=ys[k1,k2,k3];
	y13=y13s[k1,k2,k3];
	Bcounterev[{Plin_,f_}]=BCounter[k1,k2,k3,Plin,f,y,y13];
	prederiv = Table[Bcounterev[pks[[i]][[{3,9}]]],{i,1,Length[pks]}];
	Makederiv[prederiv]
	];

(*Loop derivatives*)
BLoopDerivsNR[{k1_,k2_,k3_},pks_]:=Module[
	{B222,B3211,B3212,B411,
	B222coef,B3211coef,B3212coef,B411coef,
	Bloopeval,y,y13,PTab3211,PTab3212,PTab411,prederiv},
	y=ys[k1,k2,k3];
	y13=y13s[k1,k2,k3];

	B222=B222n[k1,k2,k3];
	B3211=B3211n[k1,k2,k3];
	B3212=B3212n[k1,k2,k3];
	B411=B411n[k1,k2,k3];
	
	B222coef[f_]=B222coefs[f,y,y13];
	B3211coef[f_]=B3211coefs[f,y,y13];
	B3212coef[f_]=B3212coefs[f,y,y13];
	B411coef[f_]=B411coefs[f,y,y13];
	PTab3211[P11_] = {P11[k1],P11[k1],P11[k2],P11[k2],P11[k3],P11[k3]};
	PTab3212[P11_] = {P11[k1] P11[k2],P11[k1] P11[k3],P11[k1] P11[k2],
					  P11[k2] P11[k3],P11[k1] P11[k3],P11[k2] P11[k3]};
	PTab411[P11_] = {P11[k1] P11[k2],P11[k2] P11[k3],P11[k1] P11[k3]};
	Bloopeval[{ls_,plin_,f_}]:= (B222ev[B222, ls, B222coef[f]]
							   + B3211ev[B3211, ls, B3211coef[f], PTab3211[plin]]
							   + B3212ev[B3212, ls, B3212coef[f], PTab3212[plin]]
							   + B411ev[B411, ls, B411coef[f], PTab411[plin]]);


	(*The order of derivatives here is {lnAs, ns, h, \[CapitalOmega]m, \[Omega]b, \[Nu]m, log(b1), BAllcoefs[[2;;]]}*)
	prederiv = Table[Bloopeval[pks[[i]][[{1,3,9}]]],{i,1,Length[pks]}];
	Makederiv[prederiv]];
	

(*Full B1loop derivatives*)
B1LoopDerivsNR[triangle_,pks_]:= (BLoopDerivsNR[triangle,pks]+BTreeDerivsNR[triangle,pks]
							  + BCounterDerivsNR[triangle,pks]+BNGDerivs[triangle,pks]);
BLoopDerivsNR[triangle_,pks_]:=BLoopDerivsNR[triangle,pks]+BCounterDerivsNR[triangle,pks];
BTreeDerivsNR[triangle_,pks_]:=BTreeDerivsNR[triangle,pks]+BNGDerivs[triangle,pks];


(* ::Subsubsection:: *)
(*Covariance*)


(* Analytical covariance *)

(*Inverse covariance already decompose in combinations of u and \[Mu]1*)
Bcovlist =Import[loopTablesPath<>"Covlist.m"]/.ctP22[1]->Be1;
sb[k1_,k2_,k3_]:=Module[{sym},sym=0;
	If[k1===k2||k2===k3,sym=1];
	If[k1===k2&&k2===k3,sym=5];
	1+sym];
pisym[k1_,k2_,k3_] :=If[k1+k2-k3<10^-3,1/2,1];
GetCovs[triangs_,pks_,biasfix_,Survey_]:=Module[
	{covfast,cov,covfast2,Pback,fback,prefac},
	
	Pback = pks[[1,3]];
	fback = pks[[1,9]];
	covlist2= Bcovlist/.Survey/.f1->fback/.biasfix;
	covfast[y_,y13_,pk1_,pk2_,pk3_]=covlist2;
	
	prefac=((8 \[Pi]^4)/(dk^3 Vs) )^-1/.Survey;
	covfast2[k1_,k2_,k3_]:=prefac ( sb[k1,k2,k3]/( k1 k2 k3 pisym[k1,k2,k3]))^-1 covfast[ys[k1,k2,k3],y13s[k1,k2,k3],Pback[k1],Pback[k2],Pback[k3]];
	cov = EV2[covfast2,triangs];
	Clear[covlist2];
	cov
	];

GetCovs2[triangs_,pks_,biasfix_,Survey_,fs_]:=Module[	
	{covfast,cov,covfast2,Pback,fback,prefac},
	(*Difference is using f as argument instead of extracting from pks*)
	Pback = pks[[1,3]];
	covlist2= Bcovlist/.Survey/.f1->fs/.biasfix;
	covfast[y_,y13_,pk1_,pk2_,pk3_]=covlist2;
	
	prefac=((8 \[Pi]^4)/(dk^3 Vs) )^-1/.Survey;
	covfast2[k1_,k2_,k3_]:=prefac ( sb[k1,k2,k3]/( k1 k2 k3 pisym[k1,k2,k3]))^-1 covfast[ys[k1,k2,k3],y13s[k1,k2,k3],Pback[k1],Pback[k2],Pback[k3]];
	cov = EV2[covfast2,triangs];
	Clear[covlist2];
	cov
	];

CovBmulti[Triangs_,pks_,biasfix_,Survey_,fs_]:=Monitor[Sum[GetCovs2[Triangs,pks,biasfix[[i]],Survey[[i]],fs[[i]]],{i,1,Length[biasfix]}],i]

GetFullzCovB[Triangs_,pks_,spec_,Dg_,biasfix_,dkk_]:=CovBmulti[Triangs,pks,shiftedbfix[spec,Dg,biasfix],fixsurv[spec,dkk],spec[[-1]]]//Expand;

GetFullzCovBNS[Triangs_,pks_,spec_,Dg_,biasfix_,dkk_]:=CovBmulti[Triangs,pks,shiftedbfixNS[spec,Dg,biasfix],fixsurv[spec,dkk],spec[[-1]]]//Expand;

GetFullzCovBS1[Triangs_,pks_,spec_,Dg_,biasfix_,dkk_]:=CovBmulti[Triangs,pks,shiftedbfixS1[spec,Dg,biasfix],fixsurv[spec,dkk],spec[[-1]]]//Expand;

GetFullzCovBn[Triangs_,pks_,spec_,Dg_,biasfix_,dkk_]:=CovBmulti[Triangs,pks,shiftedbfixn[spec,Dg,biasfix],fixsurv[spec,dkk],spec[[-1]]]//Expand;



(* ::Subsubsection:: *)
(*Preliminaries for redshift integrations*)


(*List definitions for loop, tree and covariance *)
(*u = Sqrt[1-\[Mu]1^2]Sin[\[Phi]]*)

Blooplist = {1,u^2,u^4,u^6,u^8,u^10,u^12,u \[Mu]1,u^3 \[Mu]1,
			 u^5 \[Mu]1,u^7 \[Mu]1,u^9 \[Mu]1,u^11 \[Mu]1,\[Mu]1^2,u^2 \[Mu]1^2,
			 u^4 \[Mu]1^2,u^6 \[Mu]1^2,u^8 \[Mu]1^2,u^10 \[Mu]1^2,u \[Mu]1^3,u^3 \[Mu]1^3,
			 u^5 \[Mu]1^3,u^7 \[Mu]1^3,u^9 \[Mu]1^3,\[Mu]1^4,u^2 \[Mu]1^4,u^4 \[Mu]1^4,u^6 \[Mu]1^4,
			 u^8 \[Mu]1^4,u \[Mu]1^5,u^3 \[Mu]1^5,u^5 \[Mu]1^5,u^7 \[Mu]1^5,\[Mu]1^6,u^2 \[Mu]1^6,
			 u^4 \[Mu]1^6,u^6 \[Mu]1^6,u \[Mu]1^7,u^3 \[Mu]1^7,u^5 \[Mu]1^7,\[Mu]1^8,u^2 \[Mu]1^8,
			 u^4 \[Mu]1^8,u \[Mu]1^9,u^3 \[Mu]1^9,\[Mu]1^10,u^2 \[Mu]1^10,u \[Mu]1^11,\[Mu]1^12};
Bcoefstree = {1,u^2,u^4,u \[Mu]1,u^3 \[Mu]1,u^5 \[Mu]1,\[Mu]1^2,u^2 \[Mu]1^2,u^4 \[Mu]1^2,
			  u^6 \[Mu]1^2,u \[Mu]1^3,u^3 \[Mu]1^3,u^5 \[Mu]1^3,\[Mu]1^4,u^2 \[Mu]1^4,u^4 \[Mu]1^4,
			  u \[Mu]1^5,u^3 \[Mu]1^5,\[Mu]1^6,u^2 \[Mu]1^6,u \[Mu]1^7,\[Mu]1^8};
covcoefs = {1,u^2,u^4,u^6,u^8,u^10,u^12,u^14,u^16,u^18,u^20,u^22,u^24,u \[Mu]1,
			u^3 \[Mu]1,u^5 \[Mu]1,u^7 \[Mu]1,u^9 \[Mu]1,u^11 \[Mu]1,u^13 \[Mu]1,u^15 \[Mu]1,u^17 \[Mu]1,
			u^19 \[Mu]1,u^21 \[Mu]1,u^23 \[Mu]1,\[Mu]1^2,u^2 \[Mu]1^2,u^4 \[Mu]1^2,u^6 \[Mu]1^2,u^8 \[Mu]1^2,
			u^10 \[Mu]1^2,u^12 \[Mu]1^2,u^14 \[Mu]1^2,u^16 \[Mu]1^2,u^18 \[Mu]1^2,u^20 \[Mu]1^2,u^22 \[Mu]1^2,
			u^24 \[Mu]1^2,u \[Mu]1^3,u^3 \[Mu]1^3,u^5 \[Mu]1^3,u^7 \[Mu]1^3,u^9 \[Mu]1^3,u^11 \[Mu]1^3,u^13 \[Mu]1^3,
			u^15 \[Mu]1^3,u^17 \[Mu]1^3,u^19 \[Mu]1^3,u^21 \[Mu]1^3,u^23 \[Mu]1^3,\[Mu]1^4,u^2 \[Mu]1^4,u^4 \[Mu]1^4,
			u^6 \[Mu]1^4,u^8 \[Mu]1^4,u^10 \[Mu]1^4,u^12 \[Mu]1^4,u^14 \[Mu]1^4,u^16 \[Mu]1^4,u^18 \[Mu]1^4,u^20 \[Mu]1^4,
			u^22 \[Mu]1^4,u^24 \[Mu]1^4,u \[Mu]1^5,u^3 \[Mu]1^5,u^5 \[Mu]1^5,u^7 \[Mu]1^5,u^9 \[Mu]1^5,u^11 \[Mu]1^5,
			u^13 \[Mu]1^5,u^15 \[Mu]1^5,u^17 \[Mu]1^5,u^19 \[Mu]1^5,u^21 \[Mu]1^5,u^23 \[Mu]1^5,\[Mu]1^6,u^2 \[Mu]1^6,
			u^4 \[Mu]1^6,u^6 \[Mu]1^6,u^8 \[Mu]1^6,u^10 \[Mu]1^6,u^12 \[Mu]1^6,u^14 \[Mu]1^6,u^16 \[Mu]1^6,u^18 \[Mu]1^6,
			u^20 \[Mu]1^6,u^22 \[Mu]1^6,u^24 \[Mu]1^6,u \[Mu]1^7,u^3 \[Mu]1^7,u^5 \[Mu]1^7,u^7 \[Mu]1^7,u^9 \[Mu]1^7,u^11 \[Mu]1^7,
			u^13 \[Mu]1^7,u^15 \[Mu]1^7,u^17 \[Mu]1^7,u^19 \[Mu]1^7,u^21 \[Mu]1^7,u^23 \[Mu]1^7,\[Mu]1^8,u^2 \[Mu]1^8,u^4 \[Mu]1^8,
			u^6 \[Mu]1^8,u^8 \[Mu]1^8,u^10 \[Mu]1^8,u^12 \[Mu]1^8,u^14 \[Mu]1^8,u^16 \[Mu]1^8,u^18 \[Mu]1^8,u^20 \[Mu]1^8,u^22 \[Mu]1^8,
			u^24 \[Mu]1^8,u \[Mu]1^9,u^3 \[Mu]1^9,u^5 \[Mu]1^9,u^7 \[Mu]1^9,u^9 \[Mu]1^9,u^11 \[Mu]1^9,u^13 \[Mu]1^9,u^15 \[Mu]1^9,u^17 \[Mu]1^9,
			u^19 \[Mu]1^9,u^21 \[Mu]1^9,u^23 \[Mu]1^9,\[Mu]1^10,u^2 \[Mu]1^10,u^4 \[Mu]1^10,u^6 \[Mu]1^10,u^8 \[Mu]1^10,u^10 \[Mu]1^10,
			u^12 \[Mu]1^10,u^14 \[Mu]1^10,u^16 \[Mu]1^10,u^18 \[Mu]1^10,u^20 \[Mu]1^10,u^22 \[Mu]1^10,u^24 \[Mu]1^10,u \[Mu]1^11,
			u^3 \[Mu]1^11,u^5 \[Mu]1^11,u^7 \[Mu]1^11,u^9 \[Mu]1^11,u^11 \[Mu]1^11,u^13 \[Mu]1^11,u^15 \[Mu]1^11,u^17 \[Mu]1^11,
			u^19 \[Mu]1^11,u^21 \[Mu]1^11,u^23 \[Mu]1^11,\[Mu]1^12,u^2 \[Mu]1^12,u^4 \[Mu]1^12,u^6 \[Mu]1^12,u^8 \[Mu]1^12,u^10 \[Mu]1^12,
			u^12 \[Mu]1^12,u^14 \[Mu]1^12,u^16 \[Mu]1^12,u^18 \[Mu]1^12,u^20 \[Mu]1^12,u^22 \[Mu]1^12,u^24 \[Mu]1^12,u \[Mu]1^13,u^3 \[Mu]1^13,
			u^5 \[Mu]1^13,u^7 \[Mu]1^13,u^9 \[Mu]1^13,u^11 \[Mu]1^13,u^13 \[Mu]1^13,u^15 \[Mu]1^13,u^17 \[Mu]1^13,u^19 \[Mu]1^13,u^21 \[Mu]1^13,u^23 \[Mu]1^13,
			\[Mu]1^14,u^2 \[Mu]1^14,u^4 \[Mu]1^14,u^6 \[Mu]1^14,u^8 \[Mu]1^14,u^10 \[Mu]1^14,u^12 \[Mu]1^14,u^14 \[Mu]1^14,u^16 \[Mu]1^14,u^18 \[Mu]1^14,u^20 \[Mu]1^14,
			u^22 \[Mu]1^14,u \[Mu]1^15,u^3 \[Mu]1^15,u^5 \[Mu]1^15,u^7 \[Mu]1^15,u^9 \[Mu]1^15,u^11 \[Mu]1^15,u^13 \[Mu]1^15,u^15 \[Mu]1^15,u^17 \[Mu]1^15,u^19 \[Mu]1^15,
			u^21 \[Mu]1^15,\[Mu]1^16,u^2 \[Mu]1^16,u^4 \[Mu]1^16,u^6 \[Mu]1^16,u^8 \[Mu]1^16,u^10 \[Mu]1^16,u^12 \[Mu]1^16,u^14 \[Mu]1^16,u^16 \[Mu]1^16,u^18 \[Mu]1^16,
			u^20 \[Mu]1^16,u \[Mu]1^17,u^3 \[Mu]1^17,u^5 \[Mu]1^17,u^7 \[Mu]1^17,u^9 \[Mu]1^17,u^11 \[Mu]1^17,u^13 \[Mu]1^17,u^15 \[Mu]1^17,u^17 \[Mu]1^17,u^19 \[Mu]1^17,
			\[Mu]1^18,u^2 \[Mu]1^18,u^4 \[Mu]1^18,u^6 \[Mu]1^18,u^8 \[Mu]1^18,u^10 \[Mu]1^18,u^12 \[Mu]1^18,u^14 \[Mu]1^18,u^16 \[Mu]1^18,u^18 \[Mu]1^18,u \[Mu]1^19,u^3 \[Mu]1^19,
			u^5 \[Mu]1^19,u^7 \[Mu]1^19,u^9 \[Mu]1^19,u^11 \[Mu]1^19,u^13 \[Mu]1^19,u^15 \[Mu]1^19,u^17 \[Mu]1^19,\[Mu]1^20,u^2 \[Mu]1^20,u^4 \[Mu]1^20,u^6 \[Mu]1^20,
			u^8 \[Mu]1^20,u^10 \[Mu]1^20,u^12 \[Mu]1^20,u^14 \[Mu]1^20,u^16 \[Mu]1^20,u \[Mu]1^21,u^3 \[Mu]1^21,u^5 \[Mu]1^21,u^7 \[Mu]1^21,u^9 \[Mu]1^21,u^11 \[Mu]1^21,
			u^13 \[Mu]1^21,u^15 \[Mu]1^21,\[Mu]1^22,u^2 \[Mu]1^22,u^4 \[Mu]1^22,u^6 \[Mu]1^22,u^8 \[Mu]1^22,u^10 \[Mu]1^22,u^12 \[Mu]1^22,u^14 \[Mu]1^22,u \[Mu]1^23,u^3 \[Mu]1^23,
			u^5 \[Mu]1^23,u^7 \[Mu]1^23,u^9 \[Mu]1^23,u^11 \[Mu]1^23,u^13 \[Mu]1^23,\[Mu]1^24,u^2 \[Mu]1^24,u^4 \[Mu]1^24,u^6 \[Mu]1^24,u^8 \[Mu]1^24,u^10 \[Mu]1^24,u^12 \[Mu]1^24,
			u \[Mu]1^25,u^3 \[Mu]1^25,u^5 \[Mu]1^25,u^7 \[Mu]1^25,u^9 \[Mu]1^25,u^11 \[Mu]1^25,\[Mu]1^26,u^2 \[Mu]1^26,u^4 \[Mu]1^26,u^6 \[Mu]1^26,u^8 \[Mu]1^26,u^10 \[Mu]1^26,
			u \[Mu]1^27,u^3 \[Mu]1^27,u^5 \[Mu]1^27,u^7 \[Mu]1^27,u^9 \[Mu]1^27,\[Mu]1^28,u^2 \[Mu]1^28,u^4 \[Mu]1^28,u^6 \[Mu]1^28,u^8 \[Mu]1^28,u \[Mu]1^29,u^3 \[Mu]1^29,
			u^5 \[Mu]1^29,u^7 \[Mu]1^29,\[Mu]1^30,u^2 \[Mu]1^30,u^4 \[Mu]1^30,u^6 \[Mu]1^30,u \[Mu]1^31,u^3 \[Mu]1^31,u^5 \[Mu]1^31,\[Mu]1^32,u^2 \[Mu]1^32,u^4 \[Mu]1^32,
			u \[Mu]1^33,u^3 \[Mu]1^33,\[Mu]1^34,u^2 \[Mu]1^34,u \[Mu]1^35,\[Mu]1^36};

(*Angular integrations*)
(*Integration tensors for full \[Mu] dependence*)
MasterCreator[Bcoefs_,CovCoefs_]:=TensorProduct[Bcoefs,Bcoefs,CovCoefs]/.rules\[Mu]1u;

BTreemaster = MasterCreator[Bcoefstree,covcoefs];
BLoopmaster = MasterCreator[Blooplist,covcoefs];

(*Integration for monopole only*)
Blooplistmono = Blooplist/.rules\[Mu]1u;
covcoefsmono=covcoefs/.rules\[Mu]1u;
BLoopmastermono =MasterCreator[Blooplistmono,covcoefs];

(*Integration for monopole quadrupole only*)
quadlist={co[1]-(3 co[3])/35-(2 co[4])/21-co[5]/11-(12 co[6])/143-co[7]/13-co[15]/35-(2 co[16])/105-co[17]/77-(4 co[18])/429-co[19]/143-(3 co[25])/35-(2 co[26])/105-(3 co[27])/385-(4 co[28])/1001-co[29]/429-(2 co[34])/21-co[35]/77-(4 co[36])/1001-(5 co[37])/3003-co[41]/11-(4 co[42])/429-co[43]/429-(12 co[46])/143-co[47]/143-co[49]/13,
		  co[2]+(6 co[3])/7+(5 co[4])/7+(20 co[5])/33+(75 co[6])/143+(6 co[7])/13+co[15]/7+(2 co[16])/21+(5 co[17])/77+(20 co[18])/429+(5 co[19])/143+co[26]/21+(2 co[27])/77+(15 co[28])/1001+(4 co[29])/429+(5 co[35])/231+(10 co[36])/1001+(5 co[37])/1001+(5 co[42])/429+(2 co[43])/429+co[47]/143,
		  0,0,0,0,0,-co[8]-(3 co[9])/7-(5 co[10])/21-(5 co[11])/33-(15 co[12])/143-co[13]/13-(3 co[20])/7-co[21]/7-(5 co[22])/77-(5 co[23])/143-(3 co[24])/143-(5 co[30])/21-(5 co[31])/77-(25 co[32])/1001-(5 co[33])/429-(5 co[38])/33-(5 co[39])/143-(5 co[40])/429-(15 co[44])/143-(3 co[45])/143-co[48]/13,
		  0,0,0,0,0,co[14]+co[15]/7+co[16]/21+(5 co[17])/231+(5 co[18])/429+co[19]/143+(6 co[25])/7+(2 co[26])/21+(2 co[27])/77+(10 co[28])/1001+(2 co[29])/429+(5 co[34])/7+(5 co[35])/77+(15 co[36])/1001+(5 co[37])/1001+(20 co[41])/33+(20 co[42])/429+(4 co[43])/429+(75 co[46])/143+(5 co[47])/143+(6 co[49])/13,
		  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
Transformtoquad[derivs_]:=Quiet[Table[quadlist/.co[a_]->derivs[[i,a]],{i,1,Length[derivs]}]];


(* ::Subsubsection:: *)
(*Fisher matrix*)


(*Fisher with Numerical covariance*)
GetFisherBNum[derivs_,covs_]:=(2048 - Length[covs] - 2)/(2048 - 1) Transpose[derivs.Blooplistmono].Inverse[covs].(derivs.Blooplistmono);

(*Fisher with Analytical covariance*)
GetFisherB[derivs_,invcovs_,inter_]:=Monitor[Sum[derivs[[i]].inter.invcovs[[i]].Transpose[derivs[[i]]],{i,1,Length[invcovs]}],i];
GetFisherlistB[derivs_,invcovs_,inter_]:=Monitor[Table[derivs[[i]].inter.invcovs[[i]].(derivs[[i]]//Transpose),{i,1,Length[invcovs]}],i];

(*Fisher with Analytical covariance, \[OpenCurlyDoubleQuote]Transformtoquad\[CloseCurlyDoubleQuote] quad cuts off all multipoles higher thean the quadrupoles. So intergating with iter=master gives you monopole and quadtupole contributions.
*)
GetFisherBmonoquad[derivs_,invcovs_,inter_]:=Module[{derivsmonoquad},
	derivsmonoquad = Monitor[Table[Transformtoquad[derivs[[i]]],{i,1,Length[derivs]}],i];
	Monitor[Sum[derivsmonoquad[[i]].inter.invcovs[[i]].(derivsmonoquad[[i]]//Transpose),{i,1,Length[invcovs]}],i]
	];


End[]
EndPackage[]

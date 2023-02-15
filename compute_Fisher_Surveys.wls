#!/usr/bin/env wolframscript
(* ::Package:: *)

PackageDirectory = NotebookDirectory[];
AppendTo[$Path,PackageDirectory];
SetDirectory[NotebookDirectory[]]


(* Load FisherFunctions package (takes some time)*)
<<FisherFunctions`


(*tests to check it loaded well*)
Length[Bcovlist] == 319
Tfi[0.56723569] == "0.56724"


(* ::Subsection:: *)
(*BOSS Survey*)


(* ::Subsubsection:: *)
(*BOSS numbers fully binned*)


(*z,nb,V,N,Dg,fs*)
BOSSspec = {{0.05`,0.15`,0.25`,0.35`,0.45`,0.55`,0.65`,0.75`},
			{0.000289`,0.00029`,0.0003`,0.000304`,0.00027600000000000004`,0.000323`,0.00012`,0.00001`},
			{2.55`*^7,1.64`*^8,4.02`*^8,7.04`*^8,1.04`*^9,1.38`*^9,1.72`*^9,2.04`*^9},
			{7369.499999999999`,47560.`,120599.99999999999`,214016.`,287040.00000000006`,445740.`,206400.`,20400.`},
			{0.974195497032948,0.9242943699666863,0.8770024987292571,0.8325495671618781,0.7910134230680845,0.7523626109979884,0.7164923688240663,0.6832529172629023},
			{0.5505691478621921,0.605295980688119,0.6541675497977562,0.697208424265646,0.7347223934446799,0.7671789591953987,0.7951240646380944,0.8191175994924873}};
BOSScut=5;

(*more tests*)
Norm[Getzeffs[BOSSspec,BOSScut]-{0.3572735079306312`,0.5867561780711928`}]<10^-8
Norm[Getnbeff[BOSSspec,BOSScut]-
{0.00029026056499880656`,0.00025120590596841824`}]<10^-8


LOWZspec = BOSSspec[[All,;;BOSScut]];
CMASSspec = BOSSspec[[All,BOSScut+1;;]];
Dg57=0.7449704148232558;
Dg32 =0.8455802559185244;


CMASSpk = PlinsubsT["plins/CMASSpksN/",0.57];
LOWZpk = PlinsubsT["plins/LOWZpksN/",0.32];


(*BOSS numbers fully binned into 2*)
nbfixCMASSngc =nb->(435741+500872)/((5.3+9.0)hback^3 10^9);
nbfixCMASSsgc =nb->(185498+158262)/((2.0+3.3)hback^3 10^9);
nbfixLOWZngc =nb->429182/(4.7 hback^3 10^9);
nbfixLOWZsgc =nb->174819/(1.7hback^3 10^9);
VSCMASSngc =Vs->(5.3+9.0)hback^3 10^9;
VSCMASSsgc=Vs->(2.0+3.3)hback^3 10^9;
VLOWZngc =Vs-> 4.7 hback^3 10^9;
VLOWZsgc =Vs->1.7hback^3 10^9;


(*Parameters bestfit from BOSS*)
bfixBOSS={b1->1.9391,c2->1.1447,c4->-0.2895,Bb3->-0.3740,Bc1->5.4605,Bc2->-1.5439,Bc3->1.3081,Be1->1.6922,Be2->0.9091,ce2->0.5548,Bb4->0.1276,Bb6->-0.3487,Bb7->0.2234,Bb8->-0.2972,Bb9->0.0148,Bb10->0.0428,Bb11->0.0356,Bc4->-0.4785,Bc5->0.1054,Bc6->0.8713,Bc7->-0.4555,Bc8->0.4442,Bc9->-0.4152,Bc10->-0.6500,Bc11->-0.0881,Bc12->-0.3731,Bc13->-0.1642,Bc14->-0.1958,Bd1->5.4000,Bd2->-0.7210,Bd3->-0.4266,Be3->0.0737,Be4->-0.1381,Be5->6.0673,Be6->-0.0932,Be7->-0.9702,Be8->0.2617,Be9->0.2638,Be10->-0.1526,Be11->0.4270,Be12->-0.4322,Bb1->1.9391,Bb2->0.6047,Bb5->1.0142,b2->0.6047,b3->-4.8320,b4->1.0142,cct->-5.4605,cr1->-1.0607,cr2->-0.3996,ce0->1.6922,ce1->1.1865,b12->0,b13->0,b14->0,b15->0,Tst->9.6}/.normbias;
knlfixBOSS = {knl->0.7};
bfixBOSSNS = Join[{Be1->0},bfixBOSS];
coeffixBOSS = Join[bfixBOSS,fnlfix,paramfix,knlfixBOSS];
coeffixBOSSNS = Join[{Be1->0},coeffixBOSS];

(* Tab imports *)
ktabCMASS = Import["covsks/ktabCMASS.m"];
ktabLOWZ = Import["covsks/ktabLOWZ.m"];


normbias


(*Binning sizes*)
dkP = dk->0.01;
dkBCMASS=dk->0.02154;
dkBLOWZ=dk->0.02459;


(*Analytical Covariance with Full redshift*)
COVCMASSP=GetFullzCovP[ktabCMASS,CMASSpk,CMASSspec,Dg57,bfixBOSS,dkP];
COVLOWZP=GetFullzCovP[ktabLOWZ,LOWZpk,LOWZspec,Dg32,bfixBOSS,dkP];


(*check*)
Variables[COVCMASSP]=={FisherFunctions`Private`\[Mu]1}
Variables[COVLOWZP]=={FisherFunctions`Private`\[Mu]1}
Dimensions[COVCMASSP]== {21}
Dimensions[COVLOWZP]== {21}


(*Analytical Covariance with Full redshift - no shot noise*)
COVCMASSPNS=GetFullzCovPNS[ktabCMASS,CMASSpk,CMASSspec,Dg57,bfixBOSSNS,dkP];
COVLOWZPNS=GetFullzCovPNS[ktabLOWZ,LOWZpk,LOWZspec,Dg32,bfixBOSSNS,dkP];


(*check*)
Variables[COVCMASSPNS]=={FisherFunctions`Private`\[Mu]1}
Variables[COVLOWZPNS]=={FisherFunctions`Private`\[Mu]1}
Dimensions[COVCMASSPNS]== {21}
Dimensions[COVLOWZPNS]== {21}


(* Power-spectrum derivatives *)
PCMASSngcderivs = Monitor[Table[P1LoopDerivs[ktabCMASS[[i]],CMASSpk],{i,1,Length[ktabCMASS]}]/.coeffixBOSS/.nbfixCMASSngc,i];
PCMASSsgcderivs = Monitor[Table[P1LoopDerivs[ktabCMASS[[i]],CMASSpk],{i,1,Length[ktabCMASS]}]/.coeffixBOSS/.nbfixCMASSsgc,i];
PLOWZngcderivs =Monitor[ Table[P1LoopDerivs[ktabLOWZ[[i]],LOWZpk],{i,1,Length[ktabLOWZ]}]/.coeffixBOSS/.nbfixLOWZngc,i];
PLOWZsgcderivs =Monitor[ Table[P1LoopDerivs[ktabLOWZ[[i]],LOWZpk],{i,1,Length[ktabLOWZ]}]/.coeffixBOSS/.nbfixLOWZsgc,i];


(*Replicate full P result*)
(*reduce derivatives to monopole and quadrupole for the powerspectrum*)
PCMASSngcderivsmonoquad = Join[mono[PCMASSngcderivs],quad[PCMASSngcderivs]];
PCMASSsgcderivsmonoquad = Join[mono[PCMASSsgcderivs],quad[PCMASSsgcderivs]];
PLOWZngcderivsmonoquad = Join[mono[PLOWZngcderivs],quad[PLOWZngcderivs]];
PLOWZsgcderivsmonoquad = Join[mono[PLOWZsgcderivs],quad[PLOWZsgcderivs]];


(*Covariances - ERROR*)
(*covPCMASSngc=Import["covsks/covPCMASSngc.mx"];
covPCMASSsgc=Import["covsks/covPCMASSsgc.mx"];
covPLOWZngc=Import["covsks/covPLOWZngc.mx"];
covPLOWZsgc=Import["covsks/covPLOWZsgc.mx"];*)

(*Fisher - numerical*)
CMASSNGCFisherFullP =(2048 - Length[covPCMASSngc] - 2)/(2048 - 1) (PCMASSngcderivsmonoquad//Transpose).Inverse[covPCMASSngc].PCMASSngcderivsmonoquad//Quiet;
CMASSSGCFisherFullP =(2048 - Length[covPCMASSsgc] - 2)/(2048 - 1) (PCMASSsgcderivsmonoquad//Transpose).Inverse[covPCMASSsgc].PCMASSsgcderivsmonoquad//Quiet;
LOWZNGCFisherFullP =(2048 - Length[covPLOWZngc] - 2)/(2048 - 1) (PLOWZngcderivsmonoquad//Transpose).Inverse[covPLOWZngc].PLOWZngcderivsmonoquad//Quiet;
LOWZSGCFisherFullP =(2048 - Length[covPLOWZsgc] - 2)/(2048 - 1) (PLOWZsgcderivsmonoquad//Transpose).Inverse[covPLOWZsgc].PLOWZsgcderivsmonoquad//Quiet;


(*Comparisons part 2: Full redshift dependence analytical covariance*)
(*Reduce P-derivatives to monopole and quadrupole*)
PCMASSderivsmonoquad= removelg2[PCMASSngcderivs]//Expand;
PLOWZderivsmonoquad = removelg2[PLOWZngcderivs]//Expand;
(*Powerspectrum Fishers*)
FPCnew = GetFisherPnew[PCMASSderivsmonoquad,COVCMASSP];
FPLnew = GetFisherPnew[PLOWZderivsmonoquad,COVLOWZP];
(*Bispectrum Fishers*)
(*FBCnew = GetFisherB[BderivsCMASSkeffCMASSngc,COVCMASSB,BLoopmastermono];
FBLnew = GetFisherB[BderivsLOWZkeffLOWZngc,COVLOWZB,BLoopmastermono];*)



(*Independent prediction: Full redshift dependence analytical covariance*)
(*Full angular dependence (all multipoles)*)

(*Powerspectrum Fishers*)
FPCnewmu = GetFisherPnew[PCMASSngcderivs,COVCMASSP];
FPLnewmu = GetFisherPnew[PLOWZsgcderivs,COVLOWZP];

(*Bispectrum Fishers*)
(*FBCnewmu = GetFisherB[BderivsCMASSkeffCMASSngc,COVCMASSB,Loopmaster];
FBLnewmu = GetFisherB[BderivsLOWZkeffLOWZngc,COVLOWZB,Loopmaster];*)

(*Full angular dependence 2Loop w/o shot noise *)
(*Powerspectrum Fishers*)
FPCnewmubase = GetFisherPnew[PCMASSngcderivs,COVCMASSP];
FPLnewmubase = GetFisherPnew[PLOWZsgcderivs,COVLOWZP];
FPCnewmubaseNS = GetFisherPnew[PCMASSngcderivs,COVCMASSPNS];
FPLnewmubaseNS = GetFisherPnew[PLOWZsgcderivs,COVLOWZPNS];
(*FPCnewmu2L = GetFisherPnew[PCMASSngcderivs2L,COVCMASSP2L];
FPLnewmu2L = GetFisherPnew[PLOWZngcderivs2L,COVLOWZP2L];
FPCnewmu2LNS = GetFisherPnew[PCMASSngcderivs2L,COVCMASSPNS2L];
FPLnewmu2LNS = GetFisherPnew[PLOWZngcderivs2L,COVLOWZPNS2L];*)

(*no shot noise full mu*)
(*Powerspectrum Fishers*)
FPCnewNS = GetFisherPnew[PCMASSngcderivs,COVCMASSPNS];
FPLnewNS = GetFisherPnew[PLOWZsgcderivs,COVLOWZPNS];

(*Bispectrum Fishers*)
(*FBCnewNS = GetFisherB[BderivsCMASSkeffCMASSngc,COVCMASSBNS,Loopmaster];
FBLnewNS = GetFisherB[BderivsLOWZkeffLOWZngc,COVLOWZBNS,Loopmaster];*)


(* ::Subsection:: *)
(*Exports*)


Allparams = DeleteCases[AllParams,b12|b13|b14|b15];
GetFExp[Fisher_]:=SubFisher[Fisher+Priors,Allparams]
GetFExp4[Fisher_]:=SubFisher4sky[Fisher+Prior4sky,Allparams]
Exporter[Fisher_,name_]:=Export[name,GetFExp[Fisher]]
Exporter4[Fisher_,name_]:=Export[name,GetFExp4[Fisher]]


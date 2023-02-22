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


(*z,nb,V,N,Dg,fs,b1*)
BOSSspec = {{0.05`,0.15`,0.25`,0.35`,0.45`,0.55`,0.65`,0.75`},
			{0.000289`,0.00029`,0.0003`,0.000304`,0.00027600000000000004`,0.000323`,0.00012`,0.00001`},
			{2.55`*^7,1.64`*^8,4.02`*^8,7.04`*^8,1.04`*^9,1.38`*^9,1.72`*^9,2.04`*^9},
			{7369.499999999999`,47560.`,120599.99999999999`,214016.`,287040.00000000006`,445740.`,206400.`,20400.`},
			{0.973938713655033,0.9235978369817912,0.8759572952224975,0.8312353500171761,0.789497353926307,0.7506997404081879,0.7147268398976079,0.6814195172840236},
			{0.5559799005938011,0.6106395138125067,0.6593260133608193,0.7021004757219186,0.7392994196584667,0.7714187036372595,0.799023714461632,0.8226875527129273},
			{1.4798689114411347`,1.560529449372745`,1.6454017014847226`,1.7339272493130506`,1.8255939894151774`,1.9199442152509005`,2.0165768843849974`,2.115145791144549`}};
			
BOSScut=5;

(*more tests*)
Norm[Getzeffs[BOSSspec,BOSScut]-{0.3572735079306312`,0.5867561780711928`}]<10^-8
Norm[Getnbeff[BOSSspec,BOSScut]-
{0.00029026056499880656`,0.00025120590596841824`}]<10^-8


LOWZspec = BOSSspec[[All,;;BOSScut]];
CMASSspec = BOSSspec[[All,BOSScut+1;;]];
Dg57=0.7449704148232558;
Dg32 =0.8455802559185244;


(*CMASSpk = PlinsubsT["plins/CMASSpksNk/",0.57];
LOWZpk = PlinsubsT["plins/LOWZpksNk/",0.32];
Export["aux_vars/CMASSpk.m",CMASSpk]
Export["aux_vars/LOWZpk.m",LOWZpk]*)

CMASSpk=Import["aux_vars/CMASSpk.m"];
LOWpk=Import["aux_vars/LOWZpk.m"];


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


(*Binning sizes*)
dkP = dk->0.01;
dkBCMASS=dk->0.02154;
dkBLOWZ=dk->0.02459;


(* ::Subsubsection:: *)
(*Covariance*)


(*Analytical Covariance with Full redshift*)
COVCMASSP=GetFullzCovP[ktabCMASS,CMASSpk,CMASSspec,bfixBOSS,dkP];
COVLOWZP=GetFullzCovP[ktabLOWZ,LOWZpk,LOWZspec,bfixBOSS,dkP];


(*check*)
Variables[COVCMASSP]=={FisherFunctions`Private`\[Mu]1}
Variables[COVLOWZP]=={FisherFunctions`Private`\[Mu]1}
Dimensions[COVCMASSP]== {21}
Dimensions[COVLOWZP]== {21}


(*Analytical Covariance with Full redshift - no shot noise*)
COVCMASSPNS=GetFullzCovP[ktabCMASS,CMASSpk,CMASSspec,bfixBOSSNS,dkP];
COVLOWZPNS=GetFullzCovP[ktabLOWZ,LOWZpk,LOWZspec,bfixBOSSNS,dkP];


(*check*)
Variables[COVCMASSPNS]=={FisherFunctions`Private`\[Mu]1}
Variables[COVLOWZPNS]=={FisherFunctions`Private`\[Mu]1}
Dimensions[COVCMASSPNS]== {21}
Dimensions[COVLOWZPNS]== {21}


(* ::Subsubsection:: *)
(*Derivatives*)


(* Power-spectrum derivatives *)
(*CMASS*)
(*PCMASSderivs = Monitor[Table[P1LoopDerivs[ktabCMASS[[i]],CMASSpk],{i,1,Length[ktabCMASS]}]/.coeffixBOSS,i];
PCMASSngcderivs = PCMASSderivs/.nbfixCMASSngc;
PCMASSsgcderivs = PCMASSderivs/.nbfixCMASSsgc;*)
(*LOWZ*)
(*PLOWZderivs = Monitor[Table[P1LoopDerivs[ktabLOWZ[[i]],LOWZpk],{i,1,Length[ktabLOWZ]}]/.coeffixBOSS,i];
PLOWZngcderivs =PLOWZderivs/.nbfixLOWZngc;
PLOWZsgcderivs =PLOWZderivs/.nbfixLOWZsgc;*)

(*save derivatives*)
(*Export["aux_vars/PCMASSngcderivs.m",PCMASSngcderivs]
Export["aux_vars/PCMASSsgcderivs.m",PCMASSsgcderivs]
Export["aux_vars/PLOWZngcderivs.m",PLOWZngcderivs]
Export["aux_vars/PLOWZsgcderivs.m",PLOWZsgcderivs]*)

(*Import derivs*)
PCMASSngcderivs=Import["aux_vars/PCMASSngcderivs.m"];
PCMASSsgcderivs=Import["aux_vars/PCMASSsgcderivs.m"];
PLOWZngcderivs=Import["aux_vars/PLOWZngcderivs.m"];
PLOWZsgcderivs=Import["aux_vars/PLOWZsgcderivs.m"];


(*check it only depends on \[Mu]1*)
Variables[{PCMASSngcderivs,PCMASSsgcderivs,PLOWZngcderivs,PLOWZsgcderivs}]=={FisherFunctions`Private`\[Mu]1}


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


(* ::Subsubsection:: *)
(*Fisher power spectrum*)


(*Comparisons part 2: Full redshift dependence analytical covariance*)
(*Reduce P-derivatives to monopole and quadrupole*)
PCMASSderivsmonoquad= removelg2[PCMASSngcderivs]//Expand;
PLOWZderivsmonoquad = removelg2[PLOWZngcderivs]//Expand;
(*Powerspectrum Fishers*)
FPCMASS = GetFisherP[PCMASSderivsmonoquad,COVCMASSP];
FPLOWZ = GetFisherP[PLOWZderivsmonoquad,COVLOWZP];
(*Bispectrum Fishers*)
(*FBCnew = GetFisherB[BderivsCMASSkeffCMASSngc,COVCMASSB,BLoopmastermono];
FBLnew = GetFisherB[BderivsLOWZkeffLOWZngc,COVLOWZB,BLoopmastermono];*)



(*Independent prediction: Full redshift dependence analytical covariance*)
(*Full angular dependence (all multipoles)*)

(*Powerspectrum Fishers*)
FPCMASSfull = GetFisherP[PCMASSngcderivs,COVCMASSP];
FPLOWZfull = GetFisherP[PLOWZsgcderivs,COVLOWZP];

(*Bispectrum Fishers*)
(*FBCnewmu = GetFisherB[BderivsCMASSkeffCMASSngc,COVCMASSB,Loopmaster];
FBLnewmu = GetFisherB[BderivsLOWZkeffLOWZngc,COVLOWZB,Loopmaster];*)

(*no shot noise full mu*)
(*Powerspectrum Fishers*)
FPCMASSNS = GetFisherP[PCMASSngcderivs,COVCMASSPNS];
FPLOWZNS = GetFisherP[PLOWZsgcderivs,COVLOWZPNS];

(*Bispectrum Fishers*)
(*FBCnewNS = GetFisherB[BderivsCMASSkeffCMASSngc,COVCMASSBNS,Loopmaster];
FBLnewNS = GetFisherB[BderivsLOWZkeffLOWZngc,COVLOWZBNS,Loopmaster];*)


(*Check that no shot noise is better*)
FPLOWZ[[1,1]]<FPLOWZNS[[1,1]]


(* ::Subsubsection:: *)
(*Bi-spectrum derivatives*)


(*Bispectrum Fishers*)
TriCMASSkeff = Import["covsks/TriCMASSkeff.mtx"];
TriLOWZkeff = Import["covsks/TriLOWZkeff.mtx"];


(*export triangles to calculate jmats*)
(*TriBOSSkeff = Sort[Join[TriCMASSkeff,TriLOWZkeff]];
Export[ctabpath<>"BOSS_tri_eff.csv",TriBOSSkeff]*)



CMASStriangs=Import["covsks/CMASStriangs.mtx"];
LOWZtriangs=Import["covsks/LOWZtriangs.mtx"];


(* ::Text:: *)
(*Create array of derivatives for Fisher*)


(*Calculate derivatives, one at a time instead of table*)
BderivsCMASS=ConstantArray[0,Length[TriCMASSkeff]];
BderivsLOWZ=ConstantArray[0,Length[TriLOWZkeff]];


Monitor[Do[
			If[(BderivsCMASS[[j]]==0)||(Length[BderivsCMASS[[j]]]==0),
				BderivsCMASS[[j]]=B1LoopDerivs[TriCMASSkeff[[j]],CMASSpk]/.Dispatch[coeffixBOSS]];,
	{j,1,Length[TriCMASSkeff]}],j];//AbsoluteTiming
(*{791.6264653`,Null}*)


Monitor[Do[
	If[BderivsLOWZ[[j]]==0,
	BderivsLOWZ[[j]]=B1LoopDerivs[TriLOWZkeff[[j]],LOWZpk]/.Dispatch[coeffixBOSS]];,
	{j,1,Length[TriLOWZkeff]}],
	j];//AbsoluteTiming
(*{400.0957859`,Null}*)


Dimensions[BderivsCMASS]=={150,53,49}
Dimensions[BderivsLOWZ]=={62,53,49}


(* ::Text:: *)
(*Do replacements for NGC and SGC*)


BderivsCMASSbfix=BderivsCMASS/.Dispatch[coeffixBOSS];


BderivsCMASSngc=BderivsCMASSbfix/.Dispatch[nbfixCMASSngc];
BderivsCMASSsgc=BderivsCMASSbfix/.Dispatch[nbfixCMASSsgc];


BderivsLOWZbfix=BderivsLOWZ/.Dispatch[coeffixBOSS];


BderivsLOWZngc=BderivsLOWZbfix/.Dispatch[nbfixCMASSngc];
BderivsLOWZsgc=BderivsLOWZbfix/.Dispatch[nbfixCMASSsgc];


(* ::Subsection:: *)
(*Get kmax for LOWZ (this is just for validation of the formulas. We use the same kmax as in the BOSS analyses)*)


CC=1/(2\[Pi]^2) NIntegrate[k^2 CMASSpk[[1,3]][k],{k,10^-5,0.7}];
Getknl[pk_]:=Quiet[FindRoot[1/(2\[Pi]^2) NIntegrate[k^2 pk[k],{k,10^-5,kn}]-CC,{kn,0.5}]]
Get\[Sigma][ktab_, Cov_]:=Transpose@{ktab,1/(Cov//mono)}//Interpolation;
\[Sigma]CMASS = Get\[Sigma][ktabCMASS, COVCMASSP];
\[Sigma]CMASSNS = Get\[Sigma][ktabCMASS, COVCMASSPNS];
GetkmaxP[\[Sigma]1_,pk1_]:=Quiet[FindRoot[(Sqrt[\[Sigma]1[kmm]/\[Sigma]CMASS[0.22]])-(pk1[[1,3]][0.5]/CMASSpk[[1,3]][0.5])^2 pk1[[1,3]][kmm]/CMASSpk[[1,3]][0.22] (kmm/0.22)^2.4`,{kmm,0.1}]]
GetkmaxNSP[\[Sigma]1_,pk1_]:=Quiet[FindRoot[Sqrt[\[Sigma]1[kmm]/\[Sigma]CMASSNS[0.22]]-(pk1[[1,3]][0.5]/CMASSpk[[1,3]][0.5])^2 pk1[[1,3]][kmm]/CMASSpk[[1,3]][0.22] (kmm/0.22)^2.4`,{kmm,0.1}]]
GetkmaxNSP[Get\[Sigma][ktabLOWZ, COVLOWZPNS],LOWZpk]
{kmm->0.19391015781837476`}
GetkmaxNSP[\[Sigma]1_,pk1_]:=Quiet[FindRoot[Sqrt[\[Sigma]1[kmm]/\[Sigma]CMASSNS[0.22]]-pk1[[1,3]][kmm]/CMASSpk[[1,3]][0.22] (kmm/Getknl[pk1[[1,3]]][[1,2]])^((3+GetSlope[0.8kmm,1.2kmm,CMASSpk[[1,3]]])2) (0.22/0.7)^(-(3+(-1.8))2),{kmm,0.1}]]
GetkmaxNSP[Get\[Sigma][ktabLOWZ, COVLOWZPNS],LOWZpk]
{kmm->0.18513243366801113`}
GetkmaxNSP2L[\[Sigma]1_,pk1_]:=Quiet[FindRoot[Sqrt[\[Sigma]1[kmm]/\[Sigma]CMASSNS[0.22]]-pk1[[1,3]][kmm]/CMASSpk[[1,3]][0.22] (kmm/Getknl[pk1[[1,3]]][[1,2]])^((3+GetSlope[0.8kmm,1.2kmm,CMASSpk[[1,3]]])3) (0.22/0.7)^(-(3+(-1.8))2),{kmm,0.1}]]
COVLOWZPgenNS=GetFullzCovP[ktabraw2,LOWZpk,LOWZspec,bfixstochBOSSNS,dkP];
COVCMASS2PgenNS=GetFullzCovP[ktabraw2,CMASSpk,CMASSspec,bfixstochBOSSNS,dkP];



(* ::Subsection:: *)
(*Exports*)


Allparams = DeleteCases[AllParams,b12|b13|b14|b15];
GetFExp[Fisher_]:=SubFisher[Fisher+Priors,Allparams]
GetFExp4[Fisher_]:=SubFisher4sky[Fisher+Prior4sky,Allparams]
Exporter[Fisher_,name_]:=Export[name,GetFExp[Fisher]]
Exporter4[Fisher_,name_]:=Export[name,GetFExp4[Fisher]]



CombineMs

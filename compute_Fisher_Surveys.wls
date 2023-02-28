#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Subsection:: *)
(*Import FisherFunctions package*)


PackageDirectory = NotebookDirectory[];
AppendTo[$Path,PackageDirectory];
SetDirectory[NotebookDirectory[]]

(* Load FisherFunctions package (takes 1min)*)
<<FisherFunctions`


(*tests to check it loaded well*)
Length[Bcovlist] == 319(*Yaniv's covariance*)
(*Length[Bcovlist] == 88*)(*Diogo's covariance*)
Tfi[0.56723569] == "0.56724"


(* ::Subsection:: *)
(*Define paths*)


paths = Import["paths.json"]
ctabpath = "CTAB_PATH"/.paths;
FisherPath = "FISHER_PATH"/.paths;
covsPath = "COVS_KS_PATH"/.paths
auxPath = "AUX_PATH"/.paths


(* ::Subsection:: *)
(*Make k's and save them to compute jmats*)


(*Bispectrum Fishers*)
TriCMASSkeff = Import[covsPath<>"TriCMASSkeff.mtx"];
TriLOWZkeff = Import[covsPath<>"TriLOWZkeff.mtx"];

keffbase = Import[covsPath<>"baseTrikeff.mtx"];

ktabraw = Table[0.001`+i 0.005,{i,0,300}];
Basetri=makeTri[0.02,0.91,0.02];
Length[Basetri]==Length[keffbase]


(*export triangles to calculate jmats*)
(*TriBOSSkeff = Sort[Join[TriCMASSkeff,TriLOWZkeff]];
Export[ctabpath<>"BOSS_tri_eff.csv",TriBOSSkeff]*)


(*Export[ctabpath<>"base_tri_eff.csv",keffbase]*)


(* ::Subsection::Closed:: *)
(*BOSS Survey*)


(* ::Subsubsection::Closed:: *)
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

LOWZspec = BOSSspec[[All,;;BOSScut]];
CMASSspec = BOSSspec[[All,BOSScut+1;;]];

(*more tests*)
Norm[Getzeffs[BOSSspec,BOSScut]-{0.3572735079306312`,0.5867561780711928`}]<10^-8
Norm[Getnbeff[BOSSspec,BOSScut]-
{0.00029026056499880656`,0.00025120590596841824`}]<10^-8


(*CMASSpk = PlinsubsT["plins/CMASSpksNk/",0.57];
LOWZpk = PlinsubsT["plins/LOWZpksNk/",0.32];
Export["aux_vars/CMASSpk.m",CMASSpk]
Export["aux_vars/LOWZpk.m",LOWZpk]*)

CMASSpk=Import["aux_vars/CMASSpk.m"];
LOWZpk=Import["aux_vars/LOWZpk.m"];


(*BOSS numbers fully binned into 2*)
nbfixCMASSngc =nb->(435741+500872)/((5.3+9.0)hback^3 10^9);
nbfixCMASSsgc =nb->(185498+158262)/((2.0+3.3)hback^3 10^9);
nbfixLOWZngc =nb->429182/(4.7 hback^3 10^9);
nbfixLOWZsgc =nb->174819/(1.7hback^3 10^9);


(*Parameters bestfit from BOSS*)
bfixBOSS={b1->1.9391,c2->1.1447,c4->-0.2895,Bb3->-0.3740,Bc1->5.4605,Bc2->-1.5439,Bc3->1.3081,
		  Be1->1.6922,Be2->0.9091,ce2->0.5548,Bb4->0.1276,Bb6->-0.3487,Bb7->0.2234,Bb8->-0.2972,
		  Bb9->0.0148,Bb10->0.0428,Bb11->0.0356,Bc4->-0.4785,Bc5->0.1054,Bc6->0.8713,Bc7->-0.4555,
		  Bc8->0.4442,Bc9->-0.4152,Bc10->-0.6500,Bc11->-0.0881,Bc12->-0.3731,Bc13->-0.1642,Bc14->-0.1958,
		  Bd1->5.4000,Bd2->-0.7210,Bd3->-0.4266,Be3->0.0737,Be4->-0.1381,Be5->6.0673,Be6->-0.0932,Be7->-0.9702,
		  Be8->0.2617,Be9->0.2638,Be10->-0.1526,Be11->0.4270,Be12->-0.4322,Bb1->1.9391,Bb2->0.6047,Bb5->1.0142,
		  b3->-4.8320,b4->1.0142,cct->-5.4605,cr1->-1.0607,cr2->-0.3996,ce0->1.6922,ce1->1.1865,b12->0,
		  b13->0,b14->0,b15->0,Tst->9.6}/.normbias;

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


(* ::Subsubsection::Closed:: *)
(*Power-spectrum Derivatives*)


(* Power-spectrum derivatives *)
(*CMASS*)
(*PCMASSderivs = Monitor[Table[P1LoopDerivs[ktabCMASS[[i]],CMASSpk]/.coeffixBOSS,{i,1,Length[ktabCMASS]}],i];
PCMASSngcderivs = PCMASSderivs/.nbfixCMASSngc;
PCMASSsgcderivs = PCMASSderivs/.nbfixCMASSsgc;*)
(*LOWZ*)
(*PLOWZderivs = Monitor[Table[P1LoopDerivs[ktabLOWZ[[i]],LOWZpk]/.coeffixBOSS,{i,1,Length[ktabLOWZ]}],i];
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


(*checks*)
Dimensions[PCMASSsgcderivs]=={21,56}
Variables[PCMASSngcderivs]=={\[Mu]1}

Dimensions[PLOWZsgcderivs]=={21,56}
Variables[PLOWZngcderivs]=={\[Mu]1}


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
(*CMASSNGCFisherFullP =(2048 - Length[covPCMASSngc] - 2)/(2048 - 1) (PCMASSngcderivsmonoquad//Transpose).Inverse[covPCMASSngc].PCMASSngcderivsmonoquad//Quiet;
CMASSSGCFisherFullP =(2048 - Length[covPCMASSsgc] - 2)/(2048 - 1) (PCMASSsgcderivsmonoquad//Transpose).Inverse[covPCMASSsgc].PCMASSsgcderivsmonoquad//Quiet;
LOWZNGCFisherFullP =(2048 - Length[covPLOWZngc] - 2)/(2048 - 1) (PLOWZngcderivsmonoquad//Transpose).Inverse[covPLOWZngc].PLOWZngcderivsmonoquad//Quiet;
LOWZSGCFisherFullP =(2048 - Length[covPLOWZsgc] - 2)/(2048 - 1) (PLOWZsgcderivsmonoquad//Transpose).Inverse[covPLOWZsgc].PLOWZsgcderivsmonoquad//Quiet;*)


(* ::Subsubsection::Closed:: *)
(*Power spectrum Covariance*)


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


(* ::Subsubsection::Closed:: *)
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


BderivsCMASS=Import["aux_vars/BderivsCMASS1.m"];
BderivsLOWZ=Import["aux_vars/BderivsLOWZ1.m"];


dispatcher = Dispatch[coeffixBOSS];
Monitor[Do[If[(BderivsCMASS[[j]]==0)||(Length[BderivsCMASS[[j]]]==0),
				BderivsCMASS[[j]]=B1LoopDerivs[TriCMASSkeff[[j]],CMASSpk]/.dispatcher];,
	{j,1,Length[TriCMASSkeff]}],j];//AbsoluteTiming
(*{791.6264653`,Null}*)


Export["aux_vars/BderivsCMASS1.m",BderivsCMASS]


dispatcher = Dispatch[coeffixBOSS];
Monitor[Do[
	If[BderivsLOWZ[[j]]==0,
	BderivsLOWZ[[j]]=B1LoopDerivs[TriLOWZkeff[[j]],LOWZpk]/.dispatcher];,
	{j,1,Length[TriLOWZkeff]}],
	j];//AbsoluteTiming
(*{400.0957859`,Null}*)


Export["aux_vars/BderivsLOWZ1.m",BderivsLOWZ]


(*checks*)
Dimensions[BderivsCMASS]=={150,56,49}
Dimensions[BderivsLOWZ]=={62,56,49}


(* ::Text:: *)
(*Do replacements for NGC and SGC*)


BderivsCMASSngc=BderivsCMASS/.Dispatch[nbfixCMASSngc];
BderivsCMASSsgc=BderivsCMASS/.Dispatch[nbfixCMASSsgc];


BderivsLOWZngc=BderivsLOWZ/.Dispatch[nbfixCMASSngc];
BderivsLOWZsgc=BderivsLOWZ/.Dispatch[nbfixCMASSsgc];


(* ::Subsubsection::Closed:: *)
(*Bi-spectrum covariances*)


(*COVCMASSB=GetFullzCovB[CMASStriangs,CMASSpk,CMASSspec,bfixBOSS,dkBCMASS];//AbsoluteTiming
(*{189.2065861`,Null}*)
COVLOWZB=GetFullzCovB[LOWZtriangs,LOWZpk,LOWZspec,bfixBOSS,dkBLOWZ];//AbsoluteTiming
(*{145.4908871`,Null}*)
Export["aux_vars/COVCMASSB.m",COVCMASSB]
Export["aux_vars/COVLOWZB.m",COVLOWZB]*)


COVCMASSB = Import["aux_vars/COVCMASSB.m"]
COVLOWZB = Import["aux_vars/COVLOWZB.m"]


(*COVCMASSBNS=GetFullzCovB[CMASStriangs,CMASSpk,CMASSspec,bfixBOSSNS,dkBCMASS];//AbsoluteTiming
(*{189.2065861`,Null}*)
COVLOWZBNS=GetFullzCovB[LOWZtriangs,LOWZpk,LOWZspec,bfixBOSSNS,dkBLOWZ];//AbsoluteTiming
(*{145.4908871`,Null}*)
Export["aux_vars/COVCMASSBNS.m",COVCMASSBNS]
Export["aux_vars/COVLOWZBNS.m",COVLOWZBNS]*)


COVCMASSBNS = Import["aux_vars/COVCMASSBNS.m"];
COVLOWZBNS = Import["aux_vars/COVLOWZBNS.m"];


(* ::Subsubsection::Closed:: *)
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



(* ::Subsubsection::Closed:: *)
(*Fishers*)


(*Comparisons part 2: Full redshift dependence analytical covariance*)
(*Reduce P-derivatives to monopole and quadrupole*)
PCMASSderivsmonoquad= removelg2[PCMASSngcderivs]//Expand;
PLOWZderivsmonoquad = removelg2[PLOWZngcderivs]//Expand;
(*Powerspectrum Fishers*)
FPCMASS = GetFisherP[PCMASSderivsmonoquad,COVCMASSP];
FPLOWZ = GetFisherP[PLOWZderivsmonoquad,COVLOWZP];


(*Bispectrum Fishers*)
FBCMASS = GetFisherB[BderivsCMASSngc,COVCMASSB,BLoopmastermono];
FBLOWZ = GetFisherB[BderivsLOWZngc,COVLOWZB,BLoopmastermono];


(*Independent prediction: Full redshift dependence analytical covariance*)
(*Full angular dependence (all multipoles)*)

(*Powerspectrum Fishers*)
FPCMASSfull = GetFisherP[PCMASSngcderivs,COVCMASSP];
FPLOWZfull = GetFisherP[PLOWZsgcderivs,COVLOWZP];


(*Bispectrum Fishers*)
FBCMASSfull = GetFisherB[BderivsCMASSngc,COVCMASSB,BLoopmaster];
FBLOWZfull = GetFisherB[BderivsLOWZngc,COVLOWZB,BLoopmaster];


(*no shot noise full mu*)
(*Powerspectrum Fishers*)
FPCMASSNS = GetFisherP[PCMASSngcderivs,COVCMASSPNS];
FPLOWZNS = GetFisherP[PLOWZsgcderivs,COVLOWZPNS];


(*Bispectrum Fishers*)
FBCMASSNS = GetFisherB[BderivsCMASSngc,COVCMASSBNS,BLoopmaster];
FBLOWZNS = GetFisherB[BderivsLOWZngc,COVLOWZBNS,BLoopmaster];


(*Check that no shot noise is better*)
FPLOWZ[[1,1]]<FPLOWZNS[[1,1]]
FPCMASS[[1,1]]<FPCMASSNS[[1,1]]
FBLOWZ[[1,1]]<FBLOWZNS[[1,1]]
FBCMASS[[1,1]]<FBCMASSNS[[1,1]]


(* ::Subsection::Closed:: *)
(*DESI*)


(* ::Subsubsection::Closed:: *)
(*DESI specifications and derivative computations*)


Getknl[pk_]:=Quiet[FindRoot[1/(2\[Pi]^2) NIntegrate[k^2 pk[k],{k,10^-5,kn}]-1.1525490955287694`,{kn,0.5}]]
CC=1/(2\[Pi]^2) NIntegrate[k^2 CMASSpk[[1,3]][k],{k,10^-5,0.7}];


(*biases*)
bfixBOSSnonstoch={b1->1.9391`,c2->1.1447`,b3->-0.374`,b4->0.1276`,c4->-0.2895`,b6->-0.3487`,b7->0.2234`,b8->-0.2972`,b9->0.0148`,b10->0.0428`,b11->0.0356`,b12->0,b13->0,b14->0,b15->0,Bc1->5.4605`,Bc2->-1.5439`,Bc3->1.3081`,Bc4->-0.4785`,Bc5->0.1054`,Bc6->0.8713`,Bc7->-0.4555`,Bc8->0.4442`,Bc9->-0.4152`,Bc10->-0.65`,Bc11->-0.0881`,Bc12->-0.3731`,Bc13->-0.1642`,Bc14->-0.1958`};
bfixstoch={Be1->0.9999999999999999`,Be2->0.5372296418863018`,ce2->0.32785722727809946`,Be3->0.04355277154000709`,Be4->-0.0816097388015601`,Be5->3.585450892329512`,Be6->-0.05507623212386243`,Be7->-0.5733364850490485`,Be8->0.15465075050230467`,Be9->0.1558917385651814`,Be10->-0.09017846590237562`,Be11->0.25233423945160144`,Be12->-0.2554071622739629`,Bd1->3.1911121616830163`,Bd2->-0.4260725682543434`,Bd3->-0.25209786077295826`};
bfixstochNS={Be1->0,Be2->0.5372296418863018`,ce2->0.32785722727809946`,Be3->0.04355277154000709`,Be4->-0.0816097388015601`,Be5->3.585450892329512`,Be6->-0.05507623212386243`,Be7->-0.5733364850490485`,Be8->0.15465075050230467`,Be9->0.1558917385651814`,Be10->-0.09017846590237562`,Be11->0.25233423945160144`,Be12->-0.2554071622739629`,Bd1->3.1911121616830163`,Bd2->-0.4260725682543434`,Bd3->-0.25209786077295826`};

bfixDESI = Join[bfixBOSSnonstoch,bfixstoch];
bfixDESINS = Join[bfixBOSSnonstoch,bfixstochNS];


(* ::Text:: *)
(*DESI numbers fully binned*)


(*z,nb,V,N,Dg,fs,b1*)


DESIspec =  {{0.65`,0.75`,0.85`,0.95`,1.05`,1.15`,1.25`,1.35`,1.45`,1.55`,1.65`},
			 {0.00016448669201520912`,0.0010084444444444446`,0.0007375890410958904`,0.0007150243902439025`,0.0004463274336283186`,0.00038736196319018404`,0.0003585823754789272`,0.00013312727272727273`,0.00011346086956521739`,0.00007715242881072026`,0.000028682926829268292`},
			 {2.63`*^9,3.15`*^9,3.65`*^9,4.0999999999999995`*^9,4.52`*^9,4.89`*^9,5.22`*^9,5.5`*^9,5.75`*^9,5.97`*^9,6.15`*^9},
			 {432600.`,3.1766`*^6,2.6922`*^6,2.9316`*^6,2.0174`*^6,1.8942`*^6,1.8718`*^6,732200.`,652400.`,460600.`,176400.`}, 
			 {0.7147268398976079,0.6814195172840236,0.6505961593752985,0.6220672608943191,0.5956450407211369,0.5711494009737322,0.5484113115796533,0.5272744484842405,0.5075956926551383,0.4892449166420755,0.47210435196678463},
			 {0.799023714461632,0.8226875527129273,0.8429531428816396,0.8603132350398828,0.8752024918671701,0.8879970333021452,0.8990181249903924,0.9085378342658417,0.9167853192345018,0.923952996823696,0.9302022056485543},
			 {1.1752741790420782`,1.2327207816823922`,1.2911235147876785`,1.350336294490677`,1.4102358662854422`,1.4707185170253423`,1.5316970716385292`,1.5930982478190507`,1.6548603783576588`,1.716931482426688`,1.7792676481387295`}};
DESIcut=4;


(* ::Text:: *)
(*Get zeffs and define cental values*)


Getzeffs[DESIspec,DESIcut]=={0.8379757391963608`,1.226878923766816`}


D1spec = DESIspec[[All,;;DESIcut]];
D2spec = DESIspec[[All,DESIcut+1;;]];


(*DESI1pk = PlinsubsT["plins/DESIlNk/",0.84];
DESI2pk = PlinsubsT["plins/DESIhNk/",1.23];
Export["aux_vars/DESI1pk.m",DESI1pk]
Export["aux_vars/DESI2pk.m",DESI2pk]*)

DESI1pk=Import["aux_vars/DESI1pk.m"];
DESI2pk=Import["aux_vars/DESI2pk.m"];


(* ::Text:: *)
(*Shotnoise estimates (these are not entering the covariance, so they are not very important) *)


nbD1 = nb->0.0007967599137945827`;
nbD2 = nb->0.0003225432700956092`;


(* ::Text:: *)
(*Get knls*)


knlD1 = {knl->Getknl[DESI1pk[[1,3]]][[1,2]]}
knlD2 = {knl->Getknl[DESI2pk[[1,3]]][[1,2]]}


(* ::Text:: *)
(*Parameters bestfit from BOSS, rescaled to DESI. We study ELG's with b1 = 0.84/D(z) (arxiv: 1611.00036)*)


bzDESI=Transpose@{DESIspec[[1]],DESIspec[[7]]}//Interpolation;

bfixD1=shiftBiasFix[1.2840552874707236`,bfixDESI];
bfixD2=shiftBiasFix[1.517557939573598`,bfixDESI];

coeffixD1 = Join[bfixD1,fnlfix,paramfix,knlD1];
coeffixD2 = Join[bfixD2,fnlfix,paramfix,knlD2];
coeffixD1NS = Join[{Be1->0},coeffixD1];
coeffixD2NS = Join[{Be1->0},coeffixD2];


(* ::Text:: *)
(*Base tab and setting kmax *)


ktabraw = Table[0.001`+i 0.005,{i,0,300}];
(*Export["covsks/ktabBase.m",ktabraw]*)


Basetri = makeTri[0.02,0.91,0.02];
keffbase = Import["covsks/baseTrikeff.mtx"];


(*Check*)
Dimensions[Basetri]==Dimensions[keffbase]


Bins = Basetri//Flatten//DeleteDuplicates//Sort;
Setkmax[kmax_]:=Select[Bins,#<kmax&];
goodTri[kmax_] := Complement[Table[i,{i,1,Length[Basetri]}],
							Table[Position[Basetri,Complement[Bins,Setkmax[kmax]][[i]]][[All,1]],
							{i,1,Length[Complement[Bins,Setkmax[kmax]]]}]//Flatten//DeleteDuplicates//Sort];


(* ::Text:: *)
(*Binning sizes*)


dkPD = dk->0.005;
dkBD=dk->0.02;


(* ::Text:: *)
(*Powerspectrum*)


Get\[Sigma][ktab_, Cov_]:=Transpose@{ktab,1/mono[Cov]}//Interpolation;


(*COVCMASSPgen=GetFullzCovP[ktabraw,CMASSpk,CMASSspec,bfixBOSS,dkPD];
COVCMASSPgenNS=GetFullzCovP[ktabraw,CMASSpk,CMASSspec,bfixBOSSNS,dkPD];

Export["aux_vars/COVCMASSPgen.m",COVCMASSPgen]
Export["aux_vars/COVCMASSPgenNS.m",COVCMASSPgenNS]*)
COVCMASSPgen=Import["aux_vars/COVCMASSPgen.m"];
COVCMASSPgenNS=Import["aux_vars/COVCMASSPgenNS.m"];


\[Sigma]CMASS = Get\[Sigma][ktabraw, COVCMASSPgen];
\[Sigma]CMASSNS = Get\[Sigma][ktabraw, COVCMASSPgenNS];


COVD1Pkmm=GetFullzCovP[ktabraw,DESI1pk,D1spec,bfixDESI,dkPD];
COVD2Pkmm=GetFullzCovP[ktabraw,DESI2pk,D2spec,bfixDESI,dkPD];
COVD1PNSkmm=GetFullzCovP[ktabraw,DESI1pk,D1spec,bfixDESINS,dkPD];
COVD2PNSkmm=GetFullzCovP[ktabraw,DESI2pk,D2spec,bfixDESINS,dkPD];

Export["aux_vars/COVD1Pkmm.m",COVD1Pkmm]
Export["aux_vars/COVD2Pkmm.m",COVD2Pkmm]
Export["aux_vars/COVD1PNSkmm.m",COVD1PNSkmm]
Export["aux_vars/COVD2PNSkmm.m",COVD2PNSkmm]


GetSlope[ku_,kd_,pk_]:=(ku+kd)/2  (Log[pk[ku]]-Log[pk[kd]])/(ku-kd);


GetkmaxNSP[\[Sigma]1_,pk1_]:=Quiet[FindRoot[Sqrt[\[Sigma]1[kmm]/\[Sigma]CMASSNS[0.22]]-pk1[[1,3]][kmm]/CMASSpk[[1,3]][0.22] (kmm/Getknl[pk1[[1,3]]][[1,2]])^((3+GetSlope[0.8kmm,1.2kmm,CMASSpk[[1,3]]])2) (0.22/0.7)^(-(3+(-1.8))2),{kmm,0.1}]]


GetkmaxNSP2L[\[Sigma]1_,pk1_]:=Quiet[FindRoot[Sqrt[\[Sigma]1[kmm]/\[Sigma]CMASSNS[0.22]]-pk1[[1,3]][kmm]/CMASSpk[[1,3]][0.22] (kmm/Getknl[pk1[[1,3]]][[1,2]])^((3+GetSlope[0.8kmm,1.2kmm,CMASSpk[[1,3]]])3) (0.22/0.7)^(-(3+(-1.8))2),{kmm,0.1}]]
GetkmaxNSP1L[\[Sigma]1_,pk1_]:=Quiet[FindRoot[Sqrt[\[Sigma]1[kmm]/\[Sigma]CMASSNS[0.22]]-pk1[[1,3]][kmm]/CMASSpk[[1,3]][0.22] (kmm/Getknl[pk1[[1,3]]][[1,2]])^((3+GetSlope[0.8kmm,1.2kmm,CMASSpk[[1,3]]])1) (0.22/0.7)^(-(3+(-1.8))2),{kmm,0.05}]]


kmaxD1 = GetkmaxNSP[Get\[Sigma][ktabraw, COVD1PNSkmm],DESI1pk][[1,2]]
kmaxD2 =GetkmaxNSP[Get\[Sigma][ktabraw, COVD2PNSkmm],DESI2pk][[1,2]]


kmaxD12L = GetkmaxNSP2L[Get\[Sigma][ktabraw, COVD1PNSkmm],DESI1pk][[1,2]]
kmaxD22L =GetkmaxNSP2L[Get\[Sigma][ktabraw, COVD2PNSkmm],DESI2pk][[1,2]]


kmaxD1Tree = GetkmaxNSP1L[Get\[Sigma][ktabraw, COVD1PNSkmm],DESI1pk][[1,2]]
kmaxD2Tree =GetkmaxNSP1L[Get\[Sigma][ktabraw, COVD2PNSkmm],DESI2pk][[1,2]]


(* ::Text:: *)
(*Get triangles and ktabs*)


ktabDESI1 = Select[ktabraw,#<kmaxD1&];
ktabDESI2 = Select[ktabraw,#<kmaxD2&];


ktabDESI12L = Select[ktabraw,#<kmaxD12L&];
ktabDESI22L = Select[ktabraw,#<kmaxD22L&];


keffDESI1 = keffbase[[goodTri[kmaxD1]]];
keffDESI2 = keffbase[[goodTri[kmaxD2]]];
DESItriangs1 = Basetri[[goodTri[kmaxD1]]];
DESItriangs2 = Basetri[[goodTri[kmaxD2]]];


keffDESI1T = keffbase[[goodTri[kmaxD1Tree]]];
keffDESI2T = keffbase[[goodTri[kmaxD2Tree]]];
DESItriangs1T = Basetri[[goodTri[kmaxD1Tree]]];
DESItriangs2T = Basetri[[goodTri[kmaxD2Tree]]];


(* ::Text:: *)
(*Analytical Covariance*)


COVD1P=GetFullzCovP[ktabDESI1,DESI1pk,D1spec,bfixDESI,dkPD];
COVD2P=GetFullzCovP[ktabDESI2,DESI2pk,D2spec,bfixDESI,dkPD];


COVD1P2L=GetFullzCovP[ktabDESI12L,DESI1pk,D1spec,bfixDESI,dkPD];
COVD2P2L=GetFullzCovP[ktabDESI22L,DESI2pk,D2spec,bfixDESI,dkPD];


COVD1B=GetFullzCovB[DESItriangs1,DESI1pk,D1spec,bfixDESI,dkBD];
COVD2B=GetFullzCovB[DESItriangs2,DESI2pk,D2spec,bfixDESI,dkBD];
COVD1BTree=GetFullzCovB[DESItriangs1T,DESI1pk,D1spec,bfixDESI,dkBD];
COVD2BTree=GetFullzCovB[DESItriangs2T,DESI2pk,D2spec,bfixDESI,dkBD];

Export["aux_vars/COVD1B.m",COVD1B]
Export["aux_vars/COVD2B.m",COVD2B]
Export["aux_vars/COVD1BTree.m",COVD1BTree]
Export["aux_vars/COVD2BTree.m",COVD2BTree]

(*COVD1B=Import["aux_vars/COVD1B.m"];
COVD2B=Import["aux_vars/COVD2B.m"];
COVD1BTree=Import["aux_vars/COVD1BTree.m"];
COVD2BTree=Import["aux_vars/COVD2BTree.m"];*)


(* ::Text:: *)
(*Analytical Covariance with Full redshift- No shot noise*)


COVD1PNS=GetFullzCovP[ktabDESI1,DESI1pk,D1spec,bfixDESINS,dkPD];
COVD2PNS=GetFullzCovP[ktabDESI2,DESI2pk,D2spec,bfixDESINS,dkPD];


COVD1BNS=GetFullzCovB[DESItriangs1,DESI1pk,D1spec,bfixDESINS,dkBD];
COVD2BNS=GetFullzCovB[DESItriangs2,DESI2pk,D2spec,bfixDESINS,dkBD];


Export["aux_vars/COVD1BNS.m",COVD1BNS]
Export["aux_vars/COVD2BNS.m",COVD2BNS]


(* ::Text:: *)
(*Powerspectrum Derivatives*)


PD1derivs = Table[P1LoopDerivs[ktabDESI1[[i]],DESI1pk]/.coeffixD1/.nbD1,{i,1,Length[ktabDESI1]}];
PD2derivs = Table[P1LoopDerivs[ktabDESI2[[i]],DESI2pk]/.coeffixD2/.nbD2,{i,1,Length[ktabDESI2]}];


PD1derivs2L = Table[P1LoopDerivs[ktabDESI12L[[i]],DESI1pk]/.coeffixD1/.nbD1,{i,1,Length[ktabDESI12L]}];
PD2derivs2L = Table[P1LoopDerivs[ktabDESI22L[[i]],DESI2pk]/.coeffixD2/.nbD2,{i,1,Length[ktabDESI12L]}];


(* ::Text:: *)
(*Bispectrum Derivatives*)


BD1derivsT=Monitor[Table[BTreeDerivs[keffDESI1T[[si]],DESI1pk]/.Dispatch[coeffixD1]/.nbD1,{si,1,Length[keffDESI1T]}],si];//AbsoluteTiming
BD2derivsT=Monitor[Table[BTreeDerivs[keffDESI2T[[si]],DESI2pk]/.Dispatch[coeffixD2]/.nbD2,{si,1,Length[keffDESI2T]}],si];//AbsoluteTiming


(*Calculate derivatives, one at a time instead of table*)
BD1derivs=ConstantArray[0,Length[keffDESI1]];
BD2derivs=ConstantArray[0,Length[keffDESI2]];

BD1derivsimp=Import["aux_vars/BD1derivs.m"];
BD2derivsimp=Import["aux_vars/BD2derivs.m"];


BD1derivs[[1;;Length[BD1derivsimp]]]=BD1derivsimp;


dispatcher = Dispatch[Join[coeffixD1,{nbD1}]];

Monitor[Do[If[(BD1derivs[[j]]==0)||(Length[BD1derivs[[j]]]==0),
			BD1derivs[[j]]=B1LoopDerivs[keffDESI1[[j]],DESI1pk]/.dispatcher];,
	{j,1,Length[keffDESI1]}],j];//AbsoluteTiming
(*Export["aux_vars/BD1derivs.m",BD1derivs]*)


BD2derivs[[1;;Length[BD2derivsimp]]]=BD2derivsimp;


dispatcher =Dispatch[Join[coeffixD2,{nbD2}]];
Monitor[Do[
	If[BD2derivs[[j]]==0,
	BD2derivs[[j]]=B1LoopDerivs[keffDESI2[[j]],DESI2pk]/.dispatcher];,
	{j,1,Length[keffDESI2]}],
	j];//AbsoluteTiming
(*{400.0957859`,Null}*)
(*Export["aux_vars/BD2derivs.m",BD2derivs]*)


(*checks*)
Dimensions[BD1derivs]
Dimensions[BD2derivs]


(* ::Subsubsection::Closed:: *)
(*Fishers *)


(* ::Text:: *)
(*Reduce P-derivatives to monopole and quadrupole*)


PD1monoquad= removelg2[PD1derivs]//Expand;
PD2monoquad= removelg2[PD2derivs]//Expand;


(* ::Text:: *)
(*Power spectrum Fishers mono+quad*)


FPD1= GetFisherP[PD1monoquad,COVD1P];
FPD2 = GetFisherP[PD2monoquad,COVD2P];


(* ::Text:: *)
(*Bispectrum Fishers - mono*)


FBDESI1 = GetFisherB[BD1derivs,COVD1B,BLoopmastermono];
FBDESI2 = GetFisherB[BD2derivs,COVD2B,BLoopmastermono];


(* ::Text:: *)
(*Powerspectrum Fisher full \[Mu]*)


FPDESI1full= GetFisherP[PD1derivs,COVD1P];
FPDESI2full = GetFisherP[PD2derivs,COVD2P];


(* ::Text:: *)
(*Powerspectrum Fisher full \[Mu] 2 loop*)


FPD1full2L= GetFisherP[PD1derivs2L,COVD1P2L];
FPD2full2L = GetFisherP[PD2derivs2L,COVD2P2L];


(* ::Text:: *)
(*Bispectrum Fishers full \[Mu]*)


FBD1full = GetFisherB[BD1derivs,COVD1B,BLoopmaster];
FBD2full = GetFisherB[BD2derivs,COVD2B,BLoopmaster];


(* ::Text:: *)
(*Bispectrum Fishers full \[Mu] Tree*)


FBD1fullT = GetFisherB[BD1derivsT,COVD1BTree,BLoopmaster];
FBD2fullT = GetFisherB[BD2derivsT,COVD2BTree,BLoopmaster];


(* ::Text:: *)
(*Powerspectrum Fishers full NS*)


FPD1fullNS= GetFisherP[PD1derivs,COVD1PNS];
FPD2fullNS = GetFisherP[PD2derivs,COVD2PNS];


(* ::Text:: *)
(*Bispectrum Fishers full NS*)


FBD1fullNS = GetFisherB[BD1derivs,COVD1BNS,BLoopmaster];
FBD2fullNS = GetFisherB[BD2derivs,COVD2BNS,BLoopmaster];


(* ::Subsection:: *)
(*MegaMapper - optimistic*)


(* ::Subsubsection:: *)
(*MMo specifications and derivative computations*)


(* ::Text:: *)
(*MMo numbers fully binned (section 1.2 in arxiv 2209.04322, then arxiv 1903.09208 Table 1 )*)


(* ::Input:: *)
(*(*z,nb,V,N,Dg,fs,b1*)*)


MMospec = {{2,2.5,3,3.5,4,4.5,5},
			10^-4 {25,12,6,3,1.5,0.8,0.4},
			10^9 {26.986863,27.63312497,27.20703913,26.25383626,25.06712437,23.80296403,22.54324749},
			{6.74671575`*^7,3.3159749963999998`*^7,1.6324223477999998`*^7,7.876150878`*^6,3.7600686555000003`*^6,1.9042371224000002`*^6,901729.8996000001`},
			{0.42020110428342006,0.36272799579672305,0.3188745725930273,0.2844001952601862,0.25662379231412885,0.23378307001365894,0.2146779141872147},
			{0.9466075281915435,0.9606657166004888,0.9687825324004987,0.973748370874903,0.9769321736957639,0.9790517015436108,0.9805057101340311},
			{2.5,3.3,4.1,4.9,5.8,6.6,7.4}};
MMocut=4;


(*biases*)
bfixBOSSnonstoch={b1->1.9391`,c2->1.1447`,b3->-0.374`,b4->0.1276`,c4->-0.2895`,b6->-0.3487`,b7->0.2234`,b8->-0.2972`,b9->0.0148`,b10->0.0428`,b11->0.0356`,b12->0,b13->0,b14->0,b15->0,Bc1->5.4605`,Bc2->-1.5439`,Bc3->1.3081`,Bc4->-0.4785`,Bc5->0.1054`,Bc6->0.8713`,Bc7->-0.4555`,Bc8->0.4442`,Bc9->-0.4152`,Bc10->-0.65`,Bc11->-0.0881`,Bc12->-0.3731`,Bc13->-0.1642`,Bc14->-0.1958`};
bfixstoch={Be1->0.9999999999999999`,Be2->0.5372296418863018`,ce2->0.32785722727809946`,Be3->0.04355277154000709`,Be4->-0.0816097388015601`,Be5->3.585450892329512`,Be6->-0.05507623212386243`,Be7->-0.5733364850490485`,Be8->0.15465075050230467`,Be9->0.1558917385651814`,Be10->-0.09017846590237562`,Be11->0.25233423945160144`,Be12->-0.2554071622739629`,Bd1->3.1911121616830163`,Bd2->-0.4260725682543434`,Bd3->-0.25209786077295826`};
bfixstochNS={Be1->0,Be2->0.5372296418863018`,ce2->0.32785722727809946`,Be3->0.04355277154000709`,Be4->-0.0816097388015601`,Be5->3.585450892329512`,Be6->-0.05507623212386243`,Be7->-0.5733364850490485`,Be8->0.15465075050230467`,Be9->0.1558917385651814`,Be10->-0.09017846590237562`,Be11->0.25233423945160144`,Be12->-0.2554071622739629`,Bd1->3.1911121616830163`,Bd2->-0.4260725682543434`,Bd3->-0.25209786077295826`};

bfixMM = Join[bfixBOSSnonstoch,bfixstoch];
bfixMMNS = Join[bfixBOSSnonstoch,bfixstochNS];


(* ::Text:: *)
(*Get zeffs and define cental values*)


(* ::Input:: *)
(**)


MMo1spec = MMospec[[All,;;MMocut]];
MMo2spec = MMospec[[All,MMocut+1;;]];


Getzeffs[MMospec,MMocut]
Dg236=0.37721157585659637;
Dg428 =0.2433118256462636;


(*MMo1pk = PlinsubsT["plins/MMo1Nk/",2.36];
MMo2pk = PlinsubsT["plins/MMo2Nk/",4.28];
Export[auxPath<>"MMo1pk.m",MMo1pk]
Export[auxPath<>"MMo2pk.m",MMo2pk]*)
MMo1pk=Import[auxPath<>"MMo1pk.m"];
MMo2pk=Import[auxPath<>"MMo2pk.m"];


(* ::Text:: *)
(*Shotnoise estimates (these are not entering the covariance, so they are not very important) *)


Getnbeff[MMospec,MMocut]
nbMMo1 = nb->0.0017673778507420201`;
nbMMo2 = nb->0.00011459250315671164`;


(* ::Text:: *)
(*Get knls*)


Getknl[pk_]:=Quiet[FindRoot[1/(2\[Pi]^2) NIntegrate[k^2 pk[k],{k,10^-5,kn}]-1.1525490955287694`,{kn,0.5}]]
knlMMo1 = {knl->Getknl[MMo1pk[[1,3]]][[1,2]]}
knlMMo2 = {knl->Getknl[MMo2pk[[1,3]]][[1,2]]}


(* ::Text:: *)
(*Parameters bestfit from BOSS, rescaled to MMo. (arxiv 1903.09208 Table 1)*)


bzMMo=Transpose@{MMospec[[1]],MMospec[[7]]}//Interpolation;
bzMMo[Getzeffs[MMospec,MMocut]]


bfixMMo1=shiftBiasFix[3.073186555054316`,bfixMM];
bfixMMo2=shiftBiasFix[6.257622960678745`,bfixMM];


coeffixMMo1 = Join[bfixMMo1,fnlfix,paramfix,knlMMo1];
coeffixMMo2 = Join[bfixMMo2,fnlfix,paramfix,knlMMo2];


(* ::Text:: *)
(*Base tab and setting kmax *)


Bins = Basetri//Flatten//DeleteDuplicates//Sort;
Setkmax[kmax_]:=Select[Bins,#<kmax&];
goodTri[kmax_] := Complement[Table[i,{i,1,Length[Basetri]}],Table[Position[Basetri,Complement[Bins,Setkmax[kmax]][[i]]][[All,1]],
							{i,1,Length[Complement[Bins,Setkmax[kmax]]]}]//Flatten//DeleteDuplicates//Sort];


(* ::Text:: *)
(*Binning sizes*)


dkPMM = dk->0.005;
dkBMM=dk->0.02;


(* ::Text:: *)
(*Powerspectrum kmax*)


Get\[Sigma][ktab_, Cov_]:=Transpose@{ktab,1/(Cov//mono)}//Interpolation;


COVCMASSPgenNS=Import["aux_vars/COVCMASSPgenNS.m"];
CMASSpk=Import["aux_vars/CMASSpk.m"];


(*check*)
Variables[COVCMASSPgenNS[[1]]]=={\[Mu]1}


\[Sigma]CMASSNS = Get\[Sigma][ktabraw, COVCMASSPgenNS];
COVMMo1PNSkmm=GetFullzCovP[ktabraw,MMo1pk,MMo1spec,bfixMMNS,dkPMM];
COVMMo2PNSkmm=GetFullzCovP[ktabraw,MMo2pk,MMo2spec,bfixMMNS,dkPMM];


GetSlope[ku_,kd_,pk_]:=(ku+kd)/2  (Log[pk[ku]]-Log[pk[kd]])/(ku-kd);
GetkmaxNSP[\[Sigma]1_,pk1_]:=Quiet[FindRoot[Sqrt[\[Sigma]1[kmm]/\[Sigma]CMASSNS[0.22]]-pk1[[1,3]][kmm]/CMASSpk[[1,3]][0.22] (kmm/Getknl[pk1[[1,3]]][[1,2]])^((3+GetSlope[0.8kmm,1.2kmm,CMASSpk[[1,3]]])2) (0.22/0.7)^(-(3+(-1.8))2),{kmm,0.1}]]
GetkmaxNSP2L[\[Sigma]1_,pk1_]:=Quiet[FindRoot[Sqrt[\[Sigma]1[kmm]/\[Sigma]CMASSNS[0.22]]-pk1[[1,3]][kmm]/CMASSpk[[1,3]][0.22] (kmm/Getknl[pk1[[1,3]]][[1,2]])^((3+GetSlope[0.8kmm,1.2kmm,CMASSpk[[1,3]]])3) (0.22/0.7)^(-(3+(-1.8))2),{kmm,0.1}]]
GetkmaxNSP1L[\[Sigma]1_,pk1_]:=Quiet[FindRoot[Sqrt[\[Sigma]1[kmm]/\[Sigma]CMASSNS[0.22]]-pk1[[1,3]][kmm]/CMASSpk[[1,3]][0.22] (kmm/Getknl[pk1[[1,3]]][[1,2]])^((3+GetSlope[0.8kmm,1.2kmm,CMASSpk[[1,3]]])1) (0.22/0.7)^(-(3+(-1.8))2),{kmm,0.1}]]


kmaxMMo1 = GetkmaxNSP[Get\[Sigma][ktabraw, COVMMo1PNSkmm],MMo1pk][[1,2]]
kmaxMMo2 =GetkmaxNSP[Get\[Sigma][ktabraw, COVMMo2PNSkmm],MMo2pk][[1,2]]
kmaxMMo12L= GetkmaxNSP2L[Get\[Sigma][ktabraw, COVMMo1PNSkmm],MMo1pk][[1,2]]
kmaxMMo22L =GetkmaxNSP2L[Get\[Sigma][ktabraw, COVMMo2PNSkmm],MMo2pk][[1,2]]
kmaxMMo1Tree = GetkmaxNSP1L[Get\[Sigma][ktabraw, COVMMo1PNSkmm],MMo1pk][[1,2]]
kmaxMMo2Tree =GetkmaxNSP1L[Get\[Sigma][ktabraw, COVMMo2PNSkmm],MMo2pk][[1,2]]


(* ::Text:: *)
(*Get triangles and ktabs*)


ktabMMo1 = Select[ktabraw,#<kmaxMMo1&];
ktabMMo2 = Select[ktabraw,#<kmaxMMo2&];


ktabMMo12L = Select[ktabraw,#<kmaxMMo12L&];
ktabMMo22L = Select[ktabraw,#<1&];


keffMMo1 = keffbase[[goodTri[kmaxMMo1]]];
keffMMo2 = keffbase[[goodTri[kmaxMMo2]]];
MMotriangs1 = Basetri[[goodTri[kmaxMMo1]]];
MMotriangs2 = Basetri[[goodTri[kmaxMMo2]]];


keffMMo1T = keffbase[[goodTri[kmaxMMo1Tree]]];
keffMMo2T = keffbase[[goodTri[kmaxMMo2Tree]]];
MMotriangs1T = Basetri[[goodTri[kmaxMMo1Tree]]];
MMotriangs2T = Basetri[[goodTri[kmaxMMo2Tree]]];


(* ::Text:: *)
(*Analytical Covariance with Full redshift*)


COVMMo1P=GetFullzCovP[ktabMMo1,MMo1pk,MMo1spec,bfixMM,dkPMM];
COVMMo2P=GetFullzCovP[ktabMMo2,MMo2pk,MMo2spec,bfixMM,dkPMM];


COVMMo1P2L=GetFullzCovP[ktabMMo12L,MMo1pk,MMo1spec,bfixMM,dkPMM];
COVMMo2P2L=GetFullzCovP[ktabMMo22L,MMo2pk,MMo2spec,bfixMM,dkPMM];


(*COVMMo1B=GetFullzCovB[MMotriangs1,MMo1pk,MMo1spec,bfixMM,dkBMM];
COVMMo2B=GetFullzCovB[MMotriangs2,MMo2pk,MMo2spec,bfixMM,dkBMM];
Export[auxPath<>"COVMMo1B.m",COVMMo1B]
Export[auxPath<>"COVMMo2B.m",COVMMo2B]*)
COVMMo1B=Import[auxPath<>"COVMMo1B.m"];
COVMMo2B=Import[auxPath<>"COVMMo2B.m"];

(*COVMMo1BTree=GetFullzCovB[MMotriangs1T,MMo1pk,MMo1spec,bfixMM,dkBMM];
COVMMo2BTree=GetFullzCovB[MMotriangs2T,MMo2pk,MMo2spec,bfixMM,dkBMM];
Export[auxPath<>"COVMMo1BTree.m",COVMMo1BTree]
Export[auxPath<>"COVMMo2BTree.m",COVMMo2BTree]*)
COVMMo1BTree=Import[auxPath<>"COVMMo1BTree.m"];
COVMMo2BTree=Import[auxPath<>"COVMMo2BTree.m"];


(* ::Text:: *)
(*Analytical Covariance with Full redshift- No shot noise*)


COVMMo1PNS=GetFullzCovP[ktabMMo1,MMo1pk,MMo1spec,bfixMMNS,dkPMM];
COVMMo2PNS=GetFullzCovP[ktabMMo2,MMo2pk,MMo2spec,bfixMMNS,dkPMM];


COVMMo1PNS2L=GetFullzCovP[ktabMMo12L,MMo1pk,MMo1spec,bfixMMNS,dkPMM];
COVMMo2PNS2L=GetFullzCovP[ktabMMo22L,MMo2pk,MMo2spec,bfixMMNS,dkPMM];


(*COVMMo1BNS=GetFullzCovB[MMotriangs1,MMo1pk,MMo1spec,bfixMMNS,dkBMM];
COVMMo2BNS=GetFullzCovB[MMotriangs2,MMo2pk,MMo2spec,bfixMMNS,dkBMM];
Export[auxPath<>"COVMMo1BNS.m",COVMMo1BNS]
Export[auxPath<>"COVMMo2BNS.m",COVMMo2BNS]*)
COVMMo1BNS=Import[auxPath<>"COVMMo1BNS.m"];
COVMMo2BNS=Import[auxPath<>"COVMMo2BNS.m"];


(* ::Text:: *)
(*Powerspectrum Derivatives*)


PMMo1derivs =Monitor[Table[P1LoopDerivs[ktabMMo1[[i]],MMo1pk],{i,1,Length[ktabMMo1]}],i]/.coeffixMMo1/.nbMMo1;
PMMo2derivs = Monitor[Table[P1LoopDerivs[ktabMMo2[[i]],MMo2pk],{i,1,Length[ktabMMo2]}],i]/.coeffixMMo2/.nbMMo2;


PMMo1derivs2L =Monitor[Table[P1LoopDerivs[ktabMMo12L[[i]],MMo1pk],{i,1,Length[ktabMMo12L]}],i]/.coeffixMMo1/.nbMMo1;
PMMo2derivs2L = Monitor[Table[P1LoopDerivs[ktabMMo22L[[i]],MMo2pk],{i,1,Length[ktabMMo22L]}],i]/.coeffixMMo2/.nbMMo2;


(* ::Text:: *)
(*Bispectrum Derivatives*)


BMMo1derivsT=Monitor[Table[BTreeDerivs[keffMMo1T[[si]],MMo1pk]/.Dispatch[coeffixMMo1]/.nbMMo1,{si,1,Length[keffMMo1T]}],si];//AbsoluteTiming
BMMo2derivsT=Monitor[Table[BTreeDerivs[keffMMo2T[[si]],MMo2pk]/.Dispatch[coeffixMMo2]/.nbMMo2,{si,1,Length[keffMMo2T]}],si];//AbsoluteTiming


Length[keffMMo1]==805
Length[keffMMo2]==8878


(*Calculate derivatives, one at a time instead of table*)
(*BMMo1derivs=ConstantArray[0,Length[keffMMo1]];
BMMo2derivs=ConstantArray[0,Length[keffMMo2]];*)


BMMo1derivsimp=Import["aux_vars/BMMo1derivs.m"];
(*BMMo2derivsimp=Import["aux_vars/BMMo2derivs.m"];*)
(*BMMo1derivs[[1;;Length[BMMo1derivsimp]]]=BMMo1derivsimp;*)


BMMo1derivs[1]=ConstantArray[0,300];
firstIndex[1]=1;
finalIndex[1]=firstIndex[1]+Length[BMMo1derivs[1]]-1;

BMMo1derivs[2]=ConstantArray[0,300];
firstIndex[2]=finalIndex[1]+1;
finalIndex[2]=firstIndex[2]+Length[BMMo1derivs[2]]-1;

BMMo1derivs[3]=ConstantArray[0,Length[keffMMo1]-Length[BMMo1derivs[2]]-Length[BMMo1derivs[1]]];
firstIndex[3]=finalIndex[2]+1;
finalIndex[3]=firstIndex[3]+Length[BMMo1derivs[3]]-1;

Do[keffMMo1part[part]=keffMMo1[[firstIndex[part];;finalIndex[part]]],{part,1,3}]


(*BMMo1derivs[1]=BMMo1derivsimp[[1;;300]];*)
ClearAll[BMMo1derivsimp]
BMMo1derivs[1]=.


(*Column@Thread@{Names["FisherFunctions`*"], ByteCount[#] & /@ ToExpression /@ Names["FisherFunctions`*"]}*)


dispatcher = Dispatch[Join[coeffixMMo1,{nbMMo1}]];

Monitor[
Do[
	Bderivslist = BMMo1derivs[part];
	triList=keffMMo1part[part];
	Do[
		If[(Bderivslist[[j]]==0)||(Length[Bderivslist[[j]]]==0),
			Bderivslist[[j]]=B1LoopDerivs[triList[[j]],MMo1pk]/.dispatcher];,
	{j,1,Length[triList]}];
		
	Export["aux_vars/BMMo1derivs"<>ToString[part]<>".m",Bderivslist];
	(*BMMo1derivs[part]=Bderivslist;*),
{part,2,3}],
Row[{ProgressIndicator[j,{1,Length[triList]}],j," - Part:",part}," "]]//AbsoluteTiming


(*Export["aux_vars/BMMo1derivs.m",BMMo1derivs]*)


BMMo2derivs[1]=ConstantArray[0,500];
firstIndex2[1]=1;
finalIndex2[1]=firstIndex2[1]+Length[BMMo2derivs[1]]-1;
BMMo2derivs[2]=ConstantArray[0,1000];
firstIndex2[2]=finalIndex2[1]+1;
finalIndex2[2]=firstIndex2[2]+Length[BMMo2derivs[2]]-1;
BMMo2derivs[3]=ConstantArray[0,Length[keffMMo1]-Length[BMMo2derivs[2]]-Length[BMMo2derivs[1]]];
firstIndex2[3]=finalIndex2[2]+1;
finalIndex2[3]=firstIndex2[3]+Length[BMMo2derivs[3]]-1;


(*BD2derivs[[1;;Length[BD2derivsimp]]]=BD2derivsimp;*)


dispatcher = Dispatch[Join[coeffixMMo2,{nbMMo2}]];

Monitor[
Do[Bderivslist = BMMo2derivs[part];
	Do[If[(Bderivslist[[j]]==0)||(Length[Bderivslist[[j]]]==0),
			Bderivslist[[j]]=B1LoopDerivs[keffMMo2[[j]],MMo2pk]/.dispatcher];,
		{j,firstIndex2[part],finalIndex2[part]}];
		
	Export["aux_vars/BMMo2derivs"<>ToString[part]<>".m",Bderivslist];
	BMMo2derivs[part]=Bderivslist;
	,{part,1,3}],
	Row[{ProgressIndicator[j,{1,Length[keffMMo2]}],j," - Part:",part}," "]]//AbsoluteTiming


(*{400.0957859`,Null}*)
(*Export["aux_vars/BMMo2derivs.m",BMMo2derivs]*)


(*checks*)
Dimensions[BMMo1derivs]
Dimensions[BMMo2derivs]


(* ::Input:: *)
(*(*BMMo1derivs=Monitor[Table[B1LoopDerivs[keffMMo1[[si]],MMo1pk]/.Dispatch[coeffixMMo1]/.nbMMo1,{si,1,Length[keffMMo1]}],si];//AbsoluteTiming*)
(*BMMo2derivs=Monitor[Table[B1LoopDerivs[keffMMo2[[si]],MMo2pk]/.Dispatch[coeffixMMo2]/.nbMMo2,{si,1,Length[keffMMo2]}],si];//AbsoluteTiming*)
(*Export[auxPath<>"BMMo1derivs.m",BMMo1derivs]*)
(*Export[auxPath<>"BMMo2derivs.m",BMMo2derivs]*)*)


(* ::Input:: *)
(*BMMo1derivs=Import[auxPath<>"BMMo1derivs.m"];*)
(*BMMo2derivs=Import[auxPath<>"BMMo2derivs.m"];*)


(* ::Subsubsection::Closed:: *)
(*Fishers *)


(* ::Text:: *)
(*Reduce P-derivatives to monopole and quadrupole*)


(* ::Input:: *)
(*PMMo1monoquad= removelg2[PMMo1derivs]//Expand;*)
(*PMMo2monoquad= removelg2[PMMo2derivs]//Expand;*)


(* ::Text:: *)
(*Powerspectrum Fishers monoquad*)


(* ::Input:: *)
(*FPMMo1= GetFisherP[PMMo1monoquad,COVMMo1P];*)
(*FPMMo2 = GetFisherP[PMMo2monoquad,COVMMo2P];*)


(* ::Text:: *)
(*Bispectrum Fishers - mono*)


(* ::Input:: *)
(*FBMMo1 = GetFisherB[BMMo1derivs,COVMMo1B,Loopmastermono];*)
(*FBMMo2 = GetFisherB[BMMo2derivs,COVMMo2B,Loopmastermono];*)


(* ::Text:: *)
(*Powerspectrum Fishers full \[Mu]*)


(* ::Input:: *)
(*FPMMo1mu= GetFisherP[PMMo1derivs,COVMMo1P];*)
(*FPMMo2mu = GetFisherP[PMMo2derivs,COVMMo2P];*)


(* ::Input:: *)
(*FPMMo1mu2L= GetFisherP[PMMo1derivs2L,COVMMo1P2L];*)
(*FPMMo2mu2L = GetFisherP[PMMo2derivs2L,COVMMo2P2L];*)


(* ::Text:: *)
(*Bispectrum Fishers full \[Mu]*)


(* ::Input:: *)
(*FBMMo1mu = GetFisherB[BMMo1derivs,COVMMo1B,Loopmaster];*)
(*FBMMo2mu = GetFisherB[BMMo2derivs,COVMMo2B,Loopmaster];*)


(* ::Text:: *)
(*Bispectrum Fishers full \[Mu] Tree*)


(* ::Input:: *)
(*FBMMo1muT = GetFisherB[BMMo1derivsT,COVMMo1BTree,Loopmaster];*)
(*FBMMo2muT = GetFisherB[BMMo2derivsT,COVMMo2BTree,Loopmaster];*)


(* ::Text:: *)
(*Powerspectrum Fishers full NS*)


(* ::Input:: *)
(*FPMMo1muNS= GetFisherP[PMMo1derivs,COVMMo1PNS];*)
(*FPMMo2muNS = GetFisherP[PMMo2derivs,COVMMo2PNS];*)


(* ::Text:: *)
(*Bispectrum Fishers full NS*)


(* ::Input:: *)
(*FBMMo1muNS = GetFisherB[BMMo1derivs,COVMMo1BNS,Loopmaster];*)
(*FBMMo2muNS = GetFisherB[BMMo2derivs,COVMMo2BNS,Loopmaster];*)


(* ::Subsection:: *)
(*Exports*)


AllParams


Allparams = DeleteCases[AllParams,b12|b13|b14|b15];
GetFExp[Fisher_]:=SubFisher[Fisher+Prior,Allparams]
GetFExp4[Fisher_]:=SubFisher4sky[Fisher+Prior4sky,Allparams]
Exporter[Fisher_,name_]:=Export[name,GetFExp[Fisher]]
Exporter4[Fisher_,name_]:=Export[name,GetFExp4[Fisher]]



(* ::Subsubsection::Closed:: *)
(*BOSS Exports*)


(* ::Text:: *)
(*Full P results w/o 4 sky correlation*)


(* ::Input:: *)
(*Exporter4[CombineMs[{CMASSNGCFisherFullP,CMASSSGCFisherFullP,LOWZNGCFisherFullP,LOWZSGCFisherFullP}],"Fisher1/BOSS/Pfull4sky.mtx"]*)
(*Exporter[CMASSNGCFisherFullP+CMASSSGCFisherFullP+LOWZNGCFisherFullP+LOWZSGCFisherFullP,"Fisher1/BOSS/Pfullnoskycorr.mtx"]*)


(* ::Text:: *)
(*Full PB cross w/o 4 sky correlation*)


(* ::Input:: *)
(*Exporter4[CombineMs[{CMASSNGCFisherFullPB,CMASSSGCFisherFullPB,LOWZNGCFisherFullPB,LOWZSGCFisherFullPB}],"Fisher1/BOSS/PBcross4sky.mtx"]*)
(*Exporter[CMASSNGCFisherFullPB+CMASSSGCFisherFullPB+LOWZNGCFisherFullPB+LOWZSGCFisherFullPB,"Fisher1/BOSS/PBcrossnoskycorr.mtx"]*)


(* ::Text:: *)
(*Numerical covariance Diagonals  w/o 4 sky correlation,  w/o B*)


(* ::Input:: *)
(*Exporter4[CombineMs[{PFisherCMASSngcNumD,PFisherCMASSsgcNumD,PFisherLOWZngcNumD,PFisherLOWZsgcNumD}],"Fisher1/BOSS/PNumD4sky.mtx"]*)
(*Exporter[PFisherCMASSngcNumD+PFisherCMASSsgcNumD+PFisherLOWZngcNumD+PFisherLOWZsgcNumD,"Fisher1/BOSS/PNumDnoskycorr.mtx"]*)


(* ::Input:: *)
(*Exporter4[CombineMs[{PFisherCMASSngcNumD+BFisherCMASSngcNumD,PFisherCMASSsgcNumD+BFisherCMASSsgcNumD,PFisherLOWZngcNumD+BFisherLOWZngcNumD,PFisherLOWZsgcNumD+BFisherLOWZsgcNumD}],"Fisher1/BOSS/PBNumD4sky.mtx"]*)
(*Exporter[PFisherCMASSngcNumD+BFisherCMASSngcNumD+PFisherCMASSsgcNumD+BFisherCMASSsgcNumD+PFisherLOWZngcNumD+BFisherLOWZngcNumD+PFisherLOWZsgcNumD+BFisherLOWZsgcNumD,"Fisher1/BOSS/PBNumDnoskycorr.mtx"]*)


(* ::Text:: *)
(*Analytical covariance with full redshift dependence,  w/o B*)


(* ::Input:: *)
(*Exporter[FPCnew+FPLnew,"Fisher1/BOSS/Panafullz.mtx"]*)
(*Exporter[FPCnew+FPLnew+FBCnew+FBLnew,"Fisher1/BOSS/PBanafullz.mtx"]*)


(* ::Text:: *)
(*Analytical covariance with full redshift dependence,  w/o B no shot noise (full mu)*)


(* ::Input:: *)
(*Exporter[FPCNS+FPLNS,"Fisher1/BOSS/PanamufullzNS.mtx"]*)
(*Exporter[FPCNS+FPLNS+FBCNS+FBLNS,"Fisher1/BOSS/PBanamufullzNS.mtx"]*)


(* ::Text:: *)
(*Analytical covariance with full redshift dependence full mu dependence,  w/o B*)


(* ::Input:: *)
(*Exporter[FPCnewmu+FPLnewmu,"Fisher1/BOSS/Panamufullz.mtx"]*)
(*Exporter[FPCnewmu+FPLnewmu+FBCnewmu+FBLnewmu,"Fisher1/BOSS/PBanamufullz.mtx"]*)


(* ::Subsubsection:: *)
(*DESI Exports*)


(* ::Text:: *)
(*Analytical covariance with full redshift dependence,  w/o B*)


(* ::Input:: *)
(*Exporter[FPD1+FPD2,FisherPath<>"DESI/Panafullz.mtx"]*)
(*Exporter[FPD1+FPD2+FBDESI1+FBDESI2,FisherPath<>"DESI/PBanafullz.mtx"]*)


(* ::Text:: *)
(*Analytical covariance with full redshift dependence, all multipoles,  w/o B*)


(* ::Input:: *)
(*Exporter[FPDESI1full+FPDESI2full,FisherPath<>"DESI/Panamufullz.mtx"]*)
(*Exporter[FPDESI1full+FPDESI2full+FBD1full+FBD2full,FisherPath<>"DESI/PBanamufullz.mtx"]*)


(* ::Text:: *)
(*Analytical covariance with full redshift dependence, all multipoles, 2loop,  w/o B*)


(* ::Input:: *)
(*Exporter[FPD1full2L+FPD2full2L,FisherPath<>"DESI/Panamufullz2L.mtx"]*)
(*Exporter[FPD1full2L+FPD2full2L+FBD1full+FBD2full,FisherPath<>"DESI/PBanamufullz2L.mtx"]*)


(* ::Text:: *)
(*Analytical covariance with full redshift dependence, all multipoles, Tree B w/o B*)


(* ::Input:: *)
(*Exporter[FPDESI1full+FPDESI2full+FBD1fullT+FBD2fullT,FisherPath<>"DESI/PBanamufullzT.mtx"]*)


(* ::Text:: *)
(*Analytical covariance with full redshift dependence, all multipoles, NS,  w/o B*)


(* ::Input:: *)
(*Exporter[FPD1fullNS+FPD2fullNS,FisherPath<>"DESI/PanamufullzNS.mtx"]*)
(*Exporter[FPD1fullNS+FPD2fullNS+FBD1fullNS+FBD2fullNS,FisherPath<>"DESI/PBanamufullzNS.mtx"]*)

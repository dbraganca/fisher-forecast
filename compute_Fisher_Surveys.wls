#!/usr/bin/env wolframscript
(* ::Package:: *)

PackageDirectory = NotebookDirectory[];
AppendTo[$Path,PackageDirectory];
SetDirectory[NotebookDirectory[]]


(* Load FisherFunctions package (takes 1min)*)
<<FisherFunctions`


(*tests to check it loaded well*)
(*Length[Bcovlist] == 319*)(*Yaniv's covariance*)
Length[Bcovlist] == 88
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
LOWZpk=Import["aux_vars/LOWZpk.m"];


(*pktest=LOWpk[[1,4]]*)


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
		  b2->0.6047,b3->-4.8320,b4->1.0142,cct->-5.4605,cr1->-1.0607,cr2->-0.3996,ce0->1.6922,ce1->1.1865,b12->0,
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


(* ::Subsubsection:: *)
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
Monitor[Do[
			If[(BderivsCMASS[[j]]==0)||(Length[BderivsCMASS[[j]]]==0),
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


COVCMASSB=GetFullzCovB[CMASStriangs,CMASSpk,CMASSspec,bfixBOSS,dkBCMASS];//AbsoluteTiming
(*{189.2065861`,Null}*)


COVLOWZB=GetFullzCovB[LOWZtriangs,LOWZpk,LOWZspec,bfixBOSS,dkBLOWZ];//AbsoluteTiming
(*{145.4908871`,Null}*)


Export["aux_vars/COVCMASSB.m",COVCMASSB]
Export["aux_vars/COVLOWZB.m",COVLOWZB]


COVCMASSBNS=GetFullzCovB[CMASStriangs,CMASSpk,CMASSspec,bfixBOSSNS,dkBCMASS];//AbsoluteTiming
(*{189.2065861`,Null}*)


COVLOWZBNS=GetFullzCovB[LOWZtriangs,LOWZpk,LOWZspec,bfixBOSSNS,dkBLOWZ];//AbsoluteTiming
(*{145.4908871`,Null}*)


Export["aux_vars/COVCMASSBNS.m",COVCMASSBNS]
Export["aux_vars/COVLOWZBNS.m",COVLOWZBNS]


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


(* ::Subsection:: *)
(*DESI*)


(* ::Subsubsection:: *)
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


Dg84=0.6535721061961522;
Dg122 =0.5528260185623496;


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
Export["covsks/ktabBase.m",ktabraw]


makeTriYaniv[kmin_,kmax_,dk_]:=Module[{kpairs,kpairstriangle,Alltriangs,ktab},
ktab = Table[kmin+i*dk,{i,0,Floor[(kmax-kmin)/dk]}];
kpairs = Flatten[Table[Table[Table[{ktab[[i]],ktab[[j]],ktab[[k]]},{i,1,j}],{j,1,k}],{k,1,Length[ktab]}],2];
kpairstriangle =Monitor[Table[If[kpairs[[i,1]]<kpairs[[i,3]]-kpairs[[i,2]],0,kpairs[[i]]],{i,1,Length[kpairs]}],i]//DeleteDuplicates;
Alltriangs=DeleteCases[kpairstriangle,0];
Table[Alltriangs[[i]],{i,1,Length[Alltriangs]}]//N;
Table[Alltriangs[[i]],{i,1,Length[Alltriangs]}]//Sort
];

makeTriDiogo[kmin_,kmax_,dk_]:=Module[{kpairs,kpairstriangle,triangles,ktab},
ktab = Table[k,{k,kmin,kmax,dk}];
triangles = Select[Subsets[ktab, {3}], 
    (#[[1]] + #[[2]] >= #[[3]]) &];
 Sort[triangles]
];



Basetri = makeTri[0.02,0.91,0.02];
BasetriDiogo = makeTriDiogo[0.02,0.91,0.02];
BasetriYaniv = makeTriYaniv[0.02,0.91,0.02];


Dimensions[Basetri]
Dimensions[BasetriDiogo]
Dimensions[BasetriYaniv]


keffbase = Import["covsks/baseTrikeff.mtx"];


Dimensions[Basetri]
Dimensions[keffbase]


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


COVCMASSPgen=GetFullzCovP[ktabraw,CMASSpk,CMASSspec,bfixBOSS,dkPD];
COVCMASSPgenNS=GetFullzCovP[ktabraw,CMASSpk,CMASSspec,bfixBOSSNS,dkPD];


\[Sigma]CMASS = Get\[Sigma][ktabraw, COVCMASSPgen];
\[Sigma]CMASSNS = Get\[Sigma][ktabraw, COVCMASSPgenNS];


COVD1Pkmm=GetFullzCovP[ktabraw,DESI1pk,D1spec,bfixDESI,dkPD];
COVD2Pkmm=GetFullzCovP[ktabraw,DESI2pk,D2spec,bfixDESI,dkPD];
COVD1PNSkmm=GetFullzCovP[ktabraw,DESI1pk,D1spec,bfixDESINS,dkPD];
COVD2PNSkmm=GetFullzCovP[ktabraw,DESI2pk,D2spec,bfixDESINS,dkPD];


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


(*COVD1B=GetFullzCovB[DESItriangs1,DESI1pk,D1spec,bfixDESI,dkBD];
COVD2B=GetFullzCovB[DESItriangs2,DESI2pk,D2spec,bfixDESI,dkBD];
COVD1BTree=GetFullzCovB[DESItriangs1T,DESI1pk,D1spec,bfixDESI,dkBD];
COVD2BTree=GetFullzCovB[DESItriangs2T,DESI2pk,D2spec,bfixDESI,dkBD];

Export["aux_vars/COVD1B.m",COVD1B]
Export["aux_vars/COVD2B.m",COVD2B]
Export["aux_vars/COVD1BTree.m",COVD1BTree]
Export["aux_vars/COVD2BTree.m",COVD2BTree]*)

COVD1B=Import["aux_vars/COVD1B.m"];
COVD2B=Import["aux_vars/COVD2B.m"];
COVD1BTree=Import["aux_vars/COVD1BTree.m"];
COVD2BTree=Import["aux_vars/COVD2BTree.m"];


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


PD1derivs[[1]]//Variables


PD1derivs2L = Table[P1Loopderivs[ktabDESI12L[[i]],DESI1pk]/.coeffixD1/.nbD1,{i,1,Length[ktabDESI12L]}];
PD2derivs2L = Table[P1Loopderivs[ktabDESI22L[[i]],DESI2pk]/.coeffixD2/.nbD2,{i,1,Length[ktabDESI12L]}];


(* ::Text:: *)
(*Bispectrum Derivatives*)


BD1derivsT=Monitor[Table[BTreeDerivs[keffDESI1T[[si]],DESI1pk]/.Dispatch[coeffixD1]/.nbD1,{si,1,Length[keffDESI1T]}],si];//AbsoluteTiming
BD2derivsT=Monitor[Table[BTreeDerivs[keffDESI2T[[si]],DESI2pk]/.Dispatch[coeffixD2]/.nbD2,{si,1,Length[keffDESI2T]}],si];//AbsoluteTiming


(*Calculate derivatives, one at a time instead of table*)
BD1derivs=ConstantArray[0,Length[keffDESI1]];
BD2derivs=ConstantArray[0,Length[keffDESI2]];


(*BD1derivs=Import["aux_vars/BD1derivs.m"];
BD2derivs=Import["aux_vars/BD2derivs.m"];*)


dispatcher = Dispatch[Join[coeffixD1,{nbD1}]];
Monitor[Do[If[(BD1derivs[[j]]==0)||(Length[BD1derivs[[j]]]==0),
			BD1derivs[[j]]=B1LoopDerivs[keffDESI1[[j]],DESI1pk]/.dispatcher];,
	{j,1,Length[keffDESI1]}],j];//AbsoluteTiming


Export["aux_vars/BD1derivs.m",BD1derivs]


dispatcher =Dispatch[Join[coeffixD2,{nbD2}]];
Monitor[Do[
	If[BD2derivs[[j]]==0,
	BD2derivs[[j]]=B1LoopDerivs[keffDESI2[[j]],DESI2pk]/.dispatcher];,
	{j,1,Length[keffDESI2]}],
	j];//AbsoluteTiming
(*{400.0957859`,Null}*)


Export["aux_vars/BD2derivs.m",BD2derivs]


(*checks*)
Dimensions[BD1derivs]
Dimensions[BD2derivs]


(*BD1derivs=Monitor[Table[B1Loopderivs[keffDESI1[[si]],DESI1pk]/.Dispatch[coeffixD1]/.nbD1,{si,1,Length[keffDESI1]}],si];//AbsoluteTiming
BD2derivs=Monitor[Table[B1Loopderivs[keffDESI2[[si]],DESI2pk]/.Dispatch[coeffixD2]/.nbD2,{si,1,Length[keffDESI2]}],si];//AbsoluteTiming*)



(* ::Subsubsection:: *)
(*Fishers *)


(* ::Text:: *)
(*Reduce P-derivatives to monopole and quadrupole*)


PD1monoquad= removelg2[PD1derivs]//Expand;
PD2monoquad= removelg2[PD2derivs]//Expand;


(* ::Text:: *)
(*Powerspectrum Fishers monoquad*)


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
(*Exports*)


Allparams = DeleteCases[AllParams,b12|b13|b14|b15];
GetFExp[Fisher_]:=SubFisher[Fisher+Priors,Allparams]
GetFExp4[Fisher_]:=SubFisher4sky[Fisher+Prior4sky,Allparams]
Exporter[Fisher_,name_]:=Export[name,GetFExp[Fisher]]
Exporter4[Fisher_,name_]:=Export[name,GetFExp4[Fisher]]



CombineMs

(* ::Package:: *)

(* ::Section:: *)
(*Preliminaries*)


SetDirectory[NotebookDirectory[]]
(*Checks if folder exists, and if not, creates it*)
CheckDir[a_String]:=With[{dirname=FileNameJoin[{NotebookDirectory[],a}]},Switch[FileType[dirname],
None,CreateDirectory[dirname] (*create dir*),
Directory,Print["Directory already exists."] (*do nothing*),
File,Print["File with same name already exists!!"] (*error!*)]];

(*paths*)
ctabpath="../GitHub/python-integer-power-project/3. Ctabs/";
LoopTablesPath="tabs/";
CovsPath="covsks/";
AuxPath="aux_vars/";
saveLoopCoefsPath = "saveloops/";

CheckDir[ctabpath]
CheckDir[LoopTablesPath]
CheckDir[AuxPath]
CheckDir[saveLoopCoefsPath]
CheckDir[saveLoopCoefsPath<>"B222/"]
CheckDir[saveLoopCoefsPath<>"B3211/"]
CheckDir[saveLoopCoefsPath<>"B3212/"]
CheckDir[saveLoopCoefsPath<>"B411/"]

baseTriEff=Round[Import[CovsPath<>"baseTrikeff.mtx"][[All,1;;3]],0.00001];
CMASStriEff=Round[Import[CovsPath<>"TriCMASSkeff.mtx"],0.00001];
LOWZtriEff=Round[Import[CovsPath<>"TriLOWZkeff.mtx"],0.00001];


(* ::Section::Closed:: *)
(*Utility functions*)


(* ::Input:: *)
(*(*function that outputs a list of exponents and coefficients of the following form: {exp q^-2, exp k1pq^-2, exp k2mq^-2, coef}*)*)
(*getCtab[expr_]:=*)
(*Map[{-Exponent[#,q^2],-Exponent[#,k1pq^2],-Exponent[#,k2mq^2],#/.{q-> 1,k1pq->1,k2mq->1}}&,List@@expr];*)
(*getCtab[1]={{0,0,0,1}};*)


(* ::Text:: *)
(*Tabtester checks that all terms of the expression were listed in terms of the exponents of the master integral and coefficients. It passed if all is 0*)


(* ::Input:: *)
(*Tabtester[final_,tab_]:=final-Sum[q^(-2tab[[i,1]]) k1pq^(-2tab[[i,2]]) k2mq^(-2tab[[i,3]])*tab[[i,4]],{i,1,Length[tab]}]//Expand;*)


(* ::Input:: *)
(*(*function that outputs a list of exponents and coefficients of the following form: {q,Plin,qDP,q2D2P,coef}*)*)
(*getExps[list_]:=Table[{Exponent[list[[i]],q],Exponent[list[[i]],Plin[q]],Exponent[list[[i]],qDP[q]],Exponent[list[[i]],q2D2P[q]],list[[i]]/.{q-> 1,qDP[_]->1,q2D2P[_]->1,Plin[_]->1}},{i,1,Length[list]}];*)


(* ::Input:: *)
(*(*Join terms with same exponents - equivalent to function compressor*)*)
(*(*Assuming ctab like {exp q^-2, exp k1pq^-2, exp k2mq^-2, coef}*)*)
(*compressor[ctab_]:=KeyValueMap[List,GroupBy[ctab,(#[[1;;3]]&)->Last,Total]]/.{{{a_,b_,c_},d_}:>{a,b,c,d}};*)
(*compressorB411[ctab_]:=KeyValueMap[List,GroupBy[ctab,(#[[1;;6]]&)->Last,Total]]/.{{{a_,b_,c_,a1_,b1_,c1_},d_}:>{a,b,c,a1,b1,c1,d}};*)


(* ::Input:: *)
(*(*Get all bias monomials that appear in the ctab. Strategy is to go through the ctab and add monomials to list if they are not there already*)*)
(*getMonomials[vars_,ctab_]:=Module[{testi,biaslist,coefrulesi,biaslisti},*)
(*biaslist={};*)
(*Do[testi=ctab[[i,-1]];*)
(*coefrulesi=CoefficientRules[testi,vars];*)
(*biaslisti=Sort[Times@@(vars^#)&/@coefrulesi[[All,1]]];*)
(*(*Check if new biaslist contains more terms than calculated so far*)*)
(*If[biaslisti!= biaslist&&Intersection[biaslist,biaslisti]!=biaslisti,*)
(*(*Print[i];*)*)
(*biaslist=Sort[DeleteDuplicates[Join[biaslist,biaslisti]]];],{i,1,Length[ctab]}];*)
(*biaslist*)
(*];*)


(* ::Input:: *)
(*getMonomialsSingle[vars_,expr_]:=Module[{biaslist,coefrules},*)
(*coefrules=CoefficientRules[expr,vars];*)
(*biaslist=Sort[Times@@(vars^#)&/@coefrules[[All,1]]];*)
(*biaslist*)
(*]*)


(* ::Input:: *)
(*(*creates Function that obtains the coefficient of each term in biaslist in each ctab coefficient "coef" *)*)
(*biaslisterfn[vars_,biaslist_]:=Module[{tab0,biasexps,biaslister},*)
(*biasexps=Flatten[CoefficientRules[biaslist,vars][[All,All,1]],1];*)
(*tab0=Table[_,Length[vars]];*)
(*biaslister[coef_]:=(biasexps/.CoefficientRules[coef,vars])/.{tab0:>0};*)
(*biaslister*)
(*];*)


(* ::Input:: *)
(*(*function that combines 2 ctabs - corresponds to the multiplication of the two quantities*)*)
(*JoinCtabs[ctab1_,ctab2_]:=Module[{expsmat,coefsmat,jointab},*)
(*expsmat=Outer[Plus,ctab1[[All,1;;3]],ctab2[[All,1;;3]],1];*)
(*coefsmat=Outer[Times,ctab1[[All,4]],ctab2[[All,4]],1];*)
(*jointab=Flatten[MapThread[MapThread[Append,{#1,#2}]&,{expsmat,coefsmat}],1];*)
(*compressor[jointab]*)
(*];*)


(* ::Input:: *)
(*(*angle replacement to calculate monopole*)*)
(*ysub={y->(k3^2-k1^2-k2^2)/(2 k1 k2)};*)
(*(*\[Mu]toangles3=Assuming[{k1q\[Element]Reals,qq>0,kk1>0,kk2>0,y>-1,y<1},{\[Mu]1\[Rule]Cos[\[Theta]],\[Mu]2\[Rule]y Cos[\[Theta]]+Sqrt[1-y^2] Sin[\[Theta]] Sin[\[Phi]],\[Mu]3\[Rule]-(((kk1+kk2 y) Cos[\[Theta]]+kk2 Sqrt[1-y^2] Sin[\[Theta]] Sin[\[Phi]])/kk3),\[Mu]0\[Rule]x Cos[\[Theta]]+Sqrt[1-x^2] Cos[\[Beta]-\[Phi]] Sin[\[Theta]]}(*/.x\[Rule]k1q/(kk1 qq)*)//FullSimplify];*)*)


(* ::Input:: *)
(*\[Mu]toangles3={\[Mu]1->Cos[\[Theta]],\[Mu]2->y Cos[\[Theta]]+Sqrt[1-y^2] Sin[\[Theta]] Sin[\[Phi]],\[Mu]3->-(((kk1+kk2 y) Cos[\[Theta]]+kk2 Sqrt[1-y^2] Sin[\[Theta]] Sin[\[Phi]])/kk3),\[Mu]0->x Cos[\[Theta]]+Sqrt[1-x^2] (Cos[\[Beta]] Cos[\[Phi]]Sin[\[Theta]]+Sin[\[Beta]] Sin[\[Phi]]Sin[\[Theta]])};*)


(* ::Text:: *)
(*Function defining the integrals over trig functions of \[Theta] and \[Phi]*)


(* ::Input:: *)
(*(*angint[a_,b_,c_,d_]=Assuming[{a\[GreaterEqual]0,b\[GreaterEqual]0,c\[GreaterEqual]0,d\[GreaterEqual]0,a\[Element]Integers,b\[Element]Integers,c\[Element]Integers,d\[Element]Integers},1/(4 \[Pi])Integrate[\[Mu]^a(Sqrt[1-\[Mu]^2])^bSin[\[Gamma]]^cCos[\[Gamma]]^d,{\[Mu],-1,1},{\[Gamma],0,2\[Pi]}]//Simplify]*)*)


(* ::Input:: *)
(*angint[a_,b_,c_,d_]=((1+(-1)^a) (1+(-1)^c) (1+(-1)^d) Gamma[(1+a)/2] Gamma[1+b/2] Gamma[(1+c)/2] Gamma[(1+d)/2])/(16 \[Pi] Gamma[1/2 (3+a+b)] Gamma[1/2 (2+c+d)]);*)


(* ::Text:: *)
(*Use the function to do replacements, incredibly faster than the integral*)


(* ::Input:: *)
(*rules4={Cos[\[Theta]]^a_. Sin[\[Theta]]^b_. Cos[\[Phi]]^c_. Sin[\[Phi]]^d_.->angint[a,b,c,d]};*)
(*rules3={Cos[\[Theta]]^a_. Sin[\[Theta]]^b_. Cos[\[Phi]]^c_.->angint[a,b,c,0],*)
(*Cos[\[Theta]]^a_. Sin[\[Theta]]^b_. Sin[\[Phi]]^d_.->angint[a,b,0,d],*)
(*Cos[\[Theta]]^a_. Cos[\[Phi]]^c_. Sin[\[Phi]]^d_.->angint[a,0,c,d],*)
(*Sin[\[Theta]]^b_. Cos[\[Phi]]^c_. Sin[\[Phi]]^d_.->angint[0,b,c,d]};*)
(*rules2={Cos[\[Theta]]^a_. Sin[\[Theta]]^b_.->angint[a,b,0,0],*)
(*Cos[\[Theta]]^a_. Cos[\[Phi]]^c_.->angint[a,0,c,0],*)
(* Sin[\[Theta]]^b_. Cos[\[Phi]]^c_.->angint[0,b,c,0],*)
(*Cos[\[Theta]]^a_. Sin[\[Phi]]^d_.->angint[a,0,0,d],*)
(* Sin[\[Theta]]^b_. Sin[\[Phi]]^d_.->angint[0,b,0,d],*)
(*Cos[\[Phi]]^c_. Sin[\[Phi]]^d_.->angint[0,0,c,d]};*)
(*rules1={Cos[\[Theta]]^a_.->angint[a,0,0,0],Sin[\[Theta]]^b_.->angint[0,b,0,0], Cos[\[Phi]]^c_.->angint[0,0,c,0],Sin[\[Phi]]^d_.->angint[0,0,0,d]};*)


(* ::Text:: *)
(*Given a kernel, lists coefficients of \[Mu]list*)


(* ::Input:: *)
(*(*Function that obtains the coefficient of each term in \[Mu]list in expression ker*)*)
(*\[Mu]lister[ker_,\[Mu]list_]:=Module[{vars,exps},*)
(*vars=Variables[\[Mu]list];*)
(*exps=Flatten[CoefficientRules[\[Mu]list,vars][[All,All,1]],1];*)
(*(exps/.CoefficientRules[ker,vars])/.{Table[_,{i,Length[vars]}]:>0}*)
(*]*)


(* ::Text:: *)
(*Given a list of monomials involving trig functions of \[Theta] and \[Phi], gets the monopole. One function is simplified, the other is expanded*)


(* ::Input:: *)
(*monosub[monos_]:=Monitor[Table[Assuming[{y>-1,y<1},Simplify@Total[monos[[i]]/.rules4/.rules3/.rules2/.rules1]],{i,1,Length[monos]}],i*100./Length[monos]];*)
(*monosub2[monos_]:=Monitor[Table[Assuming[{y>-1,y<1},Expand@Total[monos[[i]]/.rules4/.rules3/.rules2/.rules1]],{i,1,Length[monos]}],i*100./Length[monos]];*)


(* ::Input:: *)
(*(*function to calculate monopole*)*)
(*GetMonopolenoSimp[ker_,\[Mu]list_]:=Module[{mlist,monos},*)
(*monos=MonomialList[#,{Cos[\[Theta]],Sin[\[Theta]],Cos[\[Phi]],Sin[\[Phi]]}]&/@(\[Mu]list/.\[Mu]toangles3/.{kk1->k1,kk2->k2,kk3->k3}/.ysub);*)
(*mlist=\[Mu]lister[ker,\[Mu]list];*)
(*Print["Residual is"];*)
(*Print[(ker-mlist.\[Mu]list//Expand)];*)
(*(mlist.monosub2[monos])];*)
(**)


(* ::Input:: *)
(*GetMonopoleList[coefList_,\[Mu]list_]:=Module[*)
(*{monos,monosubs,GetMonopoleSingle},*)
(*monos=MonomialList[#,{Cos[\[Theta]],Sin[\[Theta]],Cos[\[Phi]],Sin[\[Phi]]}]&/@(\[Mu]list/.\[Mu]toangles3/.{kk1->k1,kk2->k2,kk3->k3}/.ysub);*)
(*monosubs=monosub2[monos];*)
(*GetMonopoleSingle[coef_]:=\[Mu]lister[coef,\[Mu]list].monosubs;*)
(*GetMonopoleSingle/@coefList*)
(*];*)


(* ::Section::Closed:: *)
(*Halo kernels in redshift space*)


(* ::Text:: *)
(*Final output of this section are the halo kernels in redshift space*)


(* ::Subsection::Closed:: *)
(*Perturbation Theory Kernels*)


(* ::Input:: *)
(*perm12=Permutations[{q1,q2}];*)
(*perm123=Permutations[{q1,q2,q3}];*)
(*perm1234=Permutations[{q1,q2,q3,q4}];*)


(* ::Input:: *)
(*Clear[Fn,Gn,n,m];*)
(*Fn[n_,v_]:=Fn[n,v]=If[n==1,1,Sum[Gn[m,v[[1;;m]]]/((2n+3)(n-1)) ((2n+1)al[Total[v[[1;;m]]],Total[v[[m+1;;n]]]]Fn[n-m,v[[m+1;;n]]]+2be[Total[v[[1;;m]]],Total[v[[m+1;;n]]]]Gn[n-m,v[[m+1;;n]]]),{m,1,n-1}]];*)
(*Gn[n_,v_]:=Gn[n,v]= If[n==1,1,Sum[Gn[m,v[[1;;m]]]/((2n+3)(n-1)) (3al[Total[v[[1;;m]]],Total[v[[m+1;;n]]]]Fn[n-m,v[[m+1;;n]]]+2n be[Total[v[[1;;m]]],Total[v[[m+1;;n]]]]Gn[n-m,v[[m+1;;n]]]),{m,1,n-1}]];*)


(* ::Input:: *)
(*F2[q1_,q2_]=Simplify[Sum[1/Length[perm12]Fn[2,{perm12[[i,1]],perm12[[i,2]]}],{i,1,Length[perm12]}]];*)
(*F3[q1_,q2_,q3_]=Simplify[Sum[1/Length[perm123]Fn[3,{perm123[[i,1]],perm123[[i,2]],perm123[[i,3]]}],{i,1,Length[perm123]}]];*)
(*F4[q1_,q2_,q3_,q4_]=Simplify[Sum[1/Length[perm1234]Fn[4,{perm1234[[i,1]],perm1234[[i,2]],perm1234[[i,3]],perm1234[[i,4]]}],{i,1,Length[perm1234]}]]; *)
(*G2[q1_,q2_]=Simplify[Sum[1/Length[perm12]Gn[2,{perm12[[i,1]],perm12[[i,2]]}],{i,1,Length[perm12]}]];*)
(*G3[q1_,q2_,q3_]=Simplify[Sum[1/Length[perm123]Gn[3,{perm123[[i,1]],perm123[[i,2]],perm123[[i,3]]}],{i,1,Length[perm123]}]];*)
(*G4[q1_,q2_,q3_,q4_]=Simplify[Sum[1/Length[perm1234]Gn[4,{perm1234[[i,1]],perm1234[[i,2]],perm1234[[i,3]],perm1234[[i,4]]}],{i,1,Length[perm1234]}]]; *)


(* ::Input:: *)
(*alf[k1_,k2_]=(k1+k2).k1/(k1.k1);*)
(*bef[k1_,k2_]=(k1+k2).(k1+k2)(k1.k2)/2/(k1.k1)/(k2.k2) ;*)
(*alfm[k1_,k2_]=1+1/2 (mag[k1+k2]^2-mag[k1]^2-mag[k2]^2)/mag[k1]^2;*)
(*befm[k1_,k2_]=(mag[k1+k2]^2 (mag[k1+k2]^2-mag[k1]^2-mag[k2]^2))/(4mag[k1]^2 mag[k2]^2);*)


(* ::Subsection::Closed:: *)
(*Define Subscript[C, i] Kernels and bare bias coefficients*)


(* ::Text:: *)
(*initialize bare bias Subscript[c, i] values*)


Clear[c1,c2,c3,c4]


(* ::Input:: *)
(*c4 = Table[0,{i,1,15}];*)


(* ::Input:: *)
(*c4[[1]]=Subscript["\!\(\*SuperscriptBox[OverscriptBox[\(c\), \(^\)], \((4)\)]\)","\[Delta],1"];*)
(*c4[[2]]=Subscript["\!\(\*SuperscriptBox[OverscriptBox[\(c\), \(^\)], \((4)\)]\)","\[Delta],2"];*)
(*c4[[3]]=Subscript["\!\(\*SuperscriptBox[OverscriptBox[\(c\), \(^\)], \((4)\)]\)","\[Delta],3"];*)
(*c4[[4]]=Subscript["\!\(\*SuperscriptBox[OverscriptBox[\(c\), \(^\)], \((4)\)]\)","\[Delta],4"];*)
(*c4[[5]]=Subscript["\!\(\*SuperscriptBox[OverscriptBox[\(c\), \(^\)], \((4)\)]\)","\!\(\*SuperscriptBox[\(\[Delta]\), \(2\)]\),1"];*)
(*c4[[6]]=Subscript["\!\(\*SuperscriptBox[OverscriptBox[\(c\), \(^\)], \((4)\)]\)","\!\(\*SuperscriptBox[\(\[Delta]\), \(2\)]\),2"];*)
(*c4[[7]]=Subscript["\!\(\*SuperscriptBox[OverscriptBox[\(c\), \(^\)], \((4)\)]\)","\!\(\*SuperscriptBox[\(\[Delta]\), \(2\)]\),3"];*)
(*c4[[8]]=Subscript["\!\(\*SuperscriptBox[OverscriptBox[\(c\), \(^\)], \((4)\)]\)","\!\(\*SuperscriptBox[\(r\), \(2\)]\),2"];*)
(*c4[[9]]=Subscript["\!\(\*SuperscriptBox[OverscriptBox[\(c\), \(^\)], \((4)\)]\)","\!\(\*SuperscriptBox[\(r\), \(2\)]\),3"];*)
(*c4[[10]]=Subscript["\!\(\*SuperscriptBox[OverscriptBox[\(c\), \(^\)], \((4)\)]\)","\!\(\*SuperscriptBox[\(\[Delta]\), \(3\)]\),1"];*)
(*c4[[11]]=Subscript["\!\(\*SuperscriptBox[OverscriptBox[\(c\), \(^\)], \((4)\)]\)","\!\(\*SuperscriptBox[\(\[Delta]\), \(3\)]\),2"];*)
(*c4[[12]]=Subscript["\!\(\*SuperscriptBox[OverscriptBox[\(c\), \(^\)], \((4)\)]\)","\!\(\*SuperscriptBox[\(r\), \(3\)]\),2"];*)
(*c4[[13]]=Subscript["\!\(\*SuperscriptBox[OverscriptBox[\(c\), \(^\)], \((4)\)]\)","\!\(\*SuperscriptBox[\(r\), \(2\)]\)\[Delta],2"];*)
(*c4[[14]]=Subscript["\!\(\*SuperscriptBox[OverscriptBox[\(c\), \(^\)], \((4)\)]\)","\!\(\*SuperscriptBox[\(\[Delta]\), \(4\)]\),1"];*)
(*c4[[15]]=Subscript["\!\(\*SuperscriptBox[OverscriptBox[\(c\), \(^\)], \((4)\)]\)","\!\(\*SuperscriptBox[\(\[Delta]r\), \(3\)]\),1"];*)


(* ::Input:: *)
(*c3=c4//.{"\!\(\*SuperscriptBox[OverscriptBox[\(c\), \(^\)], \((4)\)]\)"-> "\!\(\*SuperscriptBox[OverscriptBox[\(c\), \(^\)], \((3)\)]\)"};*)
(*c2=c4//.{"\!\(\*SuperscriptBox[OverscriptBox[\(c\), \(^\)], \((4)\)]\)"-> "\!\(\*SuperscriptBox[OverscriptBox[\(c\), \(^\)], \((2)\)]\)"};*)
(*c1=c4//.{"\!\(\*SuperscriptBox[OverscriptBox[\(c\), \(^\)], \((4)\)]\)"-> "\!\(\*SuperscriptBox[OverscriptBox[\(c\), \(^\)], \((1)\)]\)"};*)
(*adelta = c4//.{"\!\(\*SuperscriptBox[OverscriptBox[\(c\), \(^\)], \((4)\)]\)"-> "\!\(\*SuperscriptBox[OverscriptBox[\(c\), \(~\)], \(A\)]\)"};*)


(* ::Input:: *)
(*btheta={1,2,3,4,-1,-(3/2),-2,0,0,1,4/3,0,0,-1,0};*)


(* ::Text:: *)
(*Expressions for Subscript[C, i] kernels *)


Do[
	c1[[i]][q1_]:=0;
	c2[[i]][q1_,q2_]:=0;
	c3[[i]][q1_,q2_,q3_]:=0,
	{i,1,15}
];


(* ::Input:: *)
(*c1[[1]][q1_]:=1;*)
(*c2[[1]][q1_,q2_]:=-((mag[q1]^2+mag[q2]^2-mag[q1+q2]^2)/(2 mag[q1]^2));*)
(*c2[[2]][q1_,q2_]:=((1+2 F2[q1,q2]) mag[q1]^2+mag[q2]^2-mag[q1+q2]^2)/(2 mag[q1]^2);*)
(*c2[[5]][q1_,q2_]:=1;*)
(*c3[[1]][q1_,q2_,q3_]:=1/8 (-1/mag[q2+q3]^2 2 G2[q2,q3] (2 mag[q1]^2+mag[q2]^2-mag[q1+q2]^2+mag[q3]^2-mag[q1+q3]^2)+((mag[q1]^2+mag[q2]^2-mag[q1+q2]^2) (mag[q1]^2+mag[q2]^2+2 mag[q3]^2-mag[q1+q3]^2-mag[q2+q3]^2))/(mag[q2]^2 mag[q3]^2));*)
(*c3[[2]][q1_,q2_,q3_]:=-(((mag[q1]^2+(1+2 F2[q1,q2]) mag[q2]^2-mag[q1+q2]^2) (mag[q1]^2+mag[q2]^2+2 mag[q3]^2-mag[q1+q3]^2-mag[q2+q3]^2))/(4 mag[q2]^2 mag[q3]^2));*)
(*c3[[3]][q1_,q2_,q3_]:=1/8 ((2 G2[q2,q3] (mag[q1]^2+mag[q2+q3]^2-mag[q1+q2+q3]^2))/mag[q2+q3]^2+1/(mag[q2]^2 mag[q3]^2) (mag[q1]^2 (mag[q1+q2]^2+mag[q3]^2-mag[q1+q2+q3]^2)-mag[q1+q2]^2 (mag[q1+q2]^2+mag[q3]^2-mag[q1+q2+q3]^2)+mag[q2]^2 ((1+4 F2[q1,q2]) mag[q1+q2]^2+(1+4 F2[q1,q2]+8 F3[q1,q2,q3]) mag[q3]^2-(1+4 F2[q1,q2]) mag[q1+q2+q3]^2)));c3[[5]][q1_,q2_,q3_]:=-((mag[q2]^2+mag[q3]^2-mag[q2+q3]^2)/mag[q3]^2);*)
(*c3[[6]][q1_,q2_,q3_]:=2 F2[q1,q2]+(mag[q2]^2+mag[q3]^2-mag[q2+q3]^2)/mag[q3]^2;*)
(*c3[[8]][q1_,q2_,q3_]:=1/(4 mag[q1]^2 mag[q2]^2 mag[q1+q2]^2 mag[q3]^2) (mag[q1]^4 mag[q1+q2]^2 (mag[q2]^2+mag[q3]^2-mag[q2+q3]^2)+(-mag[q2]^2 mag[q1+q2]+mag[q1+q2]^3)^2 (mag[q2]^2+mag[q3]^2-mag[q2+q3]^2)+2 mag[q1]^2 (mag[q2]^4 mag[q1+q2]^2+mag[q1+q2]^4 (-mag[q3]^2+mag[q2+q3]^2)+mag[q2]^2 ((-1+F2[q1,q2]) mag[q1+q2]^4+F2[q1,q2] (mag[q3]^2-mag[q1+q2+q3]^2)^2+mag[q1+q2]^2 ((1+2 F2[q1,q2]) mag[q3]^2-mag[q2+q3]^2-2 F2[q1,q2] mag[q1+q2+q3]^2))));*)
(*c3[[10]][q1_,q2_,q3_]:=1;*)
(**)


(* ::Input:: *)
(*c4[[1]][q1_,q2_,q3_,q4_]:=1/48 (1/(mag[q2]^2 mag[q3+q4]^2) 2 G2[q3,q4] (mag[q1]^2+mag[q2]^2-mag[q1+q2]^2) (3 mag[q1]^2+2 mag[q2]^2+5 mag[q3+q4]^2-3 mag[q1+q3+q4]^2-2 mag[q2+q3+q4]^2)+1/(mag[q2]^2 mag[q3]^2 mag[q2+q3]^2 mag[q4]^2) (-mag[q1]^6 mag[q2+q3]^2-mag[q2]^6 mag[q2+q3]^2+mag[q2]^4 mag[q2+q3]^2 (mag[q1+q2]^2-2 mag[q3]^2+mag[q2+q3]^2-2 mag[q4]^2+mag[q2+q4]^2+mag[q3+q4]^2)+mag[q1]^4 mag[q2+q3]^2 (-mag[q2]^2+mag[q1+q2]^2-4 mag[q3]^2+mag[q1+q3]^2-4 mag[q4]^2+mag[q1+q4]^2+3 mag[q3+q4]^2)+mag[q1+q2]^2 mag[q2+q3]^2 (4 mag[q3]^4+mag[q2+q3]^2 (-2 mag[q4]^2+mag[q2+q4]^2+mag[q3+q4]^2)+mag[q1+q3]^2 (-4 mag[q4]^2+mag[q1+q4]^2+3 mag[q3+q4]^2)-mag[q3]^2 (3 mag[q1+q3]^2+mag[q2+q3]^2-6 mag[q4]^2+mag[q1+q4]^2+mag[q2+q4]^2+4 mag[q3+q4]^2))+mag[q1]^2 (-mag[q2]^4 mag[q2+q3]^2+mag[q2+q3]^2 (-4 mag[q3]^4+4 mag[q1+q3]^2 mag[q4]^2+2 mag[q2+q3]^2 mag[q4]^2-mag[q1+q3]^2 mag[q1+q4]^2-mag[q2+q3]^2 mag[q2+q4]^2-3 mag[q1+q3]^2 mag[q3+q4]^2-mag[q2+q3]^2 mag[q3+q4]^2+mag[q1+q2]^2 (4 mag[q3]^2-mag[q1+q3]^2+4 mag[q4]^2-mag[q1+q4]^2-3 mag[q3+q4]^2)+mag[q3]^2 (3 mag[q1+q3]^2+mag[q2+q3]^2-6 mag[q4]^2+mag[q1+q4]^2+mag[q2+q4]^2+4 mag[q3+q4]^2))+mag[q2]^2 (mag[q2+q3]^2 (mag[q1+q3]^2+mag[q2+q3]^2-6 mag[q4]^2+mag[q1+q4]^2+mag[q2+q4]^2+4 mag[q3+q4]^2)+2 mag[q3]^2 ((-3+G2[q2,q3]) mag[q2+q3]^2+G2[q2,q3] (mag[q4]^2-mag[q2+q3+q4]^2))))+mag[q2]^2 (-4 mag[q3]^4 mag[q2+q3]^2+mag[q1+q3]^2 mag[q2+q3]^2 (4 mag[q4]^2-mag[q1+q4]^2-3 mag[q3+q4]^2)+mag[q2+q3]^4 (2 mag[q4]^2-mag[q2+q4]^2-mag[q3+q4]^2)-mag[q1+q2]^2 mag[q2+q3]^2 (-2 mag[q3]^2+mag[q2+q3]^2-2 mag[q4]^2+mag[q2+q4]^2+mag[q3+q4]^2)+mag[q3]^2 (3 mag[q1+q3]^2 mag[q2+q3]^2+(1+2 G2[q2,q3]) mag[q2+q3]^4+2 G2[q2,q3] mag[q1+q2+q3]^2 (-mag[q4]^2+mag[q2+q3+q4]^2)+mag[q2+q3]^2 (-6 mag[q4]^2+mag[q1+q4]^2+mag[q2+q4]^2+4 mag[q3+q4]^2-2 G2[q2,q3] (mag[q1+q2+q3]^2-mag[q4]^2+mag[q2+q3+q4]^2)))))-(8 G3[q2,q3,q4] (mag[q1]^2+mag[q2+q3+q4]^2-mag[q1+q2+q3+q4]^2))/mag[q2+q3+q4]^2);*)
(*c4[[2]][q1_,q2_,q3_,q4_]:=-(1/(16 mag[q2]^2))(1/(mag[q3]^2 mag[q4]^2) (-mag[q1]^6-mag[q2]^6+mag[q2]^4 (mag[q1+q2]^2-2 mag[q3]^2+mag[q2+q3]^2-2 mag[q4]^2+mag[q2+q4]^2+mag[q3+q4]^2)+mag[q1]^4 (-mag[q2]^2+mag[q1+q2]^2-4 mag[q3]^2+mag[q1+q3]^2-4 mag[q4]^2+mag[q1+q4]^2+3 mag[q3+q4]^2)-mag[q2]^2 (4 mag[q3]^4-3 mag[q3]^2 mag[q1+q3]^2-mag[q3]^2 mag[q2+q3]^2+6 mag[q3]^2 mag[q4]^2-4 mag[q1+q3]^2 mag[q4]^2-2 mag[q2+q3]^2 mag[q4]^2-mag[q3]^2 mag[q1+q4]^2+mag[q1+q3]^2 mag[q1+q4]^2-mag[q3]^2 mag[q2+q4]^2+mag[q2+q3]^2 mag[q2+q4]^2-4 mag[q3]^2 mag[q3+q4]^2+3 mag[q1+q3]^2 mag[q3+q4]^2+mag[q2+q3]^2 mag[q3+q4]^2+2 F2[q1,q2] (mag[q1+q2]^2+mag[q3]^2-mag[q1+q2+q3]^2) (mag[q1+q2]^2+mag[q3]^2+2 mag[q4]^2-mag[q1+q2+q4]^2-mag[q3+q4]^2)+mag[q1+q2]^2 (-2 mag[q3]^2+mag[q2+q3]^2-2 mag[q4]^2+mag[q2+q4]^2+mag[q3+q4]^2))+mag[q1]^2 (-mag[q2]^4-4 mag[q3]^4+3 mag[q3]^2 mag[q1+q3]^2+mag[q3]^2 mag[q2+q3]^2-6 mag[q3]^2 mag[q4]^2+4 mag[q1+q3]^2 mag[q4]^2+2 mag[q2+q3]^2 mag[q4]^2+mag[q3]^2 mag[q1+q4]^2-mag[q1+q3]^2 mag[q1+q4]^2+mag[q3]^2 mag[q2+q4]^2-mag[q2+q3]^2 mag[q2+q4]^2+4 mag[q3]^2 mag[q3+q4]^2-3 mag[q1+q3]^2 mag[q3+q4]^2-mag[q2+q3]^2 mag[q3+q4]^2+mag[q1+q2]^2 (4 mag[q3]^2-mag[q1+q3]^2+4 mag[q4]^2-mag[q1+q4]^2-3 mag[q3+q4]^2)+mag[q2]^2 (-6 mag[q3]^2+mag[q1+q3]^2+mag[q2+q3]^2-6 mag[q4]^2+mag[q1+q4]^2+mag[q2+q4]^2+4 mag[q3+q4]^2))+mag[q1+q2]^2 (4 mag[q3]^4+mag[q2+q3]^2 (-2 mag[q4]^2+mag[q2+q4]^2+mag[q3+q4]^2)+mag[q1+q3]^2 (-4 mag[q4]^2+mag[q1+q4]^2+3 mag[q3+q4]^2)-mag[q3]^2 (3 mag[q1+q3]^2+mag[q2+q3]^2-6 mag[q4]^2+mag[q1+q4]^2+mag[q2+q4]^2+4 mag[q3+q4]^2)))+1/mag[q3+q4]^2 2 G2[q3,q4] (mag[q1]^4+mag[q2]^4+mag[q1]^2 (2 mag[q2]^2-mag[q1+q2]^2+2 mag[q3+q4]^2-mag[q1+q3+q4]^2-mag[q2+q3+q4]^2)+mag[q1+q2]^2 (-2 mag[q3+q4]^2+mag[q1+q3+q4]^2+mag[q2+q3+q4]^2)+mag[q2]^2 ((-1+2 F2[q1,q2]) mag[q1+q2]^2+2 (1+F2[q1,q2]) mag[q3+q4]^2-mag[q1+q3+q4]^2-mag[q2+q3+q4]^2-2 F2[q1,q2] mag[q1+q2+q3+q4]^2)));*)
(*c4[[3]][q1_,q2_,q3_,q4_]:=1/(16 mag[q2]^2) (-1/mag[q3+q4]^2 2 G2[q3,q4] (mag[q1]^2+mag[q2]^2-mag[q1+q2]^2) (mag[q1]^2+mag[q3+q4]^2-mag[q1+q3+q4]^2)+1/(mag[q3]^2 mag[q2+q3]^2 mag[q4]^2) (-mag[q1]^6 mag[q2+q3]^2-mag[q2]^6 mag[q2+q3]^2+mag[q2]^4 mag[q2+q3]^2 (mag[q1+q2]^2-2 mag[q3]^2+mag[q2+q3]^2-2 mag[q4]^2+mag[q2+q4]^2+mag[q3+q4]^2)+mag[q1]^4 mag[q2+q3]^2 (-mag[q2]^2+mag[q1+q2]^2-4 mag[q3]^2+mag[q1+q3]^2-4 mag[q4]^2+mag[q1+q4]^2+3 mag[q3+q4]^2)+mag[q1+q2]^2 mag[q2+q3]^2 (4 mag[q3]^4+mag[q2+q3]^2 (-2 mag[q4]^2+mag[q2+q4]^2+mag[q3+q4]^2)+mag[q1+q3]^2 (-4 mag[q4]^2+mag[q1+q4]^2+3 mag[q3+q4]^2)-mag[q3]^2 (3 mag[q1+q3]^2+mag[q2+q3]^2-6 mag[q4]^2+mag[q1+q4]^2+mag[q2+q4]^2+4 mag[q3+q4]^2))+mag[q1]^2 (-mag[q2]^4 mag[q2+q3]^2+mag[q2+q3]^2 (-4 mag[q3]^4+4 mag[q1+q3]^2 mag[q4]^2+2 mag[q2+q3]^2 mag[q4]^2-mag[q1+q3]^2 mag[q1+q4]^2-mag[q2+q3]^2 mag[q2+q4]^2-3 mag[q1+q3]^2 mag[q3+q4]^2-mag[q2+q3]^2 mag[q3+q4]^2+mag[q1+q2]^2 (4 mag[q3]^2-mag[q1+q3]^2+4 mag[q4]^2-mag[q1+q4]^2-3 mag[q3+q4]^2)+mag[q3]^2 (3 mag[q1+q3]^2+mag[q2+q3]^2-6 mag[q4]^2+mag[q1+q4]^2+mag[q2+q4]^2+4 mag[q3+q4]^2))+mag[q2]^2 (mag[q2+q3]^2 (mag[q1+q3]^2+mag[q2+q3]^2-6 mag[q4]^2+mag[q1+q4]^2+mag[q2+q4]^2+4 mag[q3+q4]^2)-2 mag[q3]^2 ((3+G2[q2,q3]) mag[q2+q3]^2+G2[q2,q3] (mag[q4]^2-mag[q2+q3+q4]^2))))-mag[q2]^2 (4 mag[q3]^4 mag[q2+q3]^2-3 mag[q3]^2 mag[q1+q3]^2 mag[q2+q3]^2-mag[q3]^2 mag[q2+q3]^4+2 G2[q2,q3] mag[q3]^2 mag[q2+q3]^4+8 F3[q1,q2,q3] mag[q3]^2 mag[q2+q3]^2 mag[q1+q2+q3]^2-2 G2[q2,q3] mag[q3]^2 mag[q2+q3]^2 mag[q1+q2+q3]^2+6 mag[q3]^2 mag[q2+q3]^2 mag[q4]^2+8 F3[q1,q2,q3] mag[q3]^2 mag[q2+q3]^2 mag[q4]^2+2 G2[q2,q3] mag[q3]^2 mag[q2+q3]^2 mag[q4]^2-4 mag[q1+q3]^2 mag[q2+q3]^2 mag[q4]^2-2 mag[q2+q3]^4 mag[q4]^2-2 G2[q2,q3] mag[q3]^2 mag[q1+q2+q3]^2 mag[q4]^2-mag[q3]^2 mag[q2+q3]^2 mag[q1+q4]^2+mag[q1+q3]^2 mag[q2+q3]^2 mag[q1+q4]^2-mag[q3]^2 mag[q2+q3]^2 mag[q2+q4]^2+mag[q2+q3]^4 mag[q2+q4]^2-4 mag[q3]^2 mag[q2+q3]^2 mag[q3+q4]^2+3 mag[q1+q3]^2 mag[q2+q3]^2 mag[q3+q4]^2+mag[q2+q3]^4 mag[q3+q4]^2+4 F2[q1,q2] mag[q2+q3]^2 (mag[q1+q2]^2+mag[q3]^2-mag[q1+q2+q3]^2) (mag[q1+q2]^2+mag[q3]^2+2 mag[q4]^2-mag[q1+q2+q4]^2-mag[q3+q4]^2)+mag[q1+q2]^2 mag[q2+q3]^2 (-2 mag[q3]^2+mag[q2+q3]^2-2 mag[q4]^2+mag[q2+q4]^2+mag[q3+q4]^2)-2 G2[q2,q3] mag[q3]^2 mag[q2+q3]^2 mag[q2+q3+q4]^2+2 G2[q2,q3] mag[q3]^2 mag[q1+q2+q3]^2 mag[q2+q3+q4]^2-8 F3[q1,q2,q3] mag[q3]^2 mag[q2+q3]^2 mag[q1+q2+q3+q4]^2)));*)
(*c4[[4]][q1_,q2_,q3_,q4_]:=1/48 ((8 G3[q2,q3,q4] (mag[q1]^2+mag[q2+q3+q4]^2-mag[q1+q2+q3+q4]^2))/mag[q2+q3+q4]^2+1/(mag[q2]^2 mag[q3+q4]^2) 2 G2[q3,q4] (3 mag[q1]^4+mag[q2]^4+mag[q1]^2 (4 mag[q2]^2-3 mag[q1+q2]^2+4 mag[q3+q4]^2-3 mag[q1+q3+q4]^2-mag[q2+q3+q4]^2)+mag[q1+q2]^2 (-4 mag[q3+q4]^2+3 mag[q1+q3+q4]^2+mag[q2+q3+q4]^2)+mag[q2]^2 ((-1+6 F2[q1,q2]) mag[q1+q2]^2+(4+6 F2[q1,q2]) mag[q3+q4]^2-3 mag[q1+q3+q4]^2-mag[q2+q3+q4]^2-6 F2[q1,q2] mag[q1+q2+q3+q4]^2))+1/(mag[q2]^2 mag[q3]^2 mag[q2+q3]^2 mag[q4]^2) (mag[q1]^6 mag[q2+q3]^2+mag[q2]^6 mag[q2+q3]^2+mag[q1]^4 mag[q2+q3]^2 (mag[q2]^2-mag[q1+q2]^2+4 mag[q3]^2-mag[q1+q3]^2+4 mag[q4]^2-mag[q1+q4]^2-3 mag[q3+q4]^2)-mag[q2]^4 mag[q2+q3]^2 (mag[q1+q2]^2-2 mag[q3]^2+mag[q2+q3]^2-2 mag[q4]^2+mag[q2+q4]^2+mag[q3+q4]^2)+mag[q1+q2]^2 mag[q2+q3]^2 (-4 mag[q3]^4+mag[q1+q3]^2 (4 mag[q4]^2-mag[q1+q4]^2-3 mag[q3+q4]^2)+mag[q2+q3]^2 (2 mag[q4]^2-mag[q2+q4]^2-mag[q3+q4]^2)+mag[q3]^2 (3 mag[q1+q3]^2+mag[q2+q3]^2-6 mag[q4]^2+mag[q1+q4]^2+mag[q2+q4]^2+4 mag[q3+q4]^2))+mag[q1]^2 (mag[q2]^4 mag[q2+q3]^2+mag[q2+q3]^2 (4 mag[q3]^4-4 mag[q1+q3]^2 mag[q4]^2-2 mag[q2+q3]^2 mag[q4]^2+mag[q1+q3]^2 mag[q1+q4]^2+mag[q2+q3]^2 mag[q2+q4]^2+3 mag[q1+q3]^2 mag[q3+q4]^2+mag[q2+q3]^2 mag[q3+q4]^2+mag[q1+q2]^2 (-4 mag[q3]^2+mag[q1+q3]^2-4 mag[q4]^2+mag[q1+q4]^2+3 mag[q3+q4]^2)-mag[q3]^2 (3 mag[q1+q3]^2+mag[q2+q3]^2-6 mag[q4]^2+mag[q1+q4]^2+mag[q2+q4]^2+4 mag[q3+q4]^2))+mag[q2]^2 (-mag[q2+q3]^2 (mag[q1+q3]^2+mag[q2+q3]^2-6 mag[q4]^2+mag[q1+q4]^2+mag[q2+q4]^2+4 mag[q3+q4]^2)+mag[q3]^2 ((6+4 G2[q2,q3]) mag[q2+q3]^2+4 G2[q2,q3] (mag[q4]^2-mag[q2+q3+q4]^2))))+mag[q2]^2 (4 mag[q3]^4 mag[q2+q3]^2-3 mag[q3]^2 mag[q1+q3]^2 mag[q2+q3]^2-mag[q3]^2 mag[q2+q3]^4+4 G2[q2,q3] mag[q3]^2 mag[q2+q3]^4+24 F3[q1,q2,q3] mag[q3]^2 mag[q2+q3]^2 mag[q1+q2+q3]^2-4 G2[q2,q3] mag[q3]^2 mag[q2+q3]^2 mag[q1+q2+q3]^2+6 mag[q3]^2 mag[q2+q3]^2 mag[q4]^2+24 F3[q1,q2,q3] mag[q3]^2 mag[q2+q3]^2 mag[q4]^2+48 F4[q1,q2,q3,q4] mag[q3]^2 mag[q2+q3]^2 mag[q4]^2+4 G2[q2,q3] mag[q3]^2 mag[q2+q3]^2 mag[q4]^2-4 mag[q1+q3]^2 mag[q2+q3]^2 mag[q4]^2-2 mag[q2+q3]^4 mag[q4]^2-4 G2[q2,q3] mag[q3]^2 mag[q1+q2+q3]^2 mag[q4]^2-mag[q3]^2 mag[q2+q3]^2 mag[q1+q4]^2+mag[q1+q3]^2 mag[q2+q3]^2 mag[q1+q4]^2-mag[q3]^2 mag[q2+q3]^2 mag[q2+q4]^2+mag[q2+q3]^4 mag[q2+q4]^2-4 mag[q3]^2 mag[q2+q3]^2 mag[q3+q4]^2+3 mag[q1+q3]^2 mag[q2+q3]^2 mag[q3+q4]^2+mag[q2+q3]^4 mag[q3+q4]^2+6 F2[q1,q2] mag[q2+q3]^2 (mag[q1+q2]^2+mag[q3]^2-mag[q1+q2+q3]^2) (mag[q1+q2]^2+mag[q3]^2+2 mag[q4]^2-mag[q1+q2+q4]^2-mag[q3+q4]^2)+mag[q1+q2]^2 mag[q2+q3]^2 (-2 mag[q3]^2+mag[q2+q3]^2-2 mag[q4]^2+mag[q2+q4]^2+mag[q3+q4]^2)-4 G2[q2,q3] mag[q3]^2 mag[q2+q3]^2 mag[q2+q3+q4]^2+4 G2[q2,q3] mag[q3]^2 mag[q1+q2+q3]^2 mag[q2+q3+q4]^2-24 F3[q1,q2,q3] mag[q3]^2 mag[q2+q3]^2 mag[q1+q2+q3+q4]^2)));*)
(*c4[[5]][q1_,q2_,q3_,q4_]:=1/4 (1/(mag[q3]^2 mag[q4]^2) (mag[q2]^4+mag[q3]^4-mag[q3]^2 mag[q2+q3]^2+3 mag[q3]^2 mag[q4]^2-mag[q1+q3]^2 mag[q4]^2-2 mag[q2+q3]^2 mag[q4]^2-2 mag[q3]^2 mag[q2+q4]^2+mag[q1+q3]^2 mag[q2+q4]^2+mag[q2+q3]^2 mag[q2+q4]^2+mag[q1]^2 (mag[q2]^2+mag[q4]^2-mag[q2+q4]^2)-mag[q3]^2 mag[q3+q4]^2+mag[q2+q3]^2 mag[q3+q4]^2+mag[q2]^2 (3 mag[q3]^2-mag[q1+q3]^2-mag[q2+q3]^2+2 mag[q4]^2-mag[q2+q4]^2-mag[q3+q4]^2))-1/mag[q3+q4]^2 2 G2[q3,q4] (mag[q2]^2+mag[q3+q4]^2-mag[q2+q3+q4]^2));*)
(*c4[[6]][q1_,q2_,q3_,q4_]:=-(1/(2 mag[q3]^2 mag[q4]^2))(mag[q2]^4+mag[q3]^4+2 F2[q1,q2] mag[q3]^4-mag[q3]^2 mag[q2+q3]^2+2 F2[q2,q3] mag[q3]^2 mag[q2+q3]^2+3 mag[q3]^2 mag[q4]^2+2 F2[q1,q2] mag[q3]^2 mag[q4]^2+2 F2[q2,q3] mag[q3]^2 mag[q4]^2-mag[q1+q3]^2 mag[q4]^2-2 mag[q2+q3]^2 mag[q4]^2-2 mag[q3]^2 mag[q2+q4]^2+mag[q1+q3]^2 mag[q2+q4]^2+mag[q2+q3]^2 mag[q2+q4]^2+mag[q1]^2 (mag[q2]^2+mag[q4]^2-mag[q2+q4]^2)-mag[q3]^2 mag[q3+q4]^2-2 F2[q1,q2] mag[q3]^2 mag[q3+q4]^2+mag[q2+q3]^2 mag[q3+q4]^2-mag[q2]^2 (-3 mag[q3]^2+mag[q1+q3]^2+mag[q2+q3]^2-2 mag[q4]^2+mag[q2+q4]^2+mag[q3+q4]^2)-2 F2[q2,q3] mag[q3]^2 mag[q2+q3+q4]^2);*)
(*c4[[7]][q1_,q2_,q3_,q4_]:=1/4 (1/mag[q4]^2 4 F2[q1,q2] (mag[q3]^2+(1+F2[q3,q4]) mag[q4]^2-mag[q3+q4]^2)+1/mag[q3+q4]^2 2 G2[q3,q4] (mag[q2]^2+mag[q3+q4]^2-mag[q2+q3+q4]^2)+1/(mag[q3]^2 mag[q4]^2) (mag[q2]^4+mag[q3]^4-mag[q3]^2 mag[q2+q3]^2+4 F2[q2,q3] mag[q3]^2 mag[q2+q3]^2+3 mag[q3]^2 mag[q4]^2+4 F2[q2,q3] mag[q3]^2 mag[q4]^2+8 F3[q1,q2,q3] mag[q3]^2 mag[q4]^2-mag[q1+q3]^2 mag[q4]^2-2 mag[q2+q3]^2 mag[q4]^2-2 mag[q3]^2 mag[q2+q4]^2+mag[q1+q3]^2 mag[q2+q4]^2+mag[q2+q3]^2 mag[q2+q4]^2+mag[q1]^2 (mag[q2]^2+mag[q4]^2-mag[q2+q4]^2)-mag[q3]^2 mag[q3+q4]^2+mag[q2+q3]^2 mag[q3+q4]^2+mag[q2]^2 (3 mag[q3]^2-mag[q1+q3]^2-mag[q2+q3]^2+2 mag[q4]^2-mag[q2+q4]^2-mag[q3+q4]^2)-4 F2[q2,q3] mag[q3]^2 mag[q2+q3+q4]^2));*)
(*c4[[8]][q1_,q2_,q3_,q4_]:=(-mag[q1]^6 mag[q1+q2]^2 mag[q2+q3]^2 (mag[q2]^2+mag[q4]^2-mag[q2+q4]^2)+mag[q1+q2]^2 (-mag[q2]^8 mag[q2+q3]^2+mag[q2]^6 mag[q2+q3]^2 (2 mag[q1+q2]^2-3 mag[q3]^2+mag[q1+q3]^2+mag[q2+q3]^2-2 mag[q4]^2+mag[q2+q4]^2+mag[q3+q4]^2)-mag[q2]^4 mag[q2+q3]^2 (mag[q1+q2]^4+mag[q3]^4-mag[q1+q3]^2 mag[q4]^2-2 mag[q2+q3]^2 mag[q4]^2+mag[q1+q3]^2 mag[q2+q4]^2+mag[q2+q3]^2 mag[q2+q4]^2+mag[q2+q3]^2 mag[q3+q4]^2+2 mag[q1+q2]^2 (-3 mag[q3]^2+mag[q1+q3]^2+mag[q2+q3]^2-2 mag[q4]^2+mag[q2+q4]^2+mag[q3+q4]^2)-mag[q3]^2 (mag[q2+q3]^2-3 mag[q4]^2+2 mag[q2+q4]^2+mag[q3+q4]^2))+mag[q1+q2]^4 mag[q2+q3]^2 (-mag[q3]^4+mag[q1+q3]^2 (mag[q4]^2-mag[q2+q4]^2)+mag[q2+q3]^2 (2 mag[q4]^2-mag[q2+q4]^2-mag[q3+q4]^2)+mag[q3]^2 (mag[q2+q3]^2-3 mag[q4]^2+2 mag[q2+q4]^2+mag[q3+q4]^2))+mag[q2]^2 (mag[q1+q2]^4 mag[q2+q3]^2 (-3 mag[q3]^2+mag[q1+q3]^2+mag[q2+q3]^2-2 mag[q4]^2+mag[q2+q4]^2+mag[q3+q4]^2)+2 mag[q1+q2]^2 mag[q2+q3]^2 (mag[q3]^4+mag[q1+q3]^2 (-mag[q4]^2+mag[q2+q4]^2)+mag[q2+q3]^2 (-2 mag[q4]^2+mag[q2+q4]^2+mag[q3+q4]^2)-mag[q3]^2 (mag[q2+q3]^2-3 mag[q4]^2+2 mag[q2+q4]^2+mag[q3+q4]^2))-2 F2[q2,q3] mag[q3]^2 (mag[q2+q3]^2-mag[q1+q2+q3]^2)^2 (mag[q2+q3]^2+mag[q4]^2-mag[q2+q3+q4]^2)))+mag[q1]^4 mag[q1+q2]^2 (-3 mag[q2]^4 mag[q2+q3]^2+mag[q2+q3]^2 (-mag[q3]^4+mag[q1+q3]^2 mag[q4]^2+2 mag[q2+q3]^2 mag[q4]^2-mag[q1+q3]^2 mag[q2+q4]^2-mag[q2+q3]^2 mag[q2+q4]^2+2 mag[q1+q2]^2 (mag[q4]^2-mag[q2+q4]^2)-mag[q2+q3]^2 mag[q3+q4]^2+mag[q3]^2 (mag[q2+q3]^2-3 mag[q4]^2+2 mag[q2+q4]^2+mag[q3+q4]^2))+mag[q2]^2 (2 mag[q1+q2]^2 mag[q2+q3]^2+mag[q2+q3]^2 (mag[q1+q3]^2+mag[q2+q3]^2-4 mag[q4]^2+3 mag[q2+q4]^2+mag[q3+q4]^2)-mag[q3]^2 ((3+2 F2[q2,q3]) mag[q2+q3]^2+2 F2[q2,q3] (mag[q4]^2-mag[q2+q3+q4]^2))))-mag[q1]^2 (3 mag[q2]^6 mag[q1+q2]^2 mag[q2+q3]^2-mag[q2]^4 mag[q1+q2]^2 mag[q2+q3]^2 (4 mag[q1+q2]^2-6 mag[q3]^2+2 mag[q1+q3]^2+2 mag[q2+q3]^2-5 mag[q4]^2+3 mag[q2+q4]^2+2 mag[q3+q4]^2)+mag[q1+q2]^4 mag[q2+q3]^2 (-2 mag[q3]^4+2 mag[q1+q3]^2 mag[q4]^2+4 mag[q2+q3]^2 mag[q4]^2-2 mag[q1+q3]^2 mag[q2+q4]^2-2 mag[q2+q3]^2 mag[q2+q4]^2+mag[q1+q2]^2 (mag[q4]^2-mag[q2+q4]^2)-2 mag[q2+q3]^2 mag[q3+q4]^2+2 mag[q3]^2 (mag[q2+q3]^2-3 mag[q4]^2+2 mag[q2+q4]^2+mag[q3+q4]^2))+mag[q2]^2 (mag[q1+q2]^6 mag[q2+q3]^2+2 F2[q1,q2] mag[q2+q3]^2 (mag[q3]^2-mag[q1+q2+q3]^2)^2 (mag[q3]^2+mag[q4]^2-mag[q3+q4]^2)+2 mag[q1+q2]^4 mag[q2+q3]^2 ((-3+F2[q1,q2]) mag[q3]^2+mag[q1+q3]^2+mag[q2+q3]^2-3 mag[q4]^2+F2[q1,q2] mag[q4]^2+2 mag[q2+q4]^2+mag[q3+q4]^2-F2[q1,q2] mag[q3+q4]^2)+2 mag[q1+q2]^2 ((1+2 F2[q1,q2]) mag[q3]^4 mag[q2+q3]^2+mag[q2+q3]^2 (mag[q1+q3]^2 (-mag[q4]^2+mag[q2+q4]^2)+2 F2[q1,q2] mag[q1+q2+q3]^2 (-mag[q4]^2+mag[q3+q4]^2)+mag[q2+q3]^2 (-2 mag[q4]^2+mag[q2+q4]^2+mag[q3+q4]^2))+mag[q3]^2 ((-1+2 F2[q2,q3]) mag[q2+q3]^4+2 F2[q2,q3] mag[q1+q2+q3]^2 (-mag[q4]^2+mag[q2+q3+q4]^2)-mag[q2+q3]^2 (-3 mag[q4]^2+2 mag[q2+q4]^2+mag[q3+q4]^2+2 F2[q1,q2] (mag[q1+q2+q3]^2-mag[q4]^2+mag[q3+q4]^2)+2 F2[q2,q3] (mag[q1+q2+q3]^2-mag[q4]^2+mag[q2+q3+q4]^2)))))))/(8 mag[q1]^2 mag[q2]^2 mag[q1+q2]^2 mag[q3]^2 mag[q2+q3]^2 mag[q4]^2);*)
(*c4[[9]][q1_,q2_,q3_,q4_]:=1/16 ((2 G2[q3,q4] (mag[q1]^2+mag[q2]^2-mag[q1+q2]^2)^2 (mag[q2]^2+mag[q3+q4]^2-mag[q2+q3+q4]^2))/(mag[q1]^2 mag[q2]^2 mag[q3+q4]^2)+1/mag[q1+q2]^2 4 F2[q1,q2] (((mag[q1+q2]^2+mag[q3]^2-mag[q1+q2+q3]^2)^2 (mag[q3]^2+mag[q4]^2-mag[q3+q4]^2))/(mag[q3]^2 mag[q4]^2)+1/mag[q3+q4]^2 F2[q3,q4] (mag[q1+q2]^2+mag[q3+q4]^2-mag[q1+q2+q3+q4]^2)^2)+(mag[q1]^6 mag[q2+q3]^2 mag[q1+q2+q3]^2 (mag[q2]^2+mag[q4]^2-mag[q2+q4]^2)+mag[q1+q2+q3]^2 (mag[q2]^8 mag[q2+q3]^2-mag[q2]^6 mag[q2+q3]^2 (2 mag[q1+q2]^2-3 mag[q3]^2+mag[q1+q3]^2+mag[q2+q3]^2-2 mag[q4]^2+mag[q2+q4]^2+mag[q3+q4]^2)+mag[q1+q2]^4 mag[q2+q3]^2 (mag[q3]^4+mag[q1+q3]^2 (-mag[q4]^2+mag[q2+q4]^2)+mag[q2+q3]^2 (-2 mag[q4]^2+mag[q2+q4]^2+mag[q3+q4]^2)-mag[q3]^2 (mag[q2+q3]^2-3 mag[q4]^2+2 mag[q2+q4]^2+mag[q3+q4]^2))+mag[q2]^4 mag[q2+q3]^2 (mag[q1+q2]^4+mag[q3]^4-mag[q1+q3]^2 mag[q4]^2-2 mag[q2+q3]^2 mag[q4]^2+mag[q1+q3]^2 mag[q2+q4]^2+mag[q2+q3]^2 mag[q2+q4]^2+mag[q2+q3]^2 mag[q3+q4]^2+2 mag[q1+q2]^2 (-3 mag[q3]^2+mag[q1+q3]^2+mag[q2+q3]^2-2 mag[q4]^2+mag[q2+q4]^2+mag[q3+q4]^2)-mag[q3]^2 (mag[q2+q3]^2-3 mag[q4]^2+2 mag[q2+q4]^2+mag[q3+q4]^2))+mag[q2]^2 (-mag[q1+q2]^4 mag[q2+q3]^2 (-3 mag[q3]^2+mag[q1+q3]^2+mag[q2+q3]^2-2 mag[q4]^2+mag[q2+q4]^2+mag[q3+q4]^2)-2 mag[q1+q2]^2 mag[q2+q3]^2 (mag[q3]^4+mag[q1+q3]^2 (-mag[q4]^2+mag[q2+q4]^2)+mag[q2+q3]^2 (-2 mag[q4]^2+mag[q2+q4]^2+mag[q3+q4]^2)-mag[q3]^2 (mag[q2+q3]^2-3 mag[q4]^2+2 mag[q2+q4]^2+mag[q3+q4]^2))+4 F2[q2,q3] mag[q3]^2 (mag[q2+q3]^2-mag[q1+q2+q3]^2)^2 (mag[q2+q3]^2+mag[q4]^2-mag[q2+q3+q4]^2)))+mag[q1]^4 mag[q1+q2+q3]^2 (3 mag[q2]^4 mag[q2+q3]^2+mag[q2+q3]^2 (mag[q3]^4-mag[q1+q3]^2 mag[q4]^2-2 mag[q2+q3]^2 mag[q4]^2+mag[q1+q3]^2 mag[q2+q4]^2+mag[q2+q3]^2 mag[q2+q4]^2-2 mag[q1+q2]^2 (mag[q4]^2-mag[q2+q4]^2)+mag[q2+q3]^2 mag[q3+q4]^2-mag[q3]^2 (mag[q2+q3]^2-3 mag[q4]^2+2 mag[q2+q4]^2+mag[q3+q4]^2))-mag[q2]^2 (2 mag[q1+q2]^2 mag[q2+q3]^2+mag[q2+q3]^2 (mag[q1+q3]^2+mag[q2+q3]^2-4 mag[q4]^2+3 mag[q2+q4]^2+mag[q3+q4]^2)-mag[q3]^2 ((3+4 F2[q2,q3]) mag[q2+q3]^2+4 F2[q2,q3] (mag[q4]^2-mag[q2+q3+q4]^2))))+mag[q1]^2 (3 mag[q2]^6 mag[q2+q3]^2 mag[q1+q2+q3]^2-mag[q2]^4 mag[q2+q3]^2 mag[q1+q2+q3]^2 (4 mag[q1+q2]^2-6 mag[q3]^2+2 mag[q1+q3]^2+2 mag[q2+q3]^2-5 mag[q4]^2+3 mag[q2+q4]^2+2 mag[q3+q4]^2)+mag[q1+q2]^2 mag[q2+q3]^2 mag[q1+q2+q3]^2 (-2 mag[q3]^4+2 mag[q1+q3]^2 mag[q4]^2+4 mag[q2+q3]^2 mag[q4]^2-2 mag[q1+q3]^2 mag[q2+q4]^2-2 mag[q2+q3]^2 mag[q2+q4]^2+mag[q1+q2]^2 (mag[q4]^2-mag[q2+q4]^2)-2 mag[q2+q3]^2 mag[q3+q4]^2+2 mag[q3]^2 (mag[q2+q3]^2-3 mag[q4]^2+2 mag[q2+q4]^2+mag[q3+q4]^2))+mag[q2]^2 (mag[q1+q2]^4 mag[q2+q3]^2 mag[q1+q2+q3]^2+2 mag[q1+q2]^2 mag[q2+q3]^2 mag[q1+q2+q3]^2 (-3 mag[q3]^2+mag[q1+q3]^2+mag[q2+q3]^2-3 mag[q4]^2+2 mag[q2+q4]^2+mag[q3+q4]^2)+2 (mag[q3]^4 mag[q2+q3]^2 mag[q1+q2+q3]^2+mag[q2+q3]^2 mag[q1+q2+q3]^2 (mag[q1+q3]^2 (-mag[q4]^2+mag[q2+q4]^2)+mag[q2+q3]^2 (-2 mag[q4]^2+mag[q2+q4]^2+mag[q3+q4]^2))+mag[q3]^2 ((-1+4 F2[q2,q3]) mag[q2+q3]^4 mag[q1+q2+q3]^2+4 F2[q2,q3] mag[q1+q2+q3]^4 (-mag[q4]^2+mag[q2+q3+q4]^2)+mag[q2+q3]^2 (mag[q1+q2+q3]^2 (3 mag[q4]^2-2 mag[q2+q4]^2-mag[q3+q4]^2)-4 F2[q2,q3] mag[q1+q2+q3]^2 (mag[q1+q2+q3]^2-mag[q4]^2+mag[q2+q3+q4]^2)+4 F3[q1,q2,q3] (mag[q1+q2+q3]^2+mag[q4]^2-mag[q1+q2+q3+q4]^2)^2))))))/(mag[q1]^2 mag[q2]^2 mag[q3]^2 mag[q2+q3]^2 mag[q1+q2+q3]^2 mag[q4]^2));*)
(*c4[[10]][q1_,q2_,q3_,q4_]:=-((3 (mag[q3]^2+mag[q4]^2-mag[q3+q4]^2))/(2 mag[q4]^2));*)
(*c4[[11]][q1_,q2_,q3_,q4_]:=(3 (mag[q3]^2+(1+2 F2[q1,q2]) mag[q4]^2-mag[q3+q4]^2))/(2 mag[q4]^2);*)
(*c4[[12]][q1_,q2_,q3_,q4_]:=-((3 (mag[q1]^4 mag[q1+q2]^2 (mag[q2]^2+mag[q3]^2-mag[q2+q3]^2)+mag[q1+q2]^2 (mag[q2]^2-mag[q1+q2]^2) (mag[q3]^2-mag[q1+q3]^2) (mag[q2]^2+mag[q3]^2-mag[q2+q3]^2)+mag[q1]^2 (mag[q2]^4 mag[q1+q2]^2-mag[q1+q2]^2 (mag[q1+q2]^2-mag[q3]^2+mag[q1+q3]^2) (mag[q3]^2-mag[q2+q3]^2)+mag[q2]^2 ((-1+2 F2[q1,q2]) mag[q1+q2]^4+2 F2[q1,q2] (mag[q3]^2-mag[q1+q2+q3]^2) (mag[q4]^2-mag[q1+q2+q4]^2)+mag[q1+q2]^2 (2 (1+F2[q1,q2]) mag[q3]^2-mag[q1+q3]^2-mag[q2+q3]^2-2 F2[q1,q2] mag[q1+q2+q3]^2+2 F2[q1,q2] mag[q4]^2-2 F2[q1,q2] mag[q1+q2+q4]^2)))) (mag[q3]^2+mag[q4]^2-mag[q3+q4]^2))/(16 mag[q1]^2 mag[q2]^2 mag[q1+q2]^2 mag[q3]^2 mag[q4]^2));*)
(*c4[[13]][q1_,q2_,q3_,q4_]:=1/(8 mag[q1]^2 mag[q2]^2 mag[q1+q2]^2 mag[q3]^2 mag[q4]^2) (mag[q1]^4 mag[q1+q2]^2 mag[q3]^2 (mag[q3]^2+(1+2 F2[q3,q4]) mag[q4]^2-mag[q3+q4]^2)+mag[q1+q2]^2 (mag[q2]^2-mag[q1+q2]^2)^2 mag[q3]^2 (mag[q3]^2+(1+2 F2[q3,q4]) mag[q4]^2-mag[q3+q4]^2)+2 mag[q1]^2 (mag[q2]^4 mag[q1+q2]^2 (mag[q3]^2+mag[q4]^2-mag[q3+q4]^2)+mag[q1+q2]^2 (mag[q3]^2-mag[q2+q3]^2)^2 (mag[q3]^2+mag[q4]^2-mag[q3+q4]^2)-mag[q1+q2]^4 mag[q3]^2 (mag[q3]^2+(1+2 F2[q3,q4]) mag[q4]^2-mag[q3+q4]^2)+mag[q2]^2 (2 F2[q1,q2] mag[q1+q2]^4 mag[q4]^2+2 F2[q1,q2] (mag[q3]^2-mag[q1+q2+q3]^2)^2 mag[q4]^2+mag[q1+q2]^2 (3 mag[q3]^4-2 mag[q2+q3]^2 mag[q4]^2-4 F2[q1,q2] mag[q1+q2+q3]^2 mag[q4]^2+2 mag[q2+q3]^2 mag[q3+q4]^2+mag[q3]^2 (-2 mag[q2+q3]^2+(3+4 F2[q1,q2]+2 F2[q3,q4]) mag[q4]^2-3 mag[q3+q4]^2)))));*)
(*c4[[14]][q1_,q2_,q3_,q4_]:=1;*)
(*c4[[15]][q1_,q2_,q3_,q4_]:=-(((mag[q1]^2+mag[q2]^2-mag[q1+q2]^2) (mag[q1]^2+mag[q3]^2-mag[q1+q3]^2) (mag[q2]^2+mag[q3]^2-mag[q2+q3]^2))/(8 mag[q1]^2 mag[q2]^2 mag[q3]^2));*)


(* ::Subsection::Closed:: *)
(*Replacements, coordinates*)


(* ::Text:: *)
(*Coordinate choices. Here we take k1 along the third axis. q is the integration vector. z is the line of sight, angles are chosen such that k1.z = \[Mu]*)


(* ::Input:: *)
(*Clear[zhat];*)
(*coordchoice={k1->{0,0,kk1},k2->kk2{0,Sqrt[1-y^2],y},k3->{0,-kk2 Sqrt[1-y^2],-kk1-kk2 y},q->qq{Sqrt[1-x^2]Cos[\[Beta]],Sqrt[1-x^2] Sin[\[Beta]],x}(*q\[Rule]qq{Sqrt[1-x^2]\[Alpha],Sqrt[1-x^2] Sqrt[1-\[Alpha]^2],x}*)};zhat1=zhat->{Sqrt[1-\[Mu]^2]Cos[\[Phi]],Sqrt[1-\[Mu]^2]Sin[\[Phi]], \[Mu]};zhat2=zhat->{Sin[\[Theta]] Cos[\[Phi]],Sin[\[Theta]]Sin[\[Phi]],Cos[\[Theta]]};*)
(*ysub={y->(k3^2-k1^2-k2^2)/(2 k1 k2)};*)


(* ::Text:: *)
(*Definitions of magnitudes, which are needed to expand the kernels and PS*)


(* ::Input:: *)
(*kmq[k1] = k1mq;*)
(*kmq[k2] = k2mq;*)
(*kmq[k3] = k3mq;*)
(*kpq[k1] = k1pq;*)
(*kpq[k2] = k2pq;*)
(*kpq[k3] = k3pq;*)


(* ::Input:: *)
(*magreps ={mag[-q]->q,mag[q]->q,mag[k1]->k1,mag[k2]->k2,mag[k3]->k3,mag[-k1]->k1,mag[-k2]->k2,mag[-k1-k2]->k3,mag[k1+k2]->k3,mag[-k1+q]->k1mq,mag[k1-q]->k1mq,mag[-k1-q]->k1pq,mag[k1+q]->k1pq,mag[-k2+q]->k2mq,mag[k2-q]->k2mq,mag[-k2-q]->k2pq,mag[k2+q]->k2pq,mag[-k1-k2-q]->k3mq,mag[k1+k2+q]->k3mq,mag[-k1-k2+q]->k3pq,mag[k1+k2-q]->k3pq,mag[-k3]->k3,mag[-k1-k3]->k2,mag[k1+k3]->k2,mag[-k3+q]->k3mq,mag[k3-q]->k3mq,mag[-k3-q]->k3pq,mag[k3+q]->k3pq,mag[-k1-k3-q]->k2mq,mag[k1+k3+q]->k2mq,mag[-k1-k3+q]->k2pq,mag[k1+k3-q]->k2pq,mag[-k2-k3]->k1,mag[k2+k3]->k1,mag[-k2-k3-q]->k1mq,mag[k2+k3+q]->k1mq,mag[-k2-k3+q]->k1pq,mag[k2+k3-q]->k1pq,*)
(*mag[-q1]->q1,mag[q1]->q1,mag[-q-q1]->qpq1,mag[q-q1]->qmq1,mag[-q+q1]->qmq1,mag[q+q1]->qpq1};*)


(* ::Input:: *)
(*magrevreps  = {k1mq-> Sqrt[k1^2+q^2-2 k1 q x],k1pq-> Sqrt[k1^2+q^2+2 k1 q x],k2mq-> Sqrt[k2^2+q^2-2 k2 q (x y+Sqrt[1-x^2] Sqrt[1-y^2] Sin[\[Beta]])],k2pq-> Sqrt[k2^2+q^2+2 k2 q (x y+Sqrt[1-x^2] Sqrt[1-y^2] Sin[\[Beta]])],k3mq-> \[Sqrt](k3^2+q^2+2k1 q x+2k2 q x y+2k2 q Sqrt[1-x^2] Sqrt[1-y^2] Sin[\[Beta]]),k3pq-> \[Sqrt](k3^2+q^2-2k1 q x-2k2 q x y-2k2 q Sqrt[1-x^2] Sqrt[1-y^2] Sin[\[Beta]]),*)
(*qmq1-> Sqrt[q^2+q1^2-2\[Nu] q q1 ],qpq1-> Sqrt[q^2+q1^2+2\[Nu] q q1 ]};*)


(* ::Text:: *)
(*define line of sight angle. Note that for now we keep \[Mu]3 and intermediary also define \[Mu]0, which we will substitute later.*)
(*Try to integrate over \[Mu] and \[Phi] to get the monopole!*)


(* ::Input:: *)
(*kzsub = {Subscript[k1, z]-> k1 \[Mu]1,Subscript[k2, z]-> k2 \[Mu]2,Subscript[q, z]-> q \[Mu]0,Subscript[k3, z]->k3 \[Mu]3};*)


(* ::Input:: *)
(*\[Mu]toangles1=Assuming[{y>-1,y<1,\[Mu]>-1,\[Mu]<1,\[Phi]>0,\[Phi]<2\[Pi],x>-1,x<1},{\[Mu]1->k1.zhat/kk1,\[Mu]2->k2.zhat/kk2,\[Mu]3->k3.zhat/kk3,\[Mu]0->q.zhat/qq}/.coordchoice/.zhat1//FullSimplify];*)


(* ::Input:: *)
(*angles1to\[Mu]={\[Mu]->\[Mu]1, Sin[\[Phi]]->(\[Mu]2-y \[Mu]1)/(Sqrt[1-y^2] Sqrt[1-\[Mu]1^2])};*)


(* ::Input:: *)
(*\[Mu]toangles2=Assuming[{y>-1,y<1,\[Mu]>-1,\[Mu]<1,\[Phi]>0,\[Phi]<2\[Pi],x>-1,x<1},{\[Mu]1->k1.zhat/kk1,\[Mu]2->k2.zhat/kk2,\[Mu]3->k3.zhat/kk3,\[Mu]0->q.zhat/qq}/.coordchoice/.zhat2];*)


(* ::Input:: *)
(*angles2to\[Mu]={Cos[\[Theta]]->\[Mu]1,Sin[\[Theta]]->Sqrt[1-\[Mu]1^2],Sin[\[Phi]]->(\[Mu]2-y \[Mu]1)/(Sqrt[1-y^2] Sqrt[1-\[Mu]1^2]),Cos[\[Phi]]->Sqrt[1-(\[Mu]2-y \[Mu]1)^2/((1-y^2)(1-\[Mu]1^2))]};*)


(* ::Text:: *)
(*define symmetrization over two and three variables*)


(* ::Input:: *)
(*q2sym[f_,q1_,q2_]:=1/2 (f[q1,q2]+f[q2,q1]);*)
(*q3sym[f_,q1_,q2_,q3_]:=Total[f@@@Permutations[{q1,q2,q3}]/3!];*)
(*q4sym[f_,q1_,q2_,q3_,q4_]:=Total[f@@@Permutations[{q1,q2,q3,q4}]/4!]*)


(* ::Text:: *)
(*define final bias coefficients*)


(* ::Input:: *)
(*biases = {b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15};*)
(*tb=Thread[adelta->biases];*)
(*tbrev  =Thread[biases->adelta];*)


(* ::Text:: *)
(*Projections on matter and velocity (as a special species of biased tracer)*)


(* ::Input:: *)
(*submat=Thread[biases->{1,1,1,1,0,0,0,0,0,0,0,0,0,0,0}];*)
(*subvel=Thread[biases->{1,2,3,4,-1,-3/2,-2,0,0,1,0,0,0,0,0}];*)


(* ::Text:: *)
(*Here we only keep \[Mu]0, and leave \[Mu]1 and \[Mu]2, \[Mu]3 to the coefficients of each FFTLog integral*)


(* ::Input:: *)
(*\[Mu]listfinal  = Table[\[Mu]0^i,{i,0,6}];*)
(*\[Mu]listfinalno1 = Delete[\[Mu]listfinal,1];*)


(* ::Text:: *)
(*Replacements*)
(*Notice that we modify \[Mu]0 as well!*)
(*the replacements q-> -q, q-> k1pq and q-> k2mq  have the following coding:*)
(*k1mk2=|k1-k2+q|,*)
(*k12pq = |2k1+q|,*)
(*k22mq = |2k2-q|,*)
(*kspeck1 = |2k1+k2+q|,*)
(*kspeck2 = |k1+2k2+q|*)


(* ::Input:: *)
(*preqtok1pq= Thread[{q,k1pq, k1mq, k2pq, k2mq, k3pq, k3mq,\[Mu]0}-> {k1pqt,k12pq,qt,k3mqt,k1mk2,k2mqt,kspeck1,\[Mu]0t}];*)
(*preqtok2mq = Thread[{q,k1pq, k1mq, k2pq, k2mq, k3pq, k3mq,\[Mu]0}-> {k2mqt,k3pqt,k1mk2,k22mq,qt,k1pqt,kspeck2,\[Mu]0t}];*)
(*qtok1pq = Thread[ {k1pqt,k12pq,qt,k3mqt,k1mk2,k2mqt,kspeck1,\[Mu]0t}-> {k1pq,Sqrt[2 k1^2+2 k1pq^2-q^2],q,k3mq,Sqrt[k1^2+k1pq^2+k2^2+k2mq^2-k3^2-q^2],k2mq,Sqrt[2 k1pq^2-k2mq^2+2 k3^2],(q \[Mu]0+k1 \[Mu]1)/k1pq}];*)
(*qtok2mq = Thread[ {k2mqt,k3pqt,k1mk2,k22mq,qt,k1pqt,kspeck2,\[Mu]0t}->{k2mq,k3pq,Sqrt[k1^2+k1pq^2+k2^2+k2mq^2-k3^2-q^2],Sqrt[2 k2^2+2 k2mq^2-q^2],q,k1pq,Sqrt[-2 k1^2+k1pq^2+4 k2^2-2 k2mq^2+2 k3^2+2 q^2],(-q \[Mu]0+k2 \[Mu]2)/k2mq} ];*)
(*restrep = Thread[{k1mq,k2pq,k3mq,k3pq}-> {Sqrt[2 k1^2-k1pq^2+2 q^2],Sqrt[2 k2^2-k2mq^2+2 q^2],Sqrt[-k1^2+k1pq^2+k2^2-k2mq^2+k3^2+q^2],Sqrt[k1^2-k1pq^2-k2^2+k2mq^2+k3^2+q^2]}];*)


(* ::Input:: *)
(*qrev1={k1mq->k1pqt,k1pq->k1mqt,k2mq->k2pqt,k2pq->k2mqt,k3mq->k3pqt,k3pq->k3mqt,\[Mu]0-> \[Mu]0t};*)
(*qrev2={k1mqt->k1mq,k1pqt->k1pq,k2mqt->k2mq,k2pqt->k2pq,k3mqt->k3mq,k3pqt->k3pq,\[Mu]0t-> -\[Mu]0};*)


(* ::Text:: *)
(*Tests*)


(* ::Text:: *)
(*Tabtester checks that all terms of the expression were listed in terms of the exponents of the master integral and coefficients. It passed if all is 0*)


(* ::Input:: *)
(*Tabtester[final_,tab_]:=final-Sum[q^(-2tab[[i,1]]) k1pq^(-2tab[[i,2]]) k2mq^(-2tab[[i,3]])*tab[[i,4]],{i,1,Length[tab]}]//Expand;*)


(* ::Subsection::Closed:: *)
(*Halo Kernels (real space)*)


(* ::Input:: *)
(*Clear[K1,K2,K3,K4]*)


(* ::Input:: *)
(*K1[aa_]:=aa[[1]];*)
(*K2[q1_,q2_,aa_]:=K2[q1,q2,aa] =((Table[c2[[i]][q1,q2],{i,1,15}]/.{al->alfm,be->befm}/.{mag[-q]->mag[q]}/.magreps//Expand)/.mag[0]-> 0//Simplify).aa;*)
(*K3[q1_,q2_,q3_,aa_]:=K3[q1,q2,q3,aa]=((Table[c3[[i]][q1,q2,q3],{i,1,15}]/.{al->alfm,be->befm}/.{mag[-q]->mag[q]}/.magreps//Expand)/.mag[0]-> 0//Simplify).aa;*)
(*K4[q1_,q2_,q3_,q4_,aa_]:=K4[q1,q2,q3,q4,aa]=((Monitor[Table[q4sym[c4[[i]],q1,q2,q3,q4]/.{al->alfm,be->befm}/.magreps//Simplify,{i,1,15}],i])/.mag[0]-> 0).aa;*)


(* ::Subsection::Closed:: *)
(*Halo Kernels (real space for delta)*)


(* ::Input:: *)
(*K1d[q1_]:=K1[adelta]//.{kz-> q1.zhat}//.Dot[a_,a_]-> mag[a]^2/.{Dot[a_,zhat]-> Subscript[a, z]}/.magreps//.kzsub;*)
(*K2d[q1_,q2_]:=K2[q1,q2,adelta]/.{al->alfm,be->befm}/.Dot[a_,a_]-> mag[a]^2/.mag[-q]->mag[q]/.{kz-> (q1+q2).zhat}/.magreps//.Dot[a_+b_,c_]-> Dot[a,c]+Dot[b,c]/.{Dot[-a_,b_]-> -Dot[a,b]}/.{Dot[a_,zhat]-> Subscript[a, z]}//.kzsub;*)
(*K2dsym[q1_,q2_]:=q2sym[K2d,q1,q2];*)
(*K3d[q1_,q2_,q3_]:=K3[q1,q2,q3,adelta]/.Dot[a_,a_]-> mag[a]^2/.mag[-q]->mag[q]/.magreps/.{kz-> (q1+q2+q3).zhat}//.Dot[a_+b_,c_]-> Dot[a,c]+Dot[b,c]/.{Dot[-a_,b_]-> -Dot[a,b]}/.{Dot[a_,zhat]-> Subscript[a, z]}//.kzsub;*)
(*K3dsym[q1_,q2_,q3_] := q3sym[K3d,q1,q2,q3];*)
(*K4d[q1_,q2_,q3_,q4_]:=K4[q1,q2,q3,q4,adelta]/.Dot[a_,a_]-> mag[a]^2/.mag[-q]->mag[q]/.{kz-> (q1+q2+q3+q4).zhat}/.magreps//.Dot[a_+b_,c_]-> Dot[a,c]+Dot[b,c]/.{Dot[-a_,b_]-> -Dot[a,b]}/.{Dot[a_,zhat]-> Subscript[a, z]}//.kzsub;*)
(*K4dsym[q1_,q2_,q3_,q4_]:=q4sym[K4d,q1,q2,q3,q4];*)


(* ::Text:: *)
(*Check the <\[Delta](2)> which has to be subtracted*)


(* ::Input:: *)
(*\[Delta]2vev=K2[q1,-q1,adelta]/.tb*)


(* ::Subsection::Closed:: *)
(*Halo Kernels (redshift space)*)


(* ::Text:: *)
(*Notice that we need to normal order, which amounts to a known counterterm due to subtraction of <delta_h(2)> in the redshift space expression for delta*)


(* ::Input:: *)
(*K1r[q1_]:=K1[adelta]+f1 kz^2/(q1).(q1) K1[btheta]//.{kz-> q1.zhat}//.Dot[a_,a_]-> mag[a]^2/.{Dot[a_,zhat]-> Subscript[a, z]}/.magreps//.kzsub;*)


(* ::Input:: *)
(*K2r[q1_,q2_]:=K2[q1,q2,adelta]+f1   kz^2/(q1+q2).(q1+q2) K2[q1,q2,btheta]+ f1 kz 1/2 (q1.zhat/q1.q1+q2.zhat/q2.q2)K1[btheta]K1[adelta]+1/2 f1^2 kz^2 ((q1.zhat)(q2.zhat))/((q1.q1)(q2.q2)) K1[btheta]K1[btheta]/.{al->alfm,be->befm}/.Dot[a_,a_]-> mag[a]^2/.mag[-q]->mag[q]/.{kz-> (q1+q2).zhat}/.magreps//.Dot[a_+b_,c_]-> Dot[a,c]+Dot[b,c]/.{Dot[-a_,b_]-> -Dot[a,b]}/.{Dot[a_,zhat]-> Subscript[a, z]}//.kzsub;*)
(*K2rsym[q1_,q2_]:=q2sym[K2r,q1,q2];*)


(* ::Input:: *)
(*K3r[q1_,q2_,q3_]:=K3[q1,q2,q3,adelta]+f1 kz^2/(q1+q2+q3).(q1+q2+q3) K3[q1,q2,q3,btheta]+f1 kz q3.zhat/q3.q3 K2[q1,q2,adelta]K1[btheta]+f1 kz (q1.zhat+q2.zhat)/(q1+q2).(q1+q2) K2[q1,q2,btheta]K1[adelta]+f1^2/2 kz^2 q1.zhat/q1.q1 q2.zhat/q2.q2 K1[btheta]K1[btheta]K1[adelta]+f1^2 kz^2 (q1.zhat+q2.zhat)/(q1+q2).(q1+q2) q3.zhat/q3.q3 K2[q1,q2,btheta]K1[btheta]+f1^3/6 kz^3 q1.zhat/q1.q1 q2.zhat/q2.q2 q3.zhat/q3.q3 K1[btheta]K1[btheta]K1[btheta]/.Dot[a_,a_]-> mag[a]^2/.mag[-q]->mag[q]/.magreps/.{kz-> (q1+q2+q3).zhat}//.Dot[a_+b_,c_]-> Dot[a,c]+Dot[b,c]/.{Dot[-a_,b_]-> -Dot[a,b]}/.{Dot[a_,zhat]-> Subscript[a, z]}//.kzsub;*)
(*K3rsym[q1_,q2_,q3_] := q3sym[K3r,q1,q2,q3];*)


(* ::Text:: *)
(*Here we split the K4 kernels from the pieces of the previous ones, because of speed*)


(* ::Input:: *)
(*K4rsome[q1_,q2_,q3_,q4_]:=K4[q1,q2,q3,q4,adelta]+f1 kz^2/(q1+q2+q3+q4).(q1+q2+q3+q4) K4[q1,q2,q3,q4,btheta]/.Dot[a_,a_]-> mag[a]^2/.mag[-q]->mag[q]/.{kz-> (q1+q2+q3+q4).zhat}/.magreps//.Dot[a_+b_,c_]-> Dot[a,c]+Dot[b,c]/.{Dot[-a_,b_]-> -Dot[a,b]}/.{Dot[a_,zhat]-> Subscript[a, z]}//.kzsub;*)


(* ::Input:: *)
(*K4dist[q1_,q2_,q3_,q4_]:=f1 kz (q1.zhat+q2.zhat+q3.zhat)/(q1+q2+q3).(q1+q2+q3) K3[q1,q2,q3,btheta]K1[adelta]+*)
(*f1 kz (q1.zhat+q2.zhat)/(q1+q2).(q1+q2) K2[q1,q2,btheta]K2[q3,q4,adelta]+*)
(*f1 kz q1.zhat/q1.q1 K1[btheta]K3[q2,q3,q4,adelta]+*)
(*f1^2 kz^2 (q1.zhat+q2.zhat+q3.zhat)/(q1+q2+q3).(q1+q2+q3) q4.zhat/q4.q4 K3[q1,q2,q3,btheta]K1[btheta]+*)
(*f1^2/2 kz^2 (q1.zhat+q2.zhat)/(q1+q2).(q1+q2) (q3.zhat+q4.zhat)/(q3+q4).(q3+q4) K2[q1,q2,btheta]K2[q3,q4,btheta]+*)
(*f1^2 kz^2 (q1.zhat+q2.zhat)/(q1+q2).(q1+q2) q3.zhat/q3.q3 K2[q1,q2,btheta]K1[btheta]K1[adelta]+*)
(*f1^2/2 kz^2 q1.zhat/q1.q1 q2.zhat/q2.q2 K1[btheta]K1[btheta]K2[q3,q4,adelta]+f1^3/2 kz^3 (q1.zhat+q2.zhat)/(q1+q2).(q1+q2) q3.zhat/q3.q3 q4.zhat/q4.q4 K2[q1,q2,btheta]K1[btheta]K1[btheta]+f1^3/6 kz^3 q1.zhat/q1.q1 q2.zhat/q2.q2 q3.zhat/q3.q3 K1[btheta]K1[btheta]K1[btheta]K1[adelta]+f1^4/24 kz^4 q1.zhat/q1.q1 q2.zhat/q2.q2 q3.zhat/q3.q3 q4.zhat/q4.q4 K1[btheta]K1[btheta]K1[btheta]K1[btheta]/.Dot[a_,a_]-> mag[a]^2/.mag[-q]->mag[q]/.{kz-> (q1+q2+q3+q4).zhat}/.magreps//.Dot[a_+b_,c_]-> Dot[a,c]+Dot[b,c]/.{Dot[-a_,b_]-> -Dot[a,b]}/.{Dot[a_,zhat]-> Subscript[a, z]}//.kzsub;*)


(* ::Input:: *)
(*K4distrsym[q1_,q2_,q3_,q4_]:=q4sym[K4dist,q1,q2,q3,q4];*)


(* ::Input:: *)
(*K4rsym[q1_,q2_,q3_,q4_]:=K4distrsym[q1,q2,q3,q4]+K4rsome[q1,q2,q3,q4];*)


(* ::Section::Closed:: *)
(*PS bias*)


(* ::Subsection:: *)
(*P22*)


(* ::Subsubsection::Closed:: *)
(*Create kernel and ctab*)


(* ::Input:: *)
(*k22[k1_,q_]:=2 K2rsym[-q,q+k1]^2;*)


(* ::Input:: *)
(*P22Integrand[k1_,q_] :=k22[k1,q]*Plin[q]Plin[kpq[k1]]/.tb;*)


(* ::Input::Initialization:: *)
(* Replace \[Mu]0 and integrate \[Beta] *)
(*
ker22exp=Expand[k22[k1,q]/.tb/.magrevreps/.{\[Mu]0\[Rule]x \[Mu]1+Sqrt[1-x^2] Cos[\[Beta]-\[Phi]] Sqrt[1-\[Mu]1^2]},Cos[\[Beta]-\[Phi]]];
ruleCos\[Beta]=Cos[\[Beta]-\[Phi]]^a_.\[RuleDelayed]((1+(-1)^a) Sqrt[\[Pi]] Gamma[(1+a)/2])/Gamma[(2+a)/2];
Cos\[Beta]CoefList=CoefficientList[ker22exp,Cos[\[Beta]-\[Phi]]];
intCos\[Beta]t=Table[Cos[\[Beta]-\[Phi]]^n,{n,1,4}]/.ruleCos\[Beta];
intCos\[Beta]=Join[{2\[Pi]},intCos\[Beta]t];
intker22=1/(2\[Pi])Cos\[Beta]CoefList.intCos\[Beta];
P22temp=intker22 Plin[q]Plin[Sqrt[k1^2+q^2+2 k1 q x]];
*)
P22temp=1/(784 q^2 (k1^2+q^2+2 k1 q x)^2) (1568 b5^2 k1^4 q^2+3136 b5^2 k1^2 q^4+1568 b5^2 q^6+6272 b5^2 k1^3 q^3 x+6272 b5^2 k1 q^5 x+6272 b5^2 k1^2 q^4 x^2+32 b2^2 q^2 (7 q^2+14 k1 q x+k1^2 (5+2 x^2))^2-224 b5 f1 k1^4 q^2 \[Mu]1^2-784 b5 f1^2 k1^4 q^2 \[Mu]1^2-224 b5 f1 k1^2 q^4 \[Mu]1^2-784 b5 f1^2 k1^2 q^4 \[Mu]1^2-1568 b5 f1 k1^5 q x \[Mu]1^2-2016 b5 f1 k1^3 q^3 x \[Mu]1^2-1568 b5 f1^2 k1^3 q^3 x \[Mu]1^2-4480 b5 f1 k1^4 q^2 x^2 \[Mu]1^2+784 b5 f1^2 k1^4 q^2 x^2 \[Mu]1^2-1344 b5 f1 k1^2 q^4 x^2 \[Mu]1^2+784 b5 f1^2 k1^2 q^4 x^2 \[Mu]1^2-2688 b5 f1 k1^3 q^3 x^3 \[Mu]1^2+1568 b5 f1^2 k1^3 q^3 x^3 \[Mu]1^2+8 f1^2 k1^4 q^2 \[Mu]1^4+784 b5 f1^2 k1^4 q^2 \[Mu]1^4+56 f1^3 k1^4 q^2 \[Mu]1^4+147 f1^4 k1^4 q^2 \[Mu]1^4+784 b5 f1^2 k1^2 q^4 \[Mu]1^4+112 f1^2 k1^5 q x \[Mu]1^4-1568 b5 f1^2 k1^5 q x \[Mu]1^4+392 f1^3 k1^5 q x \[Mu]1^4+392 f1^2 k1^6 x^2 \[Mu]1^4+96 f1^2 k1^4 q^2 x^2 \[Mu]1^4-5488 b5 f1^2 k1^4 q^2 x^2 \[Mu]1^4+280 f1^3 k1^4 q^2 x^2 \[Mu]1^4-294 f1^4 k1^4 q^2 x^2 \[Mu]1^4-2352 b5 f1^2 k1^2 q^4 x^2 \[Mu]1^4+672 f1^2 k1^5 q x^3 \[Mu]1^4-392 f1^3 k1^5 q x^3 \[Mu]1^4-4704 b5 f1^2 k1^3 q^3 x^3 \[Mu]1^4+288 f1^2 k1^4 q^2 x^4 \[Mu]1^4-336 f1^3 k1^4 q^2 x^4 \[Mu]1^4+147 f1^4 k1^4 q^2 x^4 \[Mu]1^4+196 f1^4 k1^6 \[Mu]1^6-56 f1^3 k1^4 q^2 \[Mu]1^6-294 f1^4 k1^4 q^2 \[Mu]1^6-280 f1^3 k1^5 q x \[Mu]1^6+1176 f1^4 k1^5 q x \[Mu]1^6+784 f1^3 k1^6 x^2 \[Mu]1^6-196 f1^4 k1^6 x^2 \[Mu]1^6-168 f1^3 k1^4 q^2 x^2 \[Mu]1^6+1764 f1^4 k1^4 q^2 x^2 \[Mu]1^6+1848 f1^3 k1^5 q x^3 \[Mu]1^6-1176 f1^4 k1^5 q x^3 \[Mu]1^6+1008 f1^3 k1^4 q^2 x^4 \[Mu]1^6-1470 f1^4 k1^4 q^2 x^4 \[Mu]1^6-196 f1^4 k1^6 \[Mu]1^8+147 f1^4 k1^4 q^2 \[Mu]1^8-1176 f1^4 k1^5 q x \[Mu]1^8+588 f1^4 k1^6 x^2 \[Mu]1^8-1470 f1^4 k1^4 q^2 x^2 \[Mu]1^8+1960 f1^4 k1^5 q x^3 \[Mu]1^8+1715 f1^4 k1^4 q^2 x^4 \[Mu]1^8+196 b1^2 (8 q^6+32 k1 q^5 x+8 k1^3 q^3 x (3-f1 \[Mu]1^2+4 x^2 (1+f1 \[Mu]1^2))+8 k1^2 q^4 (1-f1 \[Mu]1^2+2 x^2 (3+f1 \[Mu]1^2))+k1^6 (-f1^2 \[Mu]1^2 (-1+\[Mu]1^2)+x^2 (2+4 f1 \[Mu]1^2+f1^2 \[Mu]1^2 (-1+3 \[Mu]1^2)))+4 k1^5 q x (1+f1^2 (\[Mu]1^2-2 \[Mu]1^4)+x^2 (2+4 f1 \[Mu]1^2+f1^2 \[Mu]1^2 (-1+3 \[Mu]1^2)))+2 k1^4 q^2 ((-1+f1 \[Mu]1^2)^2+x^2 (12+8 f1 \[Mu]1^2+f1^2 (2 \[Mu]1^2-6 \[Mu]1^4))+x^4 (4+8 f1 \[Mu]1^2+f1^2 (-2 \[Mu]1^2+6 \[Mu]1^4))))+16 b2 q (7 q^2+14 k1 q x+k1^2 (5+2 x^2)) (28 b5 q (k1^2+q^2+2 k1 q x)-f1 k1^2 \[Mu]1^2 (14 k1 x (1+f1 \[Mu]1^2)+q (2+12 x^2+7 f1 (1-\[Mu]1^2+x^2 (-1+3 \[Mu]1^2)))))-56 b1 (28 b5 q (k1^2+q^2+2 k1 q x) (2 q^3+4 k1 q^2 x+k1^3 (x+f1 x \[Mu]1^2)+k1^2 q (1-f1 \[Mu]1^2+2 x^2 (1+f1 \[Mu]1^2)))+4 b2 q (7 q^2+14 k1 q x+k1^2 (5+2 x^2)) (2 q^3+4 k1 q^2 x+k1^3 (x+f1 x \[Mu]1^2)+k1^2 q (1-f1 \[Mu]1^2+2 x^2 (1+f1 \[Mu]1^2)))-f1 k1^2 \[Mu]1^2 (4 k1 q^3 x (9+12 x^2+7 f1 (1+x^2 (-1+3 \[Mu]1^2)))+2 q^4 (2+12 x^2+7 f1 (1-\[Mu]1^2+x^2 (-1+3 \[Mu]1^2)))+7 k1^4 (-f1^2 \[Mu]1^2 (-1+\[Mu]1^2)+x^2 (2+4 f1 \[Mu]1^2+f1^2 \[Mu]1^2 (-1+3 \[Mu]1^2)))+k1^2 q^2 (2+72 x^2+24 x^4+7 f1^2 \[Mu]1^2 (-1+\[Mu]1^2+x^2 (7-9 \[Mu]1^2)+2 x^4 (-3+5 \[Mu]1^2))+f1 (7-9 \[Mu]1^2+2 x^4 (-7+33 \[Mu]1^2)+x^2 (7+55 \[Mu]1^2)))+k1^3 q x (8 (2+5 x^2)+7 f1^2 \[Mu]1^2 (5-7 \[Mu]1^2+x^2 (-5+11 \[Mu]1^2))+f1 (7-5 \[Mu]1^2+x^2 (-7+89 \[Mu]1^2)))))) Plin[q] Plin[Sqrt[k1^2+q^2+2 k1 q x]];


P22reduc = Expand[P22temp/.{k1^2+q^2+2 k1 q x-> k1pq^2}/.{-k1^2-q^2-2 k1 q x->- k1pq^2}/.{x-> (k1pq^2-k1^2-q^2)/(2k1 q)}/.Sqrt[k1pq^2]-> k1pq];
P22reducnew=P22reduc/.Plin[a_]->1//Expand;
P22ctab=Sort[compressor[getCtab[P22reducnew]]];
p22exps=P22ctab[[All,1;;3]];
p22coefs=P22ctab[[All,4]];


(*check*)
Tabtester[P22reducnew,P22ctab]==0


(* ::Subsubsection::Closed:: *)
(*Decompose into basis of bias coefficients*)


p22vars={b1,b2,b5,f1,\[Mu]1};


(*p22biaslist=getMonomials[p22vars,P22ctab];*)
p22biaslist={b1^2,b1 b2,b2^2,b1 b5,b2 b5,b5^2,b1 f1 \[Mu]1^2,b1^2 f1 \[Mu]1^2,
			 b2 f1 \[Mu]1^2,b1 b2 f1 \[Mu]1^2,b5 f1 \[Mu]1^2,b1 b5 f1 \[Mu]1^2,b1 f1^2 \[Mu]1^2,
			 b1^2 f1^2 \[Mu]1^2,b2 f1^2 \[Mu]1^2,b5 f1^2 \[Mu]1^2,f1^2 \[Mu]1^4,b1 f1^2 \[Mu]1^4,
			 b1^2 f1^2 \[Mu]1^4,b2 f1^2 \[Mu]1^4,b5 f1^2 \[Mu]1^4,f1^3 \[Mu]1^4,b1 f1^3 \[Mu]1^4,
			 f1^4 \[Mu]1^4,f1^3 \[Mu]1^6,b1 f1^3 \[Mu]1^6,f1^4 \[Mu]1^6,f1^4 \[Mu]1^8};


(*Function that obtains the coefficient of each term in biaslist in each ctab coefficient "coef" *)
p22biaslister=biaslisterfn[p22vars,p22biaslist];


(*apply decomposition in biases to all ctab coefficients and simplify*)
p22biascoef=Simplify[p22biaslister/@p22coefs];
p22ctabdec=MapThread[Append,{p22exps,p22biascoef}];


(*checks*)
p22biascoef.p22biaslist == p22coefs


(* ::Subsubsection::Closed:: *)
(*Export/Load*)


(*Save decomposed ctab*)
Export[LoopTablesPath<>"p22ctabdec.m",p22ctabdec]
p22ctabdec=Import[LoopTablesPath<>"p22ctabdec.m"];
p22coefsdec=p22ctabdec[[All,4]];


(*Save ctab*)
Export[ctabpath<>"P22ctab.csv",p22exps]
Export[LoopTablesPath<>"P22ctab.csv",p22exps]

(*Save bias list and decomposed ctab coefficients*)
Export[LoopTablesPath<>"P22biaslist.m",p22biaslist]
Export[LoopTablesPath<>"P22coefsdec.m",p22coefsdec]


(* ::Subsection:: *)
(*P13*)


(* ::Subsubsection::Closed:: *)
(*Create kernel and ctab*)


(* ::Text:: *)
(*Find UV limit of K3*)
(*This only depends on b1, b2, b5*)


(*
serK3k0=Normal[Series[K3rsym[k1,q,-q]/.magrevreps,{q,\[Infinity],0}]]/.tb//Simplify;
uvK3k0=1/2Integrate[serK3k0,{x,-1,1}]//Simplify;
*)
uvK3k0=1/63 (63 b10-34 b2+47 b3-42 b5+110 b6+82 b8+21 b2 f1 \[Mu]1^2+21 b5 f1 \[Mu]1^2-b1 (13+21 f1 \[Mu]1^2));

(*Create new P13 kernel*)
k31[k1_,q_]=K1r[k1](K3rsym[k1,q,-q]-uvK3k0)/.mag[0]-> 0/.tb;
k31exp=k31[k1,q]//Expand;


(* ::Text:: *)
(*Perform sequence of operations to reduce everything to k1pq. If k1mq is in the denominator, do substitution q->-q. *)
(*After that, just write k1mq as a function of k1pq.*)


(* ::Input:: *)
(*k31list=List@@k31exp;*)
(*Do[*)
(*temp=k31list[[i]];*)
(*If[Exponent[temp,k1mq]<0,*)
(*temp=temp/.qrev1/.qrev2];*)
(*temp=temp/.restrep//Expand;*)
(*k31list[[i]]=temp;*)
(*,{i,1,Length[k31list]}*)
(*];*)


(* ::Input:: *)
(*k31new=Total[k31list]//Expand;*)


(* ::Input:: *)
(*(**)
(*P31Integrand[k1_,q_] = 6Plin[k1]Plin[q]k31new;*)
(*P31simp=P31Integrand[k1,q]//Simplify;*)
(**)*)
(*P31simp=-(1/(84 k1^2 k1pq^2 q^4))(b1+f1 \[Mu]1^2) (21 b1 k1^6 k1pq^2-42 b1 k1^4 k1pq^4+21 b1 k1^2 k1pq^6+9 b1 k1^6 q^2+66 b1 k1^4 k1pq^2 q^2+24 b2 k1^4 k1pq^2 q^2-48 b6 k1^4 k1pq^2 q^2-108 b1 k1^2 k1pq^4 q^2-48 b2 k1^2 k1pq^4 q^2+96 b6 k1^2 k1pq^4 q^2+33 b1 k1pq^6 q^2+24 b2 k1pq^6 q^2-48 b6 k1pq^6 q^2-27 b1 k1^4 q^4+19 b1 k1^2 k1pq^2 q^4+16 b2 k1^2 k1pq^2 q^4-32 b6 k1^2 k1pq^2 q^4-66 b1 k1pq^4 q^4-48 b2 k1pq^4 q^4+96 b6 k1pq^4 q^4+27 b1 k1^2 q^6+42 b1 k1pq^2 q^6+24 b2 k1pq^2 q^6-48 b6 k1pq^2 q^6-9 b1 q^8+2 b3 (k1^8-25 k1pq^6 q^2+50 k1pq^4 q^4-26 k1pq^2 q^6+q^8-2 k1^6 (k1pq^2+2 q^2)+2 k1^4 (k1pq^4-11 k1pq^2 q^2+3 q^4)-k1^2 (k1pq^6-52 k1pq^4 q^2+14 k1pq^2 q^4+4 q^6))+2 b8 (15 k1^8-81 k1pq^6 q^2+162 k1pq^4 q^4-96 k1pq^2 q^6+15 q^8-30 k1^6 (k1pq^2+2 q^2)+6 k1^4 (5 k1pq^4-6 k1pq^2 q^2+15 q^4)-k1^2 (15 k1pq^6-192 k1pq^4 q^2+14 k1pq^2 q^4+60 q^6))-84 b1 f1 k1^5 k1pq^2 q \[Mu]0 \[Mu]1+84 b1 f1 k1^3 k1pq^4 q \[Mu]0 \[Mu]1+36 b1 f1 k1^5 q^3 \[Mu]0 \[Mu]1-120 b1 f1 k1^3 k1pq^2 q^3 \[Mu]0 \[Mu]1+36 b1 f1 k1 k1pq^4 q^3 \[Mu]0 \[Mu]1-72 b1 f1 k1^3 q^5 \[Mu]0 \[Mu]1-36 b1 f1 k1 k1pq^2 q^5 \[Mu]0 \[Mu]1+36 b1 f1 k1 q^7 \[Mu]0 \[Mu]1+6 f1 k1^8 \[Mu]1^2+9 f1 k1^6 k1pq^2 \[Mu]1^2-30 f1 k1^4 k1pq^4 \[Mu]1^2+15 f1 k1^2 k1pq^6 \[Mu]1^2-15 f1 k1^6 q^2 \[Mu]1^2+36 b1 f1 k1^6 q^2 \[Mu]1^2+54 f1 k1^4 k1pq^2 q^2 \[Mu]1^2-36 b1 f1 k1^4 k1pq^2 q^2 \[Mu]1^2-36 f1 k1^2 k1pq^4 q^2 \[Mu]1^2+3 f1 k1pq^6 q^2 \[Mu]1^2+9 f1 k1^4 q^4 \[Mu]1^2-72 b1 f1 k1^4 q^4 \[Mu]1^2+15 f1 k1^2 k1pq^2 q^4 \[Mu]1^2-36 b1 f1 k1^2 k1pq^2 q^4 \[Mu]1^2-6 f1 k1pq^4 q^4 \[Mu]1^2+3 f1 k1^2 q^6 \[Mu]1^2+36 b1 f1 k1^2 q^6 \[Mu]1^2+6 f1 k1pq^2 q^6 \[Mu]1^2-3 f1 q^8 \[Mu]1^2-36 f1^2 k1^6 q^2 \[Mu]0^2 \[Mu]1^2+36 f1^2 k1^4 k1pq^2 q^2 \[Mu]0^2 \[Mu]1^2+84 b1 f1^2 k1^4 k1pq^2 q^2 \[Mu]0^2 \[Mu]1^2+72 f1^2 k1^4 q^4 \[Mu]0^2 \[Mu]1^2+36 f1^2 k1^2 k1pq^2 q^4 \[Mu]0^2 \[Mu]1^2-36 f1^2 k1^2 q^6 \[Mu]0^2 \[Mu]1^2-36 f1^2 k1^7 q \[Mu]0 \[Mu]1^3-48 f1^2 k1^5 k1pq^2 q \[Mu]0 \[Mu]1^3+48 f1^2 k1^3 k1pq^4 q \[Mu]0 \[Mu]1^3+72 f1^2 k1^5 q^3 \[Mu]0 \[Mu]1^3-48 f1^2 k1^3 k1pq^2 q^3 \[Mu]0 \[Mu]1^3-36 f1^2 k1^3 q^5 \[Mu]0 \[Mu]1^3+84 f1^3 k1^4 k1pq^2 q^2 \[Mu]0^2 \[Mu]1^4) Plin[k1] Plin[q];*)


(* ::Text:: *)
(*Integrate over mu0*)


(* ::Input:: *)
(*P31\[Mu]0list=CoefficientList[P31simp,\[Mu]0];*)


(* ::Input:: *)
(*\[Mu]0explist={1,\[Mu]0,\[Mu]0^2};*)


(*check decomposition*)
P31\[Mu]0list.\[Mu]0explist == P31simp//Expand


(* ::Input:: *)
(*(*int\[Mu]0list=1/(2\[Pi])Integrate[\[Mu]0explist/.{\[Mu]0\[Rule]x \[Mu]1+Sqrt[1-x^2] Cos[\[Beta]-\[Phi]] Sqrt[1-\[Mu]1^2]},{\[Beta],0,2\[Pi]}];*)*)
(*int\[Mu]0list={1,x \[Mu]1,1/2 (1-\[Mu]1^2+x^2 (-1+3 \[Mu]1^2))};*)


(* ::Input:: *)
(*P31temp=P31\[Mu]0list.int\[Mu]0list;*)
(*P31expanded=P31temp/.{x->(k1pq^2-k1^2-q^2)/(2k1 q)}/.Plin[a_]-> 1//Expand;*)


(* ::Input:: *)
(*Clear[P31list]*)
(*P31list=List@@P31expanded;*)


(* ::Input:: *)
(*P31ctab=Sort[compressor[getCtab[P31list]]];*)


(* ::Input:: *)
(*P31exps=P31ctab[[All,1;;3]];*)
(*P13coefs=P31ctab[[All,4]];*)


(* ::Subsubsection::Closed:: *)
(*Decompose into basis of bias coefficients*)


(* ::Input:: *)
(*p13vars={b1,b2,b3,b6,b8,f1,\[Mu]1};*)


(* ::Input:: *)
(*(*p13biaslist=getMonomials[p13vars,P31ctab];*)*)
(*p13biaslist={b1^2,b1 b2,b1 b3,b1 b6,b1 b8,b1 f1 \[Mu]1^2,b1^2 f1 \[Mu]1^2,b2 f1 \[Mu]1^2,b3 f1 \[Mu]1^2,b6 f1 \[Mu]1^2,b8 f1 \[Mu]1^2,b1 f1^2 \[Mu]1^2,b1^2 f1^2 \[Mu]1^2,f1^2 \[Mu]1^4,b1 f1^2 \[Mu]1^4,b1^2 f1^2 \[Mu]1^4,f1^3 \[Mu]1^4,b1 f1^3 \[Mu]1^4,f1^3 \[Mu]1^6,b1 f1^3 \[Mu]1^6,f1^4 \[Mu]1^6,f1^4 \[Mu]1^8};*)


(* ::Input:: *)
(*(*Function that obtains the coefficient of each term in biaslist in each ctab coefficient "coef" *)*)
(*p13biaslister=biaslisterfn[p13vars,p13biaslist];*)


(* ::Input:: *)
(*(*apply decomposition in biases to all ctab coefficients and simplify*)*)
(*p13biascoef=Simplify[p13biaslister/@P13coefs];*)
(*p13ctabdec=MapThread[Append,{P31exps,p13biascoef}];*)


(*check decomposition*)
p13biascoef.p13biaslist==P13coefs


(* ::Subsubsection::Closed:: *)
(*Export/Load*)


Export[ctabpath<>"P13ctab.csv",P31exps];
Export[LoopTablesPath<>"P13ctab.csv",P31exps]
(*Export["P13coefs.mx",P13coefs]*)
(*Export["P13expanded.mx",P31expanded]*)

(*Save decomposed ctab*)
Export[LoopTablesPath<>"p13ctabdec.m",p13ctabdec]

p13ctabdec=Import[LoopTablesPath<>"p13ctabdec.m"];
p13coefsdec=p13ctabdec[[All,4]];
(*Save bias list and decomposed ctab coefficients*)
Export[LoopTablesPath<>"P13biaslist.m",p13biaslist]
Export[LoopTablesPath<>"P13coefsdec.m",p13coefsdec]


(* ::Section::Closed:: *)
(*Tables of substitutions of integrals with \[Mu]0 into integrals of \[Mu]1, \[Mu]2*)


(* ::Text:: *)
(*Final output of this section are the 4 bi-spectrum diagrams tabulated to the master integral*)


(* ::Subsection::Closed:: *)
(*Useful definitions*)


(* ::Text:: *)
(*All integrands in Bintegrand can be FFTlogged to q^(-2\[Nu]1) k1pq^(-2 \[Nu]2) k2mq^(-2 \[Nu]3) = I(\[Nu]1,\[Nu]2,\[Nu]3,k1,k2).*)
(*If such an integrand is multiplied by \[Mu]0^n, then we can expand the vector integral q...q I in scalar integrals.*)
(*We will then get some scalar coefficients and \[Mu]1 and \[Mu]2*)


(* ::Text:: *)
(*Great trick to get scalar products. Instead of vector components, we write vectors such that the scalar products come out right as a single symbol.*)
(*Notice that we need separately scalar products with q, scalar products between k1, k2, and finally scalar products k1.z, k2.z .*)
(*So we separate definitions by using the trick of getting 3 components corresponding to the 3 scalar products q.q, q.k1, q.k2.*)
(*Then we make definitions for k1, k2, and for z*)


(* ::Input:: *)
(*qv={Sqrt[q^2-1-qk2^2],1,qk2};*)
(*k1vq={0,qk1,0};*)
(*k2vq={0,0,1};*)
(*mydeltaq=IdentityMatrix[3];*)


(* ::Input:: *)
(*k1vk={Sqrt[k1^2-1],1,0};*)
(*k2vk={0,k1k2,Sqrt[k2^2-k1k2^2]};*)
(*mydeltak=IdentityMatrix[3];*)
(**)
(*{k1vk.k1vk, k2vk.k2vk, k1vk.k2vk, k1vk.mydeltak.k1vk, k1vk.mydeltak.k2vk, k2vk.mydeltak.k2vk}*)


(* ::Input:: *)
(*k1vz={0,0,k1 \[Mu]1};*)
(*k2vz={0,0,k2 \[Mu]2};*)
(*zh={0,0,1};*)
(*mydeltaz=IdentityMatrix[3];*)
(**)
(*{k1vz.zh,k2vz.zh,zh.mydeltaz.zh}*)


(* ::Text:: *)
(*In every case, we have the simple identity and traces of the deltas are fine.*)
(*Substitutions for the scalar products.*)


(* ::Input:: *)
(*kqsub={qk1->1/2 (k1pq^2-k1^2-q^2),qk2->1/2 (-k2mq^2+k2^2+q^2)};*)
(*k1k2sub={k1k2->1/2 (k3^2-k1^2-k2^2)};*)


(* ::Text:: *)
(*Function to get the product of n vectors, as in q^i... q^n and similar for z*)


(* ::Input:: *)
(*Getvtensor[v_,n_]:=TensorProduct@@Table[v,n]*)


(* ::Text:: *)
(*Functions to generate the vectors of operators in which we expand the integrals with q^i... q^n at each order.*)
(*The simplest is B3212, which only depends on k2 (for the permutation we consider), up to second order.*)
(*Finally we need non-symmetric tensors for B222 up to 6th order,  for B411 up to second order, and for B3211*)


(* ::Input:: *)
(*Tns1[k1v_,k2v_,mydelta_]:=Symmetrize/@{TensorProduct[k1v],*)
(*TensorProduct[k2v]};*)
(*Tns2[k1v_,k2v_,mydelta_]:=Symmetrize/@{TensorProduct[k1v,k1v],*)
(*TensorProduct[k1v,k2v],*)
(*TensorProduct[k2v,k2v],*)
(*mydelta};*)
(*Tns3[k1v_,k2v_,mydelta_]:=Symmetrize/@{TensorProduct[k1v,k1v,k1v],*)
(*TensorProduct[k1v,k1v,k2v],*)
(*TensorProduct[k1v,k2v,k2v],*)
(*TensorProduct[k2v,k2v,k2v],*)
(*TensorProduct[k1v,mydelta],*)
(*TensorProduct[k2v,mydelta]};*)
(*Tns4[k1v_,k2v_,mydelta_]:=Symmetrize/@{TensorProduct[k1v,k1v,k1v,k1v],*)
(*TensorProduct[k1v,k1v,k1v,k2v],*)
(*TensorProduct[k1v,k1v,k2v,k2v],*)
(*TensorProduct[k1v,k2v,k2v,k2v],*)
(*TensorProduct[k2v,k2v,k2v,k2v],*)
(*TensorProduct[k1v,k1v,mydelta],*)
(*TensorProduct[k1v,k2v,mydelta],*)
(*TensorProduct[k2v,k2v,mydelta],*)
(*TensorProduct[mydelta,mydelta]};*)
(*Tns5[k1v_,k2v_,mydelta_]:=Symmetrize/@{TensorProduct[k1v,k1v,k1v,k1v,k1v],*)
(*TensorProduct[k1v,k1v,k1v,k1v,k2v],*)
(*TensorProduct[k1v,k1v,k1v,k2v,k2v],*)
(*TensorProduct[k1v,k1v,k2v,k2v,k2v],*)
(*TensorProduct[k1v,k2v,k2v,k2v,k2v],*)
(*TensorProduct[k2v,k2v,k2v,k2v,k2v],*)
(*TensorProduct[k1v,k1v,k1v,mydelta],*)
(*TensorProduct[k1v,k1v,k2v,mydelta],*)
(*TensorProduct[k1v,k2v,k2v,mydelta],*)
(*TensorProduct[k2v,k2v,k2v,mydelta],*)
(*TensorProduct[k1v,mydelta,mydelta],*)
(*TensorProduct[k2v,mydelta,mydelta]};*)
(*Tns6[k1v_,k2v_,mydelta_]:=Symmetrize/@{TensorProduct[k1v,k1v,k1v,k1v,k1v,k1v],*)
(*TensorProduct[k1v,k1v,k1v,k1v,k1v,k2v],*)
(*TensorProduct[k1v,k1v,k1v,k1v,k2v,k2v],*)
(*TensorProduct[k1v,k1v,k1v,k2v,k2v,k2v],*)
(*TensorProduct[k1v,k1v,k2v,k2v,k2v,k2v],*)
(*TensorProduct[k1v,k2v,k2v,k2v,k2v,k2v],*)
(*TensorProduct[k2v,k2v,k2v,k2v,k2v,k2v],*)
(*TensorProduct[k1v,k1v,k1v,k1v,mydelta],*)
(*TensorProduct[k1v,k1v,k1v,k2v,mydelta],*)
(*TensorProduct[k1v,k1v,k2v,k2v,mydelta],*)
(*TensorProduct[k1v,k2v,k2v,k2v,mydelta],*)
(*TensorProduct[k2v,k2v,k2v,k2v,mydelta],*)
(*TensorProduct[k1v,k1v,mydelta,mydelta],*)
(*TensorProduct[k1v,k2v,mydelta,mydelta],*)
(*TensorProduct[k2v,k2v,mydelta,mydelta],*)
(*TensorProduct[mydelta,mydelta,mydelta]};*)


(* ::Input:: *)
(*Tv1[kv_,mydelta_]:=Symmetrize/@{TensorProduct[kv]};*)
(*Tv2[kv_,mydelta_]:=Symmetrize/@{TensorProduct[kv,kv],mydelta};*)


(* ::Subsection:: *)
(*Redshift space kernels, \[Mu]0 expansion non-symmetric in k1-k2, up to \[Mu]0^6 (for B3211, B222 and B411)*)


(* ::Subsubsection::Closed:: *)
(*First order, \[Mu]0^1*)


(* ::Text:: *)
(*Generate tensors for the different contractions we need*)


(* ::Input:: *)
(*Tq=Tns1[k1vq,k2vq,IdentityMatrix[3]];*)
(*Tk=Tns1[k1vk,k2vk,IdentityMatrix[3]];*)
(*Tz=Tns1[k1vz,k2vz,IdentityMatrix[3]];*)
(*qtensor=Getvtensor[qv/q,1];*)
(*ztensor=Getvtensor[zh,1];*)
(*av=Table[Symbol["a"<>ToString[i]],{i,Length[Tk]}]*)


(* ::Text:: *)
(*Contraction of T's with zhat's*)


(* ::Input:: *)
(*Tdotz=Tz.ztensor//Simplify*)


(* ::Text:: *)
(*Contraction with the integral*)


(* ::Input:: *)
(*Tdotint=Tq.qtensor//FullSimplify*)


(* ::Text:: *)
(*Contraction of T with itself, and inverse*)


(* ::Input:: *)
(*TdotT=Tk.Transpose[Tk]/.k1k2sub//Simplify;*)


(* ::Input:: *)
(*invT=Inverse[TdotT]//Simplify;*)


(* ::Text:: *)
(*Solution*)


(* ::Input:: *)
(*\[Mu]0sol1=Tdotz.invT.Tdotint/.kqsub//Simplify*)


(* ::Subsubsection::Closed:: *)
(*Second order, \[Mu]0^2*)


(* ::Text:: *)
(*Generate tensors for the different contractions we need*)


(* ::Input:: *)
(*Tq=Tns2[k1vq,k2vq,IdentityMatrix[3]];*)
(*Tk=Tns2[k1vk,k2vk,IdentityMatrix[3]];*)
(*Tz=Tns2[k1vz,k2vz,IdentityMatrix[3]];*)
(*qtensor=Getvtensor[qv/q,2];*)
(*ztensor=Getvtensor[zh,2];*)
(*av=Table[Symbol["a"<>ToString[i]],{i,Length[Tk]}]*)


(* ::Text:: *)
(*Contraction of T's with zhat's*)


(* ::Input:: *)
(*Tdotz=Table[Sum[Tz[[\[Alpha],i,j]]ztensor[[i,j]],{i,Length[zh]},{j,Length[zh]}],*)
(*{\[Alpha],1,Length[Tz]}]//Simplify*)


(* ::Text:: *)
(*Contraction with the integral*)


(* ::Input:: *)
(*Tdotint=Table[Sum[Tq[[\[Alpha],i,j]]qtensor[[i,j]],{i,Length[qv]},{j,Length[qv]}],*)
(*{\[Alpha],1,Length[Tq]}]//FullSimplify*)


(* ::Text:: *)
(*Contraction of T with itself, and inverse*)


(* ::Input:: *)
(*TdotT=Table[Sum[Tk[[\[Alpha],i,j]] Tk[[\[Beta],i,j]],{i,Length[qv]},{j,Length[qv]}],*)
(*{\[Alpha],1,Length[Tk]},{\[Beta],1,Length[Tk]}]//Simplify;*)


(* ::Input:: *)
(*invT=Simplify[Inverse[TdotT]]/.k1k2sub//Simplify;*)


(* ::Text:: *)
(*Solution*)


(* ::Input:: *)
(*\[Mu]0sol2=Tdotz.invT.Tdotint/.kqsub//Simplify*)


(* ::Subsubsection::Closed:: *)
(*Third order, \[Mu]0^3*)


(* ::Text:: *)
(*Generate tensors for the different contractions we need*)


(* ::Input:: *)
(*Tq=Tns3[k1vq,k2vq,IdentityMatrix[3]];*)
(*Tk=Tns3[k1vk,k2vk,IdentityMatrix[3]];*)
(*Tz=Tns3[k1vz,k2vz,IdentityMatrix[3]];*)
(*qtensor=Getvtensor[qv/q,3];*)
(*ztensor=Getvtensor[zh,3];*)
(*av=Table[Symbol["a"<>ToString[i]],{i,Length[Tk]}]*)


(* ::Text:: *)
(*Contraction of T's with zhat's*)


(* ::Input:: *)
(*Tdotz=Table[Sum[Tz[[\[Alpha],i,j,k]]ztensor[[i,j,k]],{i,Length[zh]},{j,Length[zh]},{k,Length[zh]}],*)
(*{\[Alpha],1,Length[Tz]}]//Simplify*)


(* ::Text:: *)
(*Contraction with the integral*)


(* ::Input:: *)
(*Tdotint=Table[Sum[Tq[[\[Alpha],i,j,k]]qtensor[[i,j,k]],{i,Length[qv]},{j,Length[qv]},{k,Length[qv]}],*)
(*{\[Alpha],1,Length[Tq]}]//FullSimplify;*)


(* ::Text:: *)
(*Contraction of T with itself, and inverse*)


(* ::Input:: *)
(*TdotT=Table[Sum[Tk[[\[Alpha],i,j,k]] Tk[[\[Beta],i,j,k]],{i,Length[qv]},{j,Length[qv]},{k,Length[qv]}],*)
(*{\[Alpha],1,Length[Tk]},{\[Beta],1,Length[Tk]}]//Simplify;*)


(* ::Input:: *)
(*invT=Simplify[Inverse[TdotT]]/.k1k2sub//Simplify;*)


(* ::Text:: *)
(*Solution*)


(* ::Input:: *)
(*\[Mu]0sol3=Tdotz.invT.Tdotint/.kqsub//Simplify;*)


(* ::Subsubsection::Closed:: *)
(*Fourth order, \[Mu]0^4*)


(* ::Text:: *)
(*Generate tensors for the different contractions we need*)


(* ::Input:: *)
(*Tq=Tns4[k1vq,k2vq,IdentityMatrix[3]];*)
(*Tk=Tns4[k1vk,k2vk,IdentityMatrix[3]];*)
(*Tz=Tns4[k1vz,k2vz,IdentityMatrix[3]];*)
(*qtensor=Getvtensor[qv/q,4];*)
(*ztensor=Getvtensor[zh,4];*)
(*av=Table[Symbol["a"<>ToString[i]],{i,Length[Tk]}]*)


(* ::Text:: *)
(*Contraction of T's with zhat's*)


(* ::Input:: *)
(*Tdotz=Table[Sum[Tz[[\[Alpha],i,j,k,l]]ztensor[[i,j,k,l]],{i,Length[zh]},{j,Length[zh]},{k,Length[zh]},{l,Length[zh]}],*)
(*{\[Alpha],1,Length[Tz]}]//Simplify*)


(* ::Text:: *)
(*Contraction with the integral*)


(* ::Input:: *)
(*Tdotint=Table[Sum[Tq[[\[Alpha],i,j,k,l]]qtensor[[i,j,k,l]],{i,Length[qv]},{j,Length[qv]},{k,Length[qv]},{l,Length[qv]}],*)
(*{\[Alpha],1,Length[Tq]}]/.kqsub//Simplify;*)


(* ::Text:: *)
(*Contraction of T with itself, and inverse*)


(* ::Input:: *)
(*TdotT=Table[Sum[Tk[[\[Alpha],i,j,k,l]] Tk[[\[Beta],i,j,k,l]],{i,Length[qv]},{j,Length[qv]},{k,Length[qv]},{l,Length[qv]}],*)
(*{\[Alpha],1,Length[Tk]},{\[Beta],1,Length[Tk]}]//Simplify;*)


(* ::Input:: *)
(*(*invT=Simplify[Inverse[TdotT]]/.k1k2sub//Simplify;*)*)


(* ::Text:: *)
(*Solution*)


(* ::Input:: *)
(*(*\[Mu]0sol4=Tdotz.invT.Tdotint//Simplify;*)*)
(*\[Mu]0sol4=1+(8 k1^2 k2^2 \[Mu]1^2)/(k1^4+(k2^2-k3^2)^2-2 k1^2 (k2^2+k3^2))+(16 k1^4 k2^4 \[Mu]1^4)/(k1^4+(k2^2-k3^2)^2-2 k1^2 (k2^2+k3^2))^2+(8 k1 k2 (k1^2+k2^2-k3^2) \[Mu]1 \[Mu]2)/(k1^4+(k2^2-k3^2)^2-2 k1^2 (k2^2+k3^2))+(32 k1^3 k2^3 (k1^2+k2^2-k3^2) \[Mu]1^3 \[Mu]2)/(k1^4+(k2^2-k3^2)^2-2 k1^2 (k2^2+k3^2))^2+(8 k1^2 k2^2 \[Mu]2^2)/(k1^4+(k2^2-k3^2)^2-2 k1^2 (k2^2+k3^2))+(2 k1^2 k2^2 (k1^2 k2^2+1/2 (k1^2+k2^2-k3^2)^2) \[Mu]1^2 \[Mu]2^2)/(k1^2 k2^2-1/4 (k1^2+k2^2-k3^2)^2)^2+(32 k1^3 k2^3 (k1^2+k2^2-k3^2) \[Mu]1 \[Mu]2^3)/(k1^4+(k2^2-k3^2)^2-2 k1^2 (k2^2+k3^2))^2+(16 k1^4 k2^4 \[Mu]2^4)/(k1^4+(k2^2-k3^2)^2-2 k1^2 (k2^2+k3^2))^2+((k1^2-k1pq^2+q^2)^2 (k2^2-k2mq^2+q^2)^2 ((k1^2 k2^2-1/4 (k1^2+k2^2-k3^2)^2)^2 (k1^2 k2^2+1/2 (k1^2+k2^2-k3^2)^2)-k1^2 (k1^2 k2^2-1/4 (k1^2+k2^2-k3^2)^2) (5 k1^2 k2^4+19/4 k2^2 (k1^2+k2^2-k3^2)^2) \[Mu]1^2+4 k1^4 (k1^2 k2^6+5/4 k2^4 (k1^2+k2^2-k3^2)^2) \[Mu]1^4-14 k1^5 k2^5 (k1^2+k2^2-k3^2) \[Mu]1 \[Mu]2+k1^3 k2^3 (k1^2+k2^2-k3^2)^3 \[Mu]1 \[Mu]2+5/8 k1 k2 (k1^2+k2^2-k3^2)^5 \[Mu]1 \[Mu]2+20 k1^5 k2^5 (k1^2+k2^2-k3^2) \[Mu]1^3 \[Mu]2+7 k1^3 k2^3 (k1^2+k2^2-k3^2)^3 \[Mu]1^3 \[Mu]2-k2^2 (k1^2 k2^2-1/4 (k1^2+k2^2-k3^2)^2) (5 k1^4 k2^2+19/4 k1^2 (k1^2+k2^2-k3^2)^2) \[Mu]2^2+k1^2 k2^2 (11 k1^4 k2^4+49/2 k1^2 k2^2 (k1^2+k2^2-k3^2)^2+35/16 (k1^2+k2^2-k3^2)^4) \[Mu]1^2 \[Mu]2^2+20 k1^5 k2^5 (k1^2+k2^2-k3^2) \[Mu]1 \[Mu]2^3+7 k1^3 k2^3 (k1^2+k2^2-k3^2)^3 \[Mu]1 \[Mu]2^3+4 k2^4 (k1^6 k2^2+5/4 k1^4 (k1^2+k2^2-k3^2)^2) \[Mu]2^4))/(8 (k1^2 k2^2-1/4 (k1^2+k2^2-k3^2)^2)^4 q^4)-1/(4 q^4) k1 (k1^2-k1pq^2+q^2) (k2^2-k2mq^2+q^2)^3 ((8 k1 (k1^2+k2^2-k3^2))/(k1^4+(k2^2-k3^2)^2-2 k1^2 (k2^2+k3^2))^2-(k1 (-(5/2) k1^2 k2^2 (k1^2+k2^2-k3^2)-3/8 (k1^2+k2^2-k3^2)^3) \[Mu]1^2)/(-k1^2 k2^2+1/4 (k1^2+k2^2-k3^2)^2)^3+(128 k1^3 k2^2 (k1^2+k2^2-k3^2) (k1^4+k1^2 (6 k2^2-2 k3^2)+(k2^2-k3^2)^2) \[Mu]1^4)/(k1^4+(k2^2-k3^2)^2-2 k1^2 (k2^2+k3^2))^4+(k2 (3 k1^4 k2^2+13/4 k1^2 (k1^2+k2^2-k3^2)^2) \[Mu]1 \[Mu]2)/(-k1^2 k2^2+1/4 (k1^2+k2^2-k3^2)^2)^3+(k1^2 k2 (3 k1^4 k2^4+11/2 k1^2 k2^2 (k1^2+k2^2-k3^2)^2+7/16 (k1^2+k2^2-k3^2)^4) \[Mu]1^3 \[Mu]2)/(k1^2 k2^2-1/4 (k1^2+k2^2-k3^2)^2)^4+(256 k1^3 k2^2 (k1^2+k2^2-k3^2) \[Mu]2^2)/(k1^4+(k2^2-k3^2)^2-2 k1^2 (k2^2+k3^2))^3-(4 k1 k2^2 (-(5/2) k1^4 k2^2 (k1^2+k2^2-k3^2)-7/8 k1^2 (k1^2+k2^2-k3^2)^3) \[Mu]1^2 \[Mu]2^2)/(k1^2 k2^2-1/4 (k1^2+k2^2-k3^2)^2)^4+(4 k2^3 (k1^6 k2^2+7/4 k1^4 (k1^2+k2^2-k3^2)^2) \[Mu]1 \[Mu]2^3)/(k1^2 k2^2-1/4 (k1^2+k2^2-k3^2)^2)^4+(1024 k1^5 k2^4 (k1^2+k2^2-k3^2) \[Mu]2^4)/(k1^4+(k2^2-k3^2)^2-2 k1^2 (k2^2+k3^2))^4)-1/(4 q^4) k2 (k1^2-k1pq^2+q^2)^3 (k2^2-k2mq^2+q^2) ((8 k2 (k1^2+k2^2-k3^2))/(k1^4+(k2^2-k3^2)^2-2 k1^2 (k2^2+k3^2))^2+(256 k1^2 k2^3 (k1^2+k2^2-k3^2) \[Mu]1^2)/(k1^4+(k2^2-k3^2)^2-2 k1^2 (k2^2+k3^2))^3+(1024 k1^4 k2^5 (k1^2+k2^2-k3^2) \[Mu]1^4)/(k1^4+(k2^2-k3^2)^2-2 k1^2 (k2^2+k3^2))^4+(k1 (3 k1^2 k2^4+13/4 k2^2 (k1^2+k2^2-k3^2)^2) \[Mu]1 \[Mu]2)/(-k1^2 k2^2+1/4 (k1^2+k2^2-k3^2)^2)^3+(4 k1^3 (k1^2 k2^6+7/4 k2^4 (k1^2+k2^2-k3^2)^2) \[Mu]1^3 \[Mu]2)/(k1^2 k2^2-1/4 (k1^2+k2^2-k3^2)^2)^4-(k2 (-(5/2) k1^2 k2^2 (k1^2+k2^2-k3^2)-3/8 (k1^2+k2^2-k3^2)^3) \[Mu]2^2)/(-k1^2 k2^2+1/4 (k1^2+k2^2-k3^2)^2)^3-(4 k1^2 k2 (-(5/2) k1^2 k2^4 (k1^2+k2^2-k3^2)-7/8 k2^2 (k1^2+k2^2-k3^2)^3) \[Mu]1^2 \[Mu]2^2)/(k1^2 k2^2-1/4 (k1^2+k2^2-k3^2)^2)^4+(k1 k2^2 (3 k1^4 k2^4+11/2 k1^2 k2^2 (k1^2+k2^2-k3^2)^2+7/16 (k1^2+k2^2-k3^2)^4) \[Mu]1 \[Mu]2^3)/(k1^2 k2^2-1/4 (k1^2+k2^2-k3^2)^2)^4+(128 k1^2 k2^3 (k1^2+k2^2-k3^2) (k1^4+k1^2 (6 k2^2-2 k3^2)+(k2^2-k3^2)^2) \[Mu]2^4)/(k1^4+(k2^2-k3^2)^2-2 k1^2 (k2^2+k3^2))^4)+(k2^4 (k1^2-k1pq^2+q^2)^4 (32 k1^7 k2 \[Mu]1 \[Mu]2 (1+\[Mu]2^2)+32 k1 k2 (k2^2-k3^2)^3 \[Mu]1 \[Mu]2 (1+\[Mu]2^2)+k1^8 (1+6 \[Mu]2^2+\[Mu]2^4)+(k2^2-k3^2)^4 (1+6 \[Mu]2^2+\[Mu]2^4)+32 k1^5 k2 \[Mu]1 \[Mu]2 (-3 k3^2 (1+\[Mu]2^2)+k2^2 (-1+8 \[Mu]1^2+7 \[Mu]2^2))+32 k1^3 k2 (k2^2-k3^2) \[Mu]1 \[Mu]2 (-3 k3^2 (1+\[Mu]2^2)+k2^2 (-1+8 \[Mu]1^2+7 \[Mu]2^2))+4 k1^6 (-k3^2 (1+6 \[Mu]2^2+\[Mu]2^4)+k2^2 (-1+2 \[Mu]2^2+7 \[Mu]2^4+8 \[Mu]1^2 (1+5 \[Mu]2^2)))+4 k1^2 (k2^2-k3^2)^2 (-k3^2 (1+6 \[Mu]2^2+\[Mu]2^4)+k2^2 (-1+2 \[Mu]2^2+7 \[Mu]2^4+8 \[Mu]1^2 (1+5 \[Mu]2^2)))+2 k1^4 (3 k3^4 (1+6 \[Mu]2^2+\[Mu]2^4)-2 k2^2 k3^2 (-1+10 \[Mu]2^2+15 \[Mu]2^4+16 \[Mu]1^2 (1+5 \[Mu]2^2))+k2^4 (3+64 \[Mu]1^4-14 \[Mu]2^2+35 \[Mu]2^4+32 \[Mu]1^2 (-1+7 \[Mu]2^2)))))/((k1^4+(k2^2-k3^2)^2-2 k1^2 (k2^2+k3^2))^4 q^4)+(2 k2^2 (k1^2-k1pq^2+q^2)^2 (k1^8 (1+3 \[Mu]2^2)+(k2^2-k3^2)^4 (1+3 \[Mu]2^2)+4 k1^7 k2 \[Mu]1 \[Mu]2 (5+3 \[Mu]2^2)+4 k1 k2 (k2^2-k3^2)^3 \[Mu]1 \[Mu]2 (5+3 \[Mu]2^2)+4 k1^5 k2 \[Mu]1 \[Mu]2 (-3 k3^2 (5+3 \[Mu]2^2)+k2^2 (-5+32 \[Mu]1^2+29 \[Mu]2^2))+4 k1^3 k2 (k2^2-k3^2) \[Mu]1 \[Mu]2 (-3 k3^2 (5+3 \[Mu]2^2)+k2^2 (-5+32 \[Mu]1^2+29 \[Mu]2^2))+4 k1^6 (-k3^2 (1+3 \[Mu]2^2)+k2^2 (-1+2 \[Mu]2^2+3 \[Mu]2^4+\[Mu]1^2 (5+19 \[Mu]2^2)))+4 k1^2 (k2^2-k3^2)^2 (-k3^2 (1+3 \[Mu]2^2)+k2^2 (-1+2 \[Mu]2^2+3 \[Mu]2^4+\[Mu]1^2 (5+19 \[Mu]2^2)))+2 k1^4 (3 k3^4 (1+3 \[Mu]2^2)-2 k2^2 k3^2 (-1+7 \[Mu]2^2+6 \[Mu]2^4+2 \[Mu]1^2 (5+19 \[Mu]2^2))+k2^4 (3+32 \[Mu]1^4-11 \[Mu]2^2+20 \[Mu]2^4+4 \[Mu]1^2 (-5+29 \[Mu]2^2)))))/((k1^4+(k2^2-k3^2)^2-2 k1^2 (k2^2+k3^2))^3 q^2)+(k1^4 (k2^2-k2mq^2+q^2)^4 (k1^8 (1+6 \[Mu]1^2+\[Mu]1^4)+(k2^2-k3^2)^4 (1+6 \[Mu]1^2+\[Mu]1^4)+32 k1^7 k2 \[Mu]1 (1+\[Mu]1^2) \[Mu]2+32 k1 k2 (k2^2-k3^2)^3 \[Mu]1 (1+\[Mu]1^2) \[Mu]2+32 k1^5 k2 \[Mu]1 \[Mu]2 (-3 k3^2 (1+\[Mu]1^2)+k2^2 (-1+7 \[Mu]1^2+8 \[Mu]2^2))+32 k1^3 k2 (k2^2-k3^2) \[Mu]1 \[Mu]2 (-3 k3^2 (1+\[Mu]1^2)+k2^2 (-1+7 \[Mu]1^2+8 \[Mu]2^2))+2 k1^4 (3 k3^4 (1+6 \[Mu]1^2+\[Mu]1^4)-2 k2^2 k3^2 (-1+15 \[Mu]1^4+16 \[Mu]2^2+10 \[Mu]1^2 (1+8 \[Mu]2^2))+k2^4 (3+35 \[Mu]1^4-32 \[Mu]2^2+64 \[Mu]2^4+14 \[Mu]1^2 (-1+16 \[Mu]2^2)))+4 k1^6 (-k3^2 (1+6 \[Mu]1^2+\[Mu]1^4)+k2^2 (-1+7 \[Mu]1^4+8 \[Mu]2^2+\[Mu]1^2 (2+40 \[Mu]2^2)))+4 k1^2 (k2^2-k3^2)^2 (-k3^2 (1+6 \[Mu]1^2+\[Mu]1^4)+k2^2 (-1+7 \[Mu]1^4+8 \[Mu]2^2+\[Mu]1^2 (2+40 \[Mu]2^2)))))/((k1^4+(k2^2-k3^2)^2-2 k1^2 (k2^2+k3^2))^4 q^4)-(2 (k1^2-k1pq^2+q^2) (k2^2-k2mq^2+q^2) (k1^10+(k2^2-k3^2)^5+14 k1^9 k2 \[Mu]1 \[Mu]2+14 k1 k2 (k2^2-k3^2)^4 \[Mu]1 \[Mu]2+4 k1^5 k2 \[Mu]1 \[Mu]2 (21 k3^4-26 k2^2 k3^2 (1+2 \[Mu]1^2+2 \[Mu]2^2)+19 k2^4 (-1+4 \[Mu]1^2+4 \[Mu]2^2))+8 k1^7 k2 \[Mu]1 \[Mu]2 (-7 k3^2+k2^2 (3+13 \[Mu]1^2+13 \[Mu]2^2))+8 k1^3 k2 (k2^2-k3^2)^2 \[Mu]1 \[Mu]2 (-7 k3^2+k2^2 (3+13 \[Mu]1^2+13 \[Mu]2^2))+k1^8 (-5 k3^2+k2^2 (-3+20 \[Mu]2^2+20 \[Mu]1^2 (1+2 \[Mu]2^2)))+k1^2 (k2^2-k3^2)^3 (-5 k3^2+k2^2 (-3+20 \[Mu]2^2+20 \[Mu]1^2 (1+2 \[Mu]2^2)))+2 k1^6 (5 k3^4-2 k2^2 k3^2 (-1+15 \[Mu]2^2+15 \[Mu]1^2 (1+2 \[Mu]2^2))+k2^4 (1+32 \[Mu]1^4-10 \[Mu]2^2+32 \[Mu]2^4+2 \[Mu]1^2 (-5+86 \[Mu]2^2)))+2 k1^4 (k2^2-k3^2) (5 k3^4-2 k2^2 k3^2 (-1+15 \[Mu]2^2+15 \[Mu]1^2 (1+2 \[Mu]2^2))+k2^4 (1+32 \[Mu]1^4-10 \[Mu]2^2+32 \[Mu]2^4+2 \[Mu]1^2 (-5+86 \[Mu]2^2)))))/((k1^4+(k2^2-k3^2)^2-2 k1^2 (k2^2+k3^2))^3 q^2)+(2 k1^2 (k2^2-k2mq^2+q^2)^2 (k1^8 (1+3 \[Mu]1^2)+(k2^2-k3^2)^4 (1+3 \[Mu]1^2)+4 k1^7 k2 \[Mu]1 (5+3 \[Mu]1^2) \[Mu]2+4 k1 k2 (k2^2-k3^2)^3 \[Mu]1 (5+3 \[Mu]1^2) \[Mu]2+4 k1^5 k2 \[Mu]1 \[Mu]2 (-3 k3^2 (5+3 \[Mu]1^2)+k2^2 (-5+29 \[Mu]1^2+32 \[Mu]2^2))+4 k1^3 k2 (k2^2-k3^2) \[Mu]1 \[Mu]2 (-3 k3^2 (5+3 \[Mu]1^2)+k2^2 (-5+29 \[Mu]1^2+32 \[Mu]2^2))+4 k1^6 (-k3^2 (1+3 \[Mu]1^2)+k2^2 (-1+3 \[Mu]1^4+5 \[Mu]2^2+\[Mu]1^2 (2+19 \[Mu]2^2)))+4 k1^2 (k2^2-k3^2)^2 (-k3^2 (1+3 \[Mu]1^2)+k2^2 (-1+3 \[Mu]1^4+5 \[Mu]2^2+\[Mu]1^2 (2+19 \[Mu]2^2)))+2 k1^4 (3 k3^4 (1+3 \[Mu]1^2)-2 k2^2 k3^2 (-1+6 \[Mu]1^4+10 \[Mu]2^2+\[Mu]1^2 (7+38 \[Mu]2^2))+k2^4 (3+20 \[Mu]1^4-20 \[Mu]2^2+32 \[Mu]2^4+\[Mu]1^2 (-11+116 \[Mu]2^2)))))/((k1^4+(k2^2-k3^2)^2-2 k1^2 (k2^2+k3^2))^3 q^2);*)


(* ::Subsubsection::Closed:: *)
(*Fifth order, \[Mu]0^5*)


(* ::Text:: *)
(*Generate tensors for the different contractions we need*)


(* ::Input:: *)
(*Tq=Tns5[k1vq,k2vq,IdentityMatrix[3]];*)
(*Tk=Tns5[k1vk,k2vk,IdentityMatrix[3]];*)
(*Tz=Tns5[k1vz,k2vz,IdentityMatrix[3]];*)
(*qtensor=Getvtensor[qv/q,5];*)
(*ztensor=Getvtensor[zh,5];*)
(*av=Table[Symbol["a"<>ToString[i]],{i,Length[Tk]}]*)


(* ::Text:: *)
(*Contraction of T's with zhat's*)


(* ::Input:: *)
(*Tdotz=Table[Sum[Tz[[\[Alpha],i,j,k,l,m]]ztensor[[i,j,k,l,m]],{i,Length[zh]},{j,Length[zh]},{k,Length[zh]},{l,Length[zh]},{m,Length[zh]}],*)
(*{\[Alpha],1,Length[Tz]}]//Simplify*)
(**)
(*Export[AuxPath<>"Tdotz5.m",Tdotz]*)


(* ::Text:: *)
(*Contraction with the integral and subsequent generation of ctab to speed up code*)


(* ::Input:: *)
(*Tdotint=Table[Sum[Tq[[\[Alpha],i,j,k,l,m]]qtensor[[i,j,k,l,m]],{i,Length[qv]},{j,Length[qv]},{k,Length[qv]},{l,Length[qv]},{m,Length[zh]}],*)
(*{\[Alpha],1,Length[Tq]}]/.kqsub//Simplify*)
(**)
(*Export[AuxPath<>"Tdotint5.m",Tdotint]*)


(* ::Input:: *)
(*Tdotintexp=Tdotint//Expand;*)
(*Tdotinttablist=getCtab/@Tdotintexp;*)
(*Tdotintctab=compressor/@Tdotinttablist;*)
(**)
(*Export[AuxPath<>"Tdotint5ctab.m",Tdotintctab]*)


(* ::Text:: *)
(*Contraction of T with itself*)


(* ::Input:: *)
(*TdotT=Table[Sum[Tk[[\[Alpha],i,j,k,l,m]] Tk[[\[Beta],i,j,k,l,m]],{i,Length[qv]},{j,Length[qv]},{k,Length[qv]},{l,Length[qv]},{m,Length[zh]}],*)
(*{\[Alpha],1,Length[Tk]},{\[Beta],1,Length[Tk]}]/.k1k2sub//Simplify;*)
(*Export[AuxPath<>"TdotT5.m",TdotT]*)


(* ::Subsubsection::Closed:: *)
(*Sixth order, \[Mu]0^6*)


(* ::Text:: *)
(*Generate tensors for the different contractions we need*)


(* ::Input:: *)
(*Tq=Tns6[k1vq,k2vq,IdentityMatrix[3]];*)
(*Tk=Tns6[k1vk,k2vk,IdentityMatrix[3]];*)
(*Tz=Tns6[k1vz,k2vz,IdentityMatrix[3]];*)
(*qtensor=Getvtensor[qv/q,6];*)
(*ztensor=Getvtensor[zh,6];*)
(*av=Table[Symbol["a"<>ToString[i]],{i,Length[Tk]}]*)


(* ::Text:: *)
(*Contraction of T's with zhat's*)


(* ::Input:: *)
(*Tdotz=Table[Sum[Tz[[\[Alpha],i,j,k,l,m,n]]ztensor[[i,j,k,l,m,n]],{i,Length[zh]},{j,Length[zh]},{k,Length[zh]},{l,Length[zh]},{m,Length[zh]},{n,Length[zh]}],*)
(*{\[Alpha],1,Length[Tz]}]//Simplify*)
(*Export[AuxPath<>"Tdotz6.m",Tdotz]*)


(* ::Text:: *)
(*Contraction with the integral*)


(* ::Input:: *)
(*Tdotint=Table[Sum[Tq[[\[Alpha],i,j,k,l,m,n]]qtensor[[i,j,k,l,m,n]],{i,Length[qv]},{j,Length[qv]},{k,Length[qv]},{l,Length[qv]},{m,Length[zh]},{n,Length[zh]}],*)
(*{\[Alpha],1,Length[Tq]}]/.kqsub//Simplify*)
(*Export[AuxPath<>"Tdotint6.m",Tdotint]*)


(* ::Input:: *)
(*Tdotintexp=Tdotint//Expand;*)
(*Tdotinttablist=getCtab/@Tdotintexp;*)
(*Tdotintctab=compressor/@Tdotinttablist;*)
(*Export[AuxPath<>"Tdotint6ctab.m",Tdotintctab]*)


(* ::Text:: *)
(*Contraction of T with itself*)


(* ::Input:: *)
(*TdotT=Table[Sum[Tk[[\[Alpha],i,j,k,l,m,n]] Tk[[\[Beta],i,j,k,l,m,n]],{i,Length[qv]},{j,Length[qv]},{k,Length[qv]},{l,Length[qv]},{m,Length[zh]},{n,Length[zh]}],*)
(*{\[Alpha],1,Length[Tk]},{\[Beta],1,Length[Tk]}]/.k1k2sub//Simplify;*)
(*Export[(AuxPath<>"TdotT6.m"),TdotT]*)


(* ::Subsubsection::Closed:: *)
(*Save everything up to \[Mu]0^4 (higher exponents are too slow, do numerically)*)


(*{1,\[Mu]0sol1,\[Mu]0sol2,\[Mu]0sol3,\[Mu]0sol4}>>(LoopTablesPath<>"mu0sol_nonsym_upto4.m");*)
Export[AuxPath<>"mu0sol_nonsym_upto4.m",{1,\[Mu]0sol1,\[Mu]0sol2,\[Mu]0sol3,\[Mu]0sol4}]


(* ::Subsection::Closed:: *)
(*Redshift space kernels, \[Mu]0 expansion dependent only on k1, up to \[Mu]0^2 (for B3212)*)


(* ::Subsubsection::Closed:: *)
(*First order, \[Mu]0^1*)


(* ::Text:: *)
(*Generate tensors for the different contractions we need*)


(* ::Input:: *)
(*Tq=Tv1[k1vq,IdentityMatrix[3]];*)
(*Tk=Tv1[k1vk,IdentityMatrix[3]];*)
(*Tz=Tv1[k1vz,IdentityMatrix[3]];*)
(*qtensor=Getvtensor[qv/q,1];*)
(*ztensor=Getvtensor[zh,1];*)
(*av=Table[Symbol["a"<>ToString[i]],{i,Length[Tk]}]*)


(* ::Text:: *)
(*Contraction of T's with zhat's*)


(* ::Input:: *)
(*Tdotz=Table[Sum[Tz[[\[Alpha],i]]ztensor[[i]],{i,Length[zh]}],{\[Alpha],1,Length[Tz]}]//Simplify*)


(* ::Text:: *)
(*Contraction with the integral*)


(* ::Input:: *)
(*Tdotint=Table[Sum[Tq[[\[Alpha],i]]qtensor[[i]],{i,Length[qv]}],*)
(*{\[Alpha],1,Length[Tq]}]//FullSimplify*)


(* ::Text:: *)
(*Contraction of T with itself, and inverse*)


(* ::Input:: *)
(*TdotT=Table[Sum[Tk[[\[Alpha],i]] Tk[[\[Beta],i]],{i,Length[qv]}],*)
(*{\[Alpha],1,Length[Tk]},{\[Beta],1,Length[Tk]}]/.k1k2sub//Simplify*)


(* ::Input:: *)
(*invT=Inverse[TdotT]//Simplify*)


(* ::Text:: *)
(*Solution*)


(* ::Input:: *)
(*\[Mu]0sol1=Tdotz.invT.Tdotint/.kqsub//Simplify*)


(* ::Subsubsection::Closed:: *)
(*Second order, \[Mu]0^2*)


(* ::Text:: *)
(*Generate tensors for the different contractions we need*)


(* ::Input:: *)
(*Tq=Tv2[k1vq,IdentityMatrix[3]];*)
(*Tk=Tv2[k1vk,IdentityMatrix[3]];*)
(*Tz=Tv2[k1vz,IdentityMatrix[3]];*)
(*qtensor=Getvtensor[qv/q,2];*)
(*ztensor=Getvtensor[zh,2];*)
(*av=Table[Symbol["a"<>ToString[i]],{i,Length[Tk]}]*)


(* ::Text:: *)
(*Contraction of T's with zhat's*)


(* ::Input:: *)
(*Tdotz=Table[Sum[Tz[[\[Alpha],i,j]]ztensor[[i,j]],{i,Length[zh]},{j,Length[zh]}],*)
(*{\[Alpha],1,Length[Tz]}]//Simplify*)


(* ::Text:: *)
(*Contraction with the integral*)


(* ::Input:: *)
(*Tdotint=Table[Sum[Tq[[\[Alpha],i,j]]qtensor[[i,j]],{i,Length[qv]},{j,Length[qv]}],*)
(*{\[Alpha],1,Length[Tq]}]//FullSimplify*)


(* ::Text:: *)
(*Contraction of T with itself, and inverse*)


(* ::Input:: *)
(*TdotT=Table[Sum[Tk[[\[Alpha],i,j]] Tk[[\[Beta],i,j]],{i,Length[qv]},{j,Length[qv]}],*)
(*{\[Alpha],1,Length[Tk]},{\[Beta],1,Length[Tk]}]//Simplify*)


(* ::Input:: *)
(*invT=Simplify[Inverse[TdotT]]/.k1k2sub//Simplify;*)


(* ::Text:: *)
(*Solution*)


(* ::Input:: *)
(*\[Mu]0sol2=Tdotz.invT.Tdotint/.kqsub//Simplify*)


(* ::Subsubsection:: *)
(*Save everything*)


(* ::Input:: *)
(*(*{1,\[Mu]0sol1,\[Mu]0sol2}>>LoopTablesPath<>"mu0solk1_vec_upto2.m"*)*)


Export[AuxPath<>"mu0solk1_vec_upto2.m",{1,\[Mu]0sol1,\[Mu]0sol2}]


(* ::Section:: *)
(*Loops*)


(*Formating jmat imports*)
fi[k_] := Round[k, N[10^-5]]; 
Tfi[k_] := ToString[fi[k]];

saveLoopCoefPath[triangle_,diagram_String]:=(saveLoopCoefsPath <> diagram <>"/"
<>Tfi[triangle[[1]]]<>"_"<>Tfi[triangle[[2]]]<>"_"<>Tfi[triangle[[3]]]<>"_.m");


(* ::Subsection::Closed:: *)
(*Useful rules to efficiently integrate the \[Mu]0 *)


(* ::Input:: *)
(*angint[a_,b_,c_,d_]=((1+(-1)^a) (1+(-1)^c) (1+(-1)^d) Gamma[(1+a)/2] Gamma[1+b/2] Gamma[(1+c)/2] Gamma[(1+d)/2])/(16 \[Pi] Gamma[1/2 (3+a+b)] Gamma[1/2 (2+c+d)]);*)


(* ::Input:: *)
(*rules\[Mu]04={x^a_. (Sqrt[1-x^2])^b_. Cos[\[Beta]]^c_. Sin[\[Beta]]^d_.->angint[a,b,c,d]};*)
(*rules\[Mu]03={x^a_. (Sqrt[1-x^2])^b_. Cos[\[Beta]]^c_.->angint[a,b,c,0],*)
(*x^a_. (Sqrt[1-x^2])^b_. Sin[\[Beta]]^d_.->angint[a,b,0,d],*)
(*x^a_. Cos[\[Beta]]^c_. Sin[\[Beta]]^d_.->angint[a,0,c,d],*)
(*(Sqrt[1-x^2])^b_. Cos[\[Beta]]^c_. Sin[\[Beta]]^d_.->angint[0,b,c,d]};*)
(*rules\[Mu]02={x^a_. (Sqrt[1-x^2])^b_.->angint[a,b,0,0],*)
(*x^a_. Cos[\[Beta]]^c_.->angint[a,0,c,0],*)
(* (Sqrt[1-x^2])^b_. Cos[\[Beta]]^c_.->angint[0,b,c,0],*)
(*x^a_. Sin[\[Beta]]^d_.->angint[a,0,0,d],*)
(* (Sqrt[1-x^2])^b_. Sin[\[Beta]]^d_.->angint[0,b,0,d],*)
(*Cos[\[Beta]]^c_. Sin[\[Beta]]^d_.->angint[0,0,c,d]};*)
(*rules\[Mu]01={x^a_.->angint[a,0,0,0],(Sqrt[1-x^2])^b_.->angint[0,b,0,0], Cos[\[Beta]]^c_.->angint[0,0,c,0],Sin[\[Beta]]^d_.->angint[0,0,0,d]};*)


(* ::Subsection::Closed:: *)
(*Btree with permutations*)


Clear[k211,Btree]
k211=2K1r[k1] K1r[k2]K2rsym[k1,k2]/.tb;
Btree=k211 Plin[k1]Plin[k2];


Variables[Btree]


(*The following code does the permutation*)

(*arr  = {k1,\[Mu]1,k2,\[Mu]2,k3,\[Mu]3};
arrt = {k1t,\[Mu]1t,k2t,\[Mu]2t,k3t,\[Mu]3t};
cycpt=NestList[RotateLeft,
			   {{k1t,\[Mu]1t},{k2t,\[Mu]2t},{k3t,\[Mu]3t}},2]/.{{a_,b_},{c_,d_},{e_,f_}}-> {a,b,c,d,e,f};
permreps=Thread[(arr->#)]&/@cycpt;*)
permreps = {{k1->k1t,\[Mu]1->\[Mu]1t,k2->k2t,\[Mu]2->\[Mu]2t,k3->k3t,\[Mu]3->\[Mu]3t},
			{k1->k2t,\[Mu]1->\[Mu]2t,k2->k3t,\[Mu]2->\[Mu]3t,k3->k1t,\[Mu]3->\[Mu]1t},
			{k1->k3t,\[Mu]1->\[Mu]3t,k2->k1t,\[Mu]2->\[Mu]1t,k3->k2t,\[Mu]3->\[Mu]2t}};


BtreePerms=Total[Btree/.permreps/.Thread[arrt->arr]]//Expand//Simplify;


(*check*)
Bpermcheck=(BtreePerms/.Thread[{k1,\[Mu]1,k2,\[Mu]2,k3,\[Mu]3}->{1,Sqrt[3]/2,1,Sqrt[3]/2,1,Sqrt[3]/2}])/3;
Bcheck=Btree/.Thread[{k1,\[Mu]1,k2,\[Mu]2,k3,\[Mu]3}->{1,Sqrt[3]/2,1,Sqrt[3]/2,1,Sqrt[3]/2}];
Bpermcheck==Bcheck//Simplify


Export[LoopTablesPath<>"BtreePerms.m",BtreePerms]


(* ::Subsection:: *)
(*B222*)


(* ::Subsubsection::Closed:: *)
(*Preliminaries*)


k222[k1_,k2_,q_]:=k222[k1,k2,q]=K2rsym[-q,q+k1]K2rsym[k1+q,k2-q]K2rsym[k2-q,q];
B222Integrand[k1_,k2_,q_] := B222Integrand[k1,k2,q]=8k222[k1,k2,q]*Plin[q]Plin[kmq[k2]]Plin[kpq[k1]]/.tb;


(* ::Text:: *)
(*We simply put Plin -> 1, as there are all 3 present and we don' t do redefinitions of q*)


ker222 =(B222Integrand[k1,k2,q]/.Plin[a_]-> 1);
ker222exp=ker222//Expand;
ker222expmu0=CoefficientList[ker222exp,\[Mu]0];
\[Mu]0tab222=Table[\[Mu]0^(i-1),{i,1,Length[ker222expmu0]}];

\[Mu]0solb222=Import[AuxPath<>"mu0sol_nonsym_upto4.m"];
b222temptomu4=ker222expmu0[[1;;5]].\[Mu]0solb222;
b222tomu4=Join[{b222temptomu4},ker222expmu0[[{6,7}]]];


(*checks*)
\[Mu]0tab222=={1,\[Mu]0,\[Mu]0^2,\[Mu]0^3,\[Mu]0^4,\[Mu]0^5,\[Mu]0^6}
ker222expmu0.\[Mu]0tab222==ker222exp//Expand


(* ::Subsubsection::Closed:: *)
(*Simplify first part (up to \[Mu]0^4) *)


(*b222part1exp=b222tomu4\[LeftDoubleBracket]1\[RightDoubleBracket]//Expand;
b222part1ctab=Sort[compressor[getCtab[b222part1exp]]];
Export[AuxPath<>"b222mu0to4ctab.m",b222part1ctab]*)


b222part1ctab=Import[AuxPath<>"b222mu0to4ctab.m"];
b222part1exps=b222part1ctab[[All,1;;3]];
b222part1coefs=b222part1ctab[[All,4]];


(* ::Subsubsection::Closed:: *)
(*Simplify second part (\[Mu]0^5)*)


(* ::Input:: *)
(*Tint\[Mu]5=Import[AuxPath<>"Tdotint5.m"];*)
(*Tint\[Mu]5ctab=Import[AuxPath<>"Tdotint5ctab.m"];*)
(*TT\[Mu]5=Import[AuxPath<>"TdotT5.m"];*)
(*Tz\[Mu]5=Import[AuxPath<>"Tdotz5.m"];*)


(* ::Input:: *)
(*Tint5fn[k1_,k2_,k3_]=Tint\[Mu]5;*)
(*Tint5ctabfn[k1_,k2_,k3_]=Tint\[Mu]5ctab;*)
(*TT5fn[k1_,k2_,k3_]=TT\[Mu]5;*)
(*TTinv5fn[k1_,k2_,k3_]:=Inverse[SetPrecision[TT5fn[k1,k2,k3],20]];*)
(*Tz5fn[k1_,k2_,k3_]=Tz\[Mu]5;*)
(*TzTTinv5fn[k1_,k2_,k3_]:=Tz5fn[k1,k2,k3].TTinv5fn[k1,k2,k3];*)
(**)
(*Tint5exps=Tint\[Mu]5ctab[[All,All,1;;3]];*)


(* ::Input:: *)
(*(*express \[Mu]0^5 with ctab - much faster*)*)
(*\[Mu]0ctab5[k1_,k2_,k3_]:=Module[{coefs,Tintval,newctab},*)
(*Tintval=Tint5ctabfn[k1,k2,k3];*)
(*coefs=TzTTinv5fn[k1,k2,k3]*Tintval[[All,All,4]];*)
(*newctab=Flatten[MapThread[MapThread[Append,{#1,#2}]&,{Tint5exps,coefs}],1];*)
(*compressor[newctab//Expand]];*)


(*Check against standard calculation of \[Mu]0*)
Clear[\[Mu]0sol5]
\[Mu]0sol5[k1_,k2_,k3_]:=Tz5fn[k1,k2,k3].TTinv5fn[k1,k2,k3].Tint5fn[k1,k2,k3]//Simplify;
\[Mu]0sol5check=\[Mu]0sol5[0.1,0.11,0.12];
\[Mu]0ctab5check=Sort[\[Mu]0ctab5[0.1,0.11,0.12]];
\[Mu]0sol5checkctab=Sort[compressor[getCtab[\[Mu]0sol5check//Expand]]];

\[Mu]0ctab5check[[All,1;;3]]==\[Mu]0sol5checkctab[[All,1;;3]]


(* ::Input:: *)
(*b222part2exp=b222tomu4[[2]]//Expand;*)
(*b222part2expfn[k1_,k2_,k3_]=b222part2exp;*)
(*b222part2ctabfn[k1_,k2_,k3_]=compressor[getCtab[b222part2exp]];*)


(* ::Input:: *)
(*(*Function to generate the ctab for a given set of k1, k2, k3*)*)
(*b222\[Mu]5ctab[k1_,k2_,k3_]:=Module[{b222coefctab,\[Mu]0ctab,jointab},*)
(*b222coefctab=b222part2ctabfn[k1,k2,k3];*)
(*\[Mu]0ctab=\[Mu]0ctab5[k1,k2,k3];*)
(*jointab=JoinCtabs[b222coefctab,\[Mu]0ctab]//Expand*)
(*];*)


(* ::Subsubsection::Closed:: *)
(*Simplify third part (\[Mu]0^6)*)


(* ::Input:: *)
(*Tint\[Mu]6=Import[(AuxPath<>"Tdotint6.m")];*)
(*Tint\[Mu]6ctab=Import[AuxPath<>"Tdotint6ctab.m"];*)
(*TT\[Mu]6=Import[(AuxPath<>"TdotT6.m")];*)
(*Tz\[Mu]6=Import[(AuxPath<>"Tdotz6.m")];*)
(**)
(*Tint6fn[k1_,k2_,k3_]=Tint\[Mu]6;*)
(*Tint6ctabfn[k1_,k2_,k3_]=Tint\[Mu]6ctab;*)
(*TT6fn[k1_,k2_,k3_]=TT\[Mu]6;*)
(*TTinv6fn[k1_,k2_,k3_]:=Inverse[SetPrecision[TT6fn[k1,k2,k3],20]];*)
(*Tz6fn[k1_,k2_,k3_]=Tz\[Mu]6;*)
(**)
(*TzTTinv6fn[k1_,k2_,k3_]:=Tz6fn[k1,k2,k3].TTinv6fn[k1,k2,k3];*)
(**)
(*Tint6exps=Tint\[Mu]6ctab[[All,All,1;;3]];*)


(* ::Input:: *)
(*Tint6coefs=Tint\[Mu]6ctab[[All,All,4]];*)
(*Tint6coefsfn[k1_,k2_]=Tint6coefs;*)


(* ::Input:: *)
(*(*express \[Mu]0^6 with ctab*)*)
(*\[Mu]0ctab6[k1_,k2_,k3_]:=Module[{coefs,Tintval,newctab},*)
(*coefs=TzTTinv6fn[k1,k2,k3]*Tint6coefsfn[k1,k2];*)
(*newctab=Flatten[MapThread[MapThread[Append,{#1,#2}]&,{Tint6exps,coefs}],1];*)
(*compressor[newctab]//Expand*)
(*];*)


(* ::Input:: *)
(*(**)
(*Clear[\[Mu]0sol6]*)
(*\[Mu]0sol6[k1_,k2_,k3_]:=(Tz6fn[k1,k2,k3].TTinv6fn[k1,k2,k3]).Tint6fn[k1,k2,k3]//Simplify;*)
(**)*)


(* ::Input:: *)
(*b222part3exp=b222tomu4[[3]]//Expand;*)
(*b222part3expfn[k1_,k2_,k3_]=b222part3exp;*)
(*b222part3ctabfn[k1_,k2_,k3_]=compressor[getCtab[b222part3exp]];*)


(* ::Input:: *)
(*(*Function to generate the ctab for a given set of k1, k2, k3*)*)
(*b222\[Mu]6ctab[k1_,k2_,k3_]:=Module[{b222coefctab,\[Mu]0ctab,jointab},*)
(*b222coefctab=b222part3ctabfn[k1,k2,k3];*)
(*\[Mu]0ctab=\[Mu]0ctab6[k1,k2,k3];*)
(*jointab=JoinCtabs[b222coefctab,\[Mu]0ctab]//Expand*)
(*];*)


(* ::Subsubsection::Closed:: *)
(*Join \[Mu]0^5 and \[Mu]0^6 ctabs*)


(* ::Input:: *)
(*(*Join \[Mu]5 and \[Mu]6 ctabs*)*)
(*b222\[Mu]56ctab[k1_,k2_,k3_]:=compressor[Join[b222\[Mu]5ctab[k1,k2,k3],b222\[Mu]6ctab[k1,k2,k3]]];*)


(* ::Subsubsection:: *)
(*Decompose into bias coefficients*)


(*Create example of ctab to obtain the bias monomials and variables*)
b222part23ctab=b222\[Mu]56ctab[0.11,0.12,0.13];

(*b222part1coefs//Variables*)
(*b222part23ctab//Variables*)


b222vars={b1,b2,b5,f1,\[Mu]1,\[Mu]2};


(* ::Input:: *)
(*(*Get all bias monomials that appear in the b222 kernel*)*)
(*(**)
(*b222biaslistpart1=getMonomials[b222vars,b222part1ctab];*)
(*b222biaslistpart23=getMonomials[b222vars,b222part23ctab];*)
(*b222biaslist=Union[b222biaslistpart1,b222biaslistpart23];*)
(**)*)


(* ::Input:: *)
(*b222biaslist={b1^3,b1^2 b2,b1 b2^2,b2^3,b1^2 b5,b1 b2 b5,b2^2 b5,b1 b5^2,b2 b5^2,b5^3,b1^2 f1 \[Mu]1^2,b1^3 f1 \[Mu]1^2,b1 b2 f1 \[Mu]1^2,b1^2 b2 f1 \[Mu]1^2,b2^2 f1 \[Mu]1^2,b1 b2^2 f1 \[Mu]1^2,b1 b5 f1 \[Mu]1^2,b1^2 b5 f1 \[Mu]1^2,b2 b5 f1 \[Mu]1^2,b1 b2 b5 f1 \[Mu]1^2,b5^2 f1 \[Mu]1^2,b1 b5^2 f1 \[Mu]1^2,b1^2 f1^2 \[Mu]1^2,b1^3 f1^2 \[Mu]1^2,b1 b2 f1^2 \[Mu]1^2,b1^2 b2 f1^2 \[Mu]1^2,b2^2 f1^2 \[Mu]1^2,b1 b5 f1^2 \[Mu]1^2,b1^2 b5 f1^2 \[Mu]1^2,b2 b5 f1^2 \[Mu]1^2,b5^2 f1^2 \[Mu]1^2,b1 f1^2 \[Mu]1^4,b1^2 f1^2 \[Mu]1^4,b1^3 f1^2 \[Mu]1^4,b2 f1^2 \[Mu]1^4,b1 b2 f1^2 \[Mu]1^4,b1^2 b2 f1^2 \[Mu]1^4,b2^2 f1^2 \[Mu]1^4,b5 f1^2 \[Mu]1^4,b1 b5 f1^2 \[Mu]1^4,b1^2 b5 f1^2 \[Mu]1^4,b2 b5 f1^2 \[Mu]1^4,b5^2 f1^2 \[Mu]1^4,b1 f1^3 \[Mu]1^4,b1^2 f1^3 \[Mu]1^4,b2 f1^3 \[Mu]1^4,b1 b2 f1^3 \[Mu]1^4,b5 f1^3 \[Mu]1^4,b1 b5 f1^3 \[Mu]1^4,b1 f1^4 \[Mu]1^4,b2 f1^4 \[Mu]1^4,b5 f1^4 \[Mu]1^4,b1 f1^3 \[Mu]1^6,b1^2 f1^3 \[Mu]1^6,b2 f1^3 \[Mu]1^6,b1 b2 f1^3 \[Mu]1^6,b5 f1^3 \[Mu]1^6,b1 b5 f1^3 \[Mu]1^6,b1 f1^4 \[Mu]1^6,b2 f1^4 \[Mu]1^6,b5 f1^4 \[Mu]1^6,b1 f1^4 \[Mu]1^8,b2 f1^4 \[Mu]1^8,b5 f1^4 \[Mu]1^8,b1^2 f1 \[Mu]1 \[Mu]2,b1^3 f1 \[Mu]1 \[Mu]2,b1 b2 f1 \[Mu]1 \[Mu]2,b1^2 b2 f1 \[Mu]1 \[Mu]2,b2^2 f1 \[Mu]1 \[Mu]2,b1 b2^2 f1 \[Mu]1 \[Mu]2,b1 b5 f1 \[Mu]1 \[Mu]2,b1^2 b5 f1 \[Mu]1 \[Mu]2,b2 b5 f1 \[Mu]1 \[Mu]2,b1 b2 b5 f1 \[Mu]1 \[Mu]2,b5^2 f1 \[Mu]1 \[Mu]2,b1 b5^2 f1 \[Mu]1 \[Mu]2,b1^2 f1^2 \[Mu]1 \[Mu]2,b1^3 f1^2 \[Mu]1 \[Mu]2,b1 b2 f1^2 \[Mu]1 \[Mu]2,b1^2 b2 f1^2 \[Mu]1 \[Mu]2,b2^2 f1^2 \[Mu]1 \[Mu]2,b1 b5 f1^2 \[Mu]1 \[Mu]2,b1^2 b5 f1^2 \[Mu]1 \[Mu]2,b2 b5 f1^2 \[Mu]1 \[Mu]2,b5^2 f1^2 \[Mu]1 \[Mu]2,b1 f1^2 \[Mu]1^3 \[Mu]2,b1^2 f1^2 \[Mu]1^3 \[Mu]2,b1^3 f1^2 \[Mu]1^3 \[Mu]2,b2 f1^2 \[Mu]1^3 \[Mu]2,b1 b2 f1^2 \[Mu]1^3 \[Mu]2,b1^2 b2 f1^2 \[Mu]1^3 \[Mu]2,b2^2 f1^2 \[Mu]1^3 \[Mu]2,b5 f1^2 \[Mu]1^3 \[Mu]2,b1 b5 f1^2 \[Mu]1^3 \[Mu]2,b1^2 b5 f1^2 \[Mu]1^3 \[Mu]2,b2 b5 f1^2 \[Mu]1^3 \[Mu]2,b5^2 f1^2 \[Mu]1^3 \[Mu]2,b1 f1^3 \[Mu]1^3 \[Mu]2,b1^2 f1^3 \[Mu]1^3 \[Mu]2,b1^3 f1^3 \[Mu]1^3 \[Mu]2,b2 f1^3 \[Mu]1^3 \[Mu]2,b1 b2 f1^3 \[Mu]1^3 \[Mu]2,b5 f1^3 \[Mu]1^3 \[Mu]2,b1 b5 f1^3 \[Mu]1^3 \[Mu]2,b1 f1^4 \[Mu]1^3 \[Mu]2,b1^2 f1^4 \[Mu]1^3 \[Mu]2,b2 f1^4 \[Mu]1^3 \[Mu]2,b5 f1^4 \[Mu]1^3 \[Mu]2,b1 f1^3 \[Mu]1^5 \[Mu]2,b1^2 f1^3 \[Mu]1^5 \[Mu]2,b1^3 f1^3 \[Mu]1^5 \[Mu]2,b2 f1^3 \[Mu]1^5 \[Mu]2,b1 b2 f1^3 \[Mu]1^5 \[Mu]2,b5 f1^3 \[Mu]1^5 \[Mu]2,b1 b5 f1^3 \[Mu]1^5 \[Mu]2,b1 f1^4 \[Mu]1^5 \[Mu]2,b1^2 f1^4 \[Mu]1^5 \[Mu]2,b2 f1^4 \[Mu]1^5 \[Mu]2,b5 f1^4 \[Mu]1^5 \[Mu]2,b1 f1^5 \[Mu]1^5 \[Mu]2,b1 f1^4 \[Mu]1^7 \[Mu]2,b1^2 f1^4 \[Mu]1^7 \[Mu]2,b2 f1^4 \[Mu]1^7 \[Mu]2,b5 f1^4 \[Mu]1^7 \[Mu]2,b1 f1^5 \[Mu]1^7 \[Mu]2,b1 f1^5 \[Mu]1^9 \[Mu]2,b1^2 f1 \[Mu]2^2,b1^3 f1 \[Mu]2^2,b1 b2 f1 \[Mu]2^2,b1^2 b2 f1 \[Mu]2^2,b2^2 f1 \[Mu]2^2,b1 b2^2 f1 \[Mu]2^2,b1 b5 f1 \[Mu]2^2,b1^2 b5 f1 \[Mu]2^2,b2 b5 f1 \[Mu]2^2,b1 b2 b5 f1 \[Mu]2^2,b5^2 f1 \[Mu]2^2,b1 b5^2 f1 \[Mu]2^2,b1^2 f1^2 \[Mu]2^2,b1^3 f1^2 \[Mu]2^2,b1 b2 f1^2 \[Mu]2^2,b1^2 b2 f1^2 \[Mu]2^2,b2^2 f1^2 \[Mu]2^2,b1 b5 f1^2 \[Mu]2^2,b1^2 b5 f1^2 \[Mu]2^2,b2 b5 f1^2 \[Mu]2^2,b5^2 f1^2 \[Mu]2^2,b1 f1^2 \[Mu]1^2 \[Mu]2^2,b1^2 f1^2 \[Mu]1^2 \[Mu]2^2,b1^3 f1^2 \[Mu]1^2 \[Mu]2^2,b2 f1^2 \[Mu]1^2 \[Mu]2^2,b1 b2 f1^2 \[Mu]1^2 \[Mu]2^2,b1^2 b2 f1^2 \[Mu]1^2 \[Mu]2^2,b2^2 f1^2 \[Mu]1^2 \[Mu]2^2,b5 f1^2 \[Mu]1^2 \[Mu]2^2,b1 b5 f1^2 \[Mu]1^2 \[Mu]2^2,b1^2 b5 f1^2 \[Mu]1^2 \[Mu]2^2,b2 b5 f1^2 \[Mu]1^2 \[Mu]2^2,b5^2 f1^2 \[Mu]1^2 \[Mu]2^2,b1 f1^3 \[Mu]1^2 \[Mu]2^2,b1^2 f1^3 \[Mu]1^2 \[Mu]2^2,b1^3 f1^3 \[Mu]1^2 \[Mu]2^2,b2 f1^3 \[Mu]1^2 \[Mu]2^2,b1 b2 f1^3 \[Mu]1^2 \[Mu]2^2,b5 f1^3 \[Mu]1^2 \[Mu]2^2,b1 b5 f1^3 \[Mu]1^2 \[Mu]2^2,b1 f1^4 \[Mu]1^2 \[Mu]2^2,b1^2 f1^4 \[Mu]1^2 \[Mu]2^2,b2 f1^4 \[Mu]1^2 \[Mu]2^2,b5 f1^4 \[Mu]1^2 \[Mu]2^2,f1^3 \[Mu]1^4 \[Mu]2^2,b1 f1^3 \[Mu]1^4 \[Mu]2^2,b1^2 f1^3 \[Mu]1^4 \[Mu]2^2,b1^3 f1^3 \[Mu]1^4 \[Mu]2^2,b2 f1^3 \[Mu]1^4 \[Mu]2^2,b1 b2 f1^3 \[Mu]1^4 \[Mu]2^2,b5 f1^3 \[Mu]1^4 \[Mu]2^2,b1 b5 f1^3 \[Mu]1^4 \[Mu]2^2,f1^4 \[Mu]1^4 \[Mu]2^2,b1 f1^4 \[Mu]1^4 \[Mu]2^2,b1^2 f1^4 \[Mu]1^4 \[Mu]2^2,b2 f1^4 \[Mu]1^4 \[Mu]2^2,b5 f1^4 \[Mu]1^4 \[Mu]2^2,f1^5 \[Mu]1^4 \[Mu]2^2,b1 f1^5 \[Mu]1^4 \[Mu]2^2,f1^6 \[Mu]1^4 \[Mu]2^2,f1^4 \[Mu]1^6 \[Mu]2^2,b1 f1^4 \[Mu]1^6 \[Mu]2^2,b1^2 f1^4 \[Mu]1^6 \[Mu]2^2,b2 f1^4 \[Mu]1^6 \[Mu]2^2,b5 f1^4 \[Mu]1^6 \[Mu]2^2,f1^5 \[Mu]1^6 \[Mu]2^2,b1 f1^5 \[Mu]1^6 \[Mu]2^2,f1^6 \[Mu]1^6 \[Mu]2^2,f1^5 \[Mu]1^8 \[Mu]2^2,b1 f1^5 \[Mu]1^8 \[Mu]2^2,f1^6 \[Mu]1^8 \[Mu]2^2,f1^6 \[Mu]1^10 \[Mu]2^2,b1 f1^2 \[Mu]1 \[Mu]2^3,b1^2 f1^2 \[Mu]1 \[Mu]2^3,b1^3 f1^2 \[Mu]1 \[Mu]2^3,b2 f1^2 \[Mu]1 \[Mu]2^3,b1 b2 f1^2 \[Mu]1 \[Mu]2^3,b1^2 b2 f1^2 \[Mu]1 \[Mu]2^3,b2^2 f1^2 \[Mu]1 \[Mu]2^3,b5 f1^2 \[Mu]1 \[Mu]2^3,b1 b5 f1^2 \[Mu]1 \[Mu]2^3,b1^2 b5 f1^2 \[Mu]1 \[Mu]2^3,b2 b5 f1^2 \[Mu]1 \[Mu]2^3,b5^2 f1^2 \[Mu]1 \[Mu]2^3,b1 f1^3 \[Mu]1 \[Mu]2^3,b1^2 f1^3 \[Mu]1 \[Mu]2^3,b1^3 f1^3 \[Mu]1 \[Mu]2^3,b2 f1^3 \[Mu]1 \[Mu]2^3,b1 b2 f1^3 \[Mu]1 \[Mu]2^3,b5 f1^3 \[Mu]1 \[Mu]2^3,b1 b5 f1^3 \[Mu]1 \[Mu]2^3,b1 f1^4 \[Mu]1 \[Mu]2^3,b1^2 f1^4 \[Mu]1 \[Mu]2^3,b2 f1^4 \[Mu]1 \[Mu]2^3,b5 f1^4 \[Mu]1 \[Mu]2^3,f1^3 \[Mu]1^3 \[Mu]2^3,b1 f1^3 \[Mu]1^3 \[Mu]2^3,b1^2 f1^3 \[Mu]1^3 \[Mu]2^3,b1^3 f1^3 \[Mu]1^3 \[Mu]2^3,b2 f1^3 \[Mu]1^3 \[Mu]2^3,b1 b2 f1^3 \[Mu]1^3 \[Mu]2^3,b5 f1^3 \[Mu]1^3 \[Mu]2^3,b1 b5 f1^3 \[Mu]1^3 \[Mu]2^3,f1^4 \[Mu]1^3 \[Mu]2^3,b1 f1^4 \[Mu]1^3 \[Mu]2^3,b1^2 f1^4 \[Mu]1^3 \[Mu]2^3,b2 f1^4 \[Mu]1^3 \[Mu]2^3,b5 f1^4 \[Mu]1^3 \[Mu]2^3,f1^5 \[Mu]1^3 \[Mu]2^3,b1 f1^5 \[Mu]1^3 \[Mu]2^3,f1^6 \[Mu]1^3 \[Mu]2^3,f1^4 \[Mu]1^5 \[Mu]2^3,b1 f1^4 \[Mu]1^5 \[Mu]2^3,b1^2 f1^4 \[Mu]1^5 \[Mu]2^3,b2 f1^4 \[Mu]1^5 \[Mu]2^3,b5 f1^4 \[Mu]1^5 \[Mu]2^3,f1^5 \[Mu]1^5 \[Mu]2^3,b1 f1^5 \[Mu]1^5 \[Mu]2^3,f1^6 \[Mu]1^5 \[Mu]2^3,f1^5 \[Mu]1^7 \[Mu]2^3,b1 f1^5 \[Mu]1^7 \[Mu]2^3,f1^6 \[Mu]1^7 \[Mu]2^3,f1^6 \[Mu]1^9 \[Mu]2^3,b1 f1^2 \[Mu]2^4,b1^2 f1^2 \[Mu]2^4,b1^3 f1^2 \[Mu]2^4,b2 f1^2 \[Mu]2^4,b1 b2 f1^2 \[Mu]2^4,b1^2 b2 f1^2 \[Mu]2^4,b2^2 f1^2 \[Mu]2^4,b5 f1^2 \[Mu]2^4,b1 b5 f1^2 \[Mu]2^4,b1^2 b5 f1^2 \[Mu]2^4,b2 b5 f1^2 \[Mu]2^4,b5^2 f1^2 \[Mu]2^4,b1 f1^3 \[Mu]2^4,b1^2 f1^3 \[Mu]2^4,b2 f1^3 \[Mu]2^4,b1 b2 f1^3 \[Mu]2^4,b5 f1^3 \[Mu]2^4,b1 b5 f1^3 \[Mu]2^4,b1 f1^4 \[Mu]2^4,b2 f1^4 \[Mu]2^4,b5 f1^4 \[Mu]2^4,f1^3 \[Mu]1^2 \[Mu]2^4,b1 f1^3 \[Mu]1^2 \[Mu]2^4,b1^2 f1^3 \[Mu]1^2 \[Mu]2^4,b1^3 f1^3 \[Mu]1^2 \[Mu]2^4,b2 f1^3 \[Mu]1^2 \[Mu]2^4,b1 b2 f1^3 \[Mu]1^2 \[Mu]2^4,b5 f1^3 \[Mu]1^2 \[Mu]2^4,b1 b5 f1^3 \[Mu]1^2 \[Mu]2^4,f1^4 \[Mu]1^2 \[Mu]2^4,b1 f1^4 \[Mu]1^2 \[Mu]2^4,b1^2 f1^4 \[Mu]1^2 \[Mu]2^4,b2 f1^4 \[Mu]1^2 \[Mu]2^4,b5 f1^4 \[Mu]1^2 \[Mu]2^4,f1^5 \[Mu]1^2 \[Mu]2^4,b1 f1^5 \[Mu]1^2 \[Mu]2^4,f1^6 \[Mu]1^2 \[Mu]2^4,f1^4 \[Mu]1^4 \[Mu]2^4,b1 f1^4 \[Mu]1^4 \[Mu]2^4,b1^2 f1^4 \[Mu]1^4 \[Mu]2^4,b2 f1^4 \[Mu]1^4 \[Mu]2^4,b5 f1^4 \[Mu]1^4 \[Mu]2^4,f1^5 \[Mu]1^4 \[Mu]2^4,b1 f1^5 \[Mu]1^4 \[Mu]2^4,f1^6 \[Mu]1^4 \[Mu]2^4,f1^5 \[Mu]1^6 \[Mu]2^4,b1 f1^5 \[Mu]1^6 \[Mu]2^4,f1^6 \[Mu]1^6 \[Mu]2^4,f1^6 \[Mu]1^8 \[Mu]2^4,b1 f1^3 \[Mu]1 \[Mu]2^5,b1^2 f1^3 \[Mu]1 \[Mu]2^5,b1^3 f1^3 \[Mu]1 \[Mu]2^5,b2 f1^3 \[Mu]1 \[Mu]2^5,b1 b2 f1^3 \[Mu]1 \[Mu]2^5,b5 f1^3 \[Mu]1 \[Mu]2^5,b1 b5 f1^3 \[Mu]1 \[Mu]2^5,b1 f1^4 \[Mu]1 \[Mu]2^5,b1^2 f1^4 \[Mu]1 \[Mu]2^5,b2 f1^4 \[Mu]1 \[Mu]2^5,b5 f1^4 \[Mu]1 \[Mu]2^5,b1 f1^5 \[Mu]1 \[Mu]2^5,f1^4 \[Mu]1^3 \[Mu]2^5,b1 f1^4 \[Mu]1^3 \[Mu]2^5,b1^2 f1^4 \[Mu]1^3 \[Mu]2^5,b2 f1^4 \[Mu]1^3 \[Mu]2^5,b5 f1^4 \[Mu]1^3 \[Mu]2^5,f1^5 \[Mu]1^3 \[Mu]2^5,b1 f1^5 \[Mu]1^3 \[Mu]2^5,f1^6 \[Mu]1^3 \[Mu]2^5,f1^5 \[Mu]1^5 \[Mu]2^5,b1 f1^5 \[Mu]1^5 \[Mu]2^5,f1^6 \[Mu]1^5 \[Mu]2^5,f1^6 \[Mu]1^7 \[Mu]2^5,b1 f1^3 \[Mu]2^6,b1^2 f1^3 \[Mu]2^6,b2 f1^3 \[Mu]2^6,b1 b2 f1^3 \[Mu]2^6,b5 f1^3 \[Mu]2^6,b1 b5 f1^3 \[Mu]2^6,b1 f1^4 \[Mu]2^6,b2 f1^4 \[Mu]2^6,b5 f1^4 \[Mu]2^6,f1^4 \[Mu]1^2 \[Mu]2^6,b1 f1^4 \[Mu]1^2 \[Mu]2^6,b1^2 f1^4 \[Mu]1^2 \[Mu]2^6,b2 f1^4 \[Mu]1^2 \[Mu]2^6,b5 f1^4 \[Mu]1^2 \[Mu]2^6,f1^5 \[Mu]1^2 \[Mu]2^6,b1 f1^5 \[Mu]1^2 \[Mu]2^6,f1^6 \[Mu]1^2 \[Mu]2^6,f1^5 \[Mu]1^4 \[Mu]2^6,b1 f1^5 \[Mu]1^4 \[Mu]2^6,f1^6 \[Mu]1^4 \[Mu]2^6,f1^6 \[Mu]1^6 \[Mu]2^6,b1 f1^4 \[Mu]1 \[Mu]2^7,b1^2 f1^4 \[Mu]1 \[Mu]2^7,b2 f1^4 \[Mu]1 \[Mu]2^7,b5 f1^4 \[Mu]1 \[Mu]2^7,b1 f1^5 \[Mu]1 \[Mu]2^7,f1^5 \[Mu]1^3 \[Mu]2^7,b1 f1^5 \[Mu]1^3 \[Mu]2^7,f1^6 \[Mu]1^3 \[Mu]2^7,f1^6 \[Mu]1^5 \[Mu]2^7,b1 f1^4 \[Mu]2^8,b2 f1^4 \[Mu]2^8,b5 f1^4 \[Mu]2^8,f1^5 \[Mu]1^2 \[Mu]2^8,b1 f1^5 \[Mu]1^2 \[Mu]2^8,f1^6 \[Mu]1^2 \[Mu]2^8,f1^6 \[Mu]1^4 \[Mu]2^8,b1 f1^5 \[Mu]1 \[Mu]2^9,f1^6 \[Mu]1^3 \[Mu]2^9,f1^6 \[Mu]1^2 \[Mu]2^10};*)


Export[LoopTablesPath<>"B222simpcoefslist.m",b222biaslist]


(* ::Input:: *)
(*(*Function that obtains the coefficient of each term in biaslist in each ctab coefficient "coef" *)*)
(*b222biaslister=biaslisterfn[b222vars,b222biaslist];*)


(* ::Input:: *)
(*(*apply decomposition in biases to all ctab coefficients and simplify - this takes 7 minutes*)*)
(*b222biascoef=Simplify[b222biaslister/@b222part1coefs];*)
(*b222part1ctabdec=MapThread[Append,{b222part1exps,b222biascoef}];*)


(*Save decomposed ctab*)

(*b222part1ctabdec>>(LoopTablesPath<>"b222mu0to4ctabdec.mx");*)
(*Export[AuxPath<>"b222mu0to4ctabdec.m",b222part1ctabdec]*)


(*Load decomposed ctab*)

(*b222part1ctabdec=<<(LoopTablesPath<>"b222mu0to4ctabdec.mx");*)
b222part1ctabdec = Import[AuxPath<>"b222mu0to4ctabdec.m"];


(* ::Input:: *)
(*b222\[Mu]0to4ctabdec[k1_,k2_,k3_]=b222part1ctabdec;*)


(* ::Input:: *)
(*(*apply decomposition in biases to all ctab coefficients*)*)
(*b222\[Mu]56ctabdec[k1_,k2_,k3_]:=Module[{b222ctabtemp,b222exps,b222coeflist},*)
(*b222ctabtemp=b222\[Mu]56ctab[k1,k2,k3];*)
(*b222exps=b222ctabtemp[[All,1;;3]];*)
(*b222coeflist=b222biaslister/@b222ctabtemp[[All,4]];*)
(*MapThread[Append,{b222exps,b222coeflist}]*)
(*];*)


(* ::Input:: *)
(*(*combine everything*)*)
(*b222ctabdec[k1_,k2_,k3_]:=compressor[Join[b222\[Mu]0to4ctabdec[k1,k2,k3],b222\[Mu]56ctabdec[k1,k2,k3]]];*)


(* ::Subsubsection::Closed:: *)
(*Finally, get list of exponents for b222*)


(* ::Input:: *)
(**)
(*ctabexample=b222ctabdec[0.11,0.12,0.13];*)
(*b222exps=ctabexample[[All,1;;3]]*)


(* ::InheritFromParent:: *)
(*b222exps={{-4,2,2},{-3,1,2},{-3,2,1},{-3,2,2},{-2,0,2},{-2,1,1},{-2,1,2},{-2,2,0},{-2,2,1},{-2,2,2},{-1,-1,2},{-1,0,1},{-1,0,2},{-1,1,0},{-1,1,1},{-1,1,2},{-1,2,-1},{-1,2,0},{-1,2,1},{-1,2,2},{0,-2,2},{0,-1,1},{0,-1,2},{0,0,0},{0,0,1},{0,0,2},{0,1,-1},{0,1,0},{0,1,1},{0,1,2},{0,2,-2},{0,2,-1},{0,2,0},{0,2,1},{0,2,2},{1,-3,2},{1,-2,1},{1,-2,2},{1,-1,0},{1,-1,1},{1,-1,2},{1,0,-1},{1,0,0},{1,0,1},{1,0,2},{1,1,-2},{1,1,-1},{1,1,0},{1,1,1},{1,1,2},{1,2,-3},{1,2,-2},{1,2,-1},{1,2,0},{1,2,1},{1,2,2},{2,-4,2},{2,-3,1},{2,-3,2},{2,-2,0},{2,-2,1},{2,-2,2},{2,-1,-1},{2,-1,0},{2,-1,1},{2,-1,2},{2,0,-2},{2,0,-1},{2,0,0},{2,0,1},{2,0,2},{2,1,-3},{2,1,-2},{2,1,-1},{2,1,0},{2,1,1},{2,1,2},{2,2,-4},{2,2,-3},{2,2,-2},{2,2,-1},{2,2,0},{2,2,1},{2,2,2}};*)


(* ::Subsubsection:: *)
(*Generate ctabs for python code and corresponding coefficients*)


(* ::Input:: *)
(*Export[ctabpath<>"B222ctab.csv",b222exps]*)
(*Export[LoopTablesPath<>"B222ctab.csv",b222exps]*)


Export[LoopTablesPath<>"B222simpcoefslist.m",b222biaslist]


(* ::Input:: *)
(*genb222coef[k1sub_,k2sub_,k3sub_]:=Module[{filepath},*)
(*filepath=saveLoopCoefPath[{k1sub,k2sub,k3sub},"B222"];*)
(*If[Not[FileExistsQ[filepath]],Export[filepath,b222ctabdec[k1sub,k2sub,k3sub][[All,4]]]]*)
(*];*)


(* ::Input:: *)
(*(**)
(*Monitor[Do[genb222coef@@fisherPoints\[LeftDoubleBracket]i\[RightDoubleBracket],*)
(*{i,1,Length[fisherPoints]}]//AbsoluteTiming,i];*)
(**)*)


(* ::Input:: *)
(*Monitor[Do[genb222coef@@CMASStriEff[[i]],*)
(*{i,1,Length[CMASStriEff]}]//AbsoluteTiming,i];*)


(* ::Input:: *)
(*Monitor[Do[genb222coef@@LOWZtriEff[[i]],*)
(*{i,1,Length[LOWZtriEff]}]//AbsoluteTiming,i];*)


Monitor[Do[genb222coef@@baseTriEff[[i]],
{i,1,Length[baseTriEff]}]//AbsoluteTiming,i];


(* ::Subsection:: *)
(*B321 I*)


(* ::Subsubsection:: *)
(*Create ctab*)


(* ::Text:: *)
(*IMPORTANT: We  will still do the substitution with power laws it as a convenient substitution since we only need to keep track of the argument*)


(* ::Input:: *)
(*k3211[k1_,k2_,q_]:=k3211[k1,k2,q]=K1r[k1] K2rsym[q,k2-q]K3rsym[-q,q-k2,-k1]/.tb;*)


(* ::Input:: *)
(*B3211Integrand[k1_,k2_,q_] :=B3211Integrand[k1,k2,q]= 6Plin[k1]k3211[k1,k2,q] Plin[q]Plin[kmq[k2]]/.tb*)


(* ::Input:: *)
(*ker3211temp= B3211Integrand[k1,k2,q]/.{Plin[k1]->1,Plin[q]->q^\[Nu]1,Plin[k2mq]-> k2mq^\[Nu]2}//Expand;*)


(* ::Text:: *)
(*For B321I, we have terms which involve k3pq.*)
(*We get rid of them by the following transformation.*)
(*First, we check if k3pq is in the denominator: these terms will have the shift q -> k2 - q, which also implies a modification of \[Mu]0, which we take into account as well.*)
(*Then, k3pq in the numerator can be simply expressed in terms of k1pq and k2mq.*)
(*It is important to not do the change of variable when k3pq is in the numerator because that generates a ctab with exponents that exceed the computing capabilities of the python code.*)


(* ::Input:: *)
(*qtok2mqreps={(q^a_. k2mq^b_. k3pq^c_ \[Mu]0^d_./;c<0):>k2mq^a q^b k1pq^c ((-q \[Mu]0+k2 \[Mu]2)/k2mq)^d,(q^a_. k2mq^b_. k3pq^c_/;c<0):>k2mq^a q^b k1pq^c};*)
(*(*all terms that need a change of variables have q, k2mq and k3pq because of the power spectrum exponent \[Nu]*)*)
(*k3pqreps={( k3pq^c_/;c>0):>(Sqrt[k1^2-k1pq^2-k2^2+k2mq^2+k3^2+q^2])^c};*)


(* ::Input:: *)
(*ker3211reps=ker3211temp/.Join[qtok2mqreps,k3pqreps];//AbsoluteTiming*)


(* ::Input:: *)
(*Variables[ker3211reps]*)


(* ::Input:: *)
(*ker3211exp=ker3211reps//Expand;*)


(* ::Text:: *)
(*Substitute the \[Mu]0^n transformations we computed before*)


(* ::Input:: *)
(*ker3211\[Mu]0list=CoefficientList[ker3211exp,\[Mu]0];*)


(* ::Input:: *)
(*\[Mu]0sol3211=<<(LoopTablesPath<>"mu0sol_nonsym_upto4.m");*)
(*\[Mu]0sol3211exp=\[Mu]0sol3211//Expand;*)


(* ::Input:: *)
(*(*We calculate the ctabs for \[Mu]0 and ker3211 and add the constant in the \[Mu]0 replacement by hand *)*)
(*\[Mu]0ctab3211temp=Prepend[getCtab/@\[Mu]0sol3211exp[[2;;]],{{0,0,0,1}}];*)
(*\[Mu]0ctab3211=compressor/@\[Mu]0ctab3211temp;*)


(* ::Input:: *)
(*ker3211\[Mu]0ctabtemp=getCtab/@ker3211\[Mu]0list;*)
(*ker3211\[Mu]0ctab=compressor/@ker3211\[Mu]0ctabtemp;*)


(* ::Input:: *)
(*(*Combine ctabs, join them, and compress the final table*)*)
(*b3211ctabjoin=MapThread[JoinCtabs[#1,#2]&,{ker3211\[Mu]0ctab,\[Mu]0ctab3211}];*)
(*b3211ctab=compressor[(Join@@b3211ctabjoin)//Expand];*)


(*Check*)
FreeQ[b3211ctab,k3pq]
FreeQ[b3211ctab,k2pq]
FreeQ[b3211ctab,k1mq]


(* ::Input:: *)
(*(*Create ctab setting \[Nu]'s to 0, because for this diagram they do not matter*)*)
(*b3211ctabsimp=Sort[compressor[b3211ctab/.{\[Nu]1->0,\[Nu]2->0}]];*)


(* ::Input:: *)
(*(*b3211ctabsimp>>(LoopTablesPath<>"b3211ctabsimp.mx");*)*)


Export[AuxPath<>"b3211ctabsimp.m",b3211ctabsimp]


(* ::Subsubsection:: *)
(*Decompose into a basis of bias coefficients and do permutations*)


(* ::Input:: *)
(*(*b3211ctabsimp=<<(LoopTablesPath<>"b3211ctabsimp.mx");*)*)


b3211ctabsimp=Import[AuxPath<>"b3211ctabsimp.m"];


(* ::Input:: *)
(*b3211vars={b1,b2,b3,b5,b6,b8,b10,f1,\[Mu]1,\[Mu]2};*)


(* ::Input:: *)
(*(*Get all bias monomials that appear in the b3211 kernel*)*)
(*(*b3211biaslist=getMonomials[b3211vars,b3211ctabsimp];*)*)
(*b3211biaslist={b1^3,b1^2 b10,b1^2 b2,b1 b10 b2,b1 b2^2,b1^2 b3,b1 b2 b3,b1^2 b5,b1 b10 b5,b1 b2 b5,b1 b3 b5,b1 b5^2,b1^2 b6,b1 b2 b6,b1 b5 b6,b1^2 b8,b1 b2 b8,b1 b5 b8,b1^2 f1 \[Mu]1^2,b1^3 f1 \[Mu]1^2,b1 b10 f1 \[Mu]1^2,b1 b2 f1 \[Mu]1^2,b1^2 b2 f1 \[Mu]1^2,b10 b2 f1 \[Mu]1^2,b2^2 f1 \[Mu]1^2,b1 b2^2 f1 \[Mu]1^2,b1 b3 f1 \[Mu]1^2,b2 b3 f1 \[Mu]1^2,b1 b5 f1 \[Mu]1^2,b1^2 b5 f1 \[Mu]1^2,b10 b5 f1 \[Mu]1^2,b2 b5 f1 \[Mu]1^2,b1 b2 b5 f1 \[Mu]1^2,b3 b5 f1 \[Mu]1^2,b5^2 f1 \[Mu]1^2,b1 b5^2 f1 \[Mu]1^2,b1 b6 f1 \[Mu]1^2,b2 b6 f1 \[Mu]1^2,b5 b6 f1 \[Mu]1^2,b1 b8 f1 \[Mu]1^2,b2 b8 f1 \[Mu]1^2,b5 b8 f1 \[Mu]1^2,b1^2 f1^2 \[Mu]1^2,b1^3 f1^2 \[Mu]1^2,b1 b2 f1^2 \[Mu]1^2,b1^2 b2 f1^2 \[Mu]1^2,b1 b5 f1^2 \[Mu]1^2,b1^2 b5 f1^2 \[Mu]1^2,b1 f1^2 \[Mu]1^4,b1^2 f1^2 \[Mu]1^4,b1^3 f1^2 \[Mu]1^4,b2 f1^2 \[Mu]1^4,b1 b2 f1^2 \[Mu]1^4,b1^2 b2 f1^2 \[Mu]1^4,b2^2 f1^2 \[Mu]1^4,b5 f1^2 \[Mu]1^4,b1 b5 f1^2 \[Mu]1^4,b1^2 b5 f1^2 \[Mu]1^4,b2 b5 f1^2 \[Mu]1^4,b5^2 f1^2 \[Mu]1^4,b1 f1^3 \[Mu]1^4,b1^2 f1^3 \[Mu]1^4,b2 f1^3 \[Mu]1^4,b1 b2 f1^3 \[Mu]1^4,b5 f1^3 \[Mu]1^4,b1 b5 f1^3 \[Mu]1^4,b1 f1^3 \[Mu]1^6,b1^2 f1^3 \[Mu]1^6,b2 f1^3 \[Mu]1^6,b1 b2 f1^3 \[Mu]1^6,b5 f1^3 \[Mu]1^6,b1 b5 f1^3 \[Mu]1^6,b1 f1^4 \[Mu]1^6,b2 f1^4 \[Mu]1^6,b5 f1^4 \[Mu]1^6,b1 f1^4 \[Mu]1^8,b2 f1^4 \[Mu]1^8,b5 f1^4 \[Mu]1^8,b1^2 f1 \[Mu]1 \[Mu]2,b1^3 f1 \[Mu]1 \[Mu]2,b1^2 b10 f1 \[Mu]1 \[Mu]2,b1 b2 f1 \[Mu]1 \[Mu]2,b1^2 b2 f1 \[Mu]1 \[Mu]2,b1 b2^2 f1 \[Mu]1 \[Mu]2,b1^2 b3 f1 \[Mu]1 \[Mu]2,b1 b5 f1 \[Mu]1 \[Mu]2,b1^2 b5 f1 \[Mu]1 \[Mu]2,b1 b2 b5 f1 \[Mu]1 \[Mu]2,b1 b5^2 f1 \[Mu]1 \[Mu]2,b1^2 b6 f1 \[Mu]1 \[Mu]2,b1^2 b8 f1 \[Mu]1 \[Mu]2,b1^2 f1^2 \[Mu]1 \[Mu]2,b1^3 f1^2 \[Mu]1 \[Mu]2,b1 b2 f1^2 \[Mu]1 \[Mu]2,b1^2 b2 f1^2 \[Mu]1 \[Mu]2,b1 b5 f1^2 \[Mu]1 \[Mu]2,b1^2 b5 f1^2 \[Mu]1 \[Mu]2,b1 f1^2 \[Mu]1^3 \[Mu]2,b1^2 f1^2 \[Mu]1^3 \[Mu]2,b1^3 f1^2 \[Mu]1^3 \[Mu]2,b1 b10 f1^2 \[Mu]1^3 \[Mu]2,b2 f1^2 \[Mu]1^3 \[Mu]2,b1 b2 f1^2 \[Mu]1^3 \[Mu]2,b1^2 b2 f1^2 \[Mu]1^3 \[Mu]2,b2^2 f1^2 \[Mu]1^3 \[Mu]2,b1 b3 f1^2 \[Mu]1^3 \[Mu]2,b5 f1^2 \[Mu]1^3 \[Mu]2,b1 b5 f1^2 \[Mu]1^3 \[Mu]2,b1^2 b5 f1^2 \[Mu]1^3 \[Mu]2,b2 b5 f1^2 \[Mu]1^3 \[Mu]2,b5^2 f1^2 \[Mu]1^3 \[Mu]2,b1 b6 f1^2 \[Mu]1^3 \[Mu]2,b1 b8 f1^2 \[Mu]1^3 \[Mu]2,b1 f1^3 \[Mu]1^3 \[Mu]2,b1^2 f1^3 \[Mu]1^3 \[Mu]2,b1^3 f1^3 \[Mu]1^3 \[Mu]2,b2 f1^3 \[Mu]1^3 \[Mu]2,b1 b2 f1^3 \[Mu]1^3 \[Mu]2,b5 f1^3 \[Mu]1^3 \[Mu]2,b1 b5 f1^3 \[Mu]1^3 \[Mu]2,b1 f1^3 \[Mu]1^5 \[Mu]2,b1^2 f1^3 \[Mu]1^5 \[Mu]2,b1^3 f1^3 \[Mu]1^5 \[Mu]2,b2 f1^3 \[Mu]1^5 \[Mu]2,b1 b2 f1^3 \[Mu]1^5 \[Mu]2,b5 f1^3 \[Mu]1^5 \[Mu]2,b1 b5 f1^3 \[Mu]1^5 \[Mu]2,b1 f1^4 \[Mu]1^5 \[Mu]2,b1^2 f1^4 \[Mu]1^5 \[Mu]2,b2 f1^4 \[Mu]1^5 \[Mu]2,b5 f1^4 \[Mu]1^5 \[Mu]2,b1 f1^4 \[Mu]1^7 \[Mu]2,b1^2 f1^4 \[Mu]1^7 \[Mu]2,b2 f1^4 \[Mu]1^7 \[Mu]2,b5 f1^4 \[Mu]1^7 \[Mu]2,b1 f1^5 \[Mu]1^7 \[Mu]2,b1 f1^5 \[Mu]1^9 \[Mu]2,b1^2 f1 \[Mu]2^2,b1^3 f1 \[Mu]2^2,b1 b10 f1 \[Mu]2^2,b1^2 b10 f1 \[Mu]2^2,b1 b2 f1 \[Mu]2^2,b1^2 b2 f1 \[Mu]2^2,b1 b2^2 f1 \[Mu]2^2,b1 b3 f1 \[Mu]2^2,b1^2 b3 f1 \[Mu]2^2,b1 b5 f1 \[Mu]2^2,b1^2 b5 f1 \[Mu]2^2,b1 b2 b5 f1 \[Mu]2^2,b1 b5^2 f1 \[Mu]2^2,b1 b6 f1 \[Mu]2^2,b1^2 b6 f1 \[Mu]2^2,b1 b8 f1 \[Mu]2^2,b1^2 b8 f1 \[Mu]2^2,b1^2 f1^2 \[Mu]2^2,b1^3 f1^2 \[Mu]2^2,b1 b10 f1^2 \[Mu]2^2,b1 b2 f1^2 \[Mu]2^2,b1^2 b2 f1^2 \[Mu]2^2,b1 b3 f1^2 \[Mu]2^2,b1 b5 f1^2 \[Mu]2^2,b1^2 b5 f1^2 \[Mu]2^2,b1 b6 f1^2 \[Mu]2^2,b1 b8 f1^2 \[Mu]2^2,b1 f1^2 \[Mu]1^2 \[Mu]2^2,b1^2 f1^2 \[Mu]1^2 \[Mu]2^2,b1^3 f1^2 \[Mu]1^2 \[Mu]2^2,b10 f1^2 \[Mu]1^2 \[Mu]2^2,b1 b10 f1^2 \[Mu]1^2 \[Mu]2^2,b2 f1^2 \[Mu]1^2 \[Mu]2^2,b1 b2 f1^2 \[Mu]1^2 \[Mu]2^2,b1^2 b2 f1^2 \[Mu]1^2 \[Mu]2^2,b2^2 f1^2 \[Mu]1^2 \[Mu]2^2,b3 f1^2 \[Mu]1^2 \[Mu]2^2,b1 b3 f1^2 \[Mu]1^2 \[Mu]2^2,b5 f1^2 \[Mu]1^2 \[Mu]2^2,b1 b5 f1^2 \[Mu]1^2 \[Mu]2^2,b1^2 b5 f1^2 \[Mu]1^2 \[Mu]2^2,b2 b5 f1^2 \[Mu]1^2 \[Mu]2^2,b5^2 f1^2 \[Mu]1^2 \[Mu]2^2,b6 f1^2 \[Mu]1^2 \[Mu]2^2,b1 b6 f1^2 \[Mu]1^2 \[Mu]2^2,b8 f1^2 \[Mu]1^2 \[Mu]2^2,b1 b8 f1^2 \[Mu]1^2 \[Mu]2^2,b1 f1^3 \[Mu]1^2 \[Mu]2^2,b1^2 f1^3 \[Mu]1^2 \[Mu]2^2,b1^3 f1^3 \[Mu]1^2 \[Mu]2^2,b10 f1^3 \[Mu]1^2 \[Mu]2^2,b2 f1^3 \[Mu]1^2 \[Mu]2^2,b1 b2 f1^3 \[Mu]1^2 \[Mu]2^2,b3 f1^3 \[Mu]1^2 \[Mu]2^2,b5 f1^3 \[Mu]1^2 \[Mu]2^2,b1 b5 f1^3 \[Mu]1^2 \[Mu]2^2,b6 f1^3 \[Mu]1^2 \[Mu]2^2,b8 f1^3 \[Mu]1^2 \[Mu]2^2,b1 f1^4 \[Mu]1^2 \[Mu]2^2,b1^2 f1^4 \[Mu]1^2 \[Mu]2^2,f1^3 \[Mu]1^4 \[Mu]2^2,b1 f1^3 \[Mu]1^4 \[Mu]2^2,b1^2 f1^3 \[Mu]1^4 \[Mu]2^2,b1^3 f1^3 \[Mu]1^4 \[Mu]2^2,b10 f1^3 \[Mu]1^4 \[Mu]2^2,b2 f1^3 \[Mu]1^4 \[Mu]2^2,b1 b2 f1^3 \[Mu]1^4 \[Mu]2^2,b3 f1^3 \[Mu]1^4 \[Mu]2^2,b5 f1^3 \[Mu]1^4 \[Mu]2^2,b1 b5 f1^3 \[Mu]1^4 \[Mu]2^2,b6 f1^3 \[Mu]1^4 \[Mu]2^2,b8 f1^3 \[Mu]1^4 \[Mu]2^2,f1^4 \[Mu]1^4 \[Mu]2^2,b1 f1^4 \[Mu]1^4 \[Mu]2^2,b1^2 f1^4 \[Mu]1^4 \[Mu]2^2,b2 f1^4 \[Mu]1^4 \[Mu]2^2,b5 f1^4 \[Mu]1^4 \[Mu]2^2,f1^5 \[Mu]1^4 \[Mu]2^2,b1 f1^5 \[Mu]1^4 \[Mu]2^2,f1^4 \[Mu]1^6 \[Mu]2^2,b1 f1^4 \[Mu]1^6 \[Mu]2^2,b1^2 f1^4 \[Mu]1^6 \[Mu]2^2,b2 f1^4 \[Mu]1^6 \[Mu]2^2,b5 f1^4 \[Mu]1^6 \[Mu]2^2,f1^5 \[Mu]1^6 \[Mu]2^2,b1 f1^5 \[Mu]1^6 \[Mu]2^2,f1^6 \[Mu]1^6 \[Mu]2^2,f1^5 \[Mu]1^8 \[Mu]2^2,b1 f1^5 \[Mu]1^8 \[Mu]2^2,f1^6 \[Mu]1^8 \[Mu]2^2,f1^6 \[Mu]1^10 \[Mu]2^2,b1 f1^2 \[Mu]1 \[Mu]2^3,b1^2 f1^2 \[Mu]1 \[Mu]2^3,b1^3 f1^2 \[Mu]1 \[Mu]2^3,b1 b10 f1^2 \[Mu]1 \[Mu]2^3,b1 b2 f1^2 \[Mu]1 \[Mu]2^3,b1^2 b2 f1^2 \[Mu]1 \[Mu]2^3,b1 b3 f1^2 \[Mu]1 \[Mu]2^3,b1 b5 f1^2 \[Mu]1 \[Mu]2^3,b1^2 b5 f1^2 \[Mu]1 \[Mu]2^3,b1 b6 f1^2 \[Mu]1 \[Mu]2^3,b1 b8 f1^2 \[Mu]1 \[Mu]2^3,b1 f1^3 \[Mu]1 \[Mu]2^3,b1^2 f1^3 \[Mu]1 \[Mu]2^3,b1^3 f1^3 \[Mu]1 \[Mu]2^3,b1 b2 f1^3 \[Mu]1 \[Mu]2^3,b1 b5 f1^3 \[Mu]1 \[Mu]2^3,b1 f1^4 \[Mu]1 \[Mu]2^3,b1^2 f1^4 \[Mu]1 \[Mu]2^3,f1^3 \[Mu]1^3 \[Mu]2^3,b1 f1^3 \[Mu]1^3 \[Mu]2^3,b1^2 f1^3 \[Mu]1^3 \[Mu]2^3,b1^3 f1^3 \[Mu]1^3 \[Mu]2^3,b10 f1^3 \[Mu]1^3 \[Mu]2^3,b2 f1^3 \[Mu]1^3 \[Mu]2^3,b1 b2 f1^3 \[Mu]1^3 \[Mu]2^3,b3 f1^3 \[Mu]1^3 \[Mu]2^3,b5 f1^3 \[Mu]1^3 \[Mu]2^3,b1 b5 f1^3 \[Mu]1^3 \[Mu]2^3,b6 f1^3 \[Mu]1^3 \[Mu]2^3,b8 f1^3 \[Mu]1^3 \[Mu]2^3,f1^4 \[Mu]1^3 \[Mu]2^3,b1 f1^4 \[Mu]1^3 \[Mu]2^3,b1^2 f1^4 \[Mu]1^3 \[Mu]2^3,b2 f1^4 \[Mu]1^3 \[Mu]2^3,b5 f1^4 \[Mu]1^3 \[Mu]2^3,f1^5 \[Mu]1^3 \[Mu]2^3,b1 f1^5 \[Mu]1^3 \[Mu]2^3,f1^4 \[Mu]1^5 \[Mu]2^3,b1 f1^4 \[Mu]1^5 \[Mu]2^3,b1^2 f1^4 \[Mu]1^5 \[Mu]2^3,b2 f1^4 \[Mu]1^5 \[Mu]2^3,b5 f1^4 \[Mu]1^5 \[Mu]2^3,f1^5 \[Mu]1^5 \[Mu]2^3,b1 f1^5 \[Mu]1^5 \[Mu]2^3,f1^6 \[Mu]1^5 \[Mu]2^3,f1^5 \[Mu]1^7 \[Mu]2^3,b1 f1^5 \[Mu]1^7 \[Mu]2^3,f1^6 \[Mu]1^7 \[Mu]2^3,f1^6 \[Mu]1^9 \[Mu]2^3,b1 f1^2 \[Mu]2^4,b1^2 f1^2 \[Mu]2^4,b1^3 f1^2 \[Mu]2^4,b1 b10 f1^2 \[Mu]2^4,b1 b2 f1^2 \[Mu]2^4,b1^2 b2 f1^2 \[Mu]2^4,b1 b3 f1^2 \[Mu]2^4,b1 b5 f1^2 \[Mu]2^4,b1^2 b5 f1^2 \[Mu]2^4,b1 b6 f1^2 \[Mu]2^4,b1 b8 f1^2 \[Mu]2^4,b1 f1^3 \[Mu]2^4,b1^2 f1^3 \[Mu]2^4,b1^3 f1^3 \[Mu]2^4,b1 b2 f1^3 \[Mu]2^4,b1 b5 f1^3 \[Mu]2^4,b1 f1^4 \[Mu]2^4,b1^2 f1^4 \[Mu]2^4,f1^3 \[Mu]1^2 \[Mu]2^4,b1 f1^3 \[Mu]1^2 \[Mu]2^4,b1^2 f1^3 \[Mu]1^2 \[Mu]2^4,b1^3 f1^3 \[Mu]1^2 \[Mu]2^4,b10 f1^3 \[Mu]1^2 \[Mu]2^4,b2 f1^3 \[Mu]1^2 \[Mu]2^4,b1 b2 f1^3 \[Mu]1^2 \[Mu]2^4,b3 f1^3 \[Mu]1^2 \[Mu]2^4,b5 f1^3 \[Mu]1^2 \[Mu]2^4,b1 b5 f1^3 \[Mu]1^2 \[Mu]2^4,b6 f1^3 \[Mu]1^2 \[Mu]2^4,b8 f1^3 \[Mu]1^2 \[Mu]2^4,f1^4 \[Mu]1^2 \[Mu]2^4,b1 f1^4 \[Mu]1^2 \[Mu]2^4,b1^2 f1^4 \[Mu]1^2 \[Mu]2^4,b2 f1^4 \[Mu]1^2 \[Mu]2^4,b5 f1^4 \[Mu]1^2 \[Mu]2^4,f1^5 \[Mu]1^2 \[Mu]2^4,b1 f1^5 \[Mu]1^2 \[Mu]2^4,f1^4 \[Mu]1^4 \[Mu]2^4,b1 f1^4 \[Mu]1^4 \[Mu]2^4,b1^2 f1^4 \[Mu]1^4 \[Mu]2^4,b2 f1^4 \[Mu]1^4 \[Mu]2^4,b5 f1^4 \[Mu]1^4 \[Mu]2^4,f1^5 \[Mu]1^4 \[Mu]2^4,b1 f1^5 \[Mu]1^4 \[Mu]2^4,f1^6 \[Mu]1^4 \[Mu]2^4,f1^5 \[Mu]1^6 \[Mu]2^4,b1 f1^5 \[Mu]1^6 \[Mu]2^4,f1^6 \[Mu]1^6 \[Mu]2^4,f1^6 \[Mu]1^8 \[Mu]2^4,b1 f1^3 \[Mu]1 \[Mu]2^5,b1^2 f1^3 \[Mu]1 \[Mu]2^5,b1^3 f1^3 \[Mu]1 \[Mu]2^5,b1 b2 f1^3 \[Mu]1 \[Mu]2^5,b1 b5 f1^3 \[Mu]1 \[Mu]2^5,b1 f1^4 \[Mu]1 \[Mu]2^5,b1^2 f1^4 \[Mu]1 \[Mu]2^5,b1 f1^5 \[Mu]1 \[Mu]2^5,f1^4 \[Mu]1^3 \[Mu]2^5,b1 f1^4 \[Mu]1^3 \[Mu]2^5,b1^2 f1^4 \[Mu]1^3 \[Mu]2^5,b2 f1^4 \[Mu]1^3 \[Mu]2^5,b5 f1^4 \[Mu]1^3 \[Mu]2^5,f1^5 \[Mu]1^3 \[Mu]2^5,b1 f1^5 \[Mu]1^3 \[Mu]2^5,f1^6 \[Mu]1^3 \[Mu]2^5,f1^5 \[Mu]1^5 \[Mu]2^5,b1 f1^5 \[Mu]1^5 \[Mu]2^5,f1^6 \[Mu]1^5 \[Mu]2^5,f1^6 \[Mu]1^7 \[Mu]2^5,b1 f1^3 \[Mu]2^6,b1^2 f1^3 \[Mu]2^6,b1^3 f1^3 \[Mu]2^6,b1 b2 f1^3 \[Mu]2^6,b1 b5 f1^3 \[Mu]2^6,b1 f1^4 \[Mu]2^6,b1^2 f1^4 \[Mu]2^6,f1^4 \[Mu]1^2 \[Mu]2^6,b1 f1^4 \[Mu]1^2 \[Mu]2^6,b1^2 f1^4 \[Mu]1^2 \[Mu]2^6,b2 f1^4 \[Mu]1^2 \[Mu]2^6,b5 f1^4 \[Mu]1^2 \[Mu]2^6,f1^5 \[Mu]1^2 \[Mu]2^6,b1 f1^5 \[Mu]1^2 \[Mu]2^6,f1^5 \[Mu]1^4 \[Mu]2^6,b1 f1^5 \[Mu]1^4 \[Mu]2^6,f1^6 \[Mu]1^4 \[Mu]2^6,f1^6 \[Mu]1^6 \[Mu]2^6,b1 f1^4 \[Mu]1 \[Mu]2^7,b1^2 f1^4 \[Mu]1 \[Mu]2^7,b1 f1^5 \[Mu]1 \[Mu]2^7,f1^5 \[Mu]1^3 \[Mu]2^7,b1 f1^5 \[Mu]1^3 \[Mu]2^7,f1^6 \[Mu]1^3 \[Mu]2^7,f1^6 \[Mu]1^5 \[Mu]2^7,b1 f1^4 \[Mu]2^8,b1^2 f1^4 \[Mu]2^8,f1^5 \[Mu]1^2 \[Mu]2^8,b1 f1^5 \[Mu]1^2 \[Mu]2^8,f1^6 \[Mu]1^4 \[Mu]2^8,b1 f1^5 \[Mu]1 \[Mu]2^9,f1^6 \[Mu]1^3 \[Mu]2^9};*)
(*(*b3211biaslist>>"tables_for_B1loop/B3211biaslist.mx";*)*)


(* ::Input:: *)
(*(*Function that obtains the coefficient of each term in biaslist in each ctab coefficient "coef" *)*)
(*b3211biaslister=biaslisterfn[b3211vars,b3211biaslist];*)


(* ::Input:: *)
(*b3211exps=b3211ctabsimp[[All,1;;3]];*)
(*b3211coefs=b3211ctabsimp[[All,4]];*)


(* ::Input:: *)
(*(*apply decomposition in biases to all ctab coefficients and simplify - this takes 7 minutes*)*)
(*b3211biascoef=Simplify[b3211biaslister/@b3211coefs];*)
(*b3211ctabdec=MapThread[Append,{b3211exps,b3211biascoef}];*)


(* ::Input:: *)
(*(**)
(*b3211ctabdec>>(LoopTablesPath<>"b3211ctabdec.mx");*)
(**)*)


(* ::Input:: *)
(*(*b3211ctabdec=<<"tables_for_B1loop/b3211ctabdec.mx";*)*)


(* ::Text:: *)
(*Now we just have to permute the coefficients over the 6 possible permutations of k1, k2 and k3*)


(* ::Input:: *)
(*(*The following code does the permutation*)*)
(*arr  = {k1,\[Mu]1,k2,\[Mu]2,k3,\[Mu]3};*)
(*arrt={k1t,\[Mu]1t,k2t,\[Mu]2t,k3t,\[Mu]3t};*)
(*pt = Permutations[{{k1t,\[Mu]1t},{k2t,\[Mu]2t},{k3t,\[Mu]3t}}]/.{{a_,b_},{c_,d_},{e_,f_}}-> {a,b,c,d,e,f};*)
(*permreps=Thread[(arr->#)]&/@pt;*)


(* ::Input:: *)
(*(*get list of permuted biases*)*)
(*b3211biasperm=b3211biaslist/.permreps/.Thread[arrt->arr]//Transpose;*)


(* ::Input:: *)
(*b3211coefsdec=b3211ctabdec[[All,4]];*)


(* ::Input:: *)
(*(*Get list of permuted coefficients - 13s*)*)
(*b3211permdectemp=b3211coefsdec/.permreps/.Thread[arrt->arr];//AbsoluteTiming*)


(* ::Input:: *)
(*(*transpose to get a more convenient tensor*)*)
(*b3211permdec=Transpose[b3211permdectemp,{3,1,2}];*)


(* ::Input:: *)
(*b3211biasperm>>(LoopTablesPath<>"B3211biasPerm.mx");*)


(* ::Input:: *)
(*b3211permdec>>(LoopTablesPath<>"B3211coefsPerm.mx");*)


(* ::Subsubsection::Closed:: *)
(*Exponent list*)


(* ::Input:: *)
(*(*b3211exps=b3211ctabsimp\[LeftDoubleBracket]All,1;;3\[RightDoubleBracket];*)*)
(*b3211exps={{-4,1,2},{-3,0,2},{-3,1,1},{-3,1,2},{-2,-1,2},{-2,0,1},{-2,0,2},{-2,1,0},{-2,1,1},{-2,1,2},{-1,-2,2},{-1,-1,1},{-1,-1,2},{-1,0,0},{-1,0,1},{-1,0,2},{-1,1,-1},{-1,1,0},{-1,1,1},{-1,1,2},{0,-3,2},{0,-2,1},{0,-2,2},{0,-1,0},{0,-1,1},{0,-1,2},{0,0,-1},{0,0,0},{0,0,1},{0,0,2},{0,1,-2},{0,1,-1},{0,1,0},{0,1,1},{0,1,2},{1,-4,2},{1,-3,1},{1,-3,2},{1,-2,0},{1,-2,1},{1,-2,2},{1,-1,-1},{1,-1,0},{1,-1,1},{1,-1,2},{1,0,-2},{1,0,-1},{1,0,0},{1,0,1},{1,0,2},{1,1,-2},{1,1,-1},{1,1,0},{1,1,1},{1,1,2},{2,-4,1},{2,-4,2},{2,-3,0},{2,-3,1},{2,-3,2},{2,-2,-1},{2,-2,0},{2,-2,1},{2,-2,2},{2,-1,-2},{2,-1,-1},{2,-1,0},{2,-1,1},{2,-1,2},{2,0,-3},{2,0,-2},{2,0,-1},{2,0,0},{2,0,1},{2,0,2},{2,1,-2},{2,1,-1},{2,1,0},{2,1,1},{2,1,2}};*)


(* ::Subsubsection::Closed:: *)
(*Generate ctab coefs for fisherPoints*)


(* ::Input:: *)
(*b3211permdec=<<(LoopTablesPath<>"B3211coefsPerm.mx");*)
(*b3211permdecfn[k1_,k2_,k3_]=b3211permdec;*)


(* ::Input:: *)
(*(*Export ctab*)*)
(*(*Export[ctabpath<>"B3211ctab.csv",b3211exps];*)*)


(* ::Input:: *)
(*genb3211coef[k1sub_,k2sub_,k3sub_]:=Module[{temp},*)
(*temp=b3211permdecfn[k1sub,k2sub,k3sub];*)
(*Export[LoopTablesPath<>"B3211coefs/B3211coefs_"<>ToString[k1sub//N]<>"_"<>ToString[k2sub//N]<>"_"<>ToString[k3sub//N]<>"_.mx",temp]*)
(*];*)


(* ::Input:: *)
(*Monitor[Do[genb3211coef@@CMASStriEff[[j]],*)
(*{j,1,Length[CMASStriEff]}],j]*)


(* ::Input:: *)
(*Monitor[Do[genb3211coef@@LOWZtriEff[[j]],*)
(*{j,1,Length[LOWZtriEff]}],j]*)


(* ::Input:: *)
(*Monitor[Do[genb3211coef@@fisherPoints[[j]],*)
(*{j,1,Length[fisherPoints]}],j]*)


(* ::Subsection:: *)
(*B321 II*)


(* ::Subsubsection::Closed:: *)
(*Create ctab*)


(* ::Text:: *)
(*IMPORTANT: In the FFTLog, P(q) is q^\[Nu] -- We  will still do it as a convenient substitution since we only need to keep track of the argument.*)
(*Now it's very important in this diagram to subtract the K3rsym up to zeroth order, to match the finite part of the counterterms.*)


(* ::Input:: *)
(*(*serK3k2=Normal[Series[K3rsym[k1,q,-q]/.magrevreps,{q,\[Infinity],2}]]/.tb//Simplify*)*)


(* ::Input:: *)
(*serK3k0=Normal[Series[K3rsym[k1,q,-q]/.magrevreps,{q,\[Infinity],0}]]/.tb//Simplify;*)


(* ::Text:: *)
(*Integrate over angles dividing by 1/4\[Pi]*)


(* ::Input:: *)
(*uvK3k0=1/2 Integrate[serK3k0,{x,-1,1}]//Simplify;*)


(* ::Input:: *)
(*k3212[k1_,k2_,q_]:=k3212[k1,k2,q]=K1r[k2]K2rsym[k1,k2](K3rsym[k1,q,-q]-uvK3k0)/.mag[0]-> 0*)


(* ::Input:: *)
(*B3212Integrand[k1_,k2_,q_] :=B3212Integrand[k1,k2,q] = 6Plin[k1]Plin[k2]Plin[q]k3212[k1,k2,q]/.tb*)


(* ::Input:: *)
(*ker3212temp = (B3212Integrand[k1,k2,q]/.{Plin[k1]->1,Plin[k2]->1,Plin[q]-> q^\[Nu]}//Expand);*)


(* ::Text:: *)
(*In the next loop we do 2 transformations.*)
(*If there is k1mq in the denominator, we switch the sign of q, that is described by qrevreps.*)
(*Then we only have k1pq in the denominator.*)
(*Finally, we can just express k1mq in the numerator as a function of k1pq, k1 and q*)


(* ::Input:: *)
(*qrevreps={(q^a_. k1mq^b_. \[Mu]0^d_./;b<0):>q^a k1pq^b (-\[Mu]0)^d,(q^a_. k1mq^b_./;b<0):>q^a k1pq^b};*)
(*k1mqreps={( k1mq^c_./;c>0):>(Sqrt[-k1pq^2+2k1^2+2q^2])^c};*)


(* ::Input:: *)
(*ker3212reps=ker3212temp/.Join[qrevreps,k1mqreps];//AbsoluteTiming*)


(* ::Input:: *)
(*ker3212exp=ker3212reps//Expand;*)


(* ::Input:: *)
(*ker3212\[Mu]0list=CoefficientList[ker3212exp,\[Mu]0];*)


(* ::Text:: *)
(*Substitute the \[Mu]0^n transformations we computed before and obtain ctab*)


(* ::Input:: *)
(*\[Mu]0sol3212=<<(LoopTablesPath<>"mu0solk1_vec_upto2.m");*)


(* ::Input:: *)
(*ker3212=ker3212\[Mu]0list.\[Mu]0sol3212//Expand;*)


(* ::Input:: *)
(*(*get compressed ctab*)*)
(*b3212ctab=compressor[getCtab[ker3212]];*)


(* ::Input:: *)
(*b3212ctabsimp=Sort[compressor[b3212ctab/.{\[Nu]->0}//Expand]];*)


(* ::Input:: *)
(*b3212ctabsimp>>(LoopTablesPath<>"b3212ctabsimp.mx");*)


(* ::Subsubsection::Closed:: *)
(*Decompose into a basis of bias coefficients*)


(* ::Input:: *)
(*(*b3212ctabsimp=<<(LoopTablesPath<>"b3212ctabsimp.mx");*)*)


(* ::Input:: *)
(*b3212exps=b3212ctabsimp[[All,1;;3]];*)
(*b3212coefs=b3212ctabsimp[[All,4]];*)


(* ::Input:: *)
(*b3212vars={b1,b2,b3,b5,b6,b8,f1,\[Mu]1,\[Mu]2};*)


(* ::Input:: *)
(*(*Get all bias monomials that appear in the b3212 kernel by checking each term in the ctab*)*)
(*(*b3212biaslist=getMonomials[b3212vars,b3212ctabsimp];*)*)
(*b3212biaslist={b1^3,b1^2 b2,b1 b2^2,b1^2 b3,b1 b2 b3,b1^2 b5,b1 b2 b5,b1 b3 b5,b1^2 b6,b1 b2 b6,b1 b5 b6,b1^2 b8,b1 b2 b8,b1 b5 b8,b1^2 f1 \[Mu]1^2,b1^3 f1 \[Mu]1^2,b1 b2 f1 \[Mu]1^2,b1^2 b2 f1 \[Mu]1^2,b1 b3 f1 \[Mu]1^2,b1^2 b3 f1 \[Mu]1^2,b1 b5 f1 \[Mu]1^2,b1^2 b5 f1 \[Mu]1^2,b1 b6 f1 \[Mu]1^2,b1^2 b6 f1 \[Mu]1^2,b1 b8 f1 \[Mu]1^2,b1^2 b8 f1 \[Mu]1^2,b1^2 f1^2 \[Mu]1^2,b1^3 f1^2 \[Mu]1^2,b1 b2 f1^2 \[Mu]1^2,b1^2 b2 f1^2 \[Mu]1^2,b1 b5 f1^2 \[Mu]1^2,b1^2 b5 f1^2 \[Mu]1^2,b1 f1^2 \[Mu]1^4,b1^2 f1^2 \[Mu]1^4,b1^3 f1^2 \[Mu]1^4,b1 b2 f1^2 \[Mu]1^4,b1^2 b2 f1^2 \[Mu]1^4,b1 b5 f1^2 \[Mu]1^4,b1^2 b5 f1^2 \[Mu]1^4,b1 f1^3 \[Mu]1^4,b1^2 f1^3 \[Mu]1^4,b1^3 f1^3 \[Mu]1^4,b1 b2 f1^3 \[Mu]1^4,b1 b5 f1^3 \[Mu]1^4,b1 f1^3 \[Mu]1^6,b1^2 f1^3 \[Mu]1^6,b1^3 f1^3 \[Mu]1^6,b1 b2 f1^3 \[Mu]1^6,b1 b5 f1^3 \[Mu]1^6,b1 f1^4 \[Mu]1^6,b1^2 f1^4 \[Mu]1^6,b1 f1^4 \[Mu]1^8,b1^2 f1^4 \[Mu]1^8,b1^2 f1 \[Mu]1 \[Mu]2,b1^3 f1 \[Mu]1 \[Mu]2,b1 b2 f1 \[Mu]1 \[Mu]2,b1^2 b2 f1 \[Mu]1 \[Mu]2,b1 b3 f1 \[Mu]1 \[Mu]2,b1^2 b3 f1 \[Mu]1 \[Mu]2,b1 b6 f1 \[Mu]1 \[Mu]2,b1^2 b6 f1 \[Mu]1 \[Mu]2,b1 b8 f1 \[Mu]1 \[Mu]2,b1^2 b8 f1 \[Mu]1 \[Mu]2,b1 f1^2 \[Mu]1^3 \[Mu]2,b1^2 f1^2 \[Mu]1^3 \[Mu]2,b1^3 f1^2 \[Mu]1^3 \[Mu]2,b1 b2 f1^2 \[Mu]1^3 \[Mu]2,b1 b3 f1^2 \[Mu]1^3 \[Mu]2,b1 b6 f1^2 \[Mu]1^3 \[Mu]2,b1 b8 f1^2 \[Mu]1^3 \[Mu]2,b1 f1^3 \[Mu]1^3 \[Mu]2,b1^2 f1^3 \[Mu]1^3 \[Mu]2,b1^3 f1^3 \[Mu]1^3 \[Mu]2,b1 f1^3 \[Mu]1^5 \[Mu]2,b1^2 f1^3 \[Mu]1^5 \[Mu]2,b1^3 f1^3 \[Mu]1^5 \[Mu]2,b1 f1^4 \[Mu]1^5 \[Mu]2,b1^2 f1^4 \[Mu]1^5 \[Mu]2,b1 f1^4 \[Mu]1^7 \[Mu]2,b1^2 f1^4 \[Mu]1^7 \[Mu]2,b1 f1^5 \[Mu]1^7 \[Mu]2,b1 f1^5 \[Mu]1^9 \[Mu]2,b1^2 f1 \[Mu]2^2,b1^3 f1 \[Mu]2^2,b1 b2 f1 \[Mu]2^2,b1^2 b2 f1 \[Mu]2^2,b2^2 f1 \[Mu]2^2,b1 b3 f1 \[Mu]2^2,b1^2 b3 f1 \[Mu]2^2,b2 b3 f1 \[Mu]2^2,b1 b5 f1 \[Mu]2^2,b2 b5 f1 \[Mu]2^2,b3 b5 f1 \[Mu]2^2,b1 b6 f1 \[Mu]2^2,b1^2 b6 f1 \[Mu]2^2,b2 b6 f1 \[Mu]2^2,b5 b6 f1 \[Mu]2^2,b1 b8 f1 \[Mu]2^2,b1^2 b8 f1 \[Mu]2^2,b2 b8 f1 \[Mu]2^2,b5 b8 f1 \[Mu]2^2,b1 f1^2 \[Mu]1^2 \[Mu]2^2,b1^2 f1^2 \[Mu]1^2 \[Mu]2^2,b1^3 f1^2 \[Mu]1^2 \[Mu]2^2,b2 f1^2 \[Mu]1^2 \[Mu]2^2,b1 b2 f1^2 \[Mu]1^2 \[Mu]2^2,b3 f1^2 \[Mu]1^2 \[Mu]2^2,b1 b3 f1^2 \[Mu]1^2 \[Mu]2^2,b5 f1^2 \[Mu]1^2 \[Mu]2^2,b1 b5 f1^2 \[Mu]1^2 \[Mu]2^2,b6 f1^2 \[Mu]1^2 \[Mu]2^2,b1 b6 f1^2 \[Mu]1^2 \[Mu]2^2,b8 f1^2 \[Mu]1^2 \[Mu]2^2,b1 b8 f1^2 \[Mu]1^2 \[Mu]2^2,b1 f1^3 \[Mu]1^2 \[Mu]2^2,b1^2 f1^3 \[Mu]1^2 \[Mu]2^2,b1^3 f1^3 \[Mu]1^2 \[Mu]2^2,b2 f1^3 \[Mu]1^2 \[Mu]2^2,b1 b2 f1^3 \[Mu]1^2 \[Mu]2^2,b5 f1^3 \[Mu]1^2 \[Mu]2^2,b1 b5 f1^3 \[Mu]1^2 \[Mu]2^2,f1^3 \[Mu]1^4 \[Mu]2^2,b1 f1^3 \[Mu]1^4 \[Mu]2^2,b1^2 f1^3 \[Mu]1^4 \[Mu]2^2,b1^3 f1^3 \[Mu]1^4 \[Mu]2^2,b2 f1^3 \[Mu]1^4 \[Mu]2^2,b1 b2 f1^3 \[Mu]1^4 \[Mu]2^2,b5 f1^3 \[Mu]1^4 \[Mu]2^2,b1 b5 f1^3 \[Mu]1^4 \[Mu]2^2,f1^4 \[Mu]1^4 \[Mu]2^2,b1 f1^4 \[Mu]1^4 \[Mu]2^2,b1^2 f1^4 \[Mu]1^4 \[Mu]2^2,b2 f1^4 \[Mu]1^4 \[Mu]2^2,b5 f1^4 \[Mu]1^4 \[Mu]2^2,f1^4 \[Mu]1^6 \[Mu]2^2,b1 f1^4 \[Mu]1^6 \[Mu]2^2,b1^2 f1^4 \[Mu]1^6 \[Mu]2^2,b2 f1^4 \[Mu]1^6 \[Mu]2^2,b5 f1^4 \[Mu]1^6 \[Mu]2^2,f1^5 \[Mu]1^6 \[Mu]2^2,b1 f1^5 \[Mu]1^6 \[Mu]2^2,f1^5 \[Mu]1^8 \[Mu]2^2,b1 f1^5 \[Mu]1^8 \[Mu]2^2,b1 f1^2 \[Mu]1 \[Mu]2^3,b1^2 f1^2 \[Mu]1 \[Mu]2^3,b2 f1^2 \[Mu]1 \[Mu]2^3,b1 b2 f1^2 \[Mu]1 \[Mu]2^3,b3 f1^2 \[Mu]1 \[Mu]2^3,b1 b3 f1^2 \[Mu]1 \[Mu]2^3,b6 f1^2 \[Mu]1 \[Mu]2^3,b1 b6 f1^2 \[Mu]1 \[Mu]2^3,b8 f1^2 \[Mu]1 \[Mu]2^3,b1 b8 f1^2 \[Mu]1 \[Mu]2^3,f1^3 \[Mu]1^3 \[Mu]2^3,b1 f1^3 \[Mu]1^3 \[Mu]2^3,b1^2 f1^3 \[Mu]1^3 \[Mu]2^3,b2 f1^3 \[Mu]1^3 \[Mu]2^3,b3 f1^3 \[Mu]1^3 \[Mu]2^3,b6 f1^3 \[Mu]1^3 \[Mu]2^3,b8 f1^3 \[Mu]1^3 \[Mu]2^3,f1^4 \[Mu]1^3 \[Mu]2^3,b1 f1^4 \[Mu]1^3 \[Mu]2^3,b1^2 f1^4 \[Mu]1^3 \[Mu]2^3,f1^4 \[Mu]1^5 \[Mu]2^3,b1 f1^4 \[Mu]1^5 \[Mu]2^3,b1^2 f1^4 \[Mu]1^5 \[Mu]2^3,f1^5 \[Mu]1^5 \[Mu]2^3,b1 f1^5 \[Mu]1^5 \[Mu]2^3,f1^5 \[Mu]1^7 \[Mu]2^3,b1 f1^5 \[Mu]1^7 \[Mu]2^3,f1^6 \[Mu]1^7 \[Mu]2^3,f1^6 \[Mu]1^9 \[Mu]2^3,b1 f1^2 \[Mu]2^4,b1^2 f1^2 \[Mu]2^4,b2 f1^2 \[Mu]2^4,b1 b2 f1^2 \[Mu]2^4,b3 f1^2 \[Mu]2^4,b1 b3 f1^2 \[Mu]2^4,b6 f1^2 \[Mu]2^4,b1 b6 f1^2 \[Mu]2^4,b8 f1^2 \[Mu]2^4,b1 b8 f1^2 \[Mu]2^4,f1^3 \[Mu]1^2 \[Mu]2^4,b1 f1^3 \[Mu]1^2 \[Mu]2^4,b1^2 f1^3 \[Mu]1^2 \[Mu]2^4,b2 f1^3 \[Mu]1^2 \[Mu]2^4,b3 f1^3 \[Mu]1^2 \[Mu]2^4,b6 f1^3 \[Mu]1^2 \[Mu]2^4,b8 f1^3 \[Mu]1^2 \[Mu]2^4,f1^4 \[Mu]1^2 \[Mu]2^4,b1 f1^4 \[Mu]1^2 \[Mu]2^4,b1^2 f1^4 \[Mu]1^2 \[Mu]2^4,f1^4 \[Mu]1^4 \[Mu]2^4,b1 f1^4 \[Mu]1^4 \[Mu]2^4,b1^2 f1^4 \[Mu]1^4 \[Mu]2^4,f1^5 \[Mu]1^4 \[Mu]2^4,b1 f1^5 \[Mu]1^4 \[Mu]2^4,f1^5 \[Mu]1^6 \[Mu]2^4,b1 f1^5 \[Mu]1^6 \[Mu]2^4,f1^6 \[Mu]1^6 \[Mu]2^4,f1^6 \[Mu]1^8 \[Mu]2^4,b1 f1^3 \[Mu]1 \[Mu]2^5,b2 f1^3 \[Mu]1 \[Mu]2^5,b3 f1^3 \[Mu]1 \[Mu]2^5,b6 f1^3 \[Mu]1 \[Mu]2^5,b8 f1^3 \[Mu]1 \[Mu]2^5,f1^4 \[Mu]1^3 \[Mu]2^5,b1 f1^4 \[Mu]1^3 \[Mu]2^5,f1^5 \[Mu]1^3 \[Mu]2^5,b1 f1^5 \[Mu]1^3 \[Mu]2^5,f1^5 \[Mu]1^5 \[Mu]2^5,b1 f1^5 \[Mu]1^5 \[Mu]2^5,f1^6 \[Mu]1^5 \[Mu]2^5,f1^6 \[Mu]1^7 \[Mu]2^5};*)


(* ::Input:: *)
(*(*Function that obtains the coefficient of each term in biaslist in each ctab coefficient "coef" *)*)
(*b3212biaslister=biaslisterfn[b3212vars,b3212biaslist];*)


(* ::Input:: *)
(*(*apply decomposition in biases to all ctab coefficients and simplify - this takes 10s *)*)
(**)
(*b3212biascoef=Simplify[b3212biaslister/@b3212coefs];*)
(*b3212ctabdec=MapThread[Append,{b3212exps,b3212biascoef}];*)
(*b3212ctabdec>>(LoopTablesPath<>"b3212ctabdec.mx");*)


(* ::Text:: *)
(*Now we just have to permute the coefficients over the 6 possible permutations of k1, k2 and k3*)


(* ::Input:: *)
(*(*b3212ctabdec=<<"tables_for_B1loop/b3212ctabdec.mx";*)*)


(* ::Input:: *)
(*(*The following code does the permutation*)*)
(*arr  = {k1,\[Mu]1,k2,\[Mu]2,k3,\[Mu]3};*)
(*arrt={k1t,\[Mu]1t,k2t,\[Mu]2t,k3t,\[Mu]3t};*)
(*pt = Permutations[{{k1t,\[Mu]1t},{k2t,\[Mu]2t},{k3t,\[Mu]3t}}]/.{{a_,b_},{c_,d_},{e_,f_}}-> {a,b,c,d,e,f};*)
(**)
(*permreps=Thread[(arr->#)]&/@pt;*)


(* ::Input:: *)
(*permreps*)


(* ::Input:: *)
(*(*get list of permuted biases*)*)
(*b3212biasperm=b3212biaslist/.permreps/.Thread[arrt->arr]//Transpose;*)


(* ::Input:: *)
(*b3212coefsdec=b3212ctabdec[[All,4]];*)
(*(*Get list of permuted coefficients - 1s*)*)
(*b3212permdectemp=b3212coefsdec/.permreps/.Thread[arrt->arr];*)
(*(*transpose to get a more convenient tensor*)*)
(*b3212permdec=Transpose[b3212permdectemp,{3,1,2}];*)


(* ::Input:: *)
(*b3212biasperm>>(LoopTablesPath<>"B3212biasPerm.mx");*)
(*b3212permdec>>(LoopTablesPath<>"B3212coefsPerm.mx");*)


(* ::Subsubsection::Closed:: *)
(*Ctab Exponents*)


(* ::Input:: *)
(*(*b3212ctabsimp=<<(LoopTablesPath<>"b3212ctabsimp.mx");*)*)
(*(*b3212exps=b3212ctabsimp\[LeftDoubleBracket]All,1;;3\[RightDoubleBracket];*)*)


(* ::Input:: *)
(*b3212exps={{-2,1,0},{-1,0,0},{-1,1,0},{0,-1,0},{0,0,0},{0,1,0},{1,-2,0},{1,-1,0},{1,0,0},{1,1,0},{2,-2,0},{2,-1,0},{2,0,0},{2,1,0}};*)


(* ::Subsubsection::Closed:: *)
(*Generate ctab coefs for fisherPoints*)


(* ::Input:: *)
(*b3212permdec=<<(LoopTablesPath<>"B3212coefsPerm.mx");*)
(*b3212permdecfn[k1_,k2_,k3_]=b3212permdec;*)


(* ::Input:: *)
(*(*Export ctab exponents*)*)
(*Export["GitHub/python-integer-power-project/3. Ctabs/B3212ctab.csv",b3212exps];*)


(* ::Input:: *)
(*genb3212coef[k1sub_,k2sub_,k3sub_]:=Module[{temp},*)
(*temp=b3212permdecfn[k1sub,k2sub,k3sub];*)
(*Export[LoopTablesPath<>"B3212coefs/B3212coefs_"<>ToString[k1sub//N]<>"_"<>ToString[k2sub//N]<>"_"<>ToString[k3sub//N]<>"_.mx",temp]*)
(*];*)


(* ::Input:: *)
(*Monitor[Do[genb3212coef@@CMASStriEff[[i]],*)
(*{i,1,Length[fisherPoints]}],i]*)


(* ::Input:: *)
(*Monitor[Do[genb3212coef@@LOWZtriEff[[i]],*)
(*{i,1,Length[fisherPoints]}],i]*)


(* ::Input:: *)
(*Monitor[Do[genb3212coef@@fisherPoints[[i]],*)
(*{i,1,Length[fisherPoints]}],i]*)


(* ::Subsection:: *)
(*B411*)


(* ::Subsubsection:: *)
(*General case*)


(* ::Text:: *)
(*IMPORTANT: In the FFTLog, P(q) is q^\[Nu] -- We  will still do it as a convenient substitution since we only need to keep track of the argument*)


(* ::Input:: *)
(*k411[k1_,k2_,q_]:=k411[k1,k2,q]= K1r[k1]K1r[k2] K4rsym[q,-q,-k2,-k1]/.mag[0]-> 0/.tb;*)


(* ::Input:: *)
(*B411Integrand[k1_,k2_,q_] :=B411Integrand[k1,k2,q]= 12Plin[k1]Plin[k2]Plin[q]k411[k1,k2,q]/.tb;*)


(* ::Text:: *)
(*Full integrand, setting Plin[q]->q^\[Nu] and external P(k) to 1*)


(* ::Input:: *)
(*ker411temp= B411Integrand[k1,k2,q]/.{Plin[q]-> q^\[Nu],Plin[k1]-> 1,Plin[k2]-> 1}//Expand;(*this takes some time*)*)


(* ::Text:: *)
(*Short interruption to generate the monopole kernel for matter in redshift space*)


(* ::Input:: *)
(*ker411Matter=ker411temp/.{\[Nu]->0}/.submat//Expand;*)


(* ::Input:: *)
(*\[Mu]list=getMonomialsSingle[{\[Mu]0,\[Mu]1,\[Mu]2},ker411Matter];*)


(* ::Input:: *)
(*ker411MatterMono=GetMonopolenoSimp[ker411Matter,\[Mu]list];*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]];*)
(*Export["ker411MatterMono.m",ker411MatterMono]*)


(* ::Text:: *)
(*Next we do a sequence of transformations.*)
(*First we map all terms that have k3mq in the denominator to something with k3pq in the denominator by q-> -q (qrev1 and qrev2)*)
(*Next, for terms that have k3pq*k1mq in the denominator we do the shift q-> k1pq.*)
(*For all other terms with k3pq in the denominator, we shift by q-> k2mq.*)
(*In a last step, we map all terms with k1mq or k2pq in the denominator with q-> -q.*)
(*At that point we only have k1pq, q and k2mq in the denominator and in the numerator we use /.restrep to replace all other momenta.*)
(*The integrated PS is only there for the ride, and the argument will be the correct one.*)


(* ::Input:: *)
(*(* Perform the transformations sequentially*)*)
(*ker411list = List@@(ker411temp);*)
(*Monitor[*)
(*Do[*)
(*temp= ker411temp[[i]];*)
(**)
(*arr = -(1/2)Exponent[temp,k3mq];*)
(*If[arr>0,temp=temp/.qrev1/.qrev2];*)
(**)
(*arr2 = Sign[{ -(1/2)Exponent[temp,k3pq], -(1/2)Exponent[temp,k1mq]}];*)
(*If[arr2[[1]]>0&&arr2[[2]]>0,temp=temp/.preqtok1pq/.qtok1pq];*)
(**)
(*arr3 =Sign[ -(1/2)Exponent[temp,k3pq]];*)
(*If[arr3>0,temp=temp/.preqtok2mq/.qtok2mq];*)
(**)
(*arr4 =Sign[{ -(1/2)Exponent[temp,k1mq], -(1/2)Exponent[temp,k2pq]}];*)
(*If[arr4[[1]]>0||arr4[[2]]>0,temp=temp/.qrev1/.qrev2];*)
(*ker411list[[i]] = temp/.restrep//Expand;*)
(*,{i,1,Length[ker411temp]}]//AbsoluteTiming*)
(*,i*100.0/Length[ker411temp]]*)


(* ::Input:: *)
(*ker411final =Total[ker411list]//Expand;*)


(* ::Input:: *)
(*ker411\[Mu]0list=CoefficientList[ker411final,\[Mu]0];*)


(* ::Text:: *)
(*Substitute the \[Mu]0^n transformations we computed before*)


(* ::Input:: *)
(*\[Mu]0sol411=<<(LoopTablesPath<>"mu0sol_nonsym_upto4.m");*)
(*\[Mu]0sol411exp=\[Mu]0sol411[[;;3]]//Expand;*)


(* ::Input:: *)
(*(*We calculate the ctabs for \[Mu]0 and ker411 and add the constant in the \[Mu]0 replacement by hand *)*)
(*\[Mu]0ctab411temp=Prepend[getCtab/@\[Mu]0sol411exp[[2;;]],{{0,0,0,1}}];*)
(*\[Mu]0ctab411=compressor/@\[Mu]0ctab411temp;*)


(* ::Input:: *)
(*ker411\[Mu]0ctabtemp=getCtab/@ker411\[Mu]0list;*)
(*ker411\[Mu]0ctab=compressor/@ker411\[Mu]0ctabtemp;*)


(* ::Input:: *)
(*(*Combine ctabs, join them, and compress the final table*)*)
(*b411ctabjoin=MapThread[JoinCtabs[#1,#2]&,{ker411\[Mu]0ctab,\[Mu]0ctab411}];*)
(*b411ctab=compressor[(Join@@b411ctabjoin)//Expand];*)


(* ::Text:: *)
(*Now, we have P[k1] P[k2] outside, and in the terms there is P(q), or P(k1pq), or P(k2mq), which are taken care of by the power \[Nu] as a placeholder.*)
(*In the tab, we set the external P[k1] P[k2] to 1.*)


(* ::Input:: *)
(*b411exps\[Nu]=b411ctab[[All,1;;3]];*)
(*b411coefstemp=b411ctab[[All,4]];*)


(* ::Text:: *)
(*Now we encode in the ctab what is the argument of the integrated Plin*)


(* ::Input:: *)
(*b411expsfull=b411exps\[Nu]/.{*)
(*{a_.,b_.,c_.+d_. \[Nu]}:>{a,0,b,0,c,1},*)
(*{a_.,b_.+d_. \[Nu],c_.}:>{a,0,b,1,c,0},*)
(*{a_.+d_. \[Nu],b_.,c_.}:>{a,1,b,0,c,0}};*)


(* ::Input:: *)
(*(*combine exponents and coefficients to generate ctab*)*)
(*b411ctabfull=MapThread[Append,{b411expsfull,b411coefstemp}];*)
(*b411ctabsimp=Sort[compressorB411[b411ctabfull]];*)
(*b411ctabsimp>>(LoopTablesPath<>"b411ctabsimp.mx");*)


(* ::Subsubsection::Closed:: *)
(*Get UV part and subtract from main ctab*)


(* ::Input:: *)
(*(*b411ctabsimp=<<(LoopTablesPath<>"b411ctabsimp.mx");*)*)


(* ::Text:: *)
(*Strategy is calculate the series separately for each term of k411, and then add everything together*)


(* ::Input:: *)
(*k4temp=K4rsym[q,-q,-k2,-k1]/.tb;*)
(*k4\[Mu]0list=CoefficientList[k4temp,\[Mu]0];*)
(*k4\[Mu]0listexp=Expand[k4\[Mu]0list];*)
(*serk4\[Mu]0list=ConstantArray[0,Length[k4\[Mu]0listexp]];*)
(*Monitor[*)
(*Do[serk4\[Mu]0list[[i]]=*)
(*Series[k4\[Mu]0listexp[[i]]/.magrevreps,{q,\[Infinity],2},Assumptions->{k1>0,k2>0,-1<x<1,0<y<1}],{i,1,Length[k4\[Mu]0listexp]}];*)
(*,i]*)


(* ::Input:: *)
(*serk4\[Mu]0list>>(LoopTablesPath<>"serk4\[Mu]0list.mx");*)
(**)


(* ::Input:: *)
(*(*Do the same for monopole matter only*)*)
(*Serker411MatterMono=Series[ker411MatterMono/.magrevreps,{q,\[Infinity],2},Assumptions->{k1>0,k2>0,-1<x<1,0<y<1}];*)


(* ::Input:: *)
(*(*serk411\[Mu]0list=<<(LoopTablesPath<>"serk4\[Mu]0list.mx");*)*)


(* ::Input:: *)
(*(*Check coefficient is zero for q term in series*)*)
(*Expand[SeriesCoefficient[serk4\[Mu]0list[[1]],-1]/.{Sqrt[(-1+x^2) (-1+y^2)]->Sqrt[1-x^2] Sqrt[1-y^2]}]*)
(*Expand[SeriesCoefficient[serk4\[Mu]0list[[2]],-1]/.{Sqrt[(-1+x^2) (-1+y^2)]->Sqrt[1-x^2] Sqrt[1-y^2]}]*)
(*Expand[SeriesCoefficient[serk4\[Mu]0list[[3]],-1]/.{Sqrt[(-1+x^2) (-1+y^2)]->Sqrt[1-x^2] Sqrt[1-y^2]}]*)


(* ::Input:: *)
(*(*Check coefficient is zero for 1/q term in series*)*)
(*Expand[SeriesCoefficient[serk4\[Mu]0list[[1]],1]/.{Sqrt[(-1+x^2) (-1+y^2)]->Sqrt[1-x^2] Sqrt[1-y^2]}]*)
(*Expand[SeriesCoefficient[serk4\[Mu]0list[[2]],1]/.{Sqrt[(-1+x^2) (-1+y^2)]->Sqrt[1-x^2] Sqrt[1-y^2]}]*)
(*Expand[SeriesCoefficient[serk4\[Mu]0list[[3]],1]/.{Sqrt[(-1+x^2) (-1+y^2)]->Sqrt[1-x^2] Sqrt[1-y^2]}]*)


(* ::Input:: *)
(*\[Mu]0sol4=<<(LoopTablesPath<>"mu0sol_nonsym_upto4.m");*)
(*\[Mu]0sol4=\[Mu]0sol4[[;;3]]//Expand;*)


(* ::Input:: *)
(*ser\[Mu]0sol4list=ConstantArray[0,Length[\[Mu]0sol4]];*)
(*Monitor[*)
(*Do[*)
(*ser\[Mu]0sol4list[[i]]=*)
(*Series[\[Mu]0sol4[[i]]/.magrevreps,{q,\[Infinity],2},Assumptions->{k1>0,k2>0,-1<x<1,0<y<1}],{i,1,Length[\[Mu]0sol4]}];*)
(*,i]*)


(* ::Input:: *)
(*(*Check coefficient is zero for 1/q term in series*)*)
(*Expand[SeriesCoefficient[ser\[Mu]0sol4list[[2]],2]/.{Sqrt[(-1+x^2) (-1+y^2)]->Sqrt[1-x^2] Sqrt[1-y^2]}]*)
(*Expand[SeriesCoefficient[ser\[Mu]0sol4list[[3]],2]/.{Sqrt[(-1+x^2) (-1+y^2)]->Sqrt[1-x^2] Sqrt[1-y^2]}]*)


(* ::Input:: *)
(*serk4=serk4\[Mu]0list.ser\[Mu]0sol4list;*)


(* ::Text:: *)
(*Let's consider the 0th order term, and integrate the various terms. Use rules for \[Mu]0 integration*)


(* ::Input:: *)
(*ser0k4temp=Expand[SeriesCoefficient[serk4,0]/.{Sqrt[(-1+x^2) (-1+y^2)]->Sqrt[1-x^2] Sqrt[1-y^2]}];*)
(*ser0k4int=ser0k4temp/.rules\[Mu]04/.rules\[Mu]03/.rules\[Mu]02/.rules\[Mu]01;*)
(*ser0k4intsimp=Simplify[ser0k4int];*)


(* ::Text:: *)
(*Check the coefficient 1/q, it should be 0*)


(* ::Input:: *)
(*ser1k4temp=Expand[SeriesCoefficient[serk4,1]/.{Sqrt[(-1+x^2) (-1+y^2)]->Sqrt[1-x^2] Sqrt[1-y^2]}];*)
(*ser1k4int=ser1k4temp/.rules\[Mu]04/.rules\[Mu]03/.rules\[Mu]02/.rules\[Mu]01;*)
(*ser1k4intsimp=Simplify[ser1k4int]*)


(* ::Input:: *)
(*uvK4k0=(ser0k4intsimp/.ysub)//Simplify;*)


(* ::Input:: *)
(*B411UVIntegrand[k1_,k2_,q_] :=B411UVIntegrand[k1,k2,q]= 12Plin[k1]Plin[k2]Plin[q] K1r[k1]K1r[k2] uvK4k0/.tb*)


(* ::Text:: *)
(*Full integrand, setting Plin[q]->q^\[Nu] and external P(k) to 1*)


(* ::Input:: *)
(*ker411uvtemp= B411UVIntegrand[k1,k2,q]/.{Plin[q]-> q^\[Nu],Plin[k1]-> 1,Plin[k2]-> 1}//Expand;*)


(* ::Input:: *)
(* (*ker411uvtemp>>(LoopTablesPath<>"B411_uvk0.mx")*)*)


(* ::Input:: *)
(* ker411uvtemp=<<(LoopTablesPath<>"B411_uvk0.mx");*)


(* ::Input:: *)
(*b411uvsubctab=compressor[getCtab[-ker411uvtemp]];*)


(* ::Input:: *)
(*b411uvexps\[Nu]=b411uvsubctab[[All,1;;3]];*)
(*b411uvcoefstemp=b411uvsubctab[[All,4]];*)


(* ::Text:: *)
(*Now we encode in the ctab what is the argument of the integrated Plin and we merge with the main b411ctab*)


(* ::Input:: *)
(*(*Load b411 calculated before*)*)
(*b411ctabsimp=<<(LoopTablesPath<>"b411ctabsimp.mx");*)


(* ::Input:: *)
(*(*written this way so that it is easier to include higher order UV subtraction if needed*)*)
(**)
(*b411uvexpsfull=b411uvexps\[Nu]/.{{a_.+d_. \[Nu],b_.,c_.}:>{a,1,b,0,c,0}};*)
(*b411uvsubctabfull=MapThread[Append,{b411uvexpsfull,b411uvcoefstemp}];*)
(*b411ctabtottemp=Join[b411ctabsimp,b411uvsubctabfull];*)
(*b411ctabtot=Sort[compressorB411[b411ctabtottemp//Expand]]//Expand;*)


(* ::Input:: *)
(*b411ctabtot>>(LoopTablesPath<>"b411ctab_subk0.mx");*)


(* ::Input:: *)
(*b411ctabtot=<<(LoopTablesPath<>"b411ctab_subk0.mx");*)


(* ::Subsubsection:: *)
(*Decompose into bias coefficients*)


(* ::Text:: *)
(*Now we want to decompose into a basis of bias coefficients. Note that b10 is present in b411 before uv subtracting, but not after.*)


(* ::Input:: *)
(*b411ctabtot=<<(LoopTablesPath<>"b411ctab_subk0.mx");*)


(* ::Input:: *)
(*b411exps=b411ctabtot[[All,1;;6]];*)
(*b411coefs=b411ctabtot[[All,7]]//Expand;*)


(* ::Input:: *)
(*b411vars={b1,b11,b12,b13,b15,b2,b3,b4,b5,b6,b7,b8,b9,f1,\[Mu]1,\[Mu]2};*)


(* ::Input:: *)
(*(*Get all bias monomials that appear in the b411 kernel*)*)
(*(*b411biaslist=getMonomials[b411vars,b411ctabtot];*)*)
(*b411biaslist={b1^3,b1^2 b11,b1^2 b12,b1^2 b13,b1^2 b15,b1^2 b2,b1^2 b3,b1^2 b4,b1^2 b5,b1^2 b6,b1^2 b7,b1^2 b8,b1^2 b9,b1^2 f1 \[Mu]1^2,b1^3 f1 \[Mu]1^2,b1 b11 f1 \[Mu]1^2,b1 b12 f1 \[Mu]1^2,b1 b13 f1 \[Mu]1^2,b1 b15 f1 \[Mu]1^2,b1 b2 f1 \[Mu]1^2,b1^2 b2 f1 \[Mu]1^2,b1 b3 f1 \[Mu]1^2,b1^2 b3 f1 \[Mu]1^2,b1 b4 f1 \[Mu]1^2,b1 b5 f1 \[Mu]1^2,b1^2 b5 f1 \[Mu]1^2,b1 b6 f1 \[Mu]1^2,b1^2 b6 f1 \[Mu]1^2,b1 b7 f1 \[Mu]1^2,b1 b8 f1 \[Mu]1^2,b1^2 b8 f1 \[Mu]1^2,b1 b9 f1 \[Mu]1^2,b1^2 f1^2 \[Mu]1^2,b1^3 f1^2 \[Mu]1^2,b1^2 b2 f1^2 \[Mu]1^2,b1^2 b5 f1^2 \[Mu]1^2,b1 f1^2 \[Mu]1^4,b1^2 f1^2 \[Mu]1^4,b1^3 f1^2 \[Mu]1^4,b1 b2 f1^2 \[Mu]1^4,b1^2 b2 f1^2 \[Mu]1^4,b1 b3 f1^2 \[Mu]1^4,b1 b5 f1^2 \[Mu]1^4,b1^2 b5 f1^2 \[Mu]1^4,b1 b6 f1^2 \[Mu]1^4,b1 b8 f1^2 \[Mu]1^4,b1 f1^3 \[Mu]1^4,b1^2 f1^3 \[Mu]1^4,b1^3 f1^3 \[Mu]1^4,b1 b2 f1^3 \[Mu]1^4,b1 b5 f1^3 \[Mu]1^4,b1 f1^3 \[Mu]1^6,b1^2 f1^3 \[Mu]1^6,b1^3 f1^3 \[Mu]1^6,b1 b2 f1^3 \[Mu]1^6,b1 b5 f1^3 \[Mu]1^6,b1 f1^4 \[Mu]1^6,b1^2 f1^4 \[Mu]1^6,b1 f1^4 \[Mu]1^8,b1^2 f1^4 \[Mu]1^8,b1^2 f1 \[Mu]1 \[Mu]2,b1^3 f1 \[Mu]1 \[Mu]2,b1^2 b2 f1 \[Mu]1 \[Mu]2,b1^2 b3 f1 \[Mu]1 \[Mu]2,b1^2 b5 f1 \[Mu]1 \[Mu]2,b1^2 b6 f1 \[Mu]1 \[Mu]2,b1^2 b8 f1 \[Mu]1 \[Mu]2,b1^2 f1^2 \[Mu]1 \[Mu]2,b1^3 f1^2 \[Mu]1 \[Mu]2,b1^2 b2 f1^2 \[Mu]1 \[Mu]2,b1^2 b5 f1^2 \[Mu]1 \[Mu]2,b1 f1^2 \[Mu]1^3 \[Mu]2,b1^2 f1^2 \[Mu]1^3 \[Mu]2,b1^3 f1^2 \[Mu]1^3 \[Mu]2,b1 b2 f1^2 \[Mu]1^3 \[Mu]2,b1^2 b2 f1^2 \[Mu]1^3 \[Mu]2,b1 b3 f1^2 \[Mu]1^3 \[Mu]2,b1 b5 f1^2 \[Mu]1^3 \[Mu]2,b1^2 b5 f1^2 \[Mu]1^3 \[Mu]2,b1 b6 f1^2 \[Mu]1^3 \[Mu]2,b1 b8 f1^2 \[Mu]1^3 \[Mu]2,b1 f1^3 \[Mu]1^3 \[Mu]2,b1^2 f1^3 \[Mu]1^3 \[Mu]2,b1^3 f1^3 \[Mu]1^3 \[Mu]2,b1 b2 f1^3 \[Mu]1^3 \[Mu]2,b1 b5 f1^3 \[Mu]1^3 \[Mu]2,b1 f1^3 \[Mu]1^5 \[Mu]2,b1^2 f1^3 \[Mu]1^5 \[Mu]2,b1^3 f1^3 \[Mu]1^5 \[Mu]2,b1 b2 f1^3 \[Mu]1^5 \[Mu]2,b1 b5 f1^3 \[Mu]1^5 \[Mu]2,b1 f1^4 \[Mu]1^5 \[Mu]2,b1^2 f1^4 \[Mu]1^5 \[Mu]2,b1 f1^4 \[Mu]1^7 \[Mu]2,b1^2 f1^4 \[Mu]1^7 \[Mu]2,b1 f1^5 \[Mu]1^7 \[Mu]2,b1 f1^5 \[Mu]1^9 \[Mu]2,b1^2 f1 \[Mu]2^2,b1^3 f1 \[Mu]2^2,b1 b11 f1 \[Mu]2^2,b1 b12 f1 \[Mu]2^2,b1 b13 f1 \[Mu]2^2,b1 b15 f1 \[Mu]2^2,b1 b2 f1 \[Mu]2^2,b1^2 b2 f1 \[Mu]2^2,b1 b3 f1 \[Mu]2^2,b1^2 b3 f1 \[Mu]2^2,b1 b4 f1 \[Mu]2^2,b1 b5 f1 \[Mu]2^2,b1^2 b5 f1 \[Mu]2^2,b1 b6 f1 \[Mu]2^2,b1^2 b6 f1 \[Mu]2^2,b1 b7 f1 \[Mu]2^2,b1 b8 f1 \[Mu]2^2,b1^2 b8 f1 \[Mu]2^2,b1 b9 f1 \[Mu]2^2,b1^2 f1^2 \[Mu]2^2,b1^3 f1^2 \[Mu]2^2,b1^2 b2 f1^2 \[Mu]2^2,b1^2 b5 f1^2 \[Mu]2^2,b1 f1^2 \[Mu]1^2 \[Mu]2^2,b1^2 f1^2 \[Mu]1^2 \[Mu]2^2,b1^3 f1^2 \[Mu]1^2 \[Mu]2^2,b11 f1^2 \[Mu]1^2 \[Mu]2^2,b12 f1^2 \[Mu]1^2 \[Mu]2^2,b13 f1^2 \[Mu]1^2 \[Mu]2^2,b15 f1^2 \[Mu]1^2 \[Mu]2^2,b2 f1^2 \[Mu]1^2 \[Mu]2^2,b1 b2 f1^2 \[Mu]1^2 \[Mu]2^2,b1^2 b2 f1^2 \[Mu]1^2 \[Mu]2^2,b3 f1^2 \[Mu]1^2 \[Mu]2^2,b1 b3 f1^2 \[Mu]1^2 \[Mu]2^2,b4 f1^2 \[Mu]1^2 \[Mu]2^2,b5 f1^2 \[Mu]1^2 \[Mu]2^2,b1 b5 f1^2 \[Mu]1^2 \[Mu]2^2,b1^2 b5 f1^2 \[Mu]1^2 \[Mu]2^2,b6 f1^2 \[Mu]1^2 \[Mu]2^2,b1 b6 f1^2 \[Mu]1^2 \[Mu]2^2,b7 f1^2 \[Mu]1^2 \[Mu]2^2,b8 f1^2 \[Mu]1^2 \[Mu]2^2,b1 b8 f1^2 \[Mu]1^2 \[Mu]2^2,b9 f1^2 \[Mu]1^2 \[Mu]2^2,b1 f1^3 \[Mu]1^2 \[Mu]2^2,b1^2 f1^3 \[Mu]1^2 \[Mu]2^2,b1^3 f1^3 \[Mu]1^2 \[Mu]2^2,b1 b2 f1^3 \[Mu]1^2 \[Mu]2^2,b1 b5 f1^3 \[Mu]1^2 \[Mu]2^2,f1^3 \[Mu]1^4 \[Mu]2^2,b1 f1^3 \[Mu]1^4 \[Mu]2^2,b1^2 f1^3 \[Mu]1^4 \[Mu]2^2,b1^3 f1^3 \[Mu]1^4 \[Mu]2^2,b2 f1^3 \[Mu]1^4 \[Mu]2^2,b1 b2 f1^3 \[Mu]1^4 \[Mu]2^2,b3 f1^3 \[Mu]1^4 \[Mu]2^2,b5 f1^3 \[Mu]1^4 \[Mu]2^2,b1 b5 f1^3 \[Mu]1^4 \[Mu]2^2,b6 f1^3 \[Mu]1^4 \[Mu]2^2,b8 f1^3 \[Mu]1^4 \[Mu]2^2,f1^4 \[Mu]1^4 \[Mu]2^2,b1 f1^4 \[Mu]1^4 \[Mu]2^2,b1^2 f1^4 \[Mu]1^4 \[Mu]2^2,b2 f1^4 \[Mu]1^4 \[Mu]2^2,b5 f1^4 \[Mu]1^4 \[Mu]2^2,f1^4 \[Mu]1^6 \[Mu]2^2,b1 f1^4 \[Mu]1^6 \[Mu]2^2,b1^2 f1^4 \[Mu]1^6 \[Mu]2^2,b2 f1^4 \[Mu]1^6 \[Mu]2^2,b5 f1^4 \[Mu]1^6 \[Mu]2^2,f1^5 \[Mu]1^6 \[Mu]2^2,b1 f1^5 \[Mu]1^6 \[Mu]2^2,f1^5 \[Mu]1^8 \[Mu]2^2,b1 f1^5 \[Mu]1^8 \[Mu]2^2,b1 f1^2 \[Mu]1 \[Mu]2^3,b1^2 f1^2 \[Mu]1 \[Mu]2^3,b1^3 f1^2 \[Mu]1 \[Mu]2^3,b1 b2 f1^2 \[Mu]1 \[Mu]2^3,b1^2 b2 f1^2 \[Mu]1 \[Mu]2^3,b1 b3 f1^2 \[Mu]1 \[Mu]2^3,b1 b5 f1^2 \[Mu]1 \[Mu]2^3,b1^2 b5 f1^2 \[Mu]1 \[Mu]2^3,b1 b6 f1^2 \[Mu]1 \[Mu]2^3,b1 b8 f1^2 \[Mu]1 \[Mu]2^3,b1 f1^3 \[Mu]1 \[Mu]2^3,b1^2 f1^3 \[Mu]1 \[Mu]2^3,b1^3 f1^3 \[Mu]1 \[Mu]2^3,b1 b2 f1^3 \[Mu]1 \[Mu]2^3,b1 b5 f1^3 \[Mu]1 \[Mu]2^3,f1^3 \[Mu]1^3 \[Mu]2^3,b1 f1^3 \[Mu]1^3 \[Mu]2^3,b1^2 f1^3 \[Mu]1^3 \[Mu]2^3,b1^3 f1^3 \[Mu]1^3 \[Mu]2^3,b2 f1^3 \[Mu]1^3 \[Mu]2^3,b1 b2 f1^3 \[Mu]1^3 \[Mu]2^3,b3 f1^3 \[Mu]1^3 \[Mu]2^3,b5 f1^3 \[Mu]1^3 \[Mu]2^3,b1 b5 f1^3 \[Mu]1^3 \[Mu]2^3,b6 f1^3 \[Mu]1^3 \[Mu]2^3,b8 f1^3 \[Mu]1^3 \[Mu]2^3,f1^4 \[Mu]1^3 \[Mu]2^3,b1 f1^4 \[Mu]1^3 \[Mu]2^3,b1^2 f1^4 \[Mu]1^3 \[Mu]2^3,b2 f1^4 \[Mu]1^3 \[Mu]2^3,b5 f1^4 \[Mu]1^3 \[Mu]2^3,f1^4 \[Mu]1^5 \[Mu]2^3,b1 f1^4 \[Mu]1^5 \[Mu]2^3,b1^2 f1^4 \[Mu]1^5 \[Mu]2^3,b2 f1^4 \[Mu]1^5 \[Mu]2^3,b5 f1^4 \[Mu]1^5 \[Mu]2^3,f1^5 \[Mu]1^5 \[Mu]2^3,b1 f1^5 \[Mu]1^5 \[Mu]2^3,f1^5 \[Mu]1^7 \[Mu]2^3,b1 f1^5 \[Mu]1^7 \[Mu]2^3,f1^6 \[Mu]1^7 \[Mu]2^3,f1^6 \[Mu]1^9 \[Mu]2^3,b1 f1^2 \[Mu]2^4,b1^2 f1^2 \[Mu]2^4,b1^3 f1^2 \[Mu]2^4,b1 b2 f1^2 \[Mu]2^4,b1^2 b2 f1^2 \[Mu]2^4,b1 b3 f1^2 \[Mu]2^4,b1 b5 f1^2 \[Mu]2^4,b1^2 b5 f1^2 \[Mu]2^4,b1 b6 f1^2 \[Mu]2^4,b1 b8 f1^2 \[Mu]2^4,b1 f1^3 \[Mu]2^4,b1^2 f1^3 \[Mu]2^4,b1^3 f1^3 \[Mu]2^4,b1 b2 f1^3 \[Mu]2^4,b1 b5 f1^3 \[Mu]2^4,f1^3 \[Mu]1^2 \[Mu]2^4,b1 f1^3 \[Mu]1^2 \[Mu]2^4,b1^2 f1^3 \[Mu]1^2 \[Mu]2^4,b1^3 f1^3 \[Mu]1^2 \[Mu]2^4,b2 f1^3 \[Mu]1^2 \[Mu]2^4,b1 b2 f1^3 \[Mu]1^2 \[Mu]2^4,b3 f1^3 \[Mu]1^2 \[Mu]2^4,b5 f1^3 \[Mu]1^2 \[Mu]2^4,b1 b5 f1^3 \[Mu]1^2 \[Mu]2^4,b6 f1^3 \[Mu]1^2 \[Mu]2^4,b8 f1^3 \[Mu]1^2 \[Mu]2^4,f1^4 \[Mu]1^2 \[Mu]2^4,b1 f1^4 \[Mu]1^2 \[Mu]2^4,b1^2 f1^4 \[Mu]1^2 \[Mu]2^4,b2 f1^4 \[Mu]1^2 \[Mu]2^4,b5 f1^4 \[Mu]1^2 \[Mu]2^4,f1^4 \[Mu]1^4 \[Mu]2^4,b1 f1^4 \[Mu]1^4 \[Mu]2^4,b1^2 f1^4 \[Mu]1^4 \[Mu]2^4,b2 f1^4 \[Mu]1^4 \[Mu]2^4,b5 f1^4 \[Mu]1^4 \[Mu]2^4,f1^5 \[Mu]1^4 \[Mu]2^4,b1 f1^5 \[Mu]1^4 \[Mu]2^4,f1^5 \[Mu]1^6 \[Mu]2^4,b1 f1^5 \[Mu]1^6 \[Mu]2^4,f1^6 \[Mu]1^6 \[Mu]2^4,f1^6 \[Mu]1^8 \[Mu]2^4,b1 f1^3 \[Mu]1 \[Mu]2^5,b1^2 f1^3 \[Mu]1 \[Mu]2^5,b1^3 f1^3 \[Mu]1 \[Mu]2^5,b1 b2 f1^3 \[Mu]1 \[Mu]2^5,b1 b5 f1^3 \[Mu]1 \[Mu]2^5,b1 f1^4 \[Mu]1 \[Mu]2^5,b1^2 f1^4 \[Mu]1 \[Mu]2^5,f1^4 \[Mu]1^3 \[Mu]2^5,b1 f1^4 \[Mu]1^3 \[Mu]2^5,b1^2 f1^4 \[Mu]1^3 \[Mu]2^5,b2 f1^4 \[Mu]1^3 \[Mu]2^5,b5 f1^4 \[Mu]1^3 \[Mu]2^5,f1^5 \[Mu]1^3 \[Mu]2^5,b1 f1^5 \[Mu]1^3 \[Mu]2^5,f1^5 \[Mu]1^5 \[Mu]2^5,b1 f1^5 \[Mu]1^5 \[Mu]2^5,f1^6 \[Mu]1^5 \[Mu]2^5,f1^6 \[Mu]1^7 \[Mu]2^5,b1 f1^3 \[Mu]2^6,b1^2 f1^3 \[Mu]2^6,b1^3 f1^3 \[Mu]2^6,b1 b2 f1^3 \[Mu]2^6,b1 b5 f1^3 \[Mu]2^6,b1 f1^4 \[Mu]2^6,b1^2 f1^4 \[Mu]2^6,f1^4 \[Mu]1^2 \[Mu]2^6,b1 f1^4 \[Mu]1^2 \[Mu]2^6,b1^2 f1^4 \[Mu]1^2 \[Mu]2^6,b2 f1^4 \[Mu]1^2 \[Mu]2^6,b5 f1^4 \[Mu]1^2 \[Mu]2^6,f1^5 \[Mu]1^2 \[Mu]2^6,b1 f1^5 \[Mu]1^2 \[Mu]2^6,f1^5 \[Mu]1^4 \[Mu]2^6,b1 f1^5 \[Mu]1^4 \[Mu]2^6,f1^6 \[Mu]1^4 \[Mu]2^6,f1^6 \[Mu]1^6 \[Mu]2^6,b1 f1^4 \[Mu]1 \[Mu]2^7,b1^2 f1^4 \[Mu]1 \[Mu]2^7,b1 f1^5 \[Mu]1 \[Mu]2^7,f1^5 \[Mu]1^3 \[Mu]2^7,b1 f1^5 \[Mu]1^3 \[Mu]2^7,f1^6 \[Mu]1^3 \[Mu]2^7,f1^6 \[Mu]1^5 \[Mu]2^7,b1 f1^4 \[Mu]2^8,b1^2 f1^4 \[Mu]2^8,f1^5 \[Mu]1^2 \[Mu]2^8,b1 f1^5 \[Mu]1^2 \[Mu]2^8,f1^6 \[Mu]1^4 \[Mu]2^8,b1 f1^5 \[Mu]1 \[Mu]2^9,f1^6 \[Mu]1^3 \[Mu]2^9};*)


(* ::Input:: *)
(*(*Function that obtains the coefficient of each term in biaslist in each ctab coefficient "coef" *)*)
(*b411biaslister=biaslisterfn[b411vars,b411biaslist];*)


(* ::Input:: *)
(*(*apply decomposition in biases to all ctab coefficients and simplify - this takes 4 minutes*)*)
(*b411biascoef=Simplify[b411biaslister/@b411coefs];//AbsoluteTiming*)


(* ::Input:: *)
(*b411ctabdec=MapThread[Append,{b411exps,b411biascoef}];*)


(* ::Input:: *)
(*(*b411ctabdec>>(LoopTablesPath<>"b411ctabdec.mx");*)*)


(* ::Input:: *)
(*(*b411ctabdec=<<(LoopTablesPath<>"b411ctabdec.mx");*)*)


(* ::Text:: *)
(*Now we just have to permute the coefficients over the 3 cyclic permutations of k1, k2 and k3*)


(* ::Input:: *)
(*(*The following code does the permutation*)*)
(*arr  = {k1,\[Mu]1,k2,\[Mu]2,k3,\[Mu]3};*)
(*arrt={k1t,\[Mu]1t,k2t,\[Mu]2t,k3t,\[Mu]3t};*)
(*cycpt=NestList[RotateLeft,{{k1t,\[Mu]1t},{k2t,\[Mu]2t},{k3t,\[Mu]3t}},2]/.{{a_,b_},{c_,d_},{e_,f_}}-> {a,b,c,d,e,f};*)
(*permreps=Thread[(arr->#)]&/@cycpt;*)


(* ::Input:: *)
(*(*get list of permuted biases*)*)
(*b411biasperm=b411biaslist/.permreps/.Thread[arrt->arr]//Transpose;*)
(*b411coefsdec=b411ctabdec[[All,7]];*)


(* ::Input:: *)
(*(*Get list of permuted coefficients - 13s*)*)
(*b411permdectemp=b411coefsdec/.permreps/.Thread[arrt->arr];//AbsoluteTiming*)


(* ::Input:: *)
(*(*transpose to get a more convenient tensor*)*)
(*b411permdec=Transpose[b411permdectemp,{3,1,2}];*)


(* ::Input:: *)
(*b411biasperm>>(LoopTablesPath<>"B411biasPerm.mx");*)


(* ::Input:: *)
(*b411permdec>>(LoopTablesPath<>"B411coefsPerm.mx");*)


(* ::Subsubsection::Closed:: *)
(*Exponents for ctab*)


(* ::Input:: *)
(*(*b411ctabtot=<<(LoopTablesPath<>"b411ctab_subk0.mx");*)*)
(*(*b411exps=b411ctabtot\[LeftDoubleBracket]All,1;;6\[RightDoubleBracket];*)*)


(* ::Input:: *)
(*b411exps={{-4,0,0,0,2,1},{-4,0,1,0,1,1},{-4,0,1,0,2,1},{-4,1,1,0,1,0},{-3,0,-1,0,2,1},{-3,0,0,0,1,1},{-3,0,0,0,2,1},{-3,0,1,0,0,1},{-3,0,1,0,1,1},{-3,0,1,0,2,1},{-3,1,0,0,1,0},{-3,1,1,0,0,0},{-3,1,1,0,1,0},{-2,0,-2,0,2,1},{-2,0,-1,0,1,1},{-2,0,-1,0,2,1},{-2,0,0,0,0,1},{-2,0,0,0,1,1},{-2,0,0,0,2,1},{-2,0,0,1,1,0},{-2,0,1,0,-1,1},{-2,0,1,0,0,1},{-2,0,1,0,1,1},{-2,0,1,0,2,1},{-2,0,1,1,1,0},{-2,0,2,1,1,0},{-2,1,-1,0,1,0},{-2,1,0,0,0,0},{-2,1,0,0,1,0},{-2,1,1,0,-1,0},{-2,1,1,0,0,0},{-2,1,1,0,1,0},{-1,0,-3,0,2,1},{-1,0,-2,0,1,1},{-1,0,-2,0,2,1},{-1,0,-1,0,0,1},{-1,0,-1,0,1,1},{-1,0,-1,0,2,1},{-1,0,-1,1,1,0},{-1,0,0,0,-1,1},{-1,0,0,0,0,1},{-1,0,0,0,1,1},{-1,0,0,0,2,1},{-1,0,0,1,0,0},{-1,0,0,1,1,0},{-1,0,1,0,-2,1},{-1,0,1,0,-1,1},{-1,0,1,0,0,1},{-1,0,1,0,1,1},{-1,0,1,0,2,1},{-1,0,1,1,0,0},{-1,0,1,1,1,0},{-1,0,2,1,0,0},{-1,0,2,1,1,0},{-1,1,-2,0,1,0},{-1,1,-1,0,0,0},{-1,1,-1,0,1,0},{-1,1,0,0,-1,0},{-1,1,0,0,0,0},{-1,1,0,0,1,0},{-1,1,1,0,-2,0},{-1,1,1,0,-1,0},{-1,1,1,0,0,0},{-1,1,1,0,1,0},{0,0,-4,0,2,1},{0,0,-3,0,1,1},{0,0,-3,0,2,1},{0,0,-2,0,0,1},{0,0,-2,0,1,1},{0,0,-2,0,2,1},{0,0,-2,1,1,0},{0,0,-1,0,-1,1},{0,0,-1,0,0,1},{0,0,-1,0,1,1},{0,0,-1,0,2,1},{0,0,-1,1,0,0},{0,0,-1,1,1,0},{0,0,0,0,-2,1},{0,0,0,0,-1,1},{0,0,0,0,0,1},{0,0,0,0,1,1},{0,0,0,0,2,1},{0,0,0,1,-1,0},{0,0,0,1,0,0},{0,0,0,1,1,0},{0,0,1,0,-2,1},{0,0,1,0,-1,1},{0,0,1,0,0,1},{0,0,1,0,1,1},{0,0,1,0,2,1},{0,0,1,1,-1,0},{0,0,1,1,0,0},{0,0,1,1,1,0},{0,0,2,1,-1,0},{0,0,2,1,0,0},{0,0,2,1,1,0},{0,1,-3,0,1,0},{0,1,-2,0,0,0},{0,1,-2,0,1,0},{0,1,-1,0,-1,0},{0,1,-1,0,0,0},{0,1,-1,0,1,0},{0,1,0,0,-2,0},{0,1,0,0,-1,0},{0,1,0,0,0,0},{0,1,0,0,1,0},{0,1,1,0,-3,0},{0,1,1,0,-2,0},{0,1,1,0,-1,0},{0,1,1,0,0,0},{0,1,1,0,1,0},{1,0,-2,0,0,1},{1,0,-2,0,1,1},{1,0,-2,0,2,1},{1,0,-2,1,0,0},{1,0,-2,1,1,0},{1,0,-1,0,-1,1},{1,0,-1,0,0,1},{1,0,-1,0,1,1},{1,0,-1,0,2,1},{1,0,-1,1,-1,0},{1,0,-1,1,0,0},{1,0,-1,1,1,0},{1,0,0,0,-2,1},{1,0,0,0,-1,1},{1,0,0,0,0,1},{1,0,0,0,1,1},{1,0,0,0,2,1},{1,0,0,1,-2,0},{1,0,0,1,-1,0},{1,0,0,1,0,0},{1,0,0,1,1,0},{1,0,1,0,-2,1},{1,0,1,0,-1,1},{1,0,1,0,0,1},{1,0,1,0,1,1},{1,0,1,0,2,1},{1,0,1,1,-2,0},{1,0,1,1,-1,0},{1,0,1,1,0,0},{1,0,1,1,1,0},{1,0,2,1,-2,0},{1,0,2,1,-1,0},{1,0,2,1,0,0},{1,0,2,1,1,0},{1,1,-3,0,0,0},{1,1,-3,0,1,0},{1,1,-2,0,-1,0},{1,1,-2,0,0,0},{1,1,-2,0,1,0},{1,1,-1,0,-2,0},{1,1,-1,0,-1,0},{1,1,-1,0,0,0},{1,1,-1,0,1,0},{1,1,0,0,-3,0},{1,1,0,0,-2,0},{1,1,0,0,-1,0},{1,1,0,0,0,0},{1,1,0,0,1,0},{1,1,1,0,-3,0},{1,1,1,0,-2,0},{1,1,1,0,-1,0},{1,1,1,0,0,0},{1,1,1,0,1,0},{2,1,-4,0,0,0},{2,1,-3,0,-1,0},{2,1,-3,0,0,0},{2,1,-3,0,1,0},{2,1,-2,0,-2,0},{2,1,-2,0,-1,0},{2,1,-2,0,0,0},{2,1,-2,0,1,0},{2,1,-1,0,-3,0},{2,1,-1,0,-2,0},{2,1,-1,0,-1,0},{2,1,-1,0,0,0},{2,1,-1,0,1,0},{2,1,0,0,-4,0},{2,1,0,0,-3,0},{2,1,0,0,-2,0},{2,1,0,0,-1,0},{2,1,0,0,0,0},{2,1,0,0,1,0},{2,1,1,0,-3,0},{2,1,1,0,-2,0},{2,1,1,0,-1,0},{2,1,1,0,0,0},{2,1,1,0,1,0}};*)


(* ::Subsubsection::Closed:: *)
(*Generate ctabs coefs for fisherpoints*)


(* ::Input:: *)
(*b411permdec=<<(LoopTablesPath<>"B411coefsPerm.mx");*)
(*b411permdecfn[k1_,k2_,k3_]=b411permdec;*)


(* ::Input:: *)
(*Export[ctabpath<>"B411ctab.csv",b411exps];*)


(* ::Input:: *)
(*genb411coef[k1sub_,k2sub_,k3sub_]:=Module[{temp},*)
(*temp=b411permdecfn[k1sub,k2sub,k3sub];*)
(*Export[LoopTablesPath<>"B411coefs/B411coefs_"<>ToString[k1sub//N]<>"_"<>ToString[k2sub//N]<>"_"<>ToString[k3sub//N]<>"_.mx",temp]*)
(*];*)


(* ::Input:: *)
(*Monitor[*)
(*Do[genb411coef@@CMASStriEff[[j]],*)
(*{j,1,Length[CMASStriEff]}]*)
(*,j]*)


(* ::Input:: *)
(*Monitor[*)
(*Do[genb411coef@@LOWZtriEff[[j]],*)
(*{j,1,Length[LOWZtriEff]}]*)
(*,j]*)


(* ::Input:: *)
(*Monitor[*)
(*Do[genb411coef@@fisherPoints[[j]],*)
(*{j,1,Length[fisherPoints]}]*)
(*,j]*)

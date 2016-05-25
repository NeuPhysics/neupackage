(* ::Package:: *)

BeginPackage["neumat`"]

imgsizeNeuMat=800;

(*BEGIN Stimulated Matter Effect*)

solNN::usage = "numerical solution of the schrodinger equation for N frequency perturbation";
pltNN::usage = "plot the numerical solution of the schrodinger equation for N frequency perturbation";
solNum::usage = "Numerical Solution to systems with any number of perturbations. \!\(\*SubscriptBox[\(\[Sum]\), \(n\)]\)\!\(\*SubscriptBox[\(A\), \(n\)]\)Sin[\!\(\*SubscriptBox[\(k\), \(n\)]\)x+\!\(\*SubscriptBox[\(\[Phi]\), \(n\)]\)], the parameters for solNum[kNList,aNList,phiNList,thetam,endpoint]";
pltNum::usage = "plot results from solNum; Parameters are given in the order pltNum[kNList,aNList,phiNList,thetam,endpoint,pltLabel,color,frameTicks,frameLabel,startpoint], where startpoint is optional and by defaul it's set to 0.";
hNList::usage = "given a list of n's, generate the Hamiltonian";
solNList::usage = "solution of a system composed of a given list of n's";
pltNList::usage = "plot the result of solNList";

Begin["`Private`"]

solNN[listOfWaveNumber_,listOfAmplitude_,listOfPhase_,thetam_,endpoint_]:=Module[{initM,listOfWaveNumberM,listOfAmplitudeM,listOfPhaseM,thetamM,endpointM,hamilPart1,hamilPart2,hamil,hamilConj,lengthM},

listOfWaveNumberM=listOfWaveNumber;
listOfAmplitudeM=listOfAmplitude;
listOfPhaseM=listOfPhase;
thetamM=thetam;
endpointM=endpoint;
lengthM=Length[listOfWaveNumber];
initM={psi1[0],psi2[0]}=={1,0};


hamilPart1 = Total@MapThread[#1*Sin[#2*x+#3]&,{listOfAmplitudeM,listOfWaveNumberM,listOfPhaseM}];
hamilPart2 = Total@MapThread[#1/#2*Cos[#2*x+#3]&,{listOfAmplitudeM,listOfWaveNumberM,listOfPhaseM}];

hamil = Sin[2thetamM]/2 (hamilPart1)Exp[I(-x-Cos[2thetamM](hamilPart2))];
hamilConj = Sin[2thetamM]/2 (hamilPart1)Exp[-I(-x-Cos[2thetamM](hamilPart2))];

NDSolve[I D[{psi1[x],psi2[x]},x]=={{0,hamil},{hamilConj,0}}.{psi1[x],psi2[x]}&&initM,{psi1,psi2},{x,0,endpointM}]

]


pltNN[listOfWaveNumber_,listOfAmplitude_,listOfPhase_,thetam_,endpoint_,color_]:=Plot[Evaluate[Abs[psi2[x]]^2/.solNN[listOfWaveNumber,listOfAmplitude,listOfPhase,thetam,endpoint][[1]]],{x,0,endpoint},ImageSize->imgsizeNeuMat,Frame->True,FrameLabel->{"x","Transition Probability"},PlotStyle->{Dotted,color},PlotLegends->Placed[Style["Numerical;",color],{Top,Center}]]



(*BEGIN A second method to solve the system numerically*)

(*hamil BEGIN*)
hNList[nNList_,kNList_,aNList_,phiNList_,thetam_,x_]:=Module[{thetamM,BNM,phiNListM,kNListM,aNListM,nNListM,length,PhiNM,gNM},
thetamM=thetam;
nNListM=nNList;
kNListM=kNList;
aNListM=aNList;
length=Length@kNListM;
phiNListM=phiNList;

gNM=nNListM.kNListM ;
BNM=(-I)^(Total[nNListM])Tan[2thetamM] gNM Times@@Table[BesselJ[nNListM[[n]],aNListM[[n]]/kNListM[[n]] Cos[2thetamM]],{n,1,length}];
PhiNM=Exp[I nNListM.phiNListM];


1/2 {{0,-BNM PhiNM Exp[I (gNM-1) x]},{-Conjugate[BNM] Conjugate[PhiNM] Exp[-I (gNM-1) x],0}}
]
(*hamil END*)

(*solNList BEGIN*)
solNList[listInput_,kNList_,aNList_,phiNList_,thetam_,endpoint_]:=Module[
{listInputM,aNListM,kNListM,phiNListM,thetamM,length,endpointM,hamilM,initM},

thetamM=thetam;
listInputM=listInput;
aNListM=aNList;
kNListM=kNList;
phiNListM=phiNList;
length=Length@listInputM;
endpointM=endpoint;
initM={psi1[0],psi2[0]}=={1,0};


hamilM=Total@Table[hNList[listInputM[[i]],kNListM,aNListM,phiNListM,thetamM,x],{i,1,length}];

NDSolve[I D[{psi1[x],psi2[x]},x]==hamilM.{psi1[x],psi2[x]}&&initM,{psi1,psi2},{x,0,endpointM}]

]
(*solNList END*)

pltNList[listInput_,kNList_,aNList_,phiNList_,thetam_,endpoint_,legends_,color_,frameTicks_,frameLabel_,startpoint_:0]:=Plot[Evaluate[Abs[psi2[x]]^2/.solNList[listInput,kNList,aNList,phiNList,thetam,endpoint]],{x,startpoint,endpoint},ImageSize->imgsizeNeuMat,Frame->True,FrameLabel->{{"Transition Probability",None},{"x",frameLabel}},PlotStyle->color,PlotLegends->Placed[Style[ToString[legends],color],{Top,Center}],FrameTicks->{{Automatic,None},{Automatic,frameTicks}},PlotPoints->endpoint]


(*solNum BEGIN*)
solNum[kNList_,aNList_,phiNList_,thetam_,endpoint_]:=Module[
{listInputM,aNListM,kNListM,phiNListM,thetamM,length,endpointM,hamilM,hamil12M,hamil21M,initM,maxstepsValueM},

thetamM=thetam;
aNListM=aNList;
kNListM=kNList;
phiNListM=phiNList;
length=Length@listInputM;
endpointM=endpoint;
initM={psi1[0],psi2[0]}=={1,0};
maxstepsValueM=11000000000;

hamil12M=Sin[2thetamM]/2 (Total@Table[aNListM[[i]]Sin[kNListM[[i]]x+phiNListM[[i]]],{i,1,Length@kNListM}])Exp[I(-x-Cos[2thetam]( 
Total@Table[aNListM[[i]]/kNListM[[i]] Cos[kNListM[[i]]x+phiNListM[[i]]],{i,1,Length@kNListM}]
 ))];
hamil21M=Sin[2thetamM]/2 (Total@Table[aNListM[[i]]Sin[kNListM[[i]]x+phiNListM[[i]]],{i,1,Length@kNListM}])Exp[-I(-x-Cos[2thetam]( 
Total@Table[aNListM[[i]]/kNListM[[i]] Cos[kNListM[[i]]x+phiNListM[[i]]],{i,1,Length@kNListM}]
 ))];

NDSolve[I D[{psi1[x],psi2[x]},x]=={{0,hamil12M},{hamil21M,0}}.{psi1[x],psi2[x]}&&initM,{psi1,psi2},{x,0,endpointM},MaxSteps->maxstepsValueM]

]
(*solNum END*)

(*pltNum BEGIN*)
(*A NOTE HERE: adjust PlotPoints in case the plot shows weird behavior*)
pltNum[kNList_,aNList_,phiNList_,thetam_,endpoint_,pltLabel_,color_,frameTicks_,frameLabel_,startpoint_:0]:=Plot[Evaluate[Abs[psi2[x]]^2/.solNum[kNList,aNList,phiNList,thetam,endpoint][[1]]],{x,startpoint,endpoint},ImageSize->imgsizeNeuMat,Frame->True,FrameLabel->{{"Transition Probability",None},{"x",frameLabel}},PlotStyle->color,PlotLegends->Placed[Style[pltLabel,color],{Top,Center}],FrameTicks->{{Automatic,None},{Automatic,frameTicks}},LabelStyle->Black,FrameTicksStyle->Larger,BaseStyle->{(*FontWeight\[Rule]"Bold",*)FontFamily->"Arial",FontSize->18},ImagePadding->{{Automatic,Automatic},{Automatic,60}}(*,PlotPoints\[Rule]endpoint-startpoint*)]
(*pltNum END*)

(*END A second method to solve the system numerically*)





End[]

(*END Stimulated Matter Effect*)





(*BEGIN Calculate Q values*)
listNGenerator::usage = "generate a list of numbers; Example, listNGeneratorTest2[6,3] gives us {{-3,3},{-3,3},{-3,3},{-3,3},{-3,3},{-3,3}}";
listNGeneratorFormat2::usage = "generate a complete list given the upper limits; As an example listNGeneratorFormat2[{3,3,3}] gives us{{-3,-2,-1,0,1,2,3},{-3,-2,-1,0,1,2,3},{-3,-2,-1,0,1,2,3}}";
bNList::usage = "Calculate the b coefficient of the Hamiltonian; Example, bNList[listOfN,listOfWaveNumber,listOfAmplitude,listOfPhase,thetam]. Output is in unit of \!\(\*SubscriptBox[\(\[Omega]\), \(m\)]\)";
widthNList::usage = "Calculate the width; Example, widthNList[list of n's,list of k's,list of a's,list of \[Phi]'s,thetam] and here lists are in the form {0.999,0.6,0.4}. Output is in unit of \!\(\*SubscriptBox[\(\[Omega]\), \(m\)]\) .";
distanceNList::usage = "calculate the EFFECTIVE distance of a given system from a resonance line, NOTE here is NOT the geometrical distance as we have calculated before in 2 frequency case";
qValue::usage = "calcualte the ratio of the EFFECTIVE distance and Effective Width (bCoefN)";
qValueOrderdList::usage = "gives ordered list of qValues and also their corresponding integer combinations; Example, qValueOrderdList[listNGenerator[6,3],kNListTest,aNListTest,phiNListTest,thetamTest]";

Begin["`Private`"]

listNGenerator[lengthOfList_,rangeOfElement_]:=Module[{},

Table[{-rangeOfElement,rangeOfElement},
{j,1,lengthOfList}
]

]

listNGeneratorFormat2[nRangeList_]:=Module[{lengthM,tableRangeM},
lengthM=Length@nRangeList;

Table[
Table[i,{i,-nRangeList[[j]],nRangeList[[j]]}],
{j,1,lengthM}
]
]



bNList[listOfN_,listOfWaveNumber_,listOfAmplitude_,listOfPhase_,thetam_]:=
-(-I)^(Total[listOfN]) Tan[2thetam]/2 (listOfN.listOfWaveNumber) Apply[Times,MapThread[BesselJ[#1,#3/#2 Cos[2thetam]]&,{listOfN,listOfWaveNumber,listOfAmplitude}]];

distanceNList[listOfN_,listOfWaveNumber0_]:=Abs[listOfN.listOfWaveNumber0-1];

widthNList[listOfN_,listOfWaveNumber_,listOfAmplitude_,listOfPhase_,thetam_]:=Abs[bNList[listOfN,listOfWaveNumber,listOfAmplitude,listOfPhase,thetam]];

qValue[listOfN_,listOfWaveNumber_,listOfAmplitude_,listOfPhase_,thetam_]:= distanceNList[listOfN,listOfWaveNumber]/Abs[bNList[listOfN,listOfWaveNumber,listOfAmplitude,listOfPhase,thetam]];

qValueOrderdList[listOfNRange_,listOfWaveNumber_,listOfAmplitude_,listOfPhase_,thetam_]:=Module[{maptolist,listOfThreadM},

listOfThreadM=Tuples[Table[Table[ni,{ni,listOfNRange[[n,1]],listOfNRange[[n,2]]}],{n,1,Length@listOfNRange}]];

SortBy[
Table[{listN,qValue[listN,listOfWaveNumber,listOfAmplitude,listOfPhase,thetam]},{listN,listOfThreadM}]//Quiet
,Last]

]
End[]
(*END Calculate Q values*)


EndPackage[]

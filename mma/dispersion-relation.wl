(* ::Package:: *)

BeginPackage["dispersion`"]

eta=DiagonalMatrix[{1,-1,-1,-1}];

kno[nIndex_,theta_,phi_]:={1,nIndex Cos[phi]Sin[theta],nIndex Sin[phi]Sin[theta],nIndex Cos[theta]};
kvect[omega_,nIndex_,theta_,phi_]:=omega kno[nIndex,theta,phi];
vvect[theta_,phi_]:={1,Cos[phi]Sin[theta],Sin[phi]Sin[theta],Cos[theta]};

vvMatrix[theta_,phi_]=Outer[Times,v[theta,phi],v[theta,phi]];

vSH[theta_,phi_]:=Sqrt[Pi]{2SH[0,0,theta,phi],Sqrt[2/3](SH[1,-1,theta,phi]-SH[1,1,theta,phi]),I Sqrt[2/3](SH[1,-1,theta,phi]+SH[1,1,theta,phi]),2Sqrt[1/3]SH[1,0,theta,phi]};

vvSHMatrix[theta_,phi_]=Outer[Times,vSH[theta,phi],vSH[theta,phi]];




(* BEGIN for discrete beams + Axial Symmetric (\[Phi] is integrated over 0 to 2\[Pi]) *)

arcCosDivide::usage = "Divide angle range [ArcCos[0.8],ArcCos[-0.8]] into n+1 beams; Takes only one argument n; Example: arcCosDivide[10]";
nMatIntPhi::usage = "Calculate the N matrix; N matrix is defined as: \[Omega] \[CapitalPi]=\[Omega] \[Eta] - \!\(\*
StyleBox[\"N\",\nFontSlant->\"Italic\"]\); Arguments are nMatIntPhi[theta_,nIndex_,thetaK_:0,phiK_:0]";
nMatIntPhi0::usage = "Calculate N matrix in z direction (\[Phi]=\[Theta]=0,\[Phi] is azimuthal angle); Arguments are nMatIntPhi0[theta_,n_]";
nMatFunNBeams::usage = "Calculate N matrix in z direction for discrete beams; Arguments are nMatFunNBeams[gList_,nIndex_,thetaList_]";
dataPltNBeamsPlt::usage = "Plot dispersion relation for discrete beams; Arguments are dataPltNBeamsPlt[gList_,thetaList_,numrange_:{-8.9,8.9},numstep_:0.09,pltRange_:{{-5,5},{-5,5}}]";

arcCosDivide[nDivide_]:=ArcCos@Table[0.8-n 0.8*2/nDivide,{n,0,nDivide}];

nMatIntPhi[theta_,nIndex_,thetaK_:0,phiK_:0]:=Module[{thetaM,phiKM,thetaKM},

thetaM=theta;
phiKM=phiK;
thetaKM=thetaK;

Assuming[nIndex\[Element]Reals(* && (-1<n<0||0<n<1)*),
Integrate[
vvMatrix[thetaM,phi]/(kno[nIndex,thetaKM,phiKM].eta.v[thetaM,phi]),
{phi,0,2Pi}]
]
];

nMatIntPhi0[theta_,n_]=nMatIntPhi[theta,n,0,0];

nMatFunNBeams[gList_,nIndex_,thetaList_]:=Module[{thetaListM,nMatM,lengthM},
thetaListM=thetaList;
lengthM=Length@gList;

nMatM=Total@Table[gList[[i]] nMatIntPhi0[thetaListM[[i]],nIndex]Sin[ thetaListM[[i]] ],{i,1,lengthM}];

nMatM[[1]]=-nMatM[[1]];

nMatM

];

dataPltNBeamsPlt[gList_,thetaList_,numrange_:{-8.9,8.9},numstep_:0.09,pltRange_:{{-5,5},{-5,5}}]:=Module[{eigenDataM,pltDataM},

eigenDataM=Table[
Outer[Times,
Eigenvalues@nMatFunNBeams[gList,nM,thetaList],{nM,1}],
{nM,First@numrange,Last@numrange,numstep}];

pltDataM=eigenDataM//Transpose;

ListPlot[pltDataM,Joined->False,Frame->True,ImageSize->Large,PlotRange->pltRange,FrameLabel->{"k","omega"}]

];

(* END for discrete beams *)



(* BEGIN Two Beams Hyperbola *)


(*** Usage BEGIN  ***)
TwoBeamsAxialSymHyperEqn::usage = "Given the parameters, calculates the ; TwoBeamsAxialSymHyperEqn[omega_,k_,theta1_,theta2_,g1_,g2_]";
TwoBeamsAxialSymHyperFun::usage = "Calculates the function of omega(k); TwoBeamsAxialSymHyperFun[k_,theta1_,theta2_,g1_,g2_]";
TwoBeamsAxialSymPrincipalAxis::usage = "Calculates the principal axis of the hyperbola; TwoBeamsAxialSymPrincipalAxis[k_,theta1_,theta2_,g1_,g2_]";
TwoBeamsAxialSymAsymptotic::usage = "TwoBeamsAxialSymAsymptotic[k_,theta1_,theta2_,g1_,g2_]"
(*** Usage END  ***)

TwoBeamsAxialSymHyperEqn[omega_,k_,theta1_,theta2_,g1_,g2_]:=4(omega-k Cos[theta1])(omega-k Cos[theta2])==g1(1-Cos[theta1]^2)(omega-k Cos[theta2]) + g2(1-Cos[theta2]^2)(omega-k Cos[theta1])

TwoBeamsAxialSymHyperFun[k_,theta1_,theta2_,g1_,g2_]=omega/.Solve[TwoBeamsAxialSymHyperEqn[omega,k,theta1,theta2,g1,g2],omega];

TwoBeamsAxialSymPrincipalAxis[k_,theta1_,theta2_,g1_,g2_]:=Module[{kcM,omegacM,dM,akomegaM,bkM,bomegaM,aomegaM,akM,slopeM},

akomegaM=-2(Cos@theta1+Cos@theta2);
akM=4Cos[theta1]Cos[theta2];
aomegaM=4;
bkM=(g1(1-Cos[theta1]^2)Cos[theta2]+g2 (1-Cos[theta2]^2)Cos[theta1])/2;
bomegaM=-(g1(1-Cos[theta1]^2)+g2(1-Cos[theta2]^2))/2;

dM=-4(Cos@theta1-Cos@theta2)^2;

kcM=-(Det[{ {bkM,akomegaM},{bomegaM,aomegaM} }])/dM;
omegacM=-(Det[{ {akM,bkM},{akomegaM,bomegaM} }])/dM;

slopeM=Tan[(ArcTan[Cos@theta1]+ArcTan[Cos@theta2])/2];

{slopeM(k+kcM)-omegacM, -(1/slopeM)(k+kcM)-omegacM}

];


TwoBeamsAxialSymAsymptotic[k_,theta1_,theta2_,g1_,g2_]:=Module[{kcM,omegacM,dM,deltaM,akomegaM,bkM,bomegaM,aomegaM,akM,cM,paAngleM,lambdaM,aobM},

akomegaM=-2(Cos@theta1+Cos@theta2);
akM=4Cos[theta1]Cos[theta2];
aomegaM=4;
bkM=(g1(1-Cos[theta1]^2)Cos[theta2]+g2 (1-Cos[theta2]^2)Cos[theta1])/2;
bomegaM=-(g1(1-Cos[theta1]^2)+g2(1-Cos[theta2]^2))/2;
cM=0;

dM=-4(Cos@theta1-Cos@theta2)^2;
deltaM = Det[{
{akM,akomegaM,bkM},
{akomegaM,aomegaM,bomegaM},
{bkM,bomegaM,cM}
}];
lambdaM={2(1+Cos[theta1]Cos[theta2])+ 2 Sqrt[(Cos[theta1]Cos[theta2]+1)^2+(Cos[theta1]-Cos[theta2])^2],2(1+Cos[theta1]Cos[theta2])- 2 Sqrt[(Cos[theta1]Cos[theta2]+1)^2+(Cos[theta1]-Cos[theta2])^2]};

paAngleM=(ArcTan[Cos@theta1]+ArcTan[Cos@theta2])/2;
aobM=Sqrt[-lambdaM[[2]]/lambdaM[[1]]];

kcM=-(Det[{ {bkM,akomegaM},{bomegaM,aomegaM} }])/dM;
omegacM=-(Det[{ {akM,bkM},{akomegaM,bomegaM} }])/dM;


{Tan[ArcTan@aobM+paAngleM](k+kcM)-omegacM, Tan[-ArcTan@aobM+paAngleM](k+kcM)-omegacM}


];


TwoBeamsAxialSymPlot[theta1_,theta2_,g1_,g2_,krange_:{-5,5,0.01}]:=Module[{pltDataHyperFunM,paDataM,g1M,g2M,theta1M,theta2M,asymM},
(*g1M=-0.5(2Pi);
g2M=0.5(2Pi);
*)
g1M=g1;
g2M=g2;

theta1M=theta1;
theta2M=theta2;

pltDataHyperFunM=Table[-{k,#}&/@TwoBeamsAxialSymHyperFun[k,theta1M,theta2M,g1M,g2M],
{k,-5,5,0.01}]//Transpose;
paDataM=Table[{k,#}&/@TwoBeamsAxialSymPrincipalAxis[k,theta1M,theta2M,g1M,g2M],{k,krange[[1]],krange[[2]],krange[[3]]}]//Transpose;
asymM=Table[{k,#}&/@TwoBeamsAxialSymAsymptotic[k,theta1M,theta2M,g1M,g2M],{k,krange[[1]],krange[[2]],krange[[3]]}]//Transpose;



Show[
ListPlot[pltDataHyperFunM,Joined->True,Frame->True,FrameLabel->{"k","\[Omega]"}],
ListPlot[paDataM,Joined->True,Frame->True,PlotStyle->LightGray],
ListPlot[asymM,Joined->True,Frame->True,PlotStyle->{Directive[Dashed,LightGray],Directive[Dashed,LightGray]}],
ImageSize->Large]

];

(* END Two Beams Hyperbola *)



(* BEGIN Arbitrary N Beams MAA Solution *)

NBeamsLR::usage = "Generates N Beams with crossing spectrum; NBeamsLR[beams] returns a nested list {{-g,theta1},{g,theta2}}";
NBeamsH::usage = "Generates N Beams with isotropic spectrum; NBeamsH[beams] returns a nested list {{g,theta1},{g,theta2}}";
NBeamsOmegaNPlt::usage="Plot \[Omega](n) for N beams given beams set up; NBeamsOmegaNPlt[nbeams_,nRange_:{-4.5,4.5}], where nbeams should be {{g1,cos@theta1},{g2,cos@theta2},...}";
NBeamsKNPlt::usage ="Plot k(n) for given spectra; NBeamsKNPlt[nbeams_,nRange_:{-4.5,4.5},kRange_:{-10,10}], where nbeams should be {{g1,cos@theta1},{g2,cos@theta2},...}";
NBeamsAxialSymOmegaKEqn::usage = "NBeamsAxialSymOmegaKEqn[omega_,k_,beams_]";
NBeamsAxialSymOmegaKPolyEqn::usage = "NBeamsAxialSymOmegaKPolyEqn[omega_,k_,beams_]";
NBeamsAxialSymOmegaFunPoly::usage = "NBeamsAxialSymOmegaFunPoly[k_,beams_]";
NBeamsOmegaKPlt::usage = "NBeamsOmegaKPlt[nbeams_,kRange_:{-4.5,4.5}]";

NBeamsLR[beams_]:=Module[{hombmM},

hombmM=Transpose[
({(2Pi)/beams,#}&/@Table[ArcCos[0.9+beam (-0.9+0.3)/(beams-1)],{beam,0,beams-1}])
];


hombmM[[1]]=Join[
Table[-1(hombmM)[[1,i]],{i,1,beams/2}],
Table[(hombmM)[[1,i]],{i,beams/2+1,beams}]
];

Transpose@hombmM

];

NBeamsH[beams_]:=Module[{hombmM},

hombmM=Transpose[
({(2Pi)/beams,#}&/@Table[ArcCos[0.9+beam (-0.9+0.3)/(beams-1)],{beam,0,beams-1}])
];


Transpose@hombmM

];





NBeamsOmegaNPlt[nbeams_,nRange_:{-4.5,4.5},kRange_:{-10,10}]:=Module[{beamsomegaKM,omeganM,beamsOmegaKCleanM,beamsOmegaKPltM,nbeamsM,beamsSingularityM,gridlinesM},

nbeamsM=nbeams;

omeganM[n_]:=(1/4)Total[
(#[[1]] (1-#[[2]]^2)/(1-n #[[2]]))&/@nbeamsM
];

gridlinesM=Join[1/#[[2]]&/@nbeams,{{1/(0.6),Purple}}];

beamsOmegaKPltM=Plot[omeganM[n],{n,nRange[[1]],nRange[[2]]},Frame->True,FrameLabel->{"n","\[Omega]"},ImageSize->Large,PlotLabel->"Beams: "<>ToString@Length@nbeamsM,PlotStyle->Gray,PlotRange->{nRange,kRange},Exclusions->gridlinesM];

Show[beamsOmegaKPltM,GridLines->{gridlinesM,None}]

];

NBeamsKNPlt[nbeams_,nRange_:{-4.5,4.5},kRange_:{-10,10}]:=Module[{beamsomegaKM,omeganM,knM,beamsOmegaKCleanM,beamsOmegaKPltM,beamsKNPltM,nbeamsM,beamsSingularityM,gridlinesM},

nbeamsM=nbeams;

omeganM[n_]:=(1/4)Total[
(#[[1]] (1-#[[2]]^2)/(1-n #[[2]]))&/@nbeamsM
];

knM[n_]:= n omeganM[n];

gridlinesM=Join[1/#[[2]]&/@nbeams,{{1/(0.6),Purple}}];

beamsKNPltM=Plot[knM[n],{n,nRange[[1]],nRange[[2]]},Frame->True,FrameLabel->{"n","k"},ImageSize->Large,PlotLabel->"Beams: "<>ToString@Length@nbeamsM,PlotStyle->Gray,PlotRange->{nRange,kRange},Exclusions->gridlinesM];

Show[beamsKNPltM,GridLines->{gridlinesM,None}]

];


NBeamsAxialSymOmegaKEqn[omega_,k_,beams_]:=4== Total[
(#[[1]](1-#[[2]]^2) )/(omega-k #[[2]])&/@beams
];
NBeamsAxialSymOmegaKPolyEqn[omega_,k_,beams_]:=4 Apply[Times,(omega-k #[[2]])&/@beams]== Total[
Table[
(beams[[i,1]](1-beams[[i,2]]^2) )Apply[Times,(omega-k #[[2]])&/@beams]/(omega-k beams[[i,2]]),
{i,1,Length@beams}]
];
NBeamsAxialSymOmegaFunPoly[k_,beams_]:=omega/.{ToRules@Reduce[NBeamsAxialSymOmegaKPolyEqn[omega,k,beams],omega]};

NBeamsOmegaKPlt[nbeams_,kRange_:{-4.5,4.5}]:=Module[{beamsomegaKM,beamsOmegaKCleanM,beamsOmegaKPltM,nbeamsM,beamsSingularityM},

nbeamsM=nbeams;

beamsomegaKM=Table[
{k,#}&/@NBeamsAxialSymOmegaFunPoly[k,nbeamsM],
{k,-kRange[[2]],-kRange[[1]],0.01}];


(*
beamsomegaKM=Table[
-{k,#}&/@NBeamsAxialSymOmegaKPoly[k,NBeamsLR[nbeamsM]],
{k,-4.5,4.5,0.01}];*)


beamsOmegaKCleanM=DeleteCases[beamsomegaKM,x_/;Length@x!=Length@nbeamsM];

Transpose@beamsOmegaKCleanM;

beamsSingularityM=Plot[(Transpose@nbeamsM)[[2]]x,{x,kRange[[1]],kRange[[2]]},PlotStyle->Directive[Orange,Dashing[0.03]]];

{beamsOmegaKPltM=ListPlot[beamsomegaKM,Joined->False,Frame->True,FrameLabel->{"k","\[Omega]"},ImageSize->Large,PlotLabel->ToString@Length@nbeamsM<>" Beams",PlotStyle->Gray,PlotRange->{kRange,Automatic},AspectRatio->1],
Show[beamsOmegaKPltM,beamsSingularityM],
beamsomegaKM
}

]


(* END Arbitrary N Beams MAA Solution *)


(* BEGIN: New Methods to Calculate DR for Discrete Beams which is similar to the continuous code *)

DBAxialSymIntFun0nTotal::usage = "Calculates I0; DBAxialSymIntFun0n[spect], where spect should be of the form {{g1,cos@theta1},{g2,cos@theta2},...}"
DBAxialSymIntFun1nTotal::usage = "Calculates I1; DBAxialSymIntFun1n[spect], where spect should be of the form {{g1,cos@theta1},{g2,cos@theta2},...}"
DBAxialSymIntFun2nTotal::usage = "Calculates I2; DBAxialSymIntFun2n[spect], where spect should be of the form {{g1,cos@theta1},{g2,cos@theta2},...}"

DBAxialSymIntFun0nTotal[spect_,n_]:=Total[
(  ( #[[1]] )/(1- n #[[2]]) )&/@spect
];
DBAxialSymIntFun1nTotal[spect_,n_]:=Total[
( #[[1]] #[[2]] )/(1- n #[[2]])&/@spect
];
DBAxialSymIntFun2nTotal[spect_,n_]:=Total[
( #[[1]] #[[2]]^2 )/(1- n #[[2]])&/@spect
];


DBAxialSymOmegaNMAA[spect_, n_]:=Module[{spectM},

spectM = spect;

( DBAxialSymIntFun0nTotal[spect,n] - DBAxialSymIntFun2nTotal[spect,n]) /4

]


DBAxialSymKNMAA[spect_, n_]:=Module[{spectM},

spectM = spect;

n * ( DBAxialSymIntFun0nTotal[spect,n] - DBAxialSymIntFun2nTotal[spect,n]) /4

]

DBAxialSymOmegaNMZA[spect_, n_]:=Module[{spectM},

spectM = spect;

{ ( DBAxialSymIntFun0nTotal[spect,n] - DBAxialSymIntFun2nTotal[spect,n] + Sqrt[ (DBAxialSymIntFun0nTotal[spect,n] + DBAxialSymIntFun2nTotal[spect,n] -2 DBAxialSymIntFun1nTotal[spect,n])(DBAxialSymIntFun0nTotal[spect,n] + DBAxialSymIntFun2nTotal[spect,n] + 2 DBAxialSymIntFun1nTotal[spect,n]) ])  / (-4)  ,
( DBAxialSymIntFun0nTotal[spect,n] - DBAxialSymIntFun2nTotal[spect,n] - Sqrt[ (DBAxialSymIntFun0nTotal[spect,n] + DBAxialSymIntFun2nTotal[spect,n] -2 DBAxialSymIntFun1nTotal[spect,n])(DBAxialSymIntFun0nTotal[spect,n] + DBAxialSymIntFun2nTotal[spect,n] + 2 DBAxialSymIntFun1nTotal[spect,n]) ])  / (-4) 
 }

]




(* END New Methods to Calculate DR for Discrete Beams which is similar to the continuous code*)



(* BEGIN Calculate instabilities for discrete beams *)

DBAxialSymOmegaNMAAEqnLHSComplex[omega_,k_,spect_]:=Module[{spectM,nM},

spectM = spect;
nM=k/omega;

omega - ( DBAxialSymIntFun0nTotal[spect,nM] - DBAxialSymIntFun2nTotal[spect,nM]) /4

]


DBAxialSymOmegaNMZApEqnLHSComplex[omega_,k_,spect_]:=Module[{spectM,nM},

spectM = spect;

nM = k/omega ;

omega - ( DBAxialSymIntFun0nTotal[spect,nM] - DBAxialSymIntFun2nTotal[spect,nM] + Sqrt[ (DBAxialSymIntFun0nTotal[spect,nM] + DBAxialSymIntFun2nTotal[spect,nM] -2 DBAxialSymIntFun1nTotal[spect,nM])(DBAxialSymIntFun0nTotal[spect,nM] + DBAxialSymIntFun2nTotal[spect,nM] + 2 DBAxialSymIntFun1nTotal[spect,nM]) ])  / (-4) 

]

DBAxialSymOmegaNMZAmEqnLHSComplex[omega_,k_,spect_]:=Module[{spectM,nM},

spectM = spect;

nM = k/omega ;

omega - ( DBAxialSymIntFun0nTotal[spect,nM] - DBAxialSymIntFun2nTotal[spect,nM] - Sqrt[ (DBAxialSymIntFun0nTotal[spect,nM] + DBAxialSymIntFun2nTotal[spect,nM] -2 DBAxialSymIntFun1nTotal[spect,nM])(DBAxialSymIntFun0nTotal[spect,nM] + DBAxialSymIntFun2nTotal[spect,nM] + 2 DBAxialSymIntFun1nTotal[spect,nM]) ])  / (-4) 

]


LSAMAARODataRawDB[spect_,range_:{-4,4,0.01},initk_:0.1*I]:=Module[{},
Table[
	{
		omegareal,k/.FindRoot[
			DBAxialSymOmegaNMAAEqnLHSComplex[omegareal,k,spect],
			{k,initk}
		]
	},
	{omegareal,range[[1]],range[[2]],range[[3]]}
]
];


LSAMZApRODataRawDB[spect_,range_:{-4,4,0.01},initk_:0.1*I]:=Module[{},
Table[
	{
		omegareal,k/.FindRoot[
			DBAxialSymOmegaNMZApEqnLHSComplex[omegareal,k,spect],
			{k,initk}
		]
	},
	{omegareal,range[[1]],range[[2]],range[[3]]}
]
];

LSAMZAmRODataRawDB[spect_,range_:{-4,4,0.01},initk_:0.1*I]:=Module[{},
Table[
	{
		omegareal,k/.FindRoot[
			DBAxialSymOmegaNMZAmEqnLHSComplex[omegareal,k,spect],
			{k,initk}
		]
	},
	{omegareal,range[[1]],range[[2]],range[[3]]}
]
];



LSAMAARODataDB[rawdata_,spect_,thresh_:0.01,thresh2_:0.01]:=Module[{},
{#[[1]],Re@#[[2]],Im@#[[2]]}&/@Select[rawdata,
Abs@DBAxialSymOmegaNMAAEqnLHSComplex[#[[1]],#[[2]],spect]<thresh&&Im@#[[2]]>thresh2&]
]



LSAMZApRODataDB[rawdata_,spect_,thresh_:0.01,thresh2_:0.01]:=Module[{},
{#[[1]],Re@#[[2]],Im@#[[2]]}&/@Select[rawdata,
Abs@DBAxialSymOmegaNMZApEqnLHSComplex[#[[1]],#[[2]],spect]<thresh&&Im@#[[2]]>thresh2&]
]

LSAMZAmRODataDB[rawdata_,spect_,thresh_:0.01,thresh2_:0.01]:=Module[{},
{#[[1]],Re@#[[2]],Im@#[[2]]}&/@Select[rawdata,
Abs@DBAxialSymOmegaNMZAmEqnLHSComplex[#[[1]],#[[2]],spect]<thresh&&Im@#[[2]]>thresh2&]
]




(* The following is for real k *)

LSAMAARKDataRawDB[spect_,range_:{-4,4,0.01},initomega_:0.1*I]:=Module[{omegaMAAFunM,dataM,baddataM},

omegaMAAFunM[k_]:=Last[
omega/.Solve[
DBAxialSymOmegaNMAAEqnLHSComplex[omega,k,spect]==0,
omega
]
];

dataM=Table[
	{
		omegaMAAFunM[kreal],kreal
	},
	{kreal,range[[1]],range[[2]],range[[3]]}
];

baddataM[entry_]:=Not[MatchQ[entry,{_?NumberQ,_?NumberQ}]];

DeleteCases[dataM,_?baddataM]

];

LSAMAARKDataDB[rawdata_,spect_,thresh_:0.01,thresh2_:0.01]:=Module[{},
{Re@#[[1]],#[[2]],Im@#[[1]]}&/@Select[rawdata,Im@#[[1]]>thresh2&]
]




LSAMZApRKDataRawDB[spect_,range_:{-4,4,0.01},initomega_:0.1*I]:=Module[{omegaMZApFunM,dataM,baddataM},

omegaMZApFunM[k_]:=Last[
omega/.NSolve[
DBAxialSymOmegaNMZApEqnLHSComplex[omega,k,spect]==0,
omega
]
];

dataM=Table[
	{
		omegaMZApFunM[kreal],kreal
	},
	{kreal,range[[1]],range[[2]],range[[3]]}
];

baddataM[entry_]:=Not[MatchQ[entry,{_?NumberQ,_?NumberQ}]];

DeleteCases[dataM,_?baddataM]

];

LSAMZApRKDataDB[rawdata_,spect_,thresh_:0.01,thresh2_:0.01]:=Module[{},
{Re@#[[1]],#[[2]],Im@#[[1]]}&/@Select[rawdata,Abs@Im@#[[1]]>thresh2&]
]


LSAMZAmRKDataRawDB[spect_,range_:{-4,4,0.01},initomega_:0.1*I]:=Module[{omegaMZAmFunM,dataM,baddataM},

omegaMZAmFunM[k_]:=Last[
omega/.NSolve[
DBAxialSymOmegaNMZAmEqnLHSComplex[omega,k,spect]==0,
omega
]
];

dataM=Table[
	{
		omegaMZAmFunM[kreal],kreal
	},
	{kreal,range[[1]],range[[2]],range[[3]]}
];

baddataM[entry_]:=Not[MatchQ[entry,{_?NumberQ,_?NumberQ}]];

DeleteCases[dataM,_?baddataM]

];

LSAMZAmRKDataDB[rawdata_,spect_,thresh_:0.01,thresh2_:0.01]:=Module[{},
{Re@#[[1]],#[[2]],Im@#[[1]]}&/@Select[rawdata,Abs@Im@#[[1]]>thresh2&]
]


(* Calculate omega k both complex case *)

LSAMAADataRawDB[spect_,omegarange_:{{-4,4,0.01},{-4,4,0.01}}]:=Module[{kMAAFunM,dataM,baddataM},

kMAAFunM[omega_]:=Last[
k/.Solve[
DBAxialSymOmegaNMAAEqnLHSComplex[omega,k,spect]==0,
k
]
];

dataM=Table[
	{
		omegareal+omegaimag*I,kMAAFunM[omegareal+omegaimag*I]
	},
	{omegareal,omegarange[[1,1]],omegarange[[1,2]],omegarange[[1,3]]},
	{omegaimag,omegarange[[2,1]],omegarange[[2,2]],omegarange[[2,3]]}
];
(*
baddataM[entry_]:=Not[MatchQ[entry,{_?NumberQ,_?NumberQ}]];

DeleteCases[dataM,_?baddataM]*)

dataM

];


LSAMZApDataRawDB[spect_,omegarange_:{{-4,4,0.01},{-4,4,0.01}}]:=Module[{kFunM,dataM,baddataM},

kFunM[omega_]:=Last[
k/.NSolve[
DBAxialSymOmegaNMZApEqnLHSComplex[omega,k,spect]==0,
k
]
];

dataM=Table[
	{
		omegareal+omegaimag*I,kFunM[omegareal+omegaimag*I]
	},
	{omegareal,omegarange[[1,1]],omegarange[[1,2]],omegarange[[1,3]]},
	{omegaimag,omegarange[[2,1]],omegarange[[2,2]],omegarange[[2,3]]}
];
(*
baddataM[entry_]:=Not[MatchQ[entry,{_?NumberQ,_?NumberQ}]];

DeleteCases[dataM,_?baddataM]*)

dataM

];

LSAMZAmDataRawDB[spect_,omegarange_:{{-4,4,0.01},{-4,4,0.01}}]:=Module[{kFunM,dataM,baddataM},

kFunM[omega_]:=Last[
k/.NSolve[
DBAxialSymOmegaNMZAmEqnLHSComplex[omega,k,spect]==0,
k
]
];

dataM=Table[
	{
		omegareal+omegaimag*I,kFunM[omegareal+omegaimag*I]
	},
	{omegareal,omegarange[[1,1]],omegarange[[1,2]],omegarange[[1,3]]},
	{omegaimag,omegarange[[2,1]],omegarange[[2,2]],omegarange[[2,3]]}
];
(*
baddataM[entry_]:=Not[MatchQ[entry,{_?NumberQ,_?NumberQ}]];

DeleteCases[dataM,_?baddataM]*)

dataM

];



(* END Calculate instabilities for discrete beams *)


(* BEGIN Continuous Useful *)


IntFun0n::usage = "Calculates the function \!\(\*SubsuperscriptBox[\(\[Integral]\), \(c1\), \(c2\)]\)\!\(\*FractionBox[\(1\), \(1 - n*u\)]\)\[DifferentialD]u=(Log[(1-ct1 n)/(1-ct2 n)])/n; This function doesn't perform the integral. The analytical expression is hard coded in. IntFun0n[n,c1,c2] takes three parameters, n, c1 which is the starting point of the integral, c2 which is the end point of the integral."; 
IntFun1n::usage = "Calculates the function \!\(\*SubsuperscriptBox[\(\[Integral]\), \(c1\), \(c2\)]\)\!\(\*FractionBox[\(u\), \(1 - n*u\)]\)\[DifferentialD]u = ((ct1-ct2) + Log[(ct1 n-1)/(ct2 n-1)]/n)/n; This function doesn't perform the integral. The analytical expression is hard coded in. IntFun1n[n,c1,c2] takes three parameters, n, c1 which is the starting point of the integral, c2 which is the end point of the integral."; 
IntFun2n::usage = "Calculates the function \!\(\*SubsuperscriptBox[\(\[Integral]\), \(c1\), \(c2\)]\)\!\(\*FractionBox[\(u^2\), \(1 - n*u\)]\)\[DifferentialD]u = (ct1-ct2) ( (ct1+ct2) n +2 )/(2 n^2)+ Log[(ct1 n-1)/(ct2 n-1)]/n^3; This function doesn't perform the integral. The analytical expression is hard coded in. IntFun2n[n,c1,c2] takes three parameters, n, c1 which is the starting point of the integral, c2 which is the end point of the integral."; 
IntFun0::usage = "IntFun0[omega_,k_,ct1_,ct2_]:= IntFun0n[k/omega,ct1,ct2]";
IntFun1::usage = "IntFun1[omega_,k_,ct1_,ct2_]:= IntFun1n[k/omega,ct1,ct2]";
IntFun2::usage = "IntFun2[omega_,k_,ct1_,ct2_]:= IntFun2n[k/omega,ct1,ct2]";
IntFun0nABS::usage = " Calculates the function (Log[Abs[(1-ct1 n)/(1-ct2 n)]])/n; It takes three parameters, n, c1 which is the starting point of the integral, c2 which is the end point of the integral.";
IntFun1nABS::usage = " Calculates the function ((ct1-ct2) + Log[Abs[(ct1 n-1)/(ct2 n-1)]]/n)/n; It takes three parameters, n, c1 which is the starting point of the integral, c2 which is the end point of the integral.";
IntFun2nABS::usage = " Calculates the function  (ct1-ct2) ( (ct1+ct2) n +2 )/(2 n^2)+ Log[Abs[(ct1 n-1)/(ct2 n-1)]]/n^3; It takes three parameters, n, c1 which is the starting point of the integral, c2 which is the end point of the integral.";
ConAxialSymOmegaNMAA::usage = " Calculates \[Omega](n) for MAA solution;  ConAxialSymOmegaNMAA[n_,spect_:{{{0.3,0.6},1},{{0.6,0.9},-1}}] takes up to two parameters: n, spect which is the spectrum used in calculation. By convention, spect should have the form { {{\!\(\*SubscriptBox[\(u\), \(1\)]\),\!\(\*SubscriptBox[\(u\), \(1\)]\)'},\!\(\*SubscriptBox[\(g\), \(1\)]\)}, {{\!\(\*SubscriptBox[\(u\), \(2\)]\),\!\(\*SubscriptBox[\(u\), \(2\)]\)'},\!\(\*SubscriptBox[\(g\), \(2\)]\)}, ... } ";
ConAxialSymOmegaNMAAABS::usage = " Calculates \[Omega](n) for MAA solution, similar to  ConAxialSymOmegaNMAA but with IntFun0nABS, IntFun1nABS, IntFun2nABS.";
ConAxialSymOmegaNMAAEqnLHS::usage = " Returns the real part and imaginary part of \[Omega] - (\!\(\*SubscriptBox[\(I\), \(0\)]\)-\!\(\*SubscriptBox[\(I\), \(2\)]\))/4 in a list {real part, imaginary part}; ConAxialSymOmegaNMAAEqnLHS[omegaR_?NumberQ,omegaI_?NumericQ,kR_?NumberQ,kI_?NumberQ,spect_,wp_:$MachinePrecision];";
ConAxialSymOmegaNMAAEqnLHSComplex::usage = " Returns \[Omega] - (\!\(\*SubscriptBox[\(I\), \(0\)]\)-\!\(\*SubscriptBox[\(I\), \(2\)]\))/4; ConAxialSymOmegaNMAAEqnLHSComplex[omega_?NumericQ,k_?NumberQ,spect_]; This calculation uses the IntFun0n, IntFun2n functions.";
ConAxialSymOmegaNMAAEqnLHSComplexN::usage = " Returns \[Omega] - (\!\(\*SubscriptBox[\(I\), \(0\)]\)-\!\(\*SubscriptBox[\(I\), \(2\)]\))/4; ConAxialSymOmegaNMAAEqnLHSComplexN[omega_?NumericQ,k_?NumberQ,spect_,wp_:$MachinePrecision]; This calculation uses NIntegral to calculate the integrals directly.";
ConAxialSymOmegaNMZA::usage = " Returns \[Omega](n) for MZA solutions in a list, {MZA+, MZA-}; ConAxialSymOmegaNMZA[n_,spect_:{{{0.3,0.6},1},{{0.6,0.9},-1}}] ";
ConAxialSymOmegaNMZASQRT::usage = " Returns the terms unders square root in MZA solutions, (\!\(\*SubscriptBox[\(I\), \(0\)]\)+\!\(\*SubscriptBox[\(I\), \(2\)]\)-2\!\(\*SubscriptBox[\(I\), \(1\)]\))(\!\(\*SubscriptBox[\(I\), \(0\)]\)+\!\(\*SubscriptBox[\(I\), \(2\)]\)+2\!\(\*SubscriptBox[\(I\), \(1\)]\)) ; ConAxialSymOmegaNMZASQRT[n_,spect_:{{{0.3,0.6},1},{{0.6,0.9},-1}}] ";
ConAxialSymOmegaNMZAABS::usage = " Returns \[Omega](n) for MZA solutions in a list, {MZA+, MZA-}, but with IntFun0nABS,IntFun1nABS,IntFun2nABS; ConAxialSymOmegaNMZA[n_,spect_:{{{0.3,0.6},1},{{0.6,0.9},-1}}] ";
ConAxialSymOmegaNMZApEqnLHS::usage = " Returns real part and imaginary part of \[Omega] - (\!\(\*SubscriptBox[\(I\), \(0\)]\)-\!\(\*SubscriptBox[\(I\), \(2\)]\)+\!\(\*SqrtBox[\(\((\*SubscriptBox[\(I\), \(0\)] + \*SubscriptBox[\(I\), \(2\)] - 2 \*SubscriptBox[\(I\), \(1\)])\) \((\*SubscriptBox[\(I\), \(0\)] + \*SubscriptBox[\(I\), \(2\)] + 2 \*SubscriptBox[\(I\), \(1\)])\)\)]\))/(-4), which is related to MZA+ solution, in a list {real part, imaginary part} ;  ConAxialSymOmegaNMZApEqnLHS[omegaR_?NumberQ,omegaI_?NumericQ,kR_?NumberQ,kI_?NumberQ,spect_]";
ConAxialSymOmegaNMZAmEqnLHS::usage = " Returns real part and imaginary part of  \[Omega] - (\!\(\*SubscriptBox[\(I\), \(0\)]\)-\!\(\*SubscriptBox[\(I\), \(2\)]\)-\!\(\*SqrtBox[\(\((\*SubscriptBox[\(I\), \(0\)] + \*SubscriptBox[\(I\), \(2\)] - 2 \*SubscriptBox[\(I\), \(1\)])\) \((\*SubscriptBox[\(I\), \(0\)] + \*SubscriptBox[\(I\), \(2\)] + 2 \*SubscriptBox[\(I\), \(1\)])\)\)]\))/(-4), which is related to MZA- solution, in a list {real part, imaginary part} ;  ConAxialSymOmegaNMZAmEqnLHS[omegaR_?NumberQ,omegaI_?NumericQ,kR_?NumberQ,kI_?NumberQ,spect_]";
ConAxialSymOmegaNMZApEqnLHSComplex::usage = " Returns the value of \[Omega] - (\!\(\*SubscriptBox[\(I\), \(0\)]\)-\!\(\*SubscriptBox[\(I\), \(2\)]\)+\!\(\*SqrtBox[\(\((\*SubscriptBox[\(I\), \(0\)] + \*SubscriptBox[\(I\), \(2\)] - 2 \*SubscriptBox[\(I\), \(1\)])\) \((\*SubscriptBox[\(I\), \(0\)] + \*SubscriptBox[\(I\), \(2\)] + 2 \*SubscriptBox[\(I\), \(1\)])\)\)]\))/(-4) as a complex number; ConAxialSymOmegaNMZApEqnLHSComplex[omega_?NumericQ,k_?NumberQ,spect_]";
ConAxialSymOmegaNMZAmEqnLHSComplex::usage = " Returns the value of \[Omega] - (\!\(\*SubscriptBox[\(I\), \(0\)]\)-\!\(\*SubscriptBox[\(I\), \(2\)]\)-\!\(\*SqrtBox[\(\((\*SubscriptBox[\(I\), \(0\)] + \*SubscriptBox[\(I\), \(2\)] - 2 \*SubscriptBox[\(I\), \(1\)])\) \((\*SubscriptBox[\(I\), \(0\)] + \*SubscriptBox[\(I\), \(2\)] + 2 \*SubscriptBox[\(I\), \(1\)])\)\)]\))/(-4) as a complex number; ConAxialSymOmegaNMZAmEqnLHSComplex[omega_?NumericQ,k_?NumberQ,spect_]";



IntFun0n[n_,ct1_,ct2_]:= (Log[(1-ct1 n)/(1-ct2 n)])/n;
IntFun1n[n_,ct1_,ct2_]:= ((ct1-ct2) + Log[(ct1 n-1)/(ct2 n-1)]/n)/n;
IntFun2n[n_,ct1_,ct2_]:= (ct1-ct2) ( (ct1+ct2) n +2 )/(2 n^2)+ Log[(ct1 n-1)/(ct2 n-1)]/n^3;

IntFun0[omega_,k_,ct1_,ct2_]:= IntFun0n[k/omega,ct1,ct2];
IntFun1[omega_,k_,ct1_,ct2_]:= IntFun1n[k/omega,ct1,ct2];
IntFun2[omega_,k_,ct1_,ct2_]:= IntFun2n[k/omega,ct1,ct2];

IntFun0nABS[n_,ct1_,ct2_]:= (Log[Abs[(1-ct1 n)/(1-ct2 n)]])/n;
IntFun1nABS[n_,ct1_,ct2_]:= ((ct1-ct2) + Log[Abs[(ct1 n-1)/(ct2 n-1)]]/n)/n;
IntFun2nABS[n_,ct1_,ct2_]:=  (ct1-ct2) ( (ct1+ct2) n +2 )/(2 n^2)+ Log[Abs[(ct1 n-1)/(ct2 n-1)]]/n^3;

(* END Continuous Useful *)


ConAxialSymOmegaNMAA[n_,spect_:{{{0.3,0.6},1},{{0.6,0.9},-1}}]:=Module[{spectM},

spectM=spect;
(*{{{0.9,0.6},-1},
{{0.6,0.3},1}};
*)

Total[
#[[2]](IntFun0n[n,#[[1,1]],#[[1,2]]]-IntFun2n[n,#[[1,1]],#[[1,2]]])/(4)&/@spectM
]

]


ConAxialSymOmegaNMAAABS[n_,spect_:{{{0.3,0.6},1},{{0.6,0.9},-1}}]:=Module[{spectM},

spectM=spect;
(*{{{0.9,0.6},-1},
{{0.6,0.3},1}};
*)

Total[
#[[2]](IntFun0nABS[n,#[[1,1]],#[[1,2]]]-IntFun2nABS[n,#[[1,1]],#[[1,2]]])/(4)&/@spectM
]

]



ConAxialSymOmegaNMAAEqnLHS[omegaR_?NumberQ,omegaI_?NumericQ,kR_?NumberQ,kI_?NumberQ,spect_,wp_:$MachinePrecision]:=Module[{eqnLHSM,spectM,IntFun0nFFM,IntFun1nFFM,IntFun2nFFM},


spectM=spect;

IntFun0nFFM=Total[
	NIntegrate[(#[[2]])/(1-(kR+I*kI) u/(omegaR+I*omegaI)),{u,#[[1,1]],#[[1,2]]},WorkingPrecision->wp]&/@spectM
];
IntFun2nFFM=Total[
	NIntegrate[(#[[2]]u^2)/(1-(kR+I*kI) u/(omegaR+I*omegaI)),{u,#[[1,1]],#[[1,2]]},WorkingPrecision->wp]&/@spectM
];


eqnLHSM=omegaR+I omegaI -(IntFun0nFFM-IntFun2nFFM)/(4);


{ComplexExpand[Re[
eqnLHSM
]],
ComplexExpand[Im[
eqnLHSM
]]
}

]

ConAxialSymOmegaNMAAEqnLHSComplex[omega_?NumericQ,k_?NumberQ,spect_]:=Module[{eqnLHSM,spectM,IntFun0nFFM,IntFun1nFFM,IntFun2nFFM},


spectM=spect;

(*
IntFun0nFFM=Total[
	NIntegrate[(#[[2]])/(1-k u/omega),{u,#[[1,1]],#[[1,2]]},WorkingPrecision->wp]&/@spectM
];
IntFun2nFFM=Total[
	NIntegrate[(#[[2]]u^2)/(1- k u/omega),{u,#[[1,1]],#[[1,2]]},WorkingPrecision->wp]&/@spectM
];*)


IntFun0nFFM=Total[#[[2]]* IntFun0n[k/omega,#[[1,1]],#[[1,2]]]&/@spectM];
IntFun2nFFM = Total[#[[2]]* IntFun2n[k/omega,#[[1,1]],#[[1,2]]]&/@spectM];


eqnLHSM = omega -(IntFun0nFFM-IntFun2nFFM)/(4);


(*eqnLHSM = omega -(Total[#[[2]]* IntFun0n[k/omega,#[[1,1]],#[[1,2]]]&/@spectM]-Total[#[[2]]* IntFun2n[k/omega,#[[1,1]],#[[1,2]]]&/@spectM])/(4);*)

eqnLHSM

]

ConAxialSymOmegaNMAAEqnLHSComplexN[omega_?NumericQ,k_?NumberQ,spect_,wp_:$MachinePrecision]:=Module[{eqnLHSM,spectM,IntFun0nFFM,IntFun1nFFM,IntFun2nFFM},


spectM=spect;


IntFun0nFFM=Total[
	NIntegrate[(#[[2]])/(1-k u/omega),{u,#[[1,1]],#[[1,2]]},WorkingPrecision->wp]&/@spectM
];
IntFun2nFFM=Total[
	NIntegrate[(#[[2]]u^2)/(1- k u/omega),{u,#[[1,1]],#[[1,2]]},WorkingPrecision->wp]&/@spectM
];


(*IntFun0nFFM=Total[#[[2]]* IntFun0n[k/omega,#[[1,1]],#[[1,2]]]&/@spectM];
IntFun2nFFM=Total[#[[2]]* IntFun2n[k/omega,#[[1,1]],#[[1,2]]]&/@spectM];
*)

eqnLHSM=omega -(IntFun0nFFM-IntFun2nFFM)/(4);


eqnLHSM

]


(* BEGIN Continuous MAA and MZA Solution: NO Feed in Spectrum *)

ConAxialSymOmegaNMZA[n_,spect_:{{{0.3,0.6},1},{{0.6,0.9},-1}}]:=Module[{spectM,i1mzaM,i2mzaM,i0mzaM},

spectM=spect;
(*{{{0.9,0.3},1}};*)
(*{{{0.9,0.6},-1},{{0.6,0.3},1}}; *)
(*spect;*)


i0mzaM=Total[#[[2]] IntFun0n[n,#[[1,1]],#[[1,2]]]&/@spectM];
i2mzaM=Total[#[[2]] IntFun2n[n,#[[1,1]],#[[1,2]]]&/@spectM];
i1mzaM=Total[#[[2]] IntFun1n[n,#[[1,1]],#[[1,2]]]&/@spectM];


{(i0mzaM-i2mzaM+Sqrt[(i0mzaM+i2mzaM-2i1mzaM)(i0mzaM+i2mzaM+2i1mzaM)])/(-4),(i0mzaM-i2mzaM-Sqrt[(i0mzaM+i2mzaM-2i1mzaM)(i0mzaM+i2mzaM+2i1mzaM)])/(-4)}
(*{(-(i0mzaM)(i2mzaM)+ i1mzaM^2+8(i0mzaM+i2mzaM))/(-4)}*)


]

ConAxialSymOmegaNMZASQRT[n_,spect_:{{{0.3,0.6},1},{{0.6,0.9},-1}}]:=Module[{spectM,i1mzaM,i2mzaM,i0mzaM},

spectM=spect;
(*{{{0.9,0.3},1}};*)
(*{{{0.9,0.6},-1},{{0.6,0.3},1}}; *)
(*spect;*)


i0mzaM=Total[#[[2]] IntFun0n[n,#[[1,1]],#[[1,2]]]&/@spectM];
i2mzaM=Total[#[[2]] IntFun2n[n,#[[1,1]],#[[1,2]]]&/@spectM];
i1mzaM=Total[#[[2]] IntFun1n[n,#[[1,1]],#[[1,2]]]&/@spectM];


(i0mzaM+i2mzaM-2i1mzaM)(i0mzaM+i2mzaM+2i1mzaM)
(*{(-(i0mzaM)(i2mzaM)+ i1mzaM^2+8(i0mzaM+i2mzaM))/(-4)}*)


]



ConAxialSymOmegaNMZAABS[n_,spect_:{{{0.3,0.6},1},{{0.6,0.9},-1}}]:=Module[{spectM,i1mzaM,i2mzaM,i0mzaM},

spectM=spect;
(*{{{0.9,0.3},1}};*)
(*{{{0.9,0.6},-1},{{0.6,0.3},1}}; *)
(*spect;*)


i0mzaM=Total[#[[2]] IntFun0nABS[n,#[[1,1]],#[[1,2]]]&/@spectM];
i2mzaM=Total[#[[2]] IntFun2nABS[n,#[[1,1]],#[[1,2]]]&/@spectM];
i1mzaM=Total[#[[2]] IntFun1nABS[n,#[[1,1]],#[[1,2]]]&/@spectM];


{(i0mzaM-i2mzaM+Sqrt[(i0mzaM+i2mzaM-2i1mzaM)(i0mzaM+i2mzaM+2i1mzaM)])/(-4),(i0mzaM-i2mzaM-Sqrt[(i0mzaM+i2mzaM-2i1mzaM)(i0mzaM+i2mzaM+2i1mzaM)])/(-4)}
(*{(-(i0mzaM)(i2mzaM)+ i1mzaM^2+8(i0mzaM+i2mzaM))/(-4)}*)


]


ConAxialSymOmegaNMZApEqnLHS[omegaR_?NumberQ,omegaI_?NumericQ,kR_?NumberQ,kI_?NumberQ,spect_]:=Module[{eqnLHSM,spectM,IntFun0nFFM,IntFun1nFFM,IntFun2nFFM},


spectM=spect;



IntFun0nFFM = Total[#[[2]] IntFun0n[(kR+I*kI)/(omegaR+I omegaI),#[[1,1]],#[[1,2]]]&/@spectM];
IntFun1nFFM = Total[#[[2]] IntFun2n[(kR+I*kI)/(omegaR+I omegaI),#[[1,1]],#[[1,2]]]&/@spectM];
IntFun2nFFM = Total[#[[2]] IntFun1n[(kR+I*kI)/(omegaR+I omegaI),#[[1,1]],#[[1,2]]]&/@spectM];

eqnLHSM=omegaR+I omegaI-(IntFun0nFFM-IntFun2nFFM+\[Sqrt]((IntFun0nFFM+IntFun2nFFM-2IntFun1nFFM)(IntFun0nFFM+IntFun2nFFM+2IntFun1nFFM)))/(-4);


(* NIntegrate method *)
(*IntFun0nFFM=Total[
	NIntegrate[(#[[2]])/(1-(kR+I*kI) u/(omegaR+I*omegaI)),{u,#[[1,1]],#[[1,2]]}]&/@spectM
];
IntFun1nFFM=Total[
	NIntegrate[(#[[2]]u)/(1-(kR+I*kI) u/(omegaR+I*omegaI)),{u,#[[1,1]],#[[1,2]]}]&/@spectM
];
IntFun2nFFM=Total[
	NIntegrate[(#[[2]]u^2)/(1-(kR+I*kI) u/(omegaR+I*omegaI)),{u,#[[1,1]],#[[1,2]]}]&/@spectM
];


eqnLHSM=omegaR+I omegaI-(IntFun0nFFM-IntFun2nFFM+\[Sqrt]((IntFun0nFFM+IntFun2nFFM-2IntFun1nFFM)(IntFun0nFFM+IntFun2nFFM+2IntFun1nFFM)))/(-4);
*)


{ComplexExpand[Re[
eqnLHSM
]],
ComplexExpand[Im[
eqnLHSM
]]
}


]




ConAxialSymOmegaNMZApEqnLHSComplex[omega_?NumericQ,k_?NumberQ,spect_]:=Module[{eqnLHSM,spectM,IntFun0nFFM,IntFun1nFFM,IntFun2nFFM},


spectM=spect;


IntFun0nFFM = Total[#[[2]] IntFun0n[k/omega,#[[1,1]],#[[1,2]]]&/@spectM];
IntFun1nFFM = Total[#[[2]] IntFun2n[k/omega,#[[1,1]],#[[1,2]]]&/@spectM];
IntFun2nFFM = Total[#[[2]] IntFun1n[k/omega,#[[1,1]],#[[1,2]]]&/@spectM];

eqnLHSM=omega-(IntFun0nFFM-IntFun2nFFM+\[Sqrt]((IntFun0nFFM+IntFun2nFFM-2IntFun1nFFM)(IntFun0nFFM+IntFun2nFFM+2IntFun1nFFM)))/(-4);


(* NIntegrate method *)
(*IntFun0nFFM=Total[
	NIntegrate[(#[[2]])/(1-(kR+I*kI) u/(omegaR+I*omegaI)),{u,#[[1,1]],#[[1,2]]}]&/@spectM
];
IntFun1nFFM=Total[
	NIntegrate[(#[[2]]u)/(1-(kR+I*kI) u/(omegaR+I*omegaI)),{u,#[[1,1]],#[[1,2]]}]&/@spectM
];
IntFun2nFFM=Total[
	NIntegrate[(#[[2]]u^2)/(1-(kR+I*kI) u/(omegaR+I*omegaI)),{u,#[[1,1]],#[[1,2]]}]&/@spectM
];


eqnLHSM=omegaR+I omegaI-(IntFun0nFFM-IntFun2nFFM+\[Sqrt]((IntFun0nFFM+IntFun2nFFM-2IntFun1nFFM)(IntFun0nFFM+IntFun2nFFM+2IntFun1nFFM)))/(-4);
*)



eqnLHSM



]





ConAxialSymOmegaNMZApEqnLHSComplexN[omega_?NumericQ,k_?NumberQ,spect_,wp_:$MachinePrecision]:=Module[{eqnLHSM,spectM,IntFun0nFFM,IntFun1nFFM,IntFun2nFFM},


spectM=spect;

(*
IntFun0nFFM = Total[#[[2]] IntFun0n[k/omega,#[[1,1]],#[[1,2]]]&/@spectM];
IntFun1nFFM = Total[#[[2]] IntFun2n[k/omega,#[[1,1]],#[[1,2]]]&/@spectM];
IntFun2nFFM = Total[#[[2]] IntFun1n[k/omega,#[[1,1]],#[[1,2]]]&/@spectM];

eqnLHSM=omega-(IntFun0nFFM-IntFun2nFFM+\[Sqrt]((IntFun0nFFM+IntFun2nFFM-2IntFun1nFFM)(IntFun0nFFM+IntFun2nFFM+2IntFun1nFFM)))/(-4);
*)

(* NIntegrate method *)
IntFun0nFFM[omegaF_?NumericQ,kF_?NumericQ]:=Total[
	NIntegrate[(#[[2]])/(1-(kF) u/(omegaF)),{u,#[[1,1]],#[[1,2]]},WorkingPrecision->wp]&/@spectM
];
IntFun1nFFM[omegaF_?NumericQ,kF_?NumericQ]:=Total[
	NIntegrate[(#[[2]]u)/(1-(kF) u/(omegaF)),{u,#[[1,1]],#[[1,2]]},WorkingPrecision->wp]&/@spectM
];
IntFun2nFFM[omegaF_?NumericQ,kF_?NumericQ]:=Total[
	NIntegrate[(#[[2]]u^2)/(1-(kF) u/(omegaF)),{u,#[[1,1]],#[[1,2]]},WorkingPrecision->wp]&/@spectM
];


eqnLHSM=omega -(IntFun0nFFM[omega,k]-IntFun2nFFM[omega,k]+\[Sqrt]((IntFun0nFFM[omega,k]+IntFun2nFFM[omega,k]-2IntFun1nFFM[omega,k])(IntFun0nFFM[omega,k]+IntFun2nFFM[omega,k]+2IntFun1nFFM[omega,k])))/(-4);




eqnLHSM



]


ConAxialSymOmegaNMZAmEqnLHS[omegaR_?NumberQ,omegaI_?NumericQ,kR_?NumberQ,kI_?NumberQ,spect_]:=Module[{eqnLHSM,spectM,IntFun0nFFM,IntFun1nFFM,IntFun2nFFM,wpM},


spectM=spect;
wpM=30;



IntFun0nFFM = Total[#[[2]] IntFun0n[(kR+I*kI)/(omegaR+I omegaI),#[[1,1]],#[[1,2]]]&/@spectM];
IntFun1nFFM = Total[#[[2]] IntFun2n[(kR+I*kI)/(omegaR+I omegaI),#[[1,1]],#[[1,2]]]&/@spectM];
IntFun2nFFM = Total[#[[2]] IntFun1n[(kR+I*kI)/(omegaR+I omegaI),#[[1,1]],#[[1,2]]]&/@spectM];

eqnLHSM=omegaR+I omegaI-(IntFun0nFFM-IntFun2nFFM+\[Sqrt]((IntFun0nFFM+IntFun2nFFM-2IntFun1nFFM)(IntFun0nFFM+IntFun2nFFM+2IntFun1nFFM)))/(-4);


(*IntFun0nFFM=Total[
	NIntegrate[(#[[2]])/(1-(kR+I*kI) u/(omegaR+I*omegaI)),{u,#[[1,1]],#[[1,2]]},WorkingPrecision\[Rule]wpM]&/@spectM
];
IntFun1nFFM=Total[
	NIntegrate[(#[[2]]u)/(1-(kR+I*kI) u/(omegaR+I*omegaI)),{u,#[[1,1]],#[[1,2]]},WorkingPrecision\[Rule]wpM]&/@spectM
];
IntFun2nFFM=Total[
	NIntegrate[(#[[2]]u^2)/(1-(kR+I*kI) u/(omegaR+I*omegaI)),{u,#[[1,1]],#[[1,2]]},WorkingPrecision\[Rule]wpM]&/@spectM
];


eqnLHSM=omegaR+I omegaI-(IntFun0nFFM-IntFun2nFFM - \[Sqrt]((IntFun0nFFM+IntFun2nFFM-2IntFun1nFFM)(IntFun0nFFM+IntFun2nFFM+2IntFun1nFFM)))/(-4);
*)


{ComplexExpand[Re[
eqnLHSM
]],
ComplexExpand[Im[
eqnLHSM
]]
}


]






ConAxialSymOmegaNMZAmEqnLHSComplex[omega_?NumericQ,k_?NumberQ,spect_]:=Module[{eqnLHSM,spectM,IntFun0nFFM,IntFun1nFFM,IntFun2nFFM},


spectM=spect;


IntFun0nFFM = Total[#[[2]] IntFun0n[k/omega,#[[1,1]],#[[1,2]]]&/@spectM];
IntFun1nFFM = Total[#[[2]] IntFun2n[k/omega,#[[1,1]],#[[1,2]]]&/@spectM];
IntFun2nFFM = Total[#[[2]] IntFun1n[k/omega,#[[1,1]],#[[1,2]]]&/@spectM];

eqnLHSM=omega-(IntFun0nFFM-IntFun2nFFM - \[Sqrt]((IntFun0nFFM+IntFun2nFFM-2IntFun1nFFM)(IntFun0nFFM+IntFun2nFFM+2IntFun1nFFM)))/(-4);


(* NIntegrate method *)
(*IntFun0nFFM=Total[
	NIntegrate[(#[[2]])/(1-(kR+I*kI) u/(omegaR+I*omegaI)),{u,#[[1,1]],#[[1,2]]}]&/@spectM
];
IntFun1nFFM=Total[
	NIntegrate[(#[[2]]u)/(1-(kR+I*kI) u/(omegaR+I*omegaI)),{u,#[[1,1]],#[[1,2]]}]&/@spectM
];
IntFun2nFFM=Total[
	NIntegrate[(#[[2]]u^2)/(1-(kR+I*kI) u/(omegaR+I*omegaI)),{u,#[[1,1]],#[[1,2]]}]&/@spectM
];


eqnLHSM=omegaR+I omegaI-(IntFun0nFFM-IntFun2nFFM+\[Sqrt]((IntFun0nFFM+IntFun2nFFM-2IntFun1nFFM)(IntFun0nFFM+IntFun2nFFM+2IntFun1nFFM)))/(-4);
*)



eqnLHSM



]


ConAxialSymOmegaNMZAmEqnLHSComplexABS[omega_?NumericQ,k_?NumberQ,spect_]:=Module[{eqnLHSM,spectM,IntFun0nFFM,IntFun1nFFM,IntFun2nFFM},


spectM=spect;


IntFun0nFFM = Total[#[[2]] IntFun0nABS[k/omega,#[[1,1]],#[[1,2]]]&/@spectM];
IntFun1nFFM = Total[#[[2]] IntFun2n[k/omega,#[[1,1]],#[[1,2]]]&/@spectM];
IntFun2nFFM = Total[#[[2]] IntFun1n[k/omega,#[[1,1]],#[[1,2]]]&/@spectM];

eqnLHSM=omega-(IntFun0nFFM-IntFun2nFFM - \[Sqrt]((IntFun0nFFM+IntFun2nFFM-2IntFun1nFFM)(IntFun0nFFM+IntFun2nFFM+2IntFun1nFFM)))/(-4);


(* NIntegrate method *)
(*IntFun0nFFM=Total[
	NIntegrate[(#[[2]])/(1-(kR+I*kI) u/(omegaR+I*omegaI)),{u,#[[1,1]],#[[1,2]]}]&/@spectM
];
IntFun1nFFM=Total[
	NIntegrate[(#[[2]]u)/(1-(kR+I*kI) u/(omegaR+I*omegaI)),{u,#[[1,1]],#[[1,2]]}]&/@spectM
];
IntFun2nFFM=Total[
	NIntegrate[(#[[2]]u^2)/(1-(kR+I*kI) u/(omegaR+I*omegaI)),{u,#[[1,1]],#[[1,2]]}]&/@spectM
];


eqnLHSM=omegaR+I omegaI-(IntFun0nFFM-IntFun2nFFM+\[Sqrt]((IntFun0nFFM+IntFun2nFFM-2IntFun1nFFM)(IntFun0nFFM+IntFun2nFFM+2IntFun1nFFM)))/(-4);
*)



eqnLHSM



]


(* END Continuous MAA and MZA Solution: Feed in Spectrum  *)


(* BEGIN Continuous Feed in Function as spectrum *)

(* Define the integral functions to be used in dispersion relation *)
IntFun0nFSpect[n_,ct1_,ct2_,spectrumFun_]:= NIntegrate[(spectrumFun@u)/(1-n u),{u,ct1,ct2},Exclusions->{1/n}];
IntFun1nFSpect[n_,ct1_,ct2_,spectrumFun_]:= NIntegrate[(spectrumFun@u)u/(1-n u),{u,ct1,ct2},Exclusions->{1/n}];
IntFun2nFSpect[n_,ct1_,ct2_,spectrumFun_]:= NIntegrate[(spectrumFun@u)u^2/(1-n u),{u,ct1,ct2},Exclusions->{1/n}];

SpectrumOmegaNMAA[n_?NumberQ,ct1_?NumberQ,ct2_?NumberQ,spectrumFun_]:=(IntFun0nFSpect[n,ct1,ct2,spectrumFun]-IntFun2nFSpect[n,ct1,ct2,spectrumFun])/(4);

SpectrumOmegaNMZA[n_?NumberQ,ct1_?NumberQ,ct2_?NumberQ,spectrumFun_]:={(IntFun0nFSpect[n,ct1,ct2,spectrumFun]-IntFun2nFSpect[n,ct1,ct2,spectrumFun]+\[Sqrt]((IntFun0nFSpect[n,ct1,ct2,spectrumFun]+IntFun2nFSpect[n,ct1,ct2,spectrumFun]-2IntFun1nFSpect[n,ct1,ct2,spectrumFun])(IntFun0nFSpect[n,ct1,ct2,spectrumFun]+IntFun2nFSpect[n,ct1,ct2,spectrumFun]+2IntFun1nFSpect[n,ct1,ct2,spectrumFun])))/(-4),
(IntFun0nFSpect[n,ct1,ct2,spectrumFun]-IntFun2nFSpect[n,ct1,ct2,spectrumFun]-\[Sqrt]((IntFun0nFSpect[n,ct1,ct2,spectrumFun]+IntFun2nFSpect[n,ct1,ct2,spectrumFun]-2IntFun1nFSpect[n,ct1,ct2,spectrumFun])(IntFun0nFSpect[n,ct1,ct2,spectrumFun]+IntFun2nFSpect[n,ct1,ct2,spectrumFun]+2IntFun1nFSpect[n,ct1,ct2,spectrumFun])))/(-4)
};

(* Numerically Find the Instabilities *)

SpectrumOmegaKMAAEqnLHS[omegaR_?NumberQ,omegaI_?NumericQ,kR_?NumberQ,kI_?NumberQ,spect_]:=Module[{eqnLHSM,spectM,IntFun0nFFM,IntFun1nFFM,IntFun2nFFM},


spectM=spect;

IntFun0nFFM=NIntegrate[(spectM@u)/(1-(kR+I*kI) u/(omegaR+I*omegaI)),{u,-1,1}];
IntFun2nFFM=NIntegrate[(spectM[u]u^2)/(1-(kR+I*kI) u/(omegaR+I*omegaI)),{u,-1,1}];


eqnLHSM[omega_,k_]:=omega-(IntFun0nFFM-IntFun2nFFM)/(4);


{ComplexExpand[Re[
eqnLHSM[omegaR+I omegaI,kR+I kI]
]],
ComplexExpand[Im[
eqnLHSM[omegaR+I omegaI,kR+I kI]
]]
}


]


SpectrumOmegaKMZApEqnLHS[omegaR_?NumberQ,omegaI_?NumericQ,kR_?NumberQ,kI_?NumberQ,spect_]:=Module[{eqnLHSM,spectM,IntFun0nFFM,IntFun1nFFM,IntFun2nFFM},


spectM=spect;

IntFun0nFFM=NIntegrate[(spectM@u)/(1-(kR+I*kI) u/(omegaR+I*omegaI)),{u,-1,1}];
IntFun1nFFM=NIntegrate[(spectM@u)u/(1-(kR+I*kI) u/(omegaR+I*omegaI)),{u,-1,1}];
IntFun2nFFM=NIntegrate[(spectM[u]u^2)/(1-(kR+I*kI) u/(omegaR+I*omegaI)),{u,-1,1}];


eqnLHSM[omega_,k_]:=omega-(IntFun0nFFM-IntFun2nFFM+\[Sqrt]((IntFun0nFFM+IntFun2nFFM-2IntFun1nFFM)(IntFun0nFFM+IntFun2nFFM+2IntFun1nFFM)))/(-4);



{ComplexExpand[Re[
eqnLHSM[omegaR+I omegaI,kR+I kI]
]],
ComplexExpand[Im[
eqnLHSM[omegaR+I omegaI,kR+I kI]
]]
}


]


SpectrumOmegaKMZAmEqnLHS[omegaR_?NumberQ,omegaI_?NumericQ,kR_?NumberQ,kI_?NumberQ,spect_]:=Module[{eqnLHSM,spectM,IntFun0nFFM,IntFun1nFFM,IntFun2nFFM},


spectM=spect;

IntFun0nFFM=NIntegrate[(spectM@u)/(1-(kR+I*kI) u/(omegaR+I*omegaI)),{u,-1,1}];
IntFun1nFFM=NIntegrate[(spectM@u)u/(1-(kR+I*kI) u/(omegaR+I*omegaI)),{u,-1,1}];
IntFun2nFFM=NIntegrate[(spectM[u]u^2)/(1-(kR+I*kI) u/(omegaR+I*omegaI)),{u,-1,1}];


eqnLHSM[omega_,k_]:=omega-(IntFun0nFFM-IntFun2nFFM-\[Sqrt]((IntFun0nFFM+IntFun2nFFM-2IntFun1nFFM)(IntFun0nFFM+IntFun2nFFM+2IntFun1nFFM)))/(-4);



{ComplexExpand[Re[
eqnLHSM[omegaR+I omegaI,kR+I kI]
]],
ComplexExpand[Im[
eqnLHSM[omegaR+I omegaI,kR+I kI]
]]
}


]


(* END Continuous Feed in Function as spectrum *)


EndPackage[]




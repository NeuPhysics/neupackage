(* ::Package:: *)

BeginPackage["dispersion`"]

eta=DiagonalMatrix[{1,-1,-1,-1}];

kno[nIndex_,theta_,phi_]:={1,nIndex Cos[phi]Sin[theta],nIndex Sin[phi]Sin[theta],nIndex Cos[theta]};
k[omega_,nIndex_,theta_,phi_]:=omega kno[nIndex,theta,phi];
v[theta_,phi_]:={1,Cos[phi]Sin[theta],Sin[phi]Sin[theta],Cos[theta]};

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

EndPackage[]

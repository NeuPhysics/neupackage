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



(* BEGIN Arbitrary N Beams MAA Solution *)

NBeamsLR::usage = "Generates N Beams with crossing spectrum; NBeamsLR[beams] returns a nested list {{-g,theta1},{g,theta2}}";
NBeamsH::usage = "Generates N Beams with homogeneous spectrum; NBeamsH[beams] returns a nested list {{g,theta1},{g,theta2}}";
NBeamsOmegaNPlt::usage="NBeamsOmegaNPlt[nbeams_,nRange_:{-4.5,4.5}]";
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
(#[[1]] (1-Cos[#[[2]]]^2)/(1-n Cos[#[[2]]]))&/@nbeamsM
];

gridlinesM=Join[1/Cos[#[[2]]]&/@nbeams,{{1/(0.6),Purple}}];

beamsOmegaKPltM=Plot[omeganM[n],{n,nRange[[1]],nRange[[2]]},Frame->True,FrameLabel->{"n","\[Omega]"},ImageSize->Large,PlotLabel->"Beams: "<>ToString@Length@nbeamsM,PlotStyle->Gray,PlotRange->{nRange,kRange},Exclusions->gridlinesM];

Show[beamsOmegaKPltM,GridLines->{gridlinesM,None}]

];

NBeamsKNPlt[nbeams_,nRange_:{-4.5,4.5},kRange_:{-10,10}]:=Module[{beamsomegaKM,omeganM,knM,beamsOmegaKCleanM,beamsOmegaKPltM,beamsKNPltM,nbeamsM,beamsSingularityM,gridlinesM},

nbeamsM=nbeams;

omeganM[n_]:=(1/4)Total[
(#[[1]] (1-Cos[#[[2]]]^2)/(1-n Cos[#[[2]]]))&/@nbeamsM
];

knM[n_]:= n omeganM[n];

gridlinesM=Join[1/Cos[#[[2]]]&/@nbeams,{{1/(0.6),Purple}}];

beamsKNPltM=Plot[knM[n],{n,nRange[[1]],nRange[[2]]},Frame->True,FrameLabel->{"n","k"},ImageSize->Large,PlotLabel->"Beams: "<>ToString@Length@nbeamsM,PlotStyle->Gray,PlotRange->{nRange,kRange},Exclusions->gridlinesM];

Show[beamsKNPltM,GridLines->{gridlinesM,None}]

];


NBeamsAxialSymOmegaKEqn[omega_,k_,beams_]:=4== Total[
(#[[1]](1-Cos[#[[2]]]^2) )/(omega-k Cos[#[[2]]])&/@beams
];
NBeamsAxialSymOmegaKPolyEqn[omega_,k_,beams_]:=4 Apply[Times,(omega-k Cos[#[[2]]])&/@beams]== Total[
Table[
(beams[[i,1]](1-Cos[beams[[i,2]]]^2) )Apply[Times,(omega-k Cos[#[[2]]])&/@beams]/(omega-k Cos[beams[[i,2]]]),
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

beamsSingularityM=Plot[(Transpose@Cos@nbeamsM)[[2]]x,{x,kRange[[1]],kRange[[2]]},PlotStyle->Directive[Orange,Dashing[0.03]]];

{beamsOmegaKPltM=ListPlot[beamsomegaKM,Joined->False,Frame->True,FrameLabel->{"k","\[Omega]"},ImageSize->Large,PlotLabel->ToString@Length@nbeamsM<>" Beams",PlotStyle->Gray,PlotRange->{kRange,Automatic},AspectRatio->1],
Show[beamsOmegaKPltM,beamsSingularityM],
beamsomegaKM
}

]


(* END Arbitrary N Beams MAA Solution *)


(* BEGIN Continuous Useful *)

IntFun0[omega_,k_,ct1_,ct2_]:=(omega (-Log[(1-(ct1 k)/omega)/(1-(ct2 k)/omega)]))/k;
IntFun1[omega_,k_,ct1_,ct2_]:=(omega ((-ct1+ct2) k+omega Log[(ct2 k-omega)/(ct1 k-omega)]))/k^2;
IntFun2[omega_,k_,ct1_,ct2_]:=(omega (-(ct1-ct2) k ((ct1+ct2) k+2 omega)+2 omega^2 Log[(ct2 k-omega)/(ct1 k-omega)]))/(2 k^3);
IntFun0n[n_,ct1_,ct2_]:= (-Log[(1-ct1 n)/(1-ct2 n)])/n;
IntFun1n[n_,ct1_,ct2_]:=((-ct1+ct2) +Log[(ct2 n-1)/(ct1 n-1)]/n)/n;
IntFun2n[n_,ct1_,ct2_]:= (-(ct1-ct2) ((ct1+ct2) +2/n)+2 Log[(ct2 n-1)/(ct1 n-1)]/n^2)/(2  n);


ConAxialSymOmegaNMAA[n_,spect_:{{{0.9,0.6},-1},{{0.6,0.3},1}}]:=Module[{spectM},

spectM=spect;
(*{{{0.9,0.6},-1},
{{0.6,0.3},1}};
*)

Total[
#[[2]](IntFun0n[n,#[[1,1]],#[[1,2]]]-IntFun2n[n,#[[1,1]],#[[1,2]]])/(4)&/@spectM
]

]

(* END Continuous Useful *)


(* BEGIN Continuous MAA and MZA Solution: NO CROSSING *)

ConAxialSymOmegaNMZA[n_,spect_:{{{0.9,0.6},-1},{{0.6,0.3},1}}]:=Module[{spectM,i1mzaM,i2mzaM,i0mzaM},

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


(* END Continuous MAA and MZA Solution: NO CROSSING  *)


EndPackage[]




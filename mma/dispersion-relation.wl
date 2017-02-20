(* ::Package:: *)

BeginPackage["dispersion`"]

eta=DiagonalMatrix[{1,-1,-1,-1}];

kno[nIndex_,theta_,phi_]:={1,nIndex Cos[phi]Sin[theta],nIndex Sin[phi]Sin[theta],nIndex Cos[theta]};
k[omega_,nIndex_,theta_,phi_]:=omega kno[nIndex,theta,phi];
v[theta_,phi_]:={1,Cos[phi]Sin[theta],Sin[phi]Sin[theta],Cos[theta]};

vvMatrix[theta_,phi_]=Outer[Times,v[theta,phi],v[theta,phi]];


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

pltDataM=eigenDataM;

ListPlot[pltDataM,Joined->False,Frame->True,ImageSize->Large,PlotRange->pltRange,FrameLabel->{"k","omega"}]

];

(* END for discrete beams *)

EndPackage[]

(* ::Package:: *)

BeginPackage["neutrino`"]

(*Unit Conversion*)


(*Rotation Matrices*)


(*Stimulated Matter Effect*)


(*Parameters*)

(*Parameters END*)


(*2-Frequency Matter Perturbation Transition Probability*)


solN::usage = "numerical solution of the schrodinger equation for two frequency perturbation";
pltN::usage = "plots the transition probability";
bFun::usage = "";
phase::usage = "";
h1::usage = "";
h2::usage = "";
sol0::usage = "";
plt0::usage = "";
sol0Orders::usage = "";
plt0Orders::usage = "";

Begin["`Private`"]

phi1V=0;
phi2V=0;
init={psi1[0],psi2[0]}=={1,0};
imgsize=700;

solN[k1_,k2_,a1_,a2_,thetam_,endpoint_]:=Module[{k1M,k2M,a1M,a2M,thetamM,zk1M,zk2M,n101,n201,n202,n102,n1p1,n2p2,n2p1,n1p2,n1m1,n2m1,n2m2,n1m2,h},
k1M=k1;
k2M=k2;
a1M=a1;
a2M=a2;
thetamM=thetam;
zk1M=a1M Cos[2thetamM]/k1M;
zk2M=a2M Cos[2thetamM]/k2M;
n101=Round[1/k1M];
n201=Round[(1-n101 k1M)/k2M];
n202=Round[1/k2M];
n102=Round[(1-n202 k2M)/k1M];
n1p1=n101+1;
n2p1=Round[(1-n1p1 k1M)/k2M];
n2p2=n202+1;
n1p2=Round[(1-n2p2 k2M)/k1M];
n1m1=n101-1;
n2m1=Round[(1-n1m1 k1M)/k2M];
n2m2=n202-1;
n1m2=Round[(1-n2m2 k1M)/k2M];
h=(a1M Sin[k1M x+phi1V]+a2M Sin[k2M x+phi2V])Sin[2thetamM]Exp[I(-x-Cos[2thetamM](a1M/k1M Cos[k1M x+phi1V]+a2M/k2M Cos[k2M x+phi2V]))]/2;
NDSolve[I D[{psi1[x],psi2[x]},x]=={{0,h},{Conjugate[h],0}}.{psi1[x],psi2[x]}&&init,{psi1,psi2},{x,0,endpoint}]
]

pltN[k1_,k2_,a1_,a2_,thetam_,endpoint_,color_]:=Plot[Evaluate[Abs[psi2[x]]^2/.solN[k1,k2,a1,a2,thetam,endpoint]],{x,0,endpoint},ImageSize->imgsize,Frame->True,FrameLabel->{"x","Transition Probability"},PlotStyle->{Dotted,color},PlotLegends->Placed[Style["Numerical;",color],{Top,Center}]]


bFun[n1_,n2_,zk1_,zk2_,thetam_]:=Tan[2thetam]/2 (-(-I)^(n1+n2))(n1)BesselJ[n1,zk1]BesselJ[n2,zk2]
phase[n1_,n2_,phi1_,phi2_]:=Exp[I(n1 phi1+n2 phi2)]
h1[n1_,n2p_,k1_,k2_,a1_,a2_,zk1_,zk2_,phi1_,phi2_,thetam_,x_]:=k1 bFun[n1,n2p,zk1,zk2,thetam]phase[n1,n2p,phi1,phi2]Exp[I(n1 k1+n2p k2-1)x]
h2[n2_,n1p_,k1_,k2_,a1_,a2_,zk1_,zk2_,phi1_,phi2_,thetam_,x_]:=k2 bFun[n2,n1p,zk2,zk1,thetam]phase[n1p,n2,phi1,phi2]Exp[I(n1p k1+n2 k2-1)x]

sol0[k1_,k2_,a1_,a2_,thetam_,endpoint_]:=Module[{k1M,k2M,a1M,a2M,thetamM,zk1M,zk2M,n101,n201,n202,n102,n1p1,n2p2,n2p1,n1p2,n1m1,n2m1,n2m2,n1m2,h},
k1M=k1;
k2M=k2;
a1M=a1;
a2M=a2;
thetamM=thetam;
zk1M=a1M Cos[2thetamM]/k1M;
zk2M=a2M Cos[2thetamM]/k2M;
n101=Round[1/k1M];
n201=Round[(1-n101 k1M)/k2M];
n202=Round[1/k2M];
n102=Round[(1-n202 k2M)/k1M];
n1p1=n101+1;
n2p1=Round[(1-n1p1 k1M)/k2M];
n2p2=n202+1;
n1p2=Round[(1-n2p2 k2M)/k1M];
n1m1=n101-1;
n2m1=Round[(1-n1m1 k1M)/k2M];
n2m2=n202-1;
n1m2=Round[(1-n2m2 k1M)/k2M];
h=h1[n101,n201,k1M,k2M,a1M,a2M,zk1M,zk2M,phi1V,phi2V,thetamM,x]+h2[n202,n102,k1M,k2M,a1M,a2M,zk1M,zk2M,phi1V,phi2V,thetamM,x];
NDSolve[I D[{psi1[x],psi2[x]},x]=={{0,h},{Conjugate[h],0}}.{psi1[x],psi2[x]}&&init,{psi1,psi2},{x,0,endpoint}]
]

plt0[k1_,k2_,a1_,a2_,thetam_,endpoint_,color_]:=Plot[Evaluate[Abs[psi2[x]]^2/.sol0[k1,k2,a1,a2,thetam,endpoint]],{x,0,endpoint},ImageSize->imgsize,Frame->True,FrameLabel->{"x","Transition Probability"},PlotStyle->color,PlotLegends->Placed[Style["Lowest Order;",color],{Top,Center}]]


sol0Orders[k1_,k2_,a1_,a2_,thetam_,endpoint_,order1_,order2_]:=Module[{k1M,k2M,a1M,a2M,thetamM,zk1M,zk2M,n101,n201,n202,n102,h},
k1M=k1;
k2M=k2;
a1M=a1;
a2M=a2;
thetamM = thetam;
zk1M=a1M Cos[2thetamM]/k1M;
zk2M=a2M Cos[2thetamM]/k2M;
n101=Round[1/k1M];
n201=Round[(1-n101 k1M)/k2M];
n202=Round[1/k2M];
n102=Round[(1-n202 k2M)/k1M];
h=Total@Total@Table[h1[n101+i1,n201+i2,k1M,k2M,a1M,a2M,zk1M,zk2M,phi1V,phi2V,thetamM,x],{i1,-order1,order1},{i2,-order2,order2}]+Total@Total@Table[h2[n202+i1,n102+i2,k1M,k2M,a1M,a2M,zk1M,zk2M,phi1V,phi2V,thetamM,x],{i1,-order1,order1},{i2,-order2,order2}];
NDSolve[I D[{psi1[x],psi2[x]},x]=={{0,h},{Conjugate[h],0}}.{psi1[x],psi2[x]}&&init,{psi1,psi2},{x,0,endpoint}]
]

plt0Orders[k1_,k2_,a1_,a2_,thetam_,endpoint_,order1_,order2_,legends_,color_]:=Plot[Evaluate[Abs[psi2[x]]^2/.sol0Orders[k1,k2,a1,a2,thetam,endpoint,order1,order2]],{x,0,endpoint},ImageSize->imgsize,Frame->True,FrameLabel->{"x","Transition Probability"},PlotStyle->color,PlotLegends->Placed[Style[ToString[legends],color],{Top,Center}]]

End[]


(*END 2-Frequency Matter Perturbation Transition Probability*)


(*2-Frequency Matter Perturbation Coefficients*)


bCoef::usage = "Takes six parameters [n1,n2,k1,k2,a1,a2] and returns the coefficient labeled by n1, n2. In a 2-frequency Hamiltonian";
coefDenPlt::usage = "Takes six parameters [k1,k2,a1,a2,thetam,range] and returns four density plots of the Abs@phases and Abs@coefficients as a function of n1 and n2";
coefDenPltReIm::usage = "Takes six parameters [k1,k2,a1,a2,thetam,range] and returns four density plots of the Re@coefficients and Im@coefficients as a function of n1 and n";
coefSlicings::usage = "The Abs@coefficient values through some lines on the n1-n2 plane.";
amp2Freq::usage = "Given [n1,n2,k1,k2,a1,a2], calculate the oscillation amplitude of oscillations.";
amp2FreqPlt::usage = "Given [n1,n2,a1,a2,thetam,range for k1,range for k2], displays the contour plot for the transition amplitude";
amp2FreqContourPltList::usage = "List contour plots of transition probability amplitude for a pair of integers given [list of range for n1, list of range for n2, k1, k2, amplitude of k1, amplitude of k2, thetam, list of range for k1, list of range for k2]";
amp2FreqCombinedPlt::usage = "Density plot of amplitude combined for n1Range, n2Range, given input [a1,a2,thetam,k1Range,k2Range,n1Range,n2Range]; Example, amp2FreqCombinedPltTest[a1V,a2V,thetamV,{0.01,1,0.01},{0.01,1,0.01},{-1,1},{-1,1}]";
hamiltonian0::usage = "";
plt0n1n2List::usage = "";
sol0n1n2List::usage = "";
plt0n1n2::usage = "";
sol0n1n2::usage = "";
amp2FreqCombinedPlt2::usage = "";
width2Freq::usage = "Calculate the width for given parameters [n1,n2,k1,k2,a1,a2,thetam], where k1,k2 are the point of interest or the wave vectors of the system";
distanceFun::usage = "";
qDis2Width::usage = "";




Begin["`Private`"]


phi1V=0;
phi2V=0;
init={psi1[0],psi2[0]}=={1,0};
imgsize=700;

(*All quantities are normalized. k1,a1,k2,a2 are normalized by omegam. x means omegam times the actually distance *)

bCoef[n1_,n2_,k1_,k2_,a1_,a2_,thetam_]:=Piecewise[{{
-(-I)^(n1+n2) Tan[2thetam]/2 n1 k1 BesselJ[n1,a1/k1 Cos[2thetam]]BesselJ[n2,a2/k2 Cos[2thetam]],k1!=0&&k2!=0
},{
-(-I)^(n1+n2) Tan[2thetam]/2 n1 k1 BesselJ[n1,Infinity]BesselJ[n2,a2/k2 Cos[2thetam]],k1==0||k2!=0
},{
-(-I)^(n1+n2) Tan[2thetam]/2 n1 k1 BesselJ[n1,a1/k1 Cos[2thetam]]BesselJ[n2,Infinity],k1!=0||k2==0
}}];

(*Hamiltonian if no summation is counted*)
hamiltonian0[n1_,n2_,k1_,k2_,a1_,a2_,thetam_]:=(bCoef[n1,n2,k1,k2,a1,a2,thetam]+bCoef[n2,n1,k2,k1,a2,a1,thetam])*Exp[I(n1*phi1V+n2*phi2V)]*Exp[I(n1*k1+n2*k2-1)x];

(*Effective Width on a k2-k1 plane should be*)

width2Freq[n1_,n2_,k1_,k2_,a1_,a2_,thetam_]:= Module[{k1p,k2p,k1p2,k2p2,k1p3,k2p3},
k1p=ConditionalExpression[(n2^2*k1+n2*k2+n1)/(n1^2+n2^2),n1!=0&&n2!=0]//Quiet;
k2p=ConditionalExpression[n1*k1p/n2-1/n2,n1!=0&&n2!=0]//Quiet;
k1p2=ConditionalExpression[k1,n1==0&&n2!=0]//Quiet;
k2p2=ConditionalExpression[1/n2,n1==0&&n2!=0]//Quiet;
k1p3=ConditionalExpression[1/n1,n1!=0&&n2==0]//Quiet;
k2p3=ConditionalExpression[k2,n1!=0&&n2==0]//Quiet;

Piecewise[{{
Abs[2(bCoef[n1,n2,k1p,k2p,a1,a2,thetam]+bCoef[n2,n1,k2p,k1p,a2,a1,thetam])],n1!=0&&n2!=0
},{
Abs[2(bCoef[n1,n2,k1p2,k2p2,a1,a2,thetam]+bCoef[n2,n1,k2p2,k1p2,a2,a1,thetam])],n1==0&&n2!=0
},{
Abs[2(bCoef[n1,n2,k1p3,k2p3,a1,a2,thetam]+bCoef[n2,n1,k2p3,k1p3,a2,a1,thetam])],n1!=0&&n2==0
},{
0,n1==0&&n2==0
}}]
];


distanceFun[n1_,n2_,k1Point_,k2Point_]:=Piecewise[{{
Abs[n1*k1Point+n2*k2Point-1]/Sqrt[n1^2+n2^2],Not[n1==0&&n2==0]
},{
Infinity,n1==0&&n2==0
}}]


(*Define Q value which is ratio of distance to width*)
qDis2Width[n1_,n2_,k1_,k2_,a1_,a2_,thetam_]:=Piecewise[{{
distanceFun[n1,n2,k1,k2]/width2Freq[n1,n2,k1,k2,a1,a2,thetam],width2Freq[n1,n2,k1,k2,a1,a2,thetam]!=0
},{
Infinity,width2Freq[n1,n2,k1,k2,a1,a2,thetam]==0&&distanceFun[n1,n2,k1,k2]!=0
},{
0,width2Freq[n1,n2,k1,k2,a1,a2,thetam]==0&&distanceFun[n1,n2,k1,k2]==0
}}]




sol0n1n2[n1_,n2_,k1_,k2_,a1_,a2_,thetam_,endpoint_]:=Module[{n1M,n2M,k1M,k2M,a1M,a2M,thetamM,hamil},
n1M=n1;
n2M=n2;
k1M=k1;
k2M=k2;
a1M=a1;
a2M=a2;
thetamM = thetam;

hamil = hamiltonian0[n1M,n2M,k1M,k2M,a1M,a2M,thetamM]

NDSolve[I D[{psi1[x],psi2[x]},x]=={{0,hamil},{Conjugate[hamil],0}}.{psi1[x],psi2[x]}&&init,{psi1,psi2},{x,0,endpoint}]
]

plt0n1n2[n1_,n2_,k1_,k2_,a1_,a2_,thetam_,endpoint_,legends_,color_]:=Plot[Evaluate[Abs[psi2[x]]^2/.sol0n1n2[n1,n2,k1,k2,a1,a2,thetam,endpoint]],{x,0,endpoint},ImageSize->imgsize,Frame->True,FrameLabel->{"x","Transition Probability"},PlotStyle->color,PlotLegends->Placed[Style[ToString[legends],color],{Top,Center}]]

sol0n1n2List[listInput_,k1_,k2_,a1_,a2_,thetam_,endpoint_]:=Module[{n1M,n2M,k1M,k2M,a1M,a2M,thetamM,hamil,NlistM,listInputM},
k1M=k1;
k2M=k2;
a1M=a1;
a2M=a2;
thetamM = thetam;
listInputM = listInput;

hamil = Total@Table[hamiltonian0[NlistM[[1]],NlistM[[2]],k1M,k2M,a1M,a2M,thetamM],{NlistM,listInputM}];

NDSolve[I D[{psi1[x],psi2[x]},x]=={{0,hamil},{Conjugate[hamil],0}}.{psi1[x],psi2[x]}&&init,{psi1,psi2},{x,0,endpoint}]
]

plt0n1n2List[listInput_,k1_,k2_,a1_,a2_,thetam_,endpoint_,legends_,color_]:=Plot[Evaluate[Abs[psi2[x]]^2/.sol0n1n2List[listInput,k1,k2,a1,a2,thetam,endpoint]],{x,0,endpoint},ImageSize->imgsize,Frame->True,FrameLabel->{"\!\(\*OverscriptBox[\(x\), \(^\)]\)","Transition Probability"},PlotStyle->color,PlotLegends->Placed[Style[ToString[legends],color],{Top,Center}]]


coefDenPlt[k1_,k2_,a1_,a2_,thetam_,range_]:=Module[{n1M,n2M,k1M,k2M,a1M,a2M,thetamM,lineRWA,phaseList,coeffList,coeffListRe,coeffListIm,coeffListAbs,list0,list1,list2,list3,list4,n2List1,n1List2,n2List3,n1List4,meshStyle,parList},
meshStyle=None;
k1M=k1;
k2M=k2;
a1M=a1;
a2M=a2;
thetamM=thetam;

phaseList=Flatten[Table[{n1,n2,Abs[n1 k1M+ n2 k2M-1]},{n1,-range,range},{n2,-range,range}],1];
coeffList=Flatten[Table[{n1,n2,bCoef[n1,n2,k1M,k2M,a1M,a2M,thetamM]+bCoef[n2,n1,k2M,k1M,a2M,a1M,thetamM]},{n1,-range,range},{n2,-range,range}],1];
coeffListRe=Flatten[Table[{n1,n2,Log@Re[bCoef[n1,n2,k1M,k2M,a1M,a2M,thetamM]+bCoef[n2,n1,k2M,k1M,a2M,a1M,thetamM]]},{n1,-range,range},{n2,-range,range}],1];
coeffListIm=Flatten[Table[{n1,n2,Log@Im[bCoef[n1,n2,k1M,k2M,a1M,a2M,thetamM]+bCoef[n2,n1,k2M,k1M,a2M,a1M,thetamM]]},{n1,-range,range},{n2,-range,range}],1];
coeffListAbs=Flatten[Table[{n1,n2,Log@Abs[bCoef[n1,n2,k1M,k2M,a1M,a2M,thetamM]+bCoef[n2,n1,k2M,k1M,a2M,a1M,thetamM]]},{n1,-range,range},{n2,-range,range}],1];
lineRWA=Table[{n1,Round[(1-n1 k1M)/k2M]},{n1,-range,range}];

parList=Text["Parameters:\n k1="<>ToString[k1M]<>";\n k2="<>ToString[k2M]<>";\n a1="<>ToString[a1M]<>";\n a2="<>ToString[a2M]<>";\n\n"];

Grid[
{{
Show[ListDensityPlot[phaseList,Mesh->meshStyle,FrameLabel->{"n1","n2"},ColorFunction->"AvocadoColors",PlotLegends->All,ImageSize->imgsize,PlotLabel->"Phases" ,InterpolationOrder->0],ListPlot[lineRWA,ImageSize->imgsize,PlotRange->{{-range,range},{-range,range}}]]
,
Show[ListDensityPlot[coeffListAbs,Mesh->meshStyle,FrameLabel->{"n1","n2"},ColorFunction->"AvocadoColors",PlotLegends->All,ImageSize->imgsize,PlotLabel->"Absolute Value of The Coefficient" ,InterpolationOrder->0],ListPlot[lineRWA,ImageSize->imgsize,PlotRange->{{-range,range},{-range,range}}]]
}}
]
]

coefDenPltReIm[k1_,k2_,a1_,a2_,thetam_,range_]:=Module[{n1M,n2M,k1M,k2M,a1M,a2M,thetamM,phaseList,lineRWA,coeffList,coeffListRe,coeffListIm,coeffListAbs,list0,list1,list2,list3,list4,n2List1,n1List2,n2List3,n1List4,meshStyle,parList},
meshStyle=None;
k1M=k1;
k2M=k2;
a1M=a1;
a2M=a2;
thetamM=thetam;

phaseList=Flatten[Table[{n1,n2,Abs[n1 k1M+ n2 k2M-1]},{n1,-range,range},{n2,-range,range}],1];
coeffList=Flatten[Table[{n1,n2,bCoef[n1,n2,k1M,k2M,a1M,a2M,thetamM]+bCoef[n2,n1,k2M,k1M,a2M,a1M,thetamM]},{n1,-range,range},{n2,-range,range}],1];
coeffListRe=Flatten[Table[{n1,n2,Log@Re[bCoef[n1,n2,k1M,k2M,a1M,a2M,thetamM]+bCoef[n2,n1,k2M,k1M,a2M,a1M,thetamM]]},{n1,-range,range},{n2,-range,range}],1];
coeffListIm=Flatten[Table[{n1,n2,Log@Im[bCoef[n1,n2,k1M,k2M,a1M,a2M,thetamM]+bCoef[n2,n1,k2M,k1M,a2M,a1M,thetamM]]},{n1,-range,range},{n2,-range,range}],1];
coeffListAbs=Flatten[Table[{n1,n2,Log@Abs[bCoef[n1,n2,k1M,k2M,a1M,a2M,thetamM]+bCoef[n2,n1,k2M,k1M,a2M,a1M,thetamM]]},{n1,-range,range},{n2,-range,range}],1];
lineRWA=Table[{n1,Round[(1-n1 k1M)/k2M]},{n1,-range,range}];

parList=Text["Parameters:\n k1="<>ToString[k1M]<>";\n k2="<>ToString[k2M]<>";\n a1="<>ToString[a1M]<>";\n a2="<>ToString[a2M]<>";\n\n"];

Grid[
{{
Show[ListDensityPlot[coeffListRe,Mesh->meshStyle,FrameLabel->{"n1","n2"},ColorFunction->"AvocadoColors",PlotLegends->All,ImageSize->imgsize,PlotLabel->"Real Part of The Coefficient" ,InterpolationOrder->0],ListPlot[lineRWA,ImageSize->imgsize,PlotRange->{{-range,range},{-range,range}}]],
Show[ListDensityPlot[coeffListIm,Mesh->meshStyle,FrameLabel->{"n1","n2"},ColorFunction->"AvocadoColors",PlotLegends->All,ImageSize->imgsize,PlotLabel->"Imaginary Part of The Coefficient" ,InterpolationOrder->0],ListPlot[lineRWA,ImageSize->imgsize,PlotRange->{{-range,range},{-range,range}}]]
}}
]
]


coefSlicings[k1_,k2_,a1_,a2_,thetam_,range_,range2_]:=Module[{n1M,n2M,k1M,k2M,a1M,a2M,thetamM,coeffList,coeffListRe,coeffListIm,coeffListAbs,list0,list02,list1,list2,list3,list4,n2List1,n1List2,n2List3,n1List4,meshStyle},
meshStyle=None;
k1M=k1;
k2M=k2;
a1M=a1;
a2M=a2;
thetamM=thetam;

coeffList=Flatten[Table[{n1,n2,bCoef[n1,n2,k1M,k2M,a1M,a2M,thetamM]+bCoef[n2,n1,k2M,k1M,a2M,a1M,thetamM]},{n1,-range,range},{n2,-range,range}],1];
coeffListRe=Flatten[Table[{n1,n2,Log@Re[bCoef[n1,n2,k1M,k2M,a1M,a2M,thetamM]+bCoef[n2,n1,k2M,k1M,a2M,a1M,thetamM]]},{n1,-range,range},{n2,-range,range}],1];
coeffListIm=Flatten[Table[{n1,n2,Log@Im[bCoef[n1,n2,k1M,k2M,a1M,a2M,thetamM]+bCoef[n2,n1,k2M,k1M,a2M,a1M,thetamM]]},{n1,-range,range},{n2,-range,range}],1];
coeffListAbs=Flatten[Table[{n1,n2,Log@Abs[bCoef[n1,n2,k1M,k2M,a1M,a2M,thetamM]+bCoef[n2,n1,k2M,k1M,a2M,a1M,thetamM]]},{n1,-range,range},{n2,-range,range}],1];
list0=Table[{n1,Abs[bCoef[n1,Round[(1-n1 k1M)/k2M],k1M,k2M,a1M,a2M,thetamM]+bCoef[Round[(1-n1 k1M)/k2M],n1,k2M,k1M,a2M,a1M,thetamM]]},{n1,-range,range}];
list02=Table[{n1,Abs[bCoef[n1,Round[(1-n1 k1M)/k2M],k1M,k2M,a1M,a2M,thetamM]+bCoef[Round[(1-n1 k1M)/k2M],n1,k2M,k1M,a2M,a1M,thetamM]]},{n1,-range2,range2}];
n2List1=0;
list1=Table[{n1,Abs[bCoef[n1,n2List1,k1M,k2M,a1M,a2M,thetamM]+bCoef[n2List1,n1,k2M,k1M,a2M,a1M,thetamM]]},{n1,-range,range}];
n1List2=0;
list2=Table[{n2,Abs[bCoef[n1List2,n2,k1M,k2M,a1M,a2M,thetamM]+bCoef[n2,n1List2,k2M,k1M,a2M,a1M,thetamM]]},{n2,-range,range}];
n2List3=Round[range/2];
list3=Table[{n1,Abs[bCoef[n1,n2List3,k1M,k2M,a1M,a2M,thetamM]+bCoef[n2List3,n1,k2M,k1M,a2M,a1M,thetamM]]},{n1,-range,range}];
n1List4=Round[range/2];
list4=Table[{n2,Abs[bCoef[n1List4,n2,k1M,k2M,a1M,a2M,thetamM]+bCoef[n2,n1List4,k2M,k1M,a2M,a1M,thetamM]]},{n2,-range,range}];



Grid[{{ListLogPlot[list0,ImageSize->imgsize,Frame->True,PlotLabel->"Slicings",Joined->True,PlotStyle->{Dotted,Red},PlotMarkers->{Automatic, Small},PlotLegends->Placed[{"RWA"},{Top,Center}],PlotRange->{{-range2,range2},Automatic},GridLines->{{{Round[1/k1M],Dashed},{Round[(1-Round[1/k2M] k2M)/k1M],Red}}, {1}}]
,
ListLogPlot[{list0,list1,list2,list3,list4},ImageSize->imgsize,Frame->True,PlotLabel->"Slicings",Joined->True,PlotStyle->{{Dotted,Red},Blue,Black,Gray,Magenta},PlotMarkers->{Automatic, Small},PlotLegends->Placed[{"RWA","n2="<>ToString[n2List1],"n1="<>ToString[n1List2],"n2="<>ToString[n2List3],"n1="<>ToString[n1List4]},{Top,Center}]]}}]

]



amp2Freq[n1_,n2_,k1_,k2_,a1_,a2_,thetam_]:=Module[{n1M,n2M,k1M,k2M,a1M,a2M,thetamM,fM2,gM2},
k1M=k1;
k2M=k2;
a1M=a1;
a2M=a2;
n1M=n1;
n2M=n2;
thetamM=thetam;

fM2 = Abs[2*(bCoef[n1M,n2M,k1M,k2M,a1M,a2M,thetamM] + bCoef[n2M,n1M,k2M,k1M,a2M,a1M,thetamM])]^2;
gM2 = Abs[n1M*k1M+n2M*k2M-1]^2;
fM2/(fM2+gM2)
]



amp2FreqPlt[n1_,n2_,a1_,a2_,thetam_,k1Range_,k2Range_]:=Module[{n1M,n2M,k1M,k2M,a1M,a2M,thetamM,rangeM,pltDataM,nRangeM,k1RangeM,k2RangeM},
n1M=n1;
n2M=n2;
a1M=a1;
a2M=a2;
thetamM=thetam;
k1RangeM=k1Range;
k2RangeM=k2Range;

ContourPlot[amp2Freq[n1M,n2M,k1x,k2x,a1M,a2M,thetamM],{k1x,k1RangeM[[1]],k1RangeM[[2]]},{k2x,k2RangeM[[1]],k2RangeM[[2]]},PlotLegends->Automatic,ImageSize->imgsize,PlotRange->{Automatic,Automatic,{0,1}},PlotLabel->"Transition Amplitude: "<>"(n1,n2)=("<>ToString[n1M]<>","<>ToString[n2M]<>")",FrameLabel->{"\!\(\*SubscriptBox[\(k\), \(1\)]\)","\!\(\*SubscriptBox[\(k\), \(2\)]\)"}]
]


amp2FreqContourPltList[n1n2List_,k1_,k2_,a1_,a2_,thetam_,k1Range_,k2Range_]:=Module[{n1M,n2M,k1M,k2M,a1M,a2M,thetamM,k1RangeM,k2RangeM,n1n2ListM},
k1M=k1;
k2M=k2;
a1M=a1;
a2M=a2;
thetamM=thetam;
n1n2ListM = n1n2List;
k1RangeM = k1Range;
k2RangeM = k2Range;

Grid[{Flatten@Table[ContourPlot[amp2Freq[n1n2ListTable[[1]],n1n2ListTable[[2]],k1x,k2x,a1M,a2M,thetamM],{k1x,k1RangeM[[1]],k1RangeM[[2]]},{k2x,k2RangeM[[1]],k2RangeM[[2]]},PlotLabel->"Contour Plot of Transition Amplitude: "<>"{n1,n2}="<>ToString[n1n2ListTable],PlotRange->{Automatic,Automatic,{0,1}},PlotLegends->Automatic,ImageSize->imgsize,FrameLabel->{"\!\(\*SubscriptBox[\(k\), \(1\)]\)","\!\(\*SubscriptBox[\(k\), \(2\)]\)"}],{n1n2ListTable,n1n2ListM}]}]
]



amp2FreqCombinedPlt[a1_,a2_,thetam_,k1Range_,k2Range_,n1Range_,n2Range_]:=Module[{n1M,n2M,k1M,k2M,a1M,a2M,thetamM,rangeM,pltDataRawM,pltDataM,k1RangeM,k2RangeM,n1RangeM,n2RangeM},
a1M=a1;
a2M=a2;
thetamM=thetam;
k1RangeM=k1Range;
k2RangeM=k2Range;
n1RangeM=n1Range;
n2RangeM=n2Range;

pltDataRawM[n1_,n2_]:=Flatten[Table[{k1x,k2x,amp2Freq[n1,n2,k1x,k2x,a1M,a2M,thetamM]},{k1x,k1RangeM[[1]],k1RangeM[[2]],k1RangeM[[3]]},{k2x,k2RangeM[[1]],k2RangeM[[2]],k2RangeM[[3]]}],1];
pltDataM=Total[
Flatten[Table[pltDataRawM[n1,n2],{n1,n1RangeM[[1]],n1RangeM[[2]]},{n2,n2RangeM[[1]],n2RangeM[[2]]}],1]
]/((First@Differences[n1Range]+1)*(First@Differences[n2Range]+1));

ListDensityPlot[pltDataM,ImageSize->imgsize,Mesh->All,ColorFunction->"AvocadoColors",InterpolationOrder->0,PlotLegends->Automatic,PlotRange->{Automatic,Automatic,{0,1}},FrameLabel->{"\!\(\*SubscriptBox[\(k\), \(1\)]\)","\!\(\*SubscriptBox[\(k\), \(2\)]\)"},PlotLabel->"Combined Density Plot of Transition Amplitudes;"<>"n1 range:"<>ToString[n1RangeM]<>"; n2 range:"<>ToString[n2RangeM]]
]

amp2FreqCombinedPlt2[a1_,a2_,thetam_,k1Range_,k2Range_,n1n2List_]:=Module[{n1M,n2M,k1M,k2M,a1M,a2M,thetamM,rangeM,pltDataRawM,pltDataM,k1RangeM,k2RangeM,n1n2ListM},
a1M=a1;
a2M=a2;
thetamM=thetam;
k1RangeM=k1Range;
k2RangeM=k2Range;
n1n2ListM=n1n2List;

pltDataRawM[n1_,n2_]:=Flatten[Table[{k1x,k2x,amp2Freq[n1,n2,k1x,k2x,a1M,a2M,thetamM]},{k1x,k1RangeM[[1]],k1RangeM[[2]],k1RangeM[[3]]},{k2x,k2RangeM[[1]],k2RangeM[[2]],k2RangeM[[3]]}],1];
pltDataM=Total[
Table[pltDataRawM[n1n2ListTable[[1]],n1n2ListTable[[2]]],{n1n2ListTable,n1n2ListM}]
]/(Length[n1n2ListM]);

ListDensityPlot[pltDataM,ImageSize->imgsize,Mesh->All,ColorFunction->"AvocadoColors",InterpolationOrder->0,PlotLegends->Automatic,PlotRange->{Automatic,Automatic,{0,1}},FrameLabel->{"\!\(\*SubscriptBox[\(k\), \(1\)]\)","\!\(\*SubscriptBox[\(k\), \(2\)]\)"},PlotLabel->"Combined Density Plot of Transition Amplitudes;"<>"{n1,n2}="<>ToString[n1n2ListM]]
]


End[]


(*END 2-Frequency Matter Perturbation Coefficients*)


solNN::usage = "numerical solution of the schrodinger equation for N frequency perturbation";
pltNN::usage = "plot the numerical solution of the schrodinger equation for N frequency perturbation";
bCoefN::usage = "calculate the b coefficient of the Hamiltonian";
distanceN::usage = "calculate the EFFECTIVE distance of a given system from a resonance line, NOTE here is NOT the geometrical distance as we have calculated before in 2 frequency case";
qValue::usage = "calcualte the ratio of the EFFECTIVE distance and Effective Width (bCoefN)";
qValueOrderdList::usage = "gives ordered list of qValues and also their corresponding integer combinations";

Begin["`Private`"]



solNN[listOfWaveNumber_,listOfAmplitude_,listOfPhase_,thetam_,endpoint_]:=Module[{listOfWaveNumberM,listOfAmplitudeM,listOfPhaseM,thetamM,endpointM,hamilPart1,hamilPart2,hamil,hamilConj,lengthM},

listOfWaveNumberM=listOfWaveNumber;
listOfAmplitudeM=listOfAmplitude;
listOfPhaseM=listOfPhase;
thetamM=thetam;
endpointM=endpoint;
lengthM=Length[listOfWaveNumber];

hamilPart1 = Total@MapThread[#1*Sin[#2*x+#3]&,{listOfAmplitudeM,listOfWaveNumberM,listOfPhaseM}];
hamilPart2 = Total@MapThread[#1/#2*Cos[#2*x+#3]&,{listOfAmplitudeM,listOfWaveNumberM,listOfPhaseM}];

hamil = Sin[2thetamM]/2 (hamilPart1)Exp[I(-x-Cos[2thetamM](hamilPart2))];
hamilConj = Sin[2thetamM]/2 (hamilPart1)Exp[-I(-x-Cos[2thetamM](hamilPart2))];

NDSolve[I D[{psi1[x],psi2[x]},x]=={{0,hamil},{hamilConj,0}}.{psi1[x],psi2[x]}&&init,{psi1,psi2},{x,0,endpointM}]

]

pltNN[listOfWaveNumber_,listOfAmplitude_,listOfPhase_,thetam_,endpoint_,color_]:=Plot[Evaluate[Abs[psi2[x]]^2/.solNN[listOfWaveNumber,listOfAmplitude,listOfPhase,thetam,endpoint][[1]]],{x,0,endpoint},ImageSize->imgsize,Frame->True,FrameLabel->{"x","Transition Probability"},PlotStyle->{Dotted,color},PlotLegends->Placed[Style["Numerical;",color],{Top,Center}]]


bCoefN[listOfN_,listOfWaveNumber_,listOfAmplitude_,listOfPhase_,thetam_]:=
-(-I)^(Total[listOfN]) Tan[2thetam]/2 (listOfN.listOfWaveNumber) Apply[Times,MapThread[BesselJ[#1,#3/#2 Cos[2thetam]]&,{listOfN,listOfWaveNumber,listOfAmplitude}]];

distanceN[listOfN_,listOfWaveNumber0_]:=Abs[listOfN.listOfWaveNumber0-1];

qValue[listOfN_,listOfWaveNumber_,listOfAmplitude_,listOfPhase_,thetam_]:=distanceN[listOfN,listOfWaveNumber]/Abs[bCoefN[listOfN,listOfWaveNumber,listOfAmplitude,listOfPhase,thetam]];

qValueOrderdList[listOfNRange_,listOfWaveNumber_,listOfAmplitude_,listOfPhase_,thetam_]:=Module[{maptolist,listOfThreadM},

listOfThreadM=Tuples[Table[Table[ni,{ni,listOfNRange[[n,1]],listOfNRange[[n,2]]}],{n,1,Length@listOfNRange}]];

SortBy[
Table[{listN,qValue[listN,listOfWaveNumber,listOfAmplitude,listOfPhase,thetam]},{listN,listOfThreadM}]//Quiet
,Last]

]



End[]


(*END Stimulated Matter Effect*)


EndPackage[]

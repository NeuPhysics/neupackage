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


bFun[n1_,n2_,zk1_,zk2_,thetam_]:=Sin[2thetam]/4 (-(-I)^(n1+n2))(2n1)BesselJ[n1,zk1]BesselJ[n2,zk2]
phase[n1_,n2_,phi1_,phi2_]:=Exp[I(n1 phi1+n2 phi2)]
h1[n1_,n2p_,k1_,k2_,a1_,a2_,zk1_,zk2_,phi1_,phi2_,thetam_,x_]:=k1/Cos[2thetam] bFun[n1,n2p,zk1,zk2,thetam]phase[n1,n2p,phi1,phi2]Exp[I(n1 k1+n2p k2-1)x]
h2[n2_,n1p_,k1_,k2_,a1_,a2_,zk1_,zk2_,phi1_,phi2_,thetam_,x_]:=k2/Cos[2thetam] bFun[n2,n1p,zk2,zk1,thetam]phase[n1p,n2,phi1,phi2]Exp[I(n1p k1+n2 k2-1)x]

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

Begin["`Private`"]

imgsize=700;

bCoef[n1_,n2_,k1_,k2_,a1_,a2_,thetam_]:=-(-I)^(n1+n2) Tan[2thetam]/2 n1 k1 BesselJ[n1,a1/k1 Cos[2thetam]]BesselJ[n2,a2/k2 Cos[2thetam]]

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



End[]


(*END 2-Frequency Matter Perturbation Coefficients*)


(*END Stimulated Matter Effect*)


EndPackage[]

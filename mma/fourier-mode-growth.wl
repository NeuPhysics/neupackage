(* ::Package:: *)

BeginPackage["FourierModeGrowth`"]


matRandom[n_,scale_,kd_,koff_]:=Module[{matM,diagM,midM},

matM=IdentityMatrix[n];
midM=Floor[n/2];
(*
diagM=ReplacePart[
matM,
{i_,i_}\[RuleDelayed]RandomReal[1]
];*)


Do[
matM[[midM+i,midM+i]]= scale i+kd RandomReal[1],
{i,-midM+1,midM-1}
];


(*Do[
matM[[midM+i,midM+i+1]]=k i,
{i,-midM+1,midM-1}
];

Do[
matM[[midM+i+1,midM+i]]=k i,
{i,-midM+1,midM-1}
];*)

Do[
matM[[midM+i,midM+i+1]]=koff RandomReal[1],
{i,-midM+1,midM-1}
];

Do[
matM[[midM+i+1,midM+i]]=koff RandomReal[1],
{i,-midM+1,midM-1}
];

matM

]


matRandom[n_,scale_,kd_,koff_]:=Module[{matM,diagM,midM},

matM=IdentityMatrix[n];
midM=Floor[n/2];
(*
diagM=ReplacePart[
matM,
{i_,i_}\[RuleDelayed]RandomReal[1]
];*)


Do[
matM[[midM+i,midM+i]]= scale i+kd RandomReal[1],
{i,-midM+1,midM-1}
];


(*Do[
matM[[midM+i,midM+i+1]]=k i,
{i,-midM+1,midM-1}
];

Do[
matM[[midM+i+1,midM+i]]=k i,
{i,-midM+1,midM-1}
];*)

Do[
matM[[midM+i,midM+i+1]]=koff RandomReal[1],
{i,-midM+1,midM-1}
];

Do[
matM[[midM+i+1,midM+i]]=koff RandomReal[1],
{i,-midM+1,midM-1}
];

matM

]




matRegular[dimension_,diagscale_,offdiag_]:=Module[{matM,diagM,midM,n},

n=(dimension-1)/2;

matM=IdentityMatrix[2n+1];
midM=Floor[n+1];


Do[
matM[[midM+i,midM+i]]= diagscale  i,
{i,-midM+1,midM-1}
];


Do[
matM[[midM+i,midM+i+1]]=offdiag,
{i,-midM+1,midM-2}
];

Do[
matM[[midM+i+1,midM+i]]=offdiag,
{i,-midM+1,midM-2}
];

matM

]




eigenSys[dimension_, diagscale_, offdiag_]:=Transpose@Sort[Transpose@Eigensystem[  matRegular[dimension, diagscale , offdiag] ],#1[[1]]<#2[[1]]&]

eigenVectorPlot[eigenvectors_,segment_,range_]:= ListLogPlot[
eigenvectors[[segment]],
Frame->True,ImageSize->Large,PlotLegends->Placed[("ID "<>ToString@#&/@segment),{Right,Top}],PlotRange->range,Filling->Axis,FrameLabel->{"Eigenvector Elements ID","Eigenvector Element Value"},PlotLabel->"Eigenvectors"
]

eigenVecSlopeFit[eigenvec_]:=Module[{fitEqnCoeffM,pltM},

fitEqnCoeffM=CoefficientList[
Fit[Log@Abs@eigenvec,{1,x},x],x
];

fitEqnCoeffM

]


EndPackage[]

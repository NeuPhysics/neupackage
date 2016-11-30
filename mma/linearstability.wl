(* ::Package:: *)

BeginPackage["linearstability`"]
(* This package is developed for linear stability analysis in neutrino oscillations *)


(* Four Beams Some Definitions Equations BEGIN *)

FourBeamDeltas::usage = "List of deltas which are the perturbations of density matrix for four beams model. Take one paramter [z_]";
FourBeamDensityMatrix::usage = "Generate perturbed density matrix for four beams model";
FourBeamTransMatrix::usage = "The matrix M in equation i d delta/dz = M. delta for four beam model. Takes paprameters [eta_,mk_,mu_,alpha_,theta1_,theta2_,lambda_,omegav_:0]; theta1 for neutrinos, theta2 for antineutrinos; eta = 1 for NH and -1 for IH; mu =\!\(\*SqrtBox[\(2\)]\)\!\(\*SubscriptBox[\(G\), \(F\)]\)\!\(\*SubscriptBox[\(n\), \(\[Nu]\)]\) where \!\(\*SubscriptBox[\(n\), \(\[Nu]\)]\) is for the beam";
FourBeamRotationMatrixSameAngle::usage = "A matrix that rotates the list of deltas {deltaL[z],deltaLB[z],deltaR[z],deltaRB[z]} into the symmetric and anti-symmetric modes for the case that theta1 = theta2";
FourBeamRotationMatrixDiffAngle::usage = "A matrix that rotates the list of deltas {deltaL[z],deltaLB[z],deltaR[z],deltaRB[z]} into the symmetric and anti-symmetric modes for the case that theta1 != theta2 and lambda=0";

Begin["`Private`"]

FourBeamDeltas[z_]:={deltaL[z],deltaLB[z],deltaR[z],deltaRB[z]}

FourBeamDensityMatrix[x_]:={{1, delta[x]},{Conjugate@delta[x],0}}

FourBeamTransMatrix[eta_,mk_,mu_,alpha_,theta1_,theta2_,lambda_,omegav_:0]:=Module[{lenM},

	{
		{lambda+mu-2 alpha mu-eta omegav+mu Cos[2 theta1]+2 alpha mu Sin[theta1] Sin[theta2],alpha mu-alpha mu Cos[theta1-theta2],-mu-mu Cos[2 theta1],alpha mu+alpha mu Cos[theta1+theta2]}/Sin[theta1]+mk Cot[theta1]{1,0,0,0},
		{-mu+mu Cos[theta1] Cos[theta2]+mu Sin[theta1] Sin[theta2],lambda+2 mu-alpha mu+eta omegav-alpha mu Cos[2 theta2]-2 mu Sin[theta1] Sin[theta2],-mu-mu Cos[theta1] Cos[theta2]+mu Sin[theta1] Sin[theta2],alpha mu+alpha mu Cos[2 theta2]}/Sin[theta2]+mk Cot[theta2]{0,1,0,0},
		{-mu-mu Cos[2 theta1],alpha mu+alpha mu Cos[theta1] Cos[theta2]-alpha mu Sin[theta1] Sin[theta2],lambda+mu-2 alpha mu-eta omegav+mu Cos[2 theta1]+2 alpha mu Sin[theta1] Sin[theta2],alpha mu-alpha mu Cos[theta1] Cos[theta2]-alpha mu Sin[theta1] Sin[theta2]}/Sin[theta1]-mk Cot[theta1]{0,0,1,0},
		{-mu-mu Cos[theta1] Cos[theta2]+mu Sin[theta1] Sin[theta2],alpha mu+alpha mu Cos[2 theta2],-mu+mu Cos[theta1] Cos[theta2]+mu Sin[theta1] Sin[theta2],lambda+2 mu-alpha mu+eta omegav-alpha mu Cos[2 theta2]-2 mu Sin[theta1] Sin[theta2]}/Sin[theta2]-mk Cot[theta2]{0,0,0,1}
	}

]


FourBeamRotationMatrixSameAngle[alpha_]:=1/2{
	{1,-alpha,1,-alpha},
	{1,alpha,1,alpha},
	{1,-alpha,-1,alpha},
	{1,alpha,-1,-alpha}
};

FourBeamRotationMatrixDiffAngle={
	{1,0,1,0},
	{0,1,0,1},
	{-1,0,1,0},
	{0,-1,0,1}
};



End[]
(* Four Beams Some Definitions Equations END *)


(* Solving Fourth Order Equations BEGIN *)

(* Solving Fourth Order Equations END *)


EndPackage[]

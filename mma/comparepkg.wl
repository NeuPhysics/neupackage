(* ::Package:: *)

BeginPackage["neumat`"]
(* This package is only used to compare my results with other's people's. So it contains some transformations of conventions etc. *)


(* Sajad's Equations BEGIN *)

saxi::usage = "Sajad's definition of xi[a,b], which is 1-cos(\!\(\*SubscriptBox[\(\[Theta]\), \(a\)]\)-\!\(\*SubscriptBox[\(\[Theta]\), \(b\)]\)); \[Theta]'s should be the angle defined in the same polar coordinate system.";


saxi[theta1_,theta2_]:=1-Cos[theta1-theta2];



(* Sajad's Equations END *)


EndPackage[]

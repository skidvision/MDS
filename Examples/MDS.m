(* ::Package:: *)

(************************************************************************)
(* This file was generated automatically by the Mathematica front end.  *)
(* It contains Initialization cells from a Notebook file, which         *)
(* typically will have the same name as this file except ending in      *)
(* ".nb" instead of ".m".                                               *)
(*                                                                      *)
(* This file is intended to be loaded into the Mathematica kernel using *)
(* the package loading commands Get or Needs.  Doing so is equivalent   *)
(* to using the Evaluate Initialization Cells menu command in the front *)
(* end.                                                                 *)
(*                                                                      *)
(* DO NOT EDIT THIS FILE.  This entire file is regenerated              *)
(* automatically each time the parent Notebook file is saved in the     *)
(* Mathematica front end.  Any changes you make to this file will be    *)
(* overwritten.                                                         *)
(************************************************************************)



(* :Title: Multidimensional Scaling Tools *)

(* :Context: ExpTools`MDS` *)

(* :Author: Flip Phillips
	Recent modifications of $Date: 2007-12-12 23:10:45 -0500 (Wed, 12 Dec 2007) $ by $Author: flip $ *)

(* :Summary: 
	This package provides various signal detection theory function to Mathematica.
*)

(* :Package Version: $Revision: 25 $ *)

(* :Mathematica Version: 6.0 *)

(* :Copyright: Copyright 1999-2007, Flip Phillips, All Rights Reserved.  *)

(* :History: 
$Log: MDS.nb,v $
Revision 1.3  1999/03/17 03:36:00  cvs
Fixed constant function to be less of a memory hog

Revision 1.2  1999/03/16 19:15:44  cvs
fixed external creation of symbol that had list of exported functions

Revision 1.1  1999/03/16 14:46:27  cvs
first checkin since transfer from development dir

Revision 1.3  1999/03/15 23:45:47  cvs
migrated to Packages directory

Revision 1.2  1999/03/15 16:32:57  cvs
started integration into package form

*)

(* :Keywords:
	packages, path, keywords
*)

(* :Limitations:  *)

(* :Discussion:  *)


BeginPackage["ExpTools`MDS`"]


MDS::usage="MDS.m is a package which provides a few tools for doing multidimensional scaling."


MultidimensionalScaling::usage="MultidimensionalScaling[data,dimensions] takes the dissimilarity matrix in 'data' and creates an initial configuration in 'dimensions' dimensions"


MDSModel::usage="MDSModel[data,dimensions] represents the multidimensional scaling model."


IterateMDS::usage="IterateMDS[model,times] Iterates the MDS model 'times' times."


Begin["`Private`"]


Unprotect[{MultidimensionalScaling,IterateMDS,MDSModel}];


ZeroDiagonal[m_]:=Module[{mm},
mm=(1-IdentityMatrix[Dimensions[m][[1]]])m;
mm/StandardDeviation[Flatten[mm]]]


MakeSymmetric[m_]:=Module[{ms},
ms=m^2;
Sqrt[(ms+Transpose[ms])/2]]


MakeBstar[m_]:=Module[{m2,means,mmeans,n},n=Dimensions[m][[1]];
m2=m^2;
means=Map[Mean,m2];
mmeans=Table[means,{n}];
((m2-mmeans-Transpose[mmeans])+Mean[means])/-2
]


MakeCoordinates[bstar_]:=Module[{u,\[Gamma],v},
{u,\[Gamma],v}=SingularValues[bstar,Tolerance->0];
{ColumnDrop[Transpose[u].DiagonalMatrix[Sqrt[\[Gamma]]],-1],\[Gamma]}
]


MakeB[x_,disdata_]:=Module[{dm,mD,mB},
dm=MakeDistances[x];
(* add identity to keep div/0 from happening *)
mD=dm+IdentityMatrix[Length[dm]];
mB=-2(disdata/mD);
mB+DiagonalMatrix[-Map[Apply[Plus,#]&,mB]]
]


GuttmanTransform[x_,disdata_]:=Module[{bxdot,n},
n=Length[x];
1/(2n) (MakeB[x,disdata].x)]


length[v_]:=Sqrt[v.v];
distance[{p1_,p2_}]:=length[p1-p2];


MakeDistances[x_]:=Module[{n},
n=Length[x];
Table[distance[{x[[i]],x[[j]]}],{i,1,n},{j,1,n}]]


SubtractList[l_List]:=Fold[Subtract,First[l],Rest[l]]


AddConstantsFunction[m_]:=Module[{d,parts,lists},
d=Dimensions[m][[1]];
parts=Flatten[
Table[{
{a+1,(a+b+2)},
{(a+b+2),(a+b+c+3)},
{a+1,(a+b+c+3)}},
{a,0,d-3},{b,0,d-a-3},{c,0,d-a-b-3}],2];
lists=Map[Sort[{m[[#[[1,1]],#[[1,2]]]],m[[#[[2,1]],#[[2,2]]]],m[[#[[3,1]],#[[3,2]]]]},Greater]&,parts];

Max[Append[{0},Map[SubtractList,lists]]]
]


AddConstantsFunction[m_]:=Module[{d,a,b,c},
d=Dimensions[m][[1]];
Max[Append[{0},Flatten[
Table[SubtractList[Sort[{
m[[a+1,(a+b+2)]],
m[[(a+b+2),(a+b+c+3)]],
m[[a+1,(a+b+c+3)]]},Greater]],
{a,0,d-3},{b,0,d-a-3},{c,0,d-a-b-3}],2]]]
]


TransformByConstant[m_]:=Module[{constant},
constant=AddConstantsFunction[m];
ZeroDiagonal[m+constant]
]


CalcStress[x_,distdata_]:=Module[{dists},
dists=MakeDistances[x];Sqrt[(Plus@@Flatten[(distdata-dists)^2])/(Plus@@Flatten[(distdata)^2])]
]


MDS::Dimension="Invalid number of dimensions \"`1`\"";


Format[t_MDSModel]:="MDSModel[<>]"


MDSModel[m_,dimen_]:=Module[{bstar,xmatrix,umtx,tdata,stress,\[Gamma]},
If[dimen>=Dimensions[m][[1]],Message[MDS::Dimension,dimen];Return[{}]];
umtx=ZeroDiagonal[MakeSymmetric[m]];
tdata=TransformByConstant[umtx];
bstar=MakeBstar[umtx];
{xmatrix,\[Gamma]}=MakeCoordinates[N[bstar]];xmatrix=ColumnTake[xmatrix,{1,dimen}];
stress=CalcStress[xmatrix,tdata];

{xmatrix,tdata,\[Gamma],stress}]


MultidimensionalScaling[m_,d_Integer]:=Module[{dataM},
MDSModel[If[Length[m[[1,1]]]!=0,Map[Mean,m,{2}],m],d]]


IterateMDS[model_,times_Integer:1]:=Module[{newx,news},
Do[
newx=GuttmanTransform[model[[1]],model[[2]]];
news=CalcStress[newx,model[[2]]],
{times}];
{newx,model[[2]],model[[3]],news}]


End[]


Protect[{MultidimensionalScaling,IterateMDS,MDSModel}];


EndPackage[]
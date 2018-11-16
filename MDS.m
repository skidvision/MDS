(* ::Package:: *)

(* :Title: Multidimensional Scaling Tools *)

(* :Context: ExpTools`MDS` *)

(* :Author: Flip Phillips
	Recent modifications of $Date: 21Jun2007 $ by $Author: flip $ *)

(* :Summary: 
	This package provides various signal detection theory function to Mathematica.
*)

(* :Package Version: $Revision: 9 $ *)

(* :Mathematica Version: 6.0 *)

(* :Copyright: Copyright 1999-2007, Flip Phillips, All Rights Reserved.  *)

(* :History: 
$Log: MDS.m,v $
Revision 9.0  2007/06/20
Adapted for v6.0 compatability, may work w/ 5.2 also but not tested.

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
	MDS, Scaling, Multidimensional scaling
*)

(* :Limitations:  *)

(* :Discussion:  *)


(* ::Subsection:: *)
(*Package*)


BeginPackage["MDS`"]


(* ::Subsection:: *)
(*Messages*)


MDS::usage="MDS.m is a package which provides a few tools for doing \
multidimensional scaling."


MultidimensionalScaling::usage="MultidimensionalScaling[data,dimensions] \
takes the dissimilarity matrix in 'data' and creates an initial configuration \
in 'dimensions' dimensions."


MDSModel::usage="MDSModel[data,dimensions] represents the multidimensional \
scaling model."


IterateMDS::usage="IterateMDS[model,times] Iterates the MDS model 'times' \
times."


Begin["`Private`"]


Unprotect[{MultidimensionalScaling,IterateMDS,MDSModel}];


(* ::Subsection:: *)
(*Private functions*)


ZeroDiagonal[m_]:=Module[{mm},
    mm=(1-
            IdentityMatrix[
              Dimensions[m][[1]]])m;
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
    {u,\[Gamma],v}=SingularValueDecomposition[bstar,Tolerance->0];
    {Drop[Transpose[u].Sqrt[\[Gamma]],None,-1],\[Gamma]}
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
    Table[
      distance[{x[[i]],
          x[[j]]}],{i,1,n},{j,1,n}]]


SubtractList[l_List]:=Fold[Subtract,First[l],Rest[l]]


AddConstantsFunction[m_]:=Module[{d,a,b,c},
    d=Dimensions[m][[1]];
    Max[Append[{0},Flatten[
          Table[SubtractList[Sort[{
                  m[[a+1,(a+b+2)]],
                  
                  m[[(a+b+2),(a+b+
                        c+3)]],
                  m[[a+1,(a+b+c+3)]]},
                Greater]],
            {a,0,d-3},{b,0,d-a-3},{c,0,d-a-b-3}],2]]]
    ]


TransformByConstant[m_]:=Module[{constant},
    constant=AddConstantsFunction[m];
    ZeroDiagonal[m+constant]
    ]


CalcStress[x_,distdata_]:=Module[{dists},
dists=MakeDistances[x];Sqrt[(Plus@@Flatten[(distdata-dists)^2])/(Plus@@Flatten[(distdata)^2])]
]


(* ::Subsection:: *)
(*Public functions*)


MDS::Dimension="Invalid number of dimensions \"`1`\"";


Format[t_MDSModel]:="MDSModel[<>]"


(* ::Text:: *)
(*The N[m] saves us a lot of grief herem*)


MDSModel[m_,dimen_]:=Module[{bstar,xmatrix,umtx,tdata,stress,\[Gamma]},
    If[dimen>=
        Dimensions[m][[1]],
      Message[MDS::Dimension,dimen];Return[{}]];
    umtx=ZeroDiagonal[MakeSymmetric[N[m]]];
    tdata=TransformByConstant[umtx];
    bstar=MakeBstar[umtx];
    {xmatrix,\[Gamma]}=MakeCoordinates[N[bstar]];
    xmatrix=Take[xmatrix,All,{1,dimen}];
    stress=CalcStress[xmatrix,tdata];
    
    {xmatrix,tdata,\[Gamma],stress}]


MultidimensionalScaling[m_,d_Integer]:=Module[{dataM},
    MDSModel[
      If[Length[m[[1,1]]]!=0,
        Map[Mean,m,{2}],m],d]]


IterateMDS[model_,times_Integer:1]:=Module[{newx,news},
    Do[
      newx=
        GuttmanTransform[model[[1]],
          model[[2]]];
      news=CalcStress[newx,model[[2]]],
      {times}];
    {newx,model[[2]],
      model[[3]],news}]


End[]


Protect[{MultidimensionalScaling,IterateMDS,MDSModel}];


EndPackage[]

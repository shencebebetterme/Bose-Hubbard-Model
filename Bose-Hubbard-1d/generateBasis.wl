(* ::Package:: *)

nS=ToExpression@$ScriptCommandLine[[2]];
nP=ToExpression@$ScriptCommandLine[[3]];
data=Flatten[Permutations/@(PadRight[#,nS]&/@IntegerPartitions[nP,nS]),1];
Export["data.h5",data];
Print/@Dimensions[data];

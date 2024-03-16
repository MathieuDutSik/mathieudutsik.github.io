TheGramMat:=[
[5,0,0],
[0,1,0],
[0,0,1]];

U:=rec(ListCosets:=[
[1,0  ,0  ,  0],
[1,1/4,1/2,  0],
[1,1/2,1/2,1/2],
[1,3/4,0,1/2]],
GramMat:=TheGramMat);
 

LFC:=Periodic_NonirreducibleDelaunayComputationStandardFunctions(U);

LFlagFunc:=LFC.FlagFunctions();

TheDelaneySymb:=LFlagFunc.GetDelaneySymbol();
TheFileDelaney:="SolveS14twisted.ds";
RemoveFileIfExist(TheFileDelaney);
output:=OutputTextFile(TheFileDelaney, true);
eString:=Concatenation("<1", TheDelaneySymb.eString, "\n");
WriteAll(output, eString);
CloseStream(output);

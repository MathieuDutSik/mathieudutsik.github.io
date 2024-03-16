function [TheV, test]=TD_NewtonMethod(RadiusVector, TheMat, ListDefect)
[vce, nbv]=size(TheMat);

MatrixDiff=TD_DifferentialMatrix(RadiusVector, TheMat);
for i=1:nbv
    VCT(1,i)=1;
end;
MatrixDiff(nbv,:)=VCT;

TheVect=ListDefect;
TheVect(1:nbv)=0;

TheCandidate=RadiusVector-TheVect*inv(MatrixDiff);
ListD=TD_EvaluatePosition(TheCandidate, TheMat);
SQR1=TD_SquareDefect(ListDefect);
SQR2=TD_SquareDefect(ListD);
if (SQR2<SQR1)
    test=1;
    TheV=TheCandidate;
else
    test=0;
    TheV=RadiusVector;
end;

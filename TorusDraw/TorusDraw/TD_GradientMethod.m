function [TheV, test]=TD_GradientMethod(RadiusVector, TheMat, ListDefect)
delta=100;
[vce, nbv]=size(TheMat);

MatrixDiff=TD_DifferentialMatrix(RadiusVector, TheMat);
evol(1:nbv)=0;
for i=1:nbv
    evol=evol+ListDefect(i)*MatrixDiff(i,:);
end;

test=0;
maxite=20;
nbite=0;
while(nbite<maxite)
    Vect=RadiusVector;
    nbmove=0;
    while (1==1)
        Norm=TD_SquareDefect(TD_EvaluatePosition(Vect, TheMat));
        NewVect=Vect+delta*evol;
        NewNorm=TD_SquareDefect(TD_EvaluatePosition(NewVect, TheMat));
        if (NewNorm>Norm),
            break;
        end;
        Vect=NewVect;
        nbmove=nbmove+1;
    end;
    if (nbmove>0)
        test=1;
        break;
    end;
    delta=delta/2;
    nbite=nbite+1;    
end;
TheV=Vect;

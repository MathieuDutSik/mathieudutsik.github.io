function [TheV, test]=TD_AllDirections_skull(ListRadius, TheMat, ListDefect, idx, sgn)
% this function make movements in direction idx with sign sgn
delta=0.01;
[vce, nbv]=size(TheMat);
for i=1:nbv
    if (i == idx),
        evol(i)=sgn;
    else
        evol(i)=0;
    end;
end;

test=0;
maxite=10;
nbite=0;
while(nbite<maxite)
    Vect=ListRadius;
    EVOL(1:nbv)=1;
    EVOL=EVOL+delta*evol;
    nbmove=0;
    while (1==1),
        Norm=TD_SquareDefect(TD_EvaluatePosition(Vect, TheMat));
        for i=1:nbv
            NewVect(i)=Vect(i)*EVOL(i);
        end;
        NewNorm=TD_SquareDefect(TD_EvaluatePosition(NewVect, TheMat));
        if (NewNorm>Norm),
            break;
        end;
        Vect=NewVect;
        nbmove=nbmove+1;
    end;
    if (nbmove>0),
        test=1;
        break;
    end;
    delta=delta/2;
    nbite=nbite+1;
end;
TheV=Vect;
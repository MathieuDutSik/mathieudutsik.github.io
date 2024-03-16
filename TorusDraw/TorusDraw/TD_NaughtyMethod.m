function [TheV, test]=TD_NaughtyMethod(ListRadius, TheMat, ListDefect)
% this function make movements to decrease the value of the penalty
% function. By examining the equations one sees that it should
% always work and always make an improvement.
delta=0.01;
[vce, nbv]=size(TheMat);
for i=1:nbv
    if (ListDefect(i)>0),
        evol(i)=-1;
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
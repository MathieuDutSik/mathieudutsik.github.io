function [ListRadius, test]=TD_SeveralMinimization(TheMat, minimal)
[hce, m]=size(TheMat);

ListRadius(1:m)=1;
while(1==1)
    ListDefect=TD_EvaluatePosition(ListRadius, TheMat);
    SQR=TD_SquareDefect(ListDefect)
    if (SQR< minimal),
        break;
    end;
    [TheV, test]=TD_NaughtyMethod(ListRadius, TheMat, ListDefect);
    if (test == 0)
        [TheV, test]=TD_AllDirections(ListRadius, TheMat, ListDefect);
        if (test == 0)
            [TheV, test]=TD_NewtonMethod(ListRadius, TheMat, ListDefect);
        end;
    end;
    if (test == 0)
        disp('We failed to find the solution of the equation set');
        return;
    end;
    them=min(TheV);
    ListRadius=TheV/them;
end;

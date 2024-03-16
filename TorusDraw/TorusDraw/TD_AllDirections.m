function [TheV, test]=TD_AllDirections(ListRadius, TheMat, ListDefect)


[vce, nbv]=size(TheMat);
SQRfound=nbv;
TheV=0;
test=0;
for idx=1:nbv
    for u=1:2
        if (u == 1)
            sgn=1;
        else
            sgn=-1;
        end;
        [TheVcall, testcall]=TD_AllDirections_skull(ListRadius, TheMat, ListDefect, idx, sgn);
        if (testcall == 1)
            test=1;
            ListDefect=TD_EvaluatePosition(ListRadius, TheMat);
            SQR=TD_SquareDefect(ListDefect);
            if (SQR < SQRfound)
                SQRfound=SQR;
                TheV=TheVcall;
            end;
        end;
        
        
        
        
    end;
    
    
    
end;


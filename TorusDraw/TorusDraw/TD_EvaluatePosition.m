function ListPos=TD_EvaluatePosition(ListRadius, TheMat)
[hce, m]=size(TheMat);
for i=1:m
    H=pi;
    for j=1:m
        if (TheMat(i,j)==1),
            H=H-atan(ListRadius(j)/ListRadius(i));
        end;
    end;
    ListPos(i)=H;
end;
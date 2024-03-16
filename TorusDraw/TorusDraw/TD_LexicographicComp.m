function result=TD_LexicographicComp(V1, V2)
[vce, nb]=size(V1);

if (V1(1) > V2(1))
    result=-1;
    return;
end;

if (V1(1) < V2(1))
    result=1;
    return;
end;

if (nb == 1)
    result=0;
    return;
end;

V1red=V1(2:nb);
V2red=V2(2:nb);
result=TD_LexicographicComp(V1red, V2red);
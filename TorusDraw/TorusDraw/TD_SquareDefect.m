function defect=TD_SquareDefect(theVect)
defect=0;
[hce, m]=size(theVect);
for i=1:m
    defect=defect+theVect(i)*theVect(i);
end;
function ListC=TD_ListCoord(RadiusVector, SPA, plg, maxval)
[x,m]=size(RadiusVector);
ListStatus(m)=0;
ListStatus(1)=1;

ListCoord(m,3)=0;
while(1==1)
    for iEdge=1:size(SPA)
        u=SPA(iEdge, 1);
        v=SPA(iEdge, 2);
        if ((ListStatus(u) == 1) & (ListStatus(v) == 0)),
%            HCKLK=0;
%            disp('before HCKLK='); disp(HCKLK);
            Coord=TD_FindCoordinates(ListCoord(u,:), u, v, RadiusVector, plg, maxval);
%            disp('after  HCKLK='); disp(HCKLK);
            ListCoord(v,:)=Coord;
            ListStatus(v)=1;
        end;
        if ((ListStatus(u) == 0) & (ListStatus(v) == 1)),
            Coord=TD_FindCoordinates(ListCoord(v,:), v, u, RadiusVector, plg, maxval);
            ListCoord(u,:)=Coord;
            ListStatus(u)=1;
        end;
    end;
%    disp('End of considerations m=');
%    disp(m);
    nbUndone=0;
    for i=1:m
%        disp(ListStatus);
        if (ListStatus(i) == 0),
            nbUndone=nbUndone+1;
        end;
    end;    
    if (nbUndone == 0),
        break;
    end;
end;
ListC=ListCoord;

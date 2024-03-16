 function Coordj=TD_FindCoordinates(Coordi, i, j, RadiusVector, plg, maxval)
% the vector Coord contains
% --first the x coordinate of this point
% --second the z coordinate of this number of the first adjacent vertex
% --third the angle with respect to the first adjacent vertex
% j should be adjacent to i


HCKLK=1;
% Build the angles
%disp('i=');
%disp(i);
%disp('j=');
%disp(j);
valence=0;
for k=1:maxval
    vadj=plg(i, k);
    if (vadj>0),
        valence=valence+1;
    end;
end;
Ladj=plg(i,1:valence);

ListAngle(1)=Coordi(3);
for k=2:valence
    u1=Ladj(k-1);
    u2=Ladj(k);
    angle=atan(RadiusVector(u1)/RadiusVector(i));
    angle=angle+atan(RadiusVector(u2)/RadiusVector(i));
    ListAngle(k)=ListAngle(k-1)+angle;
end;


%disp('Ladj=');
%disp(Ladj);
for k=1:valence
    if (Ladj(k) == j),
        posj=k;
    end;
end;
TheAngle=ListAngle(posj);
%dist=RadiusVector(i)+RadiusVector(j);
dist=sqrt(RadiusVector(i)*RadiusVector(i)+RadiusVector(j)*RadiusVector(j));
x=Coordi(1)+dist*cos(TheAngle);
y=Coordi(2)+dist*sin(TheAngle);

valence2=0;
for k=1:maxval
    vadj=plg(j,k);
    if (vadj>0),
        valence2=valence2+1;
    end;
end;
Ladj2=plg(j, 1:valence2);

ListAngle2(1)=0;
for k=2:valence2
    u1=Ladj2(k-1);
    u2=Ladj2(k);
    angle=atan(RadiusVector(u1)/RadiusVector(j));
    angle=angle+atan(RadiusVector(u2)/RadiusVector(j));
    ListAngle2(k)=ListAngle2(k-1)+angle;
end;
for k=1:valence2
    if (Ladj2(k) == i),
        posi=k;
    end;
end;
delta=pi+ListAngle(posj)-ListAngle2(posi);
for k=1:valence2
    ListAngle2(k)=ListAngle2(k)+delta;
end;
%disp('End of TD_FindCoordinates');
Coordj=[x,y,ListAngle2(1)];

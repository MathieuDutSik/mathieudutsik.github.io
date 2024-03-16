function test=TD_TheActualDrawing(TheBasis, RadiusVector, ListPoint, TheWindow, ListCoord, ComponentVector, Select, simiInfo, plg, ListRed)
[nbv, maxval]=size(plg);
[nbPoint, garb]=size(ListPoint);

simiVector=[simiInfo(1) simiInfo(2) 0 0];
V1=TheBasis(1,:);
V2=TheBasis(2,:);

Vect1=TD_similitude(V1, simiVector);
Vect2=TD_similitude(V2, simiVector);
HInv=inv(TheBasis);

mult=sqrt(simiInfo(1)*simiInfo(1)+simiInfo(2)*simiInfo(2));
for iVert=1:nbv
    nb=0;
    for i=1:maxval
        if (plg(iVert, i) >0)
            nb=nb+1;
        end;
    end;
    VectorValency(iVert)=nb;
end;



for iPoint=1:nbPoint
    iVert=ListPoint(iPoint, 1);
    i=ListPoint(iPoint, 2);
    j=ListPoint(iPoint, 3);
    izPos=ListCoord(iVert, 1:2);
    izPos=izPos+i*V1+j*V2;
    ifImg=TD_similitude(izPos, simiInfo);
    Coordi=[izPos(1), izPos(2), ListCoord(iVert, 3)];
    if (ComponentVector(iVert) == Select)
        for i=1:VectorValency(iVert)
            jVert=plg(iVert, i);
            Coordj=TD_FindCoordinates(Coordi, iVert, jVert, RadiusVector, plg, maxval);
            eDiff=Coordj-ListCoord(jVert, :);
            eDiff=eDiff(1:2);
            U=eDiff*HInv;
            Urnd=[round(U(1)), round(U(2))];
            jzPos=ListCoord(jVert, 1:2);
            jzPos=jzPos+Urnd(1)*V1+Urnd(2)*V2;
            jfImg=TD_similitude(jzPos, simiInfo);
            Coordj=[jzPos(1) jzPos(2) ListCoord(jVert, 3)];
            if (ListRed(3) == 1)
                line([izPos(1) jzPos(1)], [izPos(2) jzPos(2)],'Color',[0,0,0]);
            end;
            for val=1:VectorValency(jVert)
                if (plg(jVert, val) == iVert)
                    irev=val;
                    if (val < VectorValency(jVert))
                        irevNext=val+1;
                    else
                        irevNext=1;
                    end;
                end;
            end;
            kVert=plg(jVert, irevNext);
            Coordk=TD_FindCoordinates(Coordj, jVert, kVert, RadiusVector, plg, maxval);
            eDiff=Coordk-ListCoord(kVert, :);
            eDiff=eDiff(1:2);
            U=eDiff*HInv;
            Urnd=[round(U(1)), round(U(2))];
            kzPos=ListCoord(kVert, 1:2);
            kzPos=kzPos+Urnd(1)*V1+Urnd(2)*V2;
            kfImg=TD_similitude(kzPos, simiInfo);
            eVect=kzPos-izPos;
%            dist=sqrt(eVect(1)*eVect(1)+eVect(2)*eVect(2));
%            DELTA=abs(dist-RadiusVector(iVert)-RadiusVector(kVert));
            V1test=[ListPoint(iPoint, 1) i j];
            V2test=[kVert Urnd(1) Urnd(2)];
            test3=TD_LexicographicComp(V1test, V2test);
            if (ListRed(1) == 1 & test3 == 1)
                line([ifImg(1) kfImg(1)], [ifImg(2) kfImg(2)],'Color',[0,0,0]);
            end;
        end;
    end;
    if (ListRed(2) == 1 & ComponentVector(iVert) == Select),
        vce=TD_circle2(ifImg(1), ifImg(2), mult*RadiusVector(iVert));
    end;
end;
test=1;
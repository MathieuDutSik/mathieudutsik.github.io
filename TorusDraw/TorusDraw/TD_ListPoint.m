function ListPoint=TD_ListPoint(TheBasis, TheWindow, RadiusVector, simiInfo, ListCoord)
% homothety is of the form [a b c d] and correspond to the operation
% z --> (a+bi)z+  (c+di)
[nbv, uce]=size(ListCoord);
MaxRad=max(RadiusVector);
%disp('iVert='); disp(iVert)

simiVector=[simiInfo(1) simiInfo(2) 0 0];
Vect1=TD_similitude(TheBasis(1,:), simiVector);
Vect2=TD_similitude(TheBasis(2,:), simiVector);

Amat(1,:)=Vect1;
Amat(2,:)=Vect2;
Ainv=inv(Amat);

alpha11=Ainv(1,1);
alpha12=Ainv(1,2);
alpha21=Ainv(2,1);
alpha22=Ainv(2,2);

nbPosition=0;
for iVert=1:nbv
    eCoord=ListCoord(iVert, 1:2);
    fImg=TD_similitude(eCoord, simiInfo);
    Bm1=TheWindow(1)-fImg(1)-MaxRad;
    BM1=TheWindow(2)-fImg(1)+MaxRad;
    Bm2=TheWindow(3)-fImg(2)-MaxRad;
    BM2=TheWindow(4)-fImg(2)+MaxRad;
    if (alpha11>0),
        Bound1=[alpha11*Bm1, alpha11*BM1];
    else
        Bound1=[alpha11*BM1, alpha11*Bm1];
    end;
    if (alpha21>0),
        Bound2=[alpha21*Bm2, alpha21*BM2];
    else
        Bound2=[alpha21*BM2, alpha21*Bm2];
    end;
    BoundSi=Bound1+Bound2;

    if (alpha12>0),
        Bound1=[alpha12*Bm1, alpha12*BM1];
    else
        Bound1=[alpha12*BM1, alpha12*Bm1];
    end;
    if (alpha22>0),
        Bound2=[alpha22*Bm2, alpha22*BM2];
    else
        Bound2=[alpha22*BM2, alpha22*Bm2];
    end;
    BoundSj=Bound1+Bound2;
    MINi=floor(BoundSi(1));
    MAXi=ceil(BoundSi(2));
    MINj=floor(BoundSj(1));
    MAXj=ceil(BoundSj(2));
    nbcase=(MAXi+1-MINi)*(MAXj+1-MINj);
%    disp('nb case considered='); disp(nbcase);
    for i=MINi:MAXi
        for j=MINj:MAXj
            Vvect=fImg+i*Vect1+j*Vect2;
            test1=0;
            test2=0;
            if (Vvect(1)>TheWindow(1)-MaxRad & Vvect(1)<TheWindow(2)+MaxRad & Vvect(2)>TheWindow(3)-MaxRad & Vvect(2)<TheWindow(4)+MaxRad),
                test1=1;
            end;
            if (test1 == 1)
                if (Vvect(1)>TheWindow(1) & Vvect(1)<TheWindow(2) & Vvect(2)>TheWindow(3) & Vvect(2)<TheWindow(4))
                    test2=1;
                end;
            end;
            if (test1 == 1)
                nbPosition=nbPosition+1;
                ListPoint(nbPosition,:)=[iVert i j test2];
            end;
        end;
    end;
    
    
end;


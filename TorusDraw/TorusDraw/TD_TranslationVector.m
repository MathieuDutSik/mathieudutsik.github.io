function TheBasis=TD_TranslationVector(RadiusVector, ListCoord, plg, tol)
% in this procedure we find a rank 2 translation lattice
% the basis found is supposed to be optimized
% but no guarantee is given
[nbv, maxval]=size(plg);
nbvector=0;

for iVert=1:nbv
    for iAdj=1:maxval
        val=plg(iVert, iAdj);
        if (val >0),
            Coord=TD_FindCoordinates(ListCoord(iVert,:), iVert, val, RadiusVector, plg, maxval);
            eDiff=ListCoord(val,:)-Coord;
%            disp('Found');
%            disp(eDiff);
            eTrans=eDiff(1:2);
            eNorm=abs(eTrans(1))+abs(eTrans(2));
            if (eNorm> tol),
                test=1;
                for iVect=1:nbvector
                    U=eTrans-ListV(iVect,:);
                    eNorm=abs(U(1))+abs(U(2));
                    if (eNorm<tol),
                        test=0;
                    end;
                end;
                if (test == 1),
                    nbvector=nbvector+1;
                    ListV(nbvector,:)=eTrans;
                end;
            end;
        end;
    end;
end;

eScalChoosen=5;
for i=1:nbvector-1
    for j=i+1:nbvector
        U1=ListV(i,:);
        U2=ListV(j,:);
        eN1=U1(1)*U1(1)+U1(2)*U1(2);
        eN2=U2(1)*U2(1)+U2(2)*U2(2);
        eScal=U1(1)*U2(1)+U1(2)*U2(2);
        eDefect=(eScal*eScal)/(eN1*eN2);
        TheMat(1,:)=U1;
        TheMat(2,:)=U2;
        eDet=det(TheMat);
        HInv=inv(TheMat);
        test=1;
        for k=1:nbvector
            eCoord=ListV(k,:)*HInv;
            delta=abs(eCoord(1)-round(eCoord(1)))+abs(eCoord(2)-round(eCoord(2)));
            if (delta>tol),
                test=0;
            end;
%            disp(eCoord);
        end;
        if (test == 1 & eDefect<eScalChoosen),
            TheChoice=[i j];
            eScalChoosen=eDefect;
%            disp('We make another choice');
        end;
        
    end;
end;
TheBasis(1,:)=ListV(TheChoice(1), :);
TheBasis(2,:)=ListV(TheChoice(2), :);

function test=TD_FundamentalDomain(TheBasis, ePoint, simiInfo, MatUnimodular)
VNEW=MatUnimodular*TheBasis;
V1=VNEW(1,:);
V2=VNEW(2,:);


ePt1=ePoint;
ePt2=ePoint+V1;
ePt3=ePoint+V1+V2;
ePt4=ePoint+V2;

fImg1=TD_similitude(ePt1, simiInfo);
fImg2=TD_similitude(ePt2, simiInfo);
fImg3=TD_similitude(ePt3, simiInfo);
fImg4=TD_similitude(ePt4, simiInfo);

line([fImg1(1) fImg2(1)], [fImg1(2) fImg2(2)]);
line([fImg2(1) fImg3(1)], [fImg2(2) fImg3(2)]);
line([fImg3(1) fImg4(1)], [fImg3(2) fImg4(2)]);
line([fImg4(1) fImg1(1)], [fImg4(2) fImg1(2)]);
test=1;
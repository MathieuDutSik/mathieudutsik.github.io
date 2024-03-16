function fPoint=TD_similitude(ePoint, simiInfo)
x=ePoint(1);
y=ePoint(2);
a=simiInfo(1);
b=simiInfo(2);
c=simiInfo(3);
d=simiInfo(4);

fPoint=[a*x-b*y+c, a*y+b*x+d];



function MatrixDiff=TD_DifferentialMatrix(RadiusVector, TheMat)
[vc, nbv]=size(RadiusVector);
MatrixDiff(nbv,nbv)=0;

for v=1:nbv
    for u=1:nbv
        if (TheMat(v,u) == 1)
            SQV=RadiusVector(u)*RadiusVector(u)+RadiusVector(v)*RadiusVector(v);
            MatrixDiff(v, u)=MatrixDiff(v, u)-RadiusVector(v)/SQV;
            MatrixDiff(v, v)=MatrixDiff(v)+RadiusVector(u)/SQV;
        end;
    end;
end;
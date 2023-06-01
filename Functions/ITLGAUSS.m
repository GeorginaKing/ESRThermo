function [out] = ITLGAUSS(beta03, t)

global isoT measESR
mat=nan(length(measESR),length(isoT));
mat(1:length(t))=t;

s=10.^beta03(1);
Et=beta03(2);
Eu=beta03(3); 

T=isoT; 

out=[];

for j=1:length(mat(1,:));
    ok=isfinite(mat(:,j)); time=mat(:,j);
    out=[out; ThermaldecayGauss(time(ok),T(j),Et,Eu,s)];
end
out = out';

function [out] = ITLGAUSS(beta, t)

global isoT measL
mat=nan(length(measL),length(isoT));
mat(1:length(t))=t;

s=10.^beta(1);
Et=beta(2);
Eu=beta(3);
A=beta(4:end); 

T=isoT; 

out=[];

for j=1:length(mat(1,:));
    ok=isfinite(mat(:,j)); time=mat(:,j);
    out=[out; ThermaldecayGauss(time(ok),T(j),Et,Eu,s)];
end
out = out';

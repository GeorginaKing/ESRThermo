function [out] = ITLGOK(beta03, x)

global isoT measESR %initialCondition
mat=nan(length(measESR),length(isoT));
mat(1:length(x))=x;

s=10.^beta03(1); 
E=beta03(2);
b=beta03(3);
% A=beta03(4:end);

T=isoT; 

out=[];

kB=8.617343e-5; %eV/K

for j=1:length(mat(1,:));
    ok=isfinite(mat(:,j));
    time=mat(:,j); 
    
    out=[out; (1-(1-b).*(s.*exp(-E./((T(j)+273.15).*kB))).*time(ok)).^(1./(1-b))]; 
   
end
out = out';



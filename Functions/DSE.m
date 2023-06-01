function out = DSE(beta0,x)

Imax1 = beta0(1);
D01 = beta0(2);
Imax2 = beta0(3);
D02 = beta0(4);

Fun = Imax1.*(1-exp(-(x)/D01))+Imax2.*(1-exp(-(x)/D02));
Fun(isnan(Fun))=0; 

out=[Fun];
end


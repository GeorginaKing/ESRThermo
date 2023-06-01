function out = EXPLIN(beta0,x)

Imax = beta0(1);
D0 = beta0(2);
m = beta0(3);

Fun = (Imax.*(1-exp(-(x)/D0)+m.*(x))); 
Fun(isnan(Fun))=0; 


out=[Fun];
end
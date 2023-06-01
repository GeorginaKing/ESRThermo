%%Fits SAR data to a GOK function%%
%%georgina.king@geo.unibe.ch%%

function out = GOK(beta0, x);
a=beta0(1); D0=beta0(2); order=beta0(3);
% order(order>=2)=2;
Growth = a*(1-(1+(x./D0)*order).^-1/order);
%Growth = a*((1+(x./D0)*alpha).^alpha);
Growth(isnan(Growth))=0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=[Growth];


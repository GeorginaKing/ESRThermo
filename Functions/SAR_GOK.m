%% Fits SAR data to a GOK function %%

% Georgina King [georgina.king@unil.ch], August 2015


function out = SAR_GOK(beta0,t)

a=beta0(1);
D0=beta0(2);
order=beta0(3);

%%% Compute the growth matrix
Growth = a*(1-(1+(t./D0)*order).^-1/order);                                 % accumulation of ESR signal

Growth(isnan(Growth))=0;

out = Growth;
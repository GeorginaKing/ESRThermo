%% Fits SAR data to a single saturating exponential %%
%%% Solves for all aliquots at once
%%% Kinetic parameters are extracted from the fit

% Georgina King [georgina.king@unil.ch], August 2015


function out = SAR_SSE(beta0,t)

a=beta0(1); 
D0=beta0(2); 

%%% Compute the growth matrix
Growth = a.*(1-exp(-(t/D0)));                                               % accumulation of ESR signal

Growth(isnan(Growth))=0;                                                    % removes NaN data

out = Growth;
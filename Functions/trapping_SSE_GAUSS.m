function [nNf,nNtest] = trapping_SSE_GAUSS(time,temp,kparams);

Ma = 3600*24*365*1e6;
% Extract parameters
ddot = kparams.natDdot(1); %convert from Gy/s to Gy/Ma
D0 = kparams.D0(1);
s = 10.^kparams.s10_GAUSS(1);
Et = kparams.Et_GAUSS(1);
sigmaEt = kparams.sigmaEt_GAUSS(1);

% Define constants
kb = 8.617343e-5; 
magic_ratio = ddot/D0;
dEa=0.01; 
nstep = length(time);
dt = ((max(time)-min(time))/(nstep-1))*Ma;

% Create variables for GAUSS model Lambert et al (Submitted)
Ea=[(5/nstep):(5/nstep):5]'; nEa=length(Ea);
pEa = exp(-0.5*((Ea-Et)./sigmaEt).^2)/(sigmaEt*sqrt(2*pi));
npEa = sum(pEa);
T = temp+273.15;

% computes nN for the random Tt path
nN = zeros(nEa,nstep);
nNf = zeros(1,nstep);
for i = 2:nstep
	inv_tauth = s*exp(-Ea./(kb.*T(i-1)))*ones(1,nEa);
    alpha = magic_ratio+inv_tauth;
	nN(:,i) = (nN(:,i-1)+magic_ratio.*dt).*(1./(1+dt.*alpha(:,i-1)));
	nNf(i) = pEa'*nN(:,i);
end
nNf = nNf./npEa;


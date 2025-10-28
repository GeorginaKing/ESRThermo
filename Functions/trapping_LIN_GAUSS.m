function [nNf] = trapping_LIN_GAUSS(time,temp,kparams)


%%% Extract parameters
ddot = kparams.natDdot(1); 
s = 10.^kparams.s10(1);
Et = kparams.Et(1);
sigmaEt = kparams.sigmaEt(1);

%%% Define constants
Ma = 3600*24*365.25*1e6;
kb = 8.617343e-5;
magic_ratio = ddot;
nstep = length(time);
dt = ((max(time)-min(time))/(nstep-1))*Ma;

%%% Create variables for GAUSS model Lambert et al. (submitted)
Ea=(5/nstep):(5/nstep):5; 
nEa=length(Ea);
pEa = exp(-0.5*((Ea-Et)./sigmaEt).^2)/(sigmaEt*sqrt(2*pi));
npEa = sum(pEa);
T = (temp+273.15)';

%%% Computes Dose (called nN) for the random Tt path (Euler integration method)
nN = zeros(nEa,nstep);
nNf = zeros(1,nstep);

%%% Loop through the time steps
for i=2:nstep

    %%% Compute inv_tauth for the current time step (nEa x nstep)
	inv_tauth = (s*exp(-(Ea)'./(kb.*T(i-1))))*ones(1,nstep);

    xkd(:,i)=-1*inv_tauth(:,i-1).*nN(:,i-1).^(1-1);
    alpha(:,i) = magic_ratio-inv_tauth(:,i-1).*(nN(:,i-1));

	nN(:,i) = nN(:,i-1)+dt*alpha(:,i-1)./(1-dt.*xkd(:,i-1));
	nNf(i) = pEa*nN(:,i);
end

%%% Normalize nNf by the sum of the probabilities of Ea
nNf = nNf./npEa;
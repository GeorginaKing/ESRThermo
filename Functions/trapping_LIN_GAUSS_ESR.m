 function [Dose] = trapping_LIN_GAUSS_ESR(time,temp,kparams);

Ma = 3600*24*365*1e6; %s
% Extract parameters
ddot = kparams.natDdot(1); 
s = 10.^kparams.s10_GAUSS(1);
Et = kparams.Et_GAUSS(1);
sigmaEt = kparams.sigmaEt_GAUSS(1);

% Define constants
kb = 8.617343e-5; %eV
magic_ratio = ddot; %Gy s-1/Ma
dEa=0.01; 
nstep = length(time);
dt = ((max(time)-min(time))/(nstep-1))*Ma;

% Create variables for GAUSS model Lambert et al (Submitted)
Ea=[(5/nstep):(5/nstep):5]'; nEa=length(Ea);
pEa = exp(-0.5*((Ea-Et)./sigmaEt).^2)/(sigmaEt*sqrt(2*pi));
npEa = sum(pEa);
T = temp+273.15;

% computes Dose (called nN) for the random Tt path (Euler integration method)
nN = zeros(nEa,nstep);
nNf = zeros(1,nstep);
for i = 2:nstep
	inv_tauth = s*exp(-Ea./(kb.*T(i-1)))*ones(1,nEa);
    xkd(:,i)=-1*inv_tauth(:,i-1).*nN(:,i-1).^(1-1);
    alpha(:,i) = magic_ratio-inv_tauth(:,i-1).*(nN(:,i-1));
	nN(:,i) = nN(:,i-1)+dt*alpha(:,i-1)./(1-dt.*xkd(:,i-1));
	nNf(i) = pEa'*nN(:,i);
end
Dose = nNf./npEa;






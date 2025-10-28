function [nNf] = trapping_SSE_GAUSS(time,temp,kparams)

%%% Extract parameters
ddot = kparams.natDdot(1);                                                  % environmental radiation dose rate [Gy/s]
D0 = kparams.D0(1);                                                         % characteristic dose of saturation [Gy]
s = 10.^kparams.s10(1);                                                     % thermal frequency factor [s-1]
Et = kparams.Et(1);                                                         % activation energy = trap depth below the conduction band [eV]
sigmaEt = kparams.sigmaEt(1);

%%% Define constants
kb = 8.617343e-5;                                                           % Boltzmann constant [eV/K] 
magic_ratio = ddot/D0;
Ma = 1e6*365.25*24*3600;                                                    % 1 [Myr] in [s]
nstep = length(time);
dt = ((max(time)-min(time))/(nstep-1))*Ma;

%%% Create variables for detrapping GAUSS model Lambert et al. (Submitted)
Ea = (5/nstep):(5/nstep):5; 
nEa = length(Ea);
pEa = exp(-0.5*((Ea-Et)./sigmaEt).^2)/(sigmaEt*sqrt(2*pi));
npEa = sum(pEa);
T = temp+273.15;                                                            % transforms temperatures from [°C] to [K]

%%% Precompute nN for the random Tt path
nN = zeros(nEa,nstep);
nNf = zeros(1,nstep);

%%% Loop through the time steps
for i = 2:nstep

    %%% Compute inv_tauth for the current time step (nEa x nstep)
	inv_tauth = (s*exp(-(Ea)'./(kb.*T(i-1))))*ones(1,nstep);

    %%% Compute alpha for all combinations of nstep and Ea for the current time step (nEa x nstep)
    alpha = magic_ratio+inv_tauth;

    %%% Update nN for this time step using vectorized operations (nEa x nstep)
	nN(:,i) = (nN(:,i-1)+magic_ratio.*dt).*(1./(1+dt.*alpha(:,i-1)));

    %%% Compute nNf for this time step using matrix operations ([1 x nEa] * [nEa x nstep] > 1 value of nNf)
	nNf(i) = pEa*nN(:,i);
end

%%% Normalize nNf by the sum of the probabilities of Ea
nNf = nNf./npEa;
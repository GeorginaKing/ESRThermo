%% Model ESR trapping with GOK and GOK %%
%%% Assumes a General order kinetic (GOK) for trapping. GOK assumes
%%% trapping slower than SSE, due to e.g., Coulomb repulsive forces. Also,
%%% the effective trap lifetime varies with time (Guralnik et al., 2015).
%%% Thermal detrapping assumes a general order kinetic decay. A slower than
%%% exponential thermal detrapping is explained by electron retrapping or
%%% distant-dependent probabilities of detrapping.

% Georgina King [georgina.king@unil.ch], August 2015


function [nNf] = trapping_GOK_GOK(time,temp,kparams)

%%% Extract parameters
ddot = kparams.natDdot(1);
D0 = kparams.D0(1);
s = 10.^kparams.s10(1);
E = kparams.Et(1);
b = kparams.GOK_b(1);
a = kparams.GOK_a(1);


%%% Define constants
kb = 8.617343e-5;                                                           % Boltzmann constant [eV/K]
Hs = 3e15;                                                                  % Hs=s_ath value after Huntley (2006) J. Phys. D.
magic_ratio = ddot/D0;
Ma = 1e6*365.25*24*3600;                                                    % 1 [Myr] in [s]
nstep = length(time);
dt = ((max(time)-min(time))/(nstep-1))*Ma;                                  % time step in [s]

%%% Precompute nN for the random Tt path (Euler integration method)
nN=zeros(nstep);
T = temp+273.15;

%%% Loop through the time steps
for i=2:nstep

    %%% Compute inv_tauth for the current time step (1 x 1)
    inv_tauth = s*exp(-E/(kb*T(i-1))); %convert to Ma

    %%% Compute xkd and xk for all combinations of rprime for the current time step (1 x nrp)
    xkd = -a*magic_ratio*(1-nN(i-1,:)).^(a-1)-b*inv_tauth;
    xk = magic_ratio*(1-nN(i-1,:)).^a-inv_tauth.*(nN(i-1,:)).^b;
    nN(i,:) = nN(i-1,:)+dt*xk./(1-dt*xkd);
end

nNf = nN(:,1);
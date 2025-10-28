function [nNf] = trapping_LIN_GOK(time,temp,kparams)


%%% Extract parameters
ddot = kparams.natDdot(1);
s = 10.^kparams.s10(1);
E = kparams.Et(1);
b = kparams.GOK_b(1);

%%% Define constants
Ma = 3600*24*365.25*1e6;
kb = 8.617343e-5;
magic_ratio = ddot;
nstep = length(time);
dt = ((max(time)-min(time))/(nstep-1))*Ma;

%%% Computes Dose (called nN) for the random Tt path (Euler integration method)
nN=zeros(1,nstep);

for j=2:nstep
    inv_tauth=s*exp(-E/(kb*(temp(j-1)+273.15)));
    xkd=-b*inv_tauth.*(nN(j-1)).^(b-1);
    xk=magic_ratio-inv_tauth.*(nN(j-1)).^b;
    nN(j) = nN(j-1)+dt*xk./(1-dt*xkd);
end

nNf = nN;
function [nNf] = trapping_SSE_GOK(time,temp,kparams)


%%% Extract parameters
ddot = kparams.natDdot(1);
D0 = kparams.D0(1);
s = 10.^kparams.s10(1);
b = kparams.GOK_b(1);
E = kparams.Et(1);
a = 1;

%%% Define constants
Ma = 3600*24*365.25*1e6;
kb = 8.617343e-5; Hs = 3e15;                                                % s value after Huntley (2006) J. Phys. D.
magic_ratio = ddot/D0;
nstep = length(time);
dt = ((max(time)-min(time))/(nstep-1))*Ma;
T = temp+273.15;

%%% Computes nN for the random Tt path (Euler integration method)
nN=zeros(nstep);
for j=2:nstep
    inv_tauth=s*exp(-E/(kb*(T(j-1))));
    xkd=-a*magic_ratio*(1-nN(j-1,:)).^(a-1)-b*inv_tauth;
    xk=magic_ratio*(1-nN(j-1,:)).^a-inv_tauth.*nN(j-1,:);
    nN(j,:) = nN(j-1,:)+dt*xk./(1-dt*xkd);
end

nNf = nN(:,1);
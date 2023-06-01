 function [Dose] = trapping_LIN_GOK_ESR(time,temp,kparams);

% Extract parameters
ddot = kparams.natDdot(1); 
s = 10.^kparams.s10_GOK(1);
E = kparams.Et_GOK(1); 
b = kparams.GOKorder(1);

% Define constants
Ma = 3600*24*365*1e6; %s
kb = 8.617343e-5; %eV
magic_ratio = ddot; %Gy s-1/Ma
nstep = length(time);
dt = ((max(time)-min(time))/(nstep-1))*Ma;
 
% computes Dose (called nN) for the random Tt path (Euler integration method)
nN=zeros(1,nstep);
    for j=2:nstep
            inv_tauth=s*exp(-E/(kb*(temp(j-1)+273.15))); 
            xkd=-b*inv_tauth.*(nN(j-1)).^(b-1);          
            xk=magic_ratio-inv_tauth.*(nN(j-1)).^b;
            nN(j) = nN(j-1)+dt*xk./(1-dt*xkd);
    end 
Dose=nN;






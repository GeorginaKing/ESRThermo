function [nNf] = trapping_SSE_GOK(time,temp,kparams);

%%requires testing and confirmation
 
 
Ma = 3600*24*365*1e6;
% Extract parameters
ddot = kparams.natDdot(1); %convert from Gy/s to Gy/Ma
D0 = kparams.D0(1);
s = 10.^kparams.s10_GOK(1);
b = kparams.GOKorder(1);
E = kparams.Et_GOK(1); 
a = 1;
% ddot = 1.4117e-10; %convert from Gy/s to Gy/Ma
% D0 = 2.2564e+03;
% s = 10.^13.2834;
% b = 4.8453;
% E = 1.4737; 
% a = 1;
% temp = linspace(100, 0, 1000);
% time = linspace(1E6, 0, 1000);
 


% Define constants
kb = 8.617343e-5; Hs = 3e15; %s value after Huntley (2006) J. Phys. D.
magic_ratio = ddot/D0;
nstep = length(time);
dt = ((max(time)-min(time))/(nstep-1))*Ma;
T = temp+273.15;

 
 % computes nN for the random Tt path (Euler integration method)
nN=zeros(nstep);
    for j=2:nstep
            inv_tauth=s*exp(-E/(kb*(temp(j-1)+273.15))); %convert to Ma
            xkd=-a*magic_ratio*(1-nN(j-1,:)).^(a-1)-b*inv_tauth;     
            xk=magic_ratio*(1-nN(j-1,:)).^a-inv_tauth.*nN(j-1,:);
            nN(j,:) = nN(j-1,:)+dt*xk./(1-dt*xkd); 
    end

 
   nNf = nN(:,1);




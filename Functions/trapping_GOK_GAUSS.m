 function [nNf] = trapping_GOK_GAUSS(time,temp,kparams);

% Does not function

 
Ma = 3600*24*365*1e6;
% Extract parameters
ddot = kparams.natDdot(1); %convert from Gy/s to Gy/Ma
D0 = kparams.D0(1);
s = 10.^kparams.s10_GAUSS(1);
Et = kparams.Et_GAUSS(1);
sigmaEt = kparams.sigmaEt_GAUSS(1);
a = kparams.GOKgrowthorder(1);

% % Et=1.6;
% % s=10.^13.96;
% % sigmaEt = 0.0944;
% % a = 1.0134;
% % D0 = 2407;
% % ddot = 1.4117e-10;
% %  
% % temp = linspace(100, 0, 1000);
% % time = linspace(1E6, 0, 1000);

% Define constants
kb = 8.617343e-5; 
magic_ratio = ddot/D0;
nstep = length(time);
dt = ((max(time)-min(time))/(nstep-1))*Ma;

% Create variables for GAUSS model Lambert et al (Submitted)
Ea=[(5/nstep):(5/nstep):5]'; nEa=length(Ea);
pEa = exp(-0.5*((Ea-Et)./sigmaEt).^2)/(sigmaEt*sqrt(2*pi));
npEa = sum(pEa);
T = (temp+273.15)';


 % computes nN for the random Tt path (Euler integration method)
nN = zeros(nEa,nstep);
nNf = zeros(1,nstep);
    for j=2:nstep
            inv_tauth = s*exp(-Ea./(kb.*T(j-1)))*ones(1,nEa); 
            xkd=-a*magic_ratio*(1-nN(:,j-1)).^(a-1)-inv_tauth;          
            xk=magic_ratio*(1-nN(:,j-1)).^a-inv_tauth.*(nN(:,j-1));
            nN(:,j) = nN(:,j-1)+dt*xk(:,j-1)./(1-dt*xkd(:,j-1));
            nNf(j) = pEa'*nN(:,j);
            
    end

 
   nNf = nNf./npEa;

   plot(nNf,temp);



%%%% STAGE 3a, Inversion %%%%
%%%% Inverts data to determine cooling histories
%%%% Requires output of Stage 2a

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arnaud Duverger, 2015 arnaud.duverger@ens.fr %
% Frédéric Herman, 2015, revised 2018, frederic.herman@unil.ch %
% Georgina King, 2015, revised 2018, georgina.king@unil.ch %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
clearvars -except filename filenamevec NITL nSARA SARA_fittype ITL_fittype SARA_MODEL ITL_MODEL; close all; clc;

load(['./ComputeData/' filename '_' SARA_MODEL '_' ITL_MODEL '_ESRfitpar.mat']);

nt = length(records); %number of signals
niter = 2000; %number of random realisations

%% Import Kinetic data
kp = [records.params];
if SARA_fittype==3; nNnat = kp.SARADeGy(1); snNnat = kp.SARADeGy(2); snNnat(snNnat<(0.1*nNnat))=0.1*nNnat;  
elseif SARA_fittype==1 || SARA_fittype == 2; nNnat = kp.SARAnNnat(1); snNnat = kp.SARAnNnat(2); snNnat(snNnat<(0.1*nNnat))=0.1*nNnat; end

%%Model parameters (TO CHANGE)
Tmax = 200;                 % Maximum temp in C
Tsurf = 25;  TsurfErr = 10; % Average modern surface or sample temp in C, and error
tmax = 5;                   % Time in Ma
tmin = 0;                   % Time in Ma
nstep = 500;                % discretisation of the model in time
time = linspace(tmin,tmax,nstep);

%Define uncertainties (OPTIONAL)
% pererr = [snNnat-nNnat<=0.05];
% snNnat(pererr) = 0.05;
snNat = 0.5*nNnat;

misOUT = zeros(niter,1);
Time = zeros(niter,nstep);
Temp = zeros(niter,nstep);
nNmod = zeros(niter,nstep,nt);

for i = 1:niter 
    Tmin = Tsurf-(TsurfErr*rand);
	% create a random path
	tT = randpathAD([tmin Tmax],[tmax Tmin]); 
	    
    % patch to keep model strictly monotonic, added GK 27.07.2016
    for ii=1:length(tT)-1; if (tT(ii+1)==tT(ii) || tT(ii+1)<tT(ii));  tT(ii+1)=tT(ii)+1e-7; end;end
    
    % interpolates the random path on a regular grid
	temp = interp1(tT(:,1),tT(:,2),time,'linear');  
    
	nNf = zeros(nt,1); v = zeros(nstep,nt); 
    
   
    if SARA_fittype==1              %SSE
       if ITL_fittype==1
           v = trapping_SSE_GOK(time,temp,kp);
       elseif ITL_fittype==2
           v = trapping_SSE_GAUSS(time,temp,kp);
        end
   
    
    elseif SARA_fittype==2          %GOK
       if ITL_fittype==1
           v = trapping_GOK_GOK(time,temp,kp);
       elseif ITL_fittype==2                              
           v = trapping_GOK_GAUSS(time,temp,kp);
       end
    
    elseif SARA_fittype==3          %LIN
       if ITL_fittype==1
           v = trapping_LIN_GOK_ESR(time,temp,kp);
       elseif ITL_fittype==2
           v = trapping_LIN_GAUSS_ESR(time,temp,kp);
       end
    end      
    
    nNf = v(end);
    nNmod(i,:)=v;
    
    misOUT(i)=(1/2*(nNnat/snNnat)*(log(nNnat./nNf))).^2;
        
	Time(i,:) = time; % saves time 
	Temp(i,:) = temp; % save temp
end

% Sort data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[sortedmisOUT,IX] = sort(misOUT);
sortedTime = Time(IX,:);
sortedTemp = Temp(IX,:);
sortednNmod = nNmod(IX,:,:);

% Save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tt.misfit = sortedmisOUT;
Tt.time = sortedTime;
Tt.temp = sortedTemp;
Tt.nNmod = sortednNmod;
Tt.nNnat = nNnat;
Tt.snNnat = snNnat;

save(['./ComputeData/' filename '_' SARA_fittype '_' ITL_fittype '_Tt.mat'], 'Tt'); 
toc

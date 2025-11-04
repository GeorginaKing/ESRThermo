%%%% STAGE 2a, FitParameters %%%%
%%%% Fits data with selected models, requires input from Stage 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Georgina King, 2022       georgina.king@unil.ch %
% Melanie Bartz, 2022       melanie.bartz@unil.ch %
% Chloé Bouscary, 2025  	chloebouscary@gmail.com %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clearvars -except filename filenamevec nSAR nSARvec nITH nITHvec ...
    SAR_fittype SAR_MODEL ITH_fittype ITH_MODEL tStartall;
close all;

tStart=tic;                                                                 % associated with the ‘toc’ at the end, records time of execution in [min] and [s]


%%% Load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['./ComputeData/' filename '.mat']);                                   % loads data file from Stage 1
nt = length(records);                                                       % number of signals (e.g. Al-center, Ti-center etc.)

%%% Define global parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global isoT measESR                                                         % makes parameters available to the functions used when fitting the data

%%% Compute vectors for fitting the data with the different models
tvec = logspace(0,6,100);                                                   % creates a vector with 100 columns, generates 100 points between 10^0 and 10^6
tmat = tvec'*ones(1,nITH);                                                  % creates a matrix with nITH columns (nb of isothermal holding temperatures) and 100 lines, with the values of tvec

%%% Define constants
ka = 1e3*365.25*24*3600;                                                    % 1 [kyr] in [s]

%%% Fit the data with the different models %%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nt                                                                  % loops through the number of traps

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Dose Response Curve, regenerative dose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Extract data from .mat file
    labDdot = records(i).rawdata(1).labDdot;                                % laboratory dose rate [Gy/s]
    t = [records(i).rawdata(1:nSAR).t];                                     % doses = irradiation time [s]
    ESR = [records(i).rawdata(1:nSAR).ESR];                            	    % ESR signal intensity (a.u.), non normalised
    ok = isfinite(t); x = t(ok)*labDdot; y = ESR(ok);                       % removes NaN values, and converts measurement times from [s] to [Gy]
    NAT = ESR(~ok); rmNA = isfinite(NAT); NAT = NAT(rmNA);                  % extracts natural signal
    maxR = max(t)*labDdot;                                                  % maximum regenerative dose in [Gy]
    modvec = linspace(0,maxR+1000,1e3);                                     % creates dose vector for graphs

    x = reshape(x,length(x)/nSAR,nSAR);
    y = reshape(y,length(y)/nSAR,nSAR);

    ESR_out = NaN(size(y)); x_out = NaN(size(y));

    for k=1:nSAR
        ESR_y = y(:,k); ESR_x = x(:,k);
        [my,yIX] = (max(ESR_y)); ESR_y = ESR_y(1:yIX); ESR_x = x(1:yIX);
        ESR_out(1:length(ESR_y),k) = (ESR_y);
        x_out(1:length(ESR_x),k) = (ESR_x);
    end

    NAT = mean(NAT./max(ESR_out)); 
    if length(NAT)>1; sNAT = std(NAT./max(ESR_out)); else; sNAT = NAT.*0.05; end

    ok = isfinite(x_out);
    x = x_out(ok); y = ESR_out(ok);

    %%% Single saturating exponential (SSE) fit %%%
    if SAR_fittype==1 % SSE fit

        %%% nlinfit solution
        beta0 = [1 1000]';                                                  % initial guesses on parameters: scaling parameter, D0 estimation in [s]
        [beta,R,J,Cov,~] = nlinfit(x,y,@SAR_SSE,beta0);

        %%% Confidence intervals
        [Ypred,d] = nlpredci(@SAR_SSE,modvec,beta,R,'covar',Cov,'alpha',.05,'predopt','curve');
        sbeta = nlparci(beta,R,'covar',Cov,'alpha',1-.68);
        sigma = abs(beta-sbeta(:,2));                                       % calculates parameter uncertainties

        delta = (Ypred-d); dl=fliplr(Ypred+d);
        xx = [modvec fliplr(modvec)]'; delta = horzcat(delta, dl)';

        %%% Extract D0 [s]
        D0 = beta(2); sD0 = mean(sigma(2,:));

        %%% Extract De [s] from dose response curve
        if NAT <= max(Ypred)
            De_interp = (interp1(Ypred,modvec,[NAT NAT+sNAT NAT-sNAT], 'linear'));
            sDe_interp = ((De_interp(2)-De_interp(1))+(De_interp(1)-De_interp(3)))/2;
        else
            De_interp = NaN;
            sDe_interp = NaN;
        end

        %%% Calculate error on De value [Gy]
        De = abs(De_interp(1));
        sDe = abs(De_interp(1)*sqrt(((d/Ypred)^2)+((sDe_interp/De_interp(1))^2)+(((mean(labDdot)*2/100)/mean(labDdot))^2)));

        %%% Calculate n/N
        N = beta(1);
        nNnat = NAT/N;
        snNnat = nNnat*sqrt(((d/Ypred)^2)+(sNAT/NAT)^2);                    % calculates nN from SAR curve

    %%% General order kinetics (GOK) fit %%%
    elseif SAR_fittype==2 % GOK fit

        %%% nlinfit solution
        beta0 = [1 1000 2]';                                                % initial guesses on parameters: a, D0 in [s], GOK growth order
        [beta,R,J,Cov,~] = nlinfit(x,y,@SAR_GOK,beta0);

        %%% Confidence intervals
        [Ypred,d] = nlpredci(@SAR_GOK,modvec,beta,R,'covar',Cov,'alpha',.05,'predopt','curve');
        sbeta = nlparci(beta,R,'covar',Cov,'alpha',1-.68);
        sigma = abs(beta-sbeta(:,2));                                       % calculates parameter uncertainties

        delta = (Ypred-d); dl=fliplr(Ypred+d);
        xx = [modvec fliplr(modvec)]'; delta = horzcat(delta, dl)';

        %%% Extract D0 [s]
        D0 = beta(2); sD0 = mean(sigma(2,:));

        %%% Extract GOK growth order
        GOK_a = beta(3); sGOK_a = mean(sigma(3,:));

        %%% Extract De [s] from dose response curve
        De_interp = (interp1(Ypred,modvec,[NAT NAT+sNAT NAT-sNAT], 'linear'));
        sDe_interp = ((De_interp(2)-De_interp(1))+(De_interp(1)-De_interp(3)))/2;

        %%% Calculate error on De value
        De = abs(De_interp(1));
        sDe = abs(De_interp(1)*sqrt(((d/Ypred)^2)+((sDe_interp/De_interp(1))^2)+(((mean(labDdot)*2/100)/mean(labDdot))^2)));

        N =  beta(1); nNnat = NAT/N; snNnat = nNnat*sqrt(((d/Ypred)^2)+(sNAT/NAT)^2); % calculate nN from SAR curve

    %%% Linear (LIN) fit %%%
    elseif SAR_fittype==3

        %%% fit solution
        [beta,SARgof,SARoutput] = fit(x, y,'poly1');                        % SARgof = goodness of fit, fit data using a linear fit  

        %%% Confidence intervals
        [d,Ypred] = predint(beta,modvec,0.68);                              % calculates 1 sigma uncertainties on the fit
        sbeta = confint(beta,0.68);
        sigma = abs(coeffvalues(beta)-sbeta(1,:));                          % calculates parameter uncertainties

        xx = [modvec fliplr(modvec)]'; delta = vertcat(d(:,1), flipud(d(:,2)))';

        %%% Calculate De [s]
        De_interp = (interp1(Ypred,modvec,[NAT NAT+sNAT NAT-sNAT], 'linear')); % extract De from dose response curve
        sDe_interp = ((De_interp(2)-De_interp(1))+(De_interp(1)-De_interp(3)))/2;

        %%% Calculate error on De value
        De = abs(De_interp(1));
        sDe = abs(De_interp(1)*sqrt(((sigma(1)/beta(1))^2)+((sigma(2)/beta(2))^2)+((sDe_interp/De_interp(1))^2)+(((mean(labDdot)*2/100)/mean(labDdot))^2)));

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Isothermal Decay %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Extract data from .mat file
    isoT = [records(i).rawdata(nSAR+1:nSAR+nITH).T];                        % creates a vector with the different temperatures T [°C] for isothermal decay
    measESR = [records(i).rawdata(nSAR+1).t]; measESR = measESR(isfinite(measESR)); % removes the NaN data from measESR
    itht = [records(i).rawdata(nSAR+1:nSAR+nITH).t];                        % creates a matrix with all the time records of the isothermal decay [s]
    ithESR = [records(i).rawdata(nSAR+1:nSAR+nITH).ESR];                    % creates a matrix with all the ESR records of the isothermal decay
    ok = isfinite(itht); ithx = itht(ok); ithy = ithESR(ok);                % removes the NaN data from the data
    ithsort = NaN(length(measESR),length(isoT));

    %%% General order kinetics, GOK %%%
    if ITH_fittype==2 % GOK model

        %%% nlinfit solution
        beta03 = [15 1.6 4]';                                               % initial estimates of s [s-1], Et [eV], GOK order b
        [beta3,R,J,Cov,~] = nlinfit(ithx,ithy,@ITH_GOK,beta03);

        %%% Confidence intervals
        measESR = tvec;                                                     % redefine time dimension for model
        [IsoPred,delta2] = nlpredci(@ITH_GOK,tmat(:),beta3,R,'covar',Cov,'alpha',.05,'predopt','curve');
        sbeta3 = nlparci(beta3,R,'covar',Cov,'alpha',1-.68);
        sigma3 = abs(beta3-sbeta3(:,2));                                    % calculate error

        %%% Extract parameters
        s10 = beta3(1); ss10 = mean(sigma3(1,:));
        Et = beta3(2); sEt = mean(sigma3(2,:));
        GOK_b = beta3(3); sGOK_b = mean(sigma3(3,:));

    %%% Gauss model, GAUSS %%%
    elseif ITH_fittype==3 % GAUSS model

        %%% nlinfit solution
        beta03 = [10 1.4 0.5]';                                             % initial estimates of s [s-1], mu(Et) [eV], sigma(Et) [eV]
        [beta3,R,J,Cov,~] = nlinfit(ithx,ithy,@ITH_GAUSS,beta03);

        %%% Confidence intervals
        measESR = tvec;                                                     % redefine time dimension for model
        [IsoPred,delta2] = nlpredci(@ITH_GAUSS,tmat(:),beta3,R,'covar',Cov,'alpha',.05,'predopt','curve');
        sbeta3 = nlparci(beta3,R,'covar',Cov,'alpha',1-.68);
        sigma3 = [abs(beta3-sbeta3(:,1)),abs(beta3-sbeta3(:,2))];           % calculates parameter lower and upper uncertainties

        %%% Extract parameters
        s10 = beta3(1); ss10 = mean(sigma3(1,:));
        Et = beta3(2); sEt = mean(sigma3(2,:));
        sigmaEt = beta3(3); ssigmaEt = mean(sigma3(3,:));

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Save parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    records(i).params.De = [De,sDe];                                        % in [Gy]
    if SAR_fittype==3; records(i).params.D0 = [NaN,NaN];                    % for LIN model
    else records(i).params.D0 = [D0,sD0]; end                               % for GOK and GAUSS models, in [Gy]
    if SAR_fittype==2; records(i).params.GOK_a = [GOK_a,sGOK_a];            % for GOK model
    else records(i).params.GOK_a = [NaN,NaN]; end                           % for LIN and GAUSS models
    
    records(i).params.s10 = [s10,ss10];                                     % in [s-1]
    records(i).params.Et = [Et,sEt];                                        % in [eV] 
    if ITH_fittype==2; records(i).params.GOK_b = [GOK_b,sGOK_b];            % for GOK model
    else records(i).params.GOK_b = [NaN,NaN]; end                           % for GAUSS model
    if ITH_fittype==3; records(i).params.sigmaEt = [sigmaEt,ssigmaEt];      % for GAUSS model, in [eV]
    else records(i).params.sigmaEt = [NaN,NaN]; end                         % for GOK model

    AgeESR = records(i).params.De(1)./(records(i).params.natDdot(1).*ka);
    sAgeESR = AgeESR.*sqrt((records(i).params.De(2)/records(i).params.De(1))^2+((records(i).params.natDdot(2).*ka)/(records(i).params.natDdot(1).*ka))^2);
    records(i).params.AgeESR = [AgeESR,sAgeESR];

    if SAR_fittype==3; records(i).params.nNnat = [NaN,NaN];                 % for LIN model
    else records(i).params.nNnat = [nNnat,snNnat]; end                      % for GOK and GAUSS models

    records(i).SAR_model = SAR_MODEL;
    records(i).ITH_model = ITH_MODEL;

    %%% Save data for plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Dose response curve (SAR)
    records(i).plot.x       =   modvec;
    records(i).plot.y       =   Ypred;
    records(i).plot.sy      =   delta;
    records(i).plot.ESR     =   ESR;
    records(i).plot.t       =   t;

    %%% Isothermal decay (ITH)
    records(i).plot.xith    =   tmat;
    records(i).plot.yith    =   reshape(IsoPred,size(tmat));
    records(i).plot.syith   =   reshape(delta2,size(tmat));
    ithsort(1:length(ithx)) =   ithx;
    records(i).plot.tith    =   ithsort;
    ithsort(1:length(ithy)) =   ithy;
    records(i).plot.ESRith  =   ithsort;

end

save(['./ComputeData/' filename '_' SAR_MODEL '_' ITH_MODEL '_ESRfitpar.mat'],'records')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Excel file with the parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

records2 = records;
samples = strings(nt,1);
nNnat_record = zeros(nt,2);
headers = {'Sample Temp (°C)','Unc. Temp','Dr (Gy/kyr)','Unc. Dr',...
    'De (Gy)','Unc. De','D0 (Gy)','Unc. D0','GOK_a','Unc. a',...
    's10','Unc. s10','Et (eV)','Unc. Et','GOK_b','Unc. b','sigmaEt (eV)','Unc. sigmaEt',...
    'Age (ka)','Unc. age','n/N_nat','Unc. n/N_nat'};
for i=1:nt
    records2(i).params.natDdot = records(i).params.natDdot.*ka;             % transform from [Gy/s] to [Gy/kyr]
    row(i,:) = table2array(struct2table(records2(i).params));
    samples(i) = string(records2(i).id);
end
tab = array2table(row,'VariableNames',headers);
tab2 = array2table(samples,'VariableNames',{'Samples'});
tab4 = [tab2,tab];
writetable(tab4,['./ComputeData/' filename '_' SAR_MODEL '_' ITH_MODEL '_fitpar.xls'])


%%% Running time
tEnd = toc(tStart);
fprintf('Stage2a_Fitparameters took %d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));
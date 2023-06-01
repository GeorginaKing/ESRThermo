%%%% STAGE 2a, FitParameters %%%%
%%%% Fits data with selected models, requires input from Stage 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Georgina King, 2022 georgina.king@unil.ch %
% Melanie Bartz, 2022 melanie.bartz@unil.ch %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic %records time of execution
clearvars -except filename filenamevec NITL nSARA nSARAvec SARA_fittype ITL_fittype SARA_MODEL ITL_MODEL; close all; clc;
%removes variables except those needed in ESRThermo to loop through multiple files i.e. NITL etc.

load(['./ComputeData/' filename '.mat']); %loads data file from Stage1
nt = length(records);                     %number of signals
global isoT measESR                       %makes parameters available to the functions used when fitting the data

%Compute vectors for fitting the data with the different models
tvec = logspace(0,6,100);                 %Creates a matrix with 100 columns, generates 100 points between 10^0 and 10^6
tmat = tvec'*ones(1,NITL);                %NITL = Number of isothermal holding temperatures %Creates a matrix with n columns and 100 lines

kb = 8.617343*10^(-5);                    %Boltzmann constant eV/K

for i=1:nt %loop through the number of traps
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%% Dose Response Curve, regenerative dose %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    labDdot = records(i).rawdata(1).Ddot;                                           %laboratory dose rate                                      
	t = [records(i).rawdata(1:nSARA).t];                                            %Irradiation time (s)
	ESR = [records(i).rawdata(1:nSARA).ESR];                                        %ESR signal intensity (a.u.)
    ok = isfinite(t); x = t(ok)*labDdot; y = ESR(ok);                               %remove NaN values
    NAT = ESR(~ok); rmNA = isfinite(NAT); NAT=NAT(rmNA);
    maxR = max(t)*labDdot;
    modvec = linspace(0,maxR,1E3);                                                  %Create vector for graphs
                                                  
    y=reshape(y,length(y)/nSARA,nSARA); 
    x = reshape(x,length(x)/nSARA,nSARA); 
        
    ESR_out=NaN(size(y)); x_out=NaN(size(y));
    for j=1:nSARA
       ESR_y=y(:,j); ESR_x=x(:,j);
       [my yIX] = (max(ESR_y)); ESR_y=ESR_y(1:yIX); ESR_x=x(1:yIX);
       ESR_out(1:length(ESR_y),j)=(ESR_y);
       x_out(1:length(ESR_x),j)=(ESR_x);
    end
    
    sNAT = std(NAT./max(ESR_out)); NAT=mean(NAT./max(ESR_out)); 
    ESR_out = ESR_out./max(ESR_out); %normalise dose response data to max dose
    
    ok = isfinite(x_out);
    x = x_out(ok); y = ESR_out(ok);
        
   if SARA_fittype==1
        
       %%% Single saturating exponential (SSE) fit %%%
       
        %%nlinfit solution%
        beta0 = [1 1000]';  %%initial parameters, scaling parameter, D0 estimation  
        [beta,R,J,Cov,MSE] = nlinfit(x,y,@SSE,beta0);
    
        %Confidence intervals%
        [Ypred,d] = nlpredci(@SSE,modvec,beta,R,'covar',Cov,'alpha',.05,'predopt','curve');
        sbeta = nlparci(beta,R,'covar',Cov,'alpha',1-.68);
        sigma = abs(beta-sbeta(:,2)); %calculate parameter uncertainties                                            
        
        delta=(Ypred-d); dl=fliplr(Ypred+d); 
        xx = [modvec fliplr(modvec)]'; delta = horzcat(delta, dl)';
	
        %%Save params & calculate De%%
        D0 = [beta(2) sigma(2)];

        %Extract De from dose response curve 
        De_interp=(interp1(Ypred,modvec,[NAT NAT+sNAT NAT-sNAT], 'linear'));       
        sDe=((De_interp(2)-De_interp(1))+(De_interp(1)-De_interp(3)))/2;
        
        %%Calculate Error on De value
        SARADeGyErr=De_interp(1)*sqrt(((d/Ypred)^2)+((sDe/De_interp(1))^2)+(((labDdot/100)/labDdot)^2));  
        SARADeGy=[abs(De_interp(1)) SARADeGyErr]; 

        N=beta(1); SARAnNnat=NAT/N; SARAsnNnat=SARAnNnat*sqrt(((d/Ypred)^2)+(sNAT/NAT)^2);                  % calculate nN from SAAD curve        
        
    elseif  SARA_fittype==2    
        
        %%% General order kinetics (GOK) fit %%%

        beta0 = [1,1000, 2]'; %% initial parameters a, D0, GOK growth order
        [beta,R,J,Cov,MSE] = nlinfit(x,y,@GOK,beta0);
    
        %Confidence intervals%
        [Ypred,d] = nlpredci(@GOK,modvec,beta,R,'covar',Cov,'alpha',.05,'predopt','curve');
        sbeta = nlparci(beta,R,'covar',Cov,'alpha',1-.68);
        sigma = abs(beta-sbeta(:,2)); %calculate parameter uncertainties
        
        delta=(Ypred-d); dl=fliplr(Ypred+d); 
        xx = [modvec fliplr(modvec)]'; delta = horzcat(delta, dl)';
        
        %%Save params & calculate De%%
        D0 = [beta(2) sigma(2)];
        GOKgrowthorder = [beta(3) sigma(3)];
        
        % Extract De from dose response curve 
        De_interp=(interp1(Ypred,modvec,[NAT NAT+sNAT NAT-sNAT], 'linear'));
        sDe=((De_interp(2)-De_interp(1))+(De_interp(1)-De_interp(3)))/2;
        
                %%Calculate Error on De value
        SARADeGyErr=De_interp(1)*sqrt(((d/Ypred)^2)+((sDe/De_interp(1))^2)+((0.005/labDdot)^2));    
        SARADeGy=[abs(De_interp) abs(SARADeGyErr)]; 
        
        N=beta(1); SARAnNnat=NAT/N; SARAsnNnat=SARAnNnat*sqrt(((d/Ypred)^2)+(sNAT/NAT)^2);             % calculate nN from SARA curve         

        
 elseif  SARA_fittype==3         
        
        %%%linear (LIN) fit%%%

        [beta,SARAgof,SARAoutput] = fit(x, y,'poly1');                             % Fit data using a linear fit
	
        %Confidence intervals%
        [d,Ypred] = predint(beta,modvec,0.68);                                     % Calculate 1 sigma uncertainties on the fit
        sbeta = confint(beta,0.68); sigma = abs(coeffvalues(beta)-sbeta(1,:));     % calculate uncertainties
        
        xx = [modvec fliplr(modvec)]'; delta = vertcat(d(:,1), flipud(d(:,2)))';
        
        %%Calculate De%%
        De_interp=(interp1(Ypred,modvec,[NAT NAT+sNAT NAT-sNAT], 'linear'));       % Extract De from dose response curve 
        sDe=((De_interp(2)-De_interp(1))+(De_interp(1)-De_interp(3)))/2;
        
        %%Calculate Error on De value
        SARADeGyErr=De_interp(1)*sqrt(((sigma(1)/beta(1))^2)+((sigma(2)/beta(2))^2)+((sDe/De_interp(1))^2)+((0.005/labDdot)^2));
        SARADeGy=[abs(De_interp) abs(SARADeGyErr)];
        
      end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Isothermal Decay %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	isoT = [records(i).rawdata(1+nSARA:end).T]; ok=isfinite(isoT); isoT=isoT(ok);
	measESR = [records(i).rawdata(1+nSARA).t]; ok=isfinite(measESR); measESR=measESR(ok);
	itlsort = NaN(length(measESR),length(isoT));	
	itlt = [records(i).rawdata(1+nSARA:end).t];
	itlESR = [records(i).rawdata(1+nSARA:end).ESR];
	ok = isfinite(itlt); itlx = itlt(ok); itly = itlESR(ok);

    if ITL_fittype==1
          
    %%% general order kinetics, GOK %%%
    beta03 = [15 1.6 4]';                                                           %s, E, b;
	[paramsGOK,resn,jacob,Cov,MSE] = nlinfit(itlx,itly,@ITLGOK,beta03);
       
    %Confidence intervals%
	measESR = tvec;                                                                 %redefine time dim for model
	[IsoPred,delta2] = nlpredci(@ITLGOK,tmat(:),paramsGOK,resn,'covar',Cov,'alpha',.05,'predopt','curve');
	sbeta4 = nlparci(paramsGOK,resn,'covar',Cov,'alpha',1-.68);
	sigma4 = abs(paramsGOK-sbeta4(:,2));                                            %calculate error
    
    %Calculate parameters%%
 	s10_GOK = [paramsGOK(1) sigma4(1)]; muEt_GOK = [paramsGOK(2) sigma4(2)]; GOKorder = [paramsGOK(3)];  
    Tau_GOK = ((1/(10^s10_GOK(1)))*exp(paramsGOK(2)/(kb*(15+273.15)))/(3600*24*365));
    
    %save params%
    records(i).params.Et_GOK = [muEt_GOK];
 	records(i).params.s10_GOK = [s10_GOK];
    records(i).params.GOKorder = [GOKorder]; 
     
    elseif ITL_fittype==2 
    
    %%% Gauss model, GAUSS %%%
         
  	beta03 = [10 1.4 0.5]'; %s, muEt, sigmaEt; A
	[paramsGAUSS,resn,jacob,Cov,MSE] = nlinfit(itlx,itly,@ITLGAUSS,beta03);        %can manipulate matrix size in @ITLBTS_raw to allow for irregular data 
        
    %Confidence intervals%
	measESR = tvec;                                                                %redefine time dimension for model
	[IsoPred,delta2] = nlpredci(@ITLGAUSS,tmat(:),paramsGAUSS,resn,'covar',Cov,'alpha',.05,'predopt','curve');
	sbeta4 = nlparci(paramsGAUSS,resn,'covar',Cov,'alpha',1-.68);
	sigma4 = abs(paramsGAUSS-sbeta4(:,2));                                         %calculate error
    
    %Calculate parameters%%
 	s10_GAUSS = [paramsGAUSS(1) sigma4(1)]; muEt_GAUSS = [paramsGAUSS(2) sigma4(2)]; sigmaEt_GAUSS = [paramsGAUSS(3) sigma4(3)];   
    
    %save params%
    records(i).params.Et_GAUSS = [muEt_GAUSS];
 	records(i).params.s10_GAUSS = [s10_GAUSS];
    records(i).params.sigmaEt_GAUSS = [sigmaEt_GAUSS]; 
    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %%%Save parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if SARA_fittype==1 || SARA_fittype==2; 
        records(i).params.SARADeGy = [SARADeGy];
        records(i).params.SARAnNnat = [SARAnNnat,SARAsnNnat]; 
        records(i).params.D0 = [D0];
        records(i).params.cov = [Cov];
    elseif SARA_fittype==3; 
        records(i).params.SARADeGy = [SARADeGy];
        records(i).params.cov = [Cov];
    end
    if SARA_fittype==2; records(i).params.GOKgrowthorder = [GOKgrowthorder]; end
    records(i).SARA_fittype=SARA_fittype;
    records(i).ITL_fittype=ITL_fittype;
    
	%%%Save data for plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Dose response curve (SARA)
    records(i).plot.x = modvec;
    records(i).plot.y = Ypred;
    records(i).plot.sy = delta;
    records(i).plot.t = t;
	records(i).plot.ESR = ESR;
    
    %Isothermal decay 
	records(i).plot.xitl = tmat;
	records(i).plot.yitl = reshape(IsoPred,size(tmat));
	records(i).plot.syitl = reshape(delta2,size(tmat));
 	itlsort(1:length(itlx)) = itlx; 
 	records(i).plot.titl = itlsort;
 	itlsort(1:length(itly)) = itly; 
 	records(i).plot.ESRitl = itlsort;

end

save(['./ComputeData/' filename '_' SARA_MODEL '_' ITL_MODEL '_ESRfitpar.mat'],'records')
toc

    
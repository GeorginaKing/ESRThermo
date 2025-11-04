%%%% STAGE 2b, PlotLxTx %%%%
%%%% Plots data fits with selected models, requires input from Stage 2a

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Georgina King, 2022       georgina.king@unil.ch %
% Melanie Bartz, 2022       melanie.bartz@unil.ch %
% Chloé Bouscary, 2025  	chloebouscary@gmail.com %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clearvars -except filename filenamevec nSAR nSARvec nITH nITHvec ...
    SAR_fittype SAR_MODEL ITH_fittype ITH_MODEL tStartall;
close all;

tStart=tic;                                                                 % associated with the ‘toc’ at the end, records time of execution in [min] and [s]


%%% Load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['./ComputeData/' filename '_' SAR_MODEL '_' ITH_MODEL '_ESRfitpar.mat']);
nt = length(records);                                                       % number of signals


%%% Plot the figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nt
    figure(i);                                                              % one figure per signal
    set(gcf,'units','normalized','outerposition',[0.2 0.2 0.8 0.8]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Dose Response Curve (SAR) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Read in vectors and data for plotting
    labDdot = records(i).rawdata.labDdot;
    x = records(i).plot.x;
    y = records(i).plot.y;
    delta = records(i).plot.sy;
    t = records(i).plot.t*labDdot; ok = isfinite(t); t = t(ok);
    ESR = records(i).plot.ESR; ESR = ESR(ok);
    De = records(i).params.De;
    if isnan(De(1)); De(1)=0; end
    if De(1)<0; De(1)=1; end

    ESR=reshape(ESR,length(ESR)/nSAR,nSAR);
    t = reshape(t,length(t)/nSAR,nSAR);

    ESR_out=NaN(size(ESR)); x_out=ESR_out; ESR_out2=ESR_out; x_out2=ESR_out;

    for j=1:nSAR
        ESR_y=ESR(:,j); ESR_x=t(:,j);
        normFac=max(ESR_y); ESR_y=ESR_y./normFac;
        [my,yIX] = (max(ESR_y)); ESR_yout=ESR_y(1:yIX); ESR_xout=t(1:yIX);
        ESR_y2 = ESR_y(yIX+1:end); ESR_x2 = ESR_x(yIX+1:end);
        ESR_out(1:length(ESR_yout),j)=(ESR_yout);
        x_out(1:length(ESR_xout),j)=(ESR_xout);
        ESR_out2(length(ESR_yout)+1:end,j)=(ESR_y2);
        x_out2(length(ESR_xout)+1:end,j)=(ESR_x2);
    end

    ok = isfinite(x_out);
    t = x_out(ok); ESR = ESR_out(ok);
    ok = isfinite(x_out2);
    x2 = x_out2(ok); y2 = ESR_out2(ok);

    xx = [x, fliplr(x)];

    %%% Plot dose response data
    subplot(1,2,1);
    xlabel('Dose (Gy)'); ylabel('ESR intensity (a.u.)');
    axis square; box on; hold on;
    xlim([0 round(max(t)+1000,-3)]);
    if SAR_fittype==1 || SAR_fittype==2
        ylim([0 1.1]);
    end

    fill(xx,delta,[0.5 0.5 1], 'EdgeColor','blue');
    plot(x,y,'-r');

    scatter(t,ESR,[],'MarkerFaceColor','red','MarkerEdgeColor','black');
    scatter(x2,y2,[],'MarkerFaceColor','white','MarkerEdgeColor','black');

    if SAR_fittype==1 || SAR_fittype==2
        nN=records(i).params.nNnat;
        n=scatter(0,nN(1),'o','MarkerFaceColor', 'y','MarkerEdgeColor','black');
        legend(n,{'(n/N)_n_a_t'},'Location','SouthWest');
    end


    %%% Add text to figure
    text(round(max(t)+1000,-3)-1000,1.05,'A','fontweight','bold');
    text((0+500),1.05,SAR_MODEL,'fontweight','bold');                       % dose response curve model


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Isothermal Decay %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Define colors for plot
    colmat = colormap(parula);
    index = round(linspace(1,length(colmat),nITH));

    %%% Read in vectors and data for plotting
    ithx = records(i).plot.xith;
    ithy = records(i).plot.yith;
    ithdelta = records(i).plot.syith;
    itht = records(i).plot.tith;
    ithL = records(i).plot.ESRith;

    isoT = [records(i).rawdata(nSAR+1:end).T]; ok=isfinite(isoT); isoT=isoT(ok); % isothermal temperatures [°C]

    d1 = (ithy-ithdelta); d1(d1<0)=0.0001;
    d2 = (ithy+ithdelta); d2(d2<0)=0.0001;

    ithxx = [ithx; flipud(ithx)];
    ithd = [d1; flipud(d2)];

    %%% Plot isothermal decay data
    subplot(1,2,2,'xScale','log');
    xlabel('Time (s)'); ylabel('ESR intensity (a.u.)');
    axis([1 1e6 0 1.1],'square'); box on; hold on;
    set(gca,'XTick',logspace(0,6,7),'XScale','log');

    fill(ithxx,ithd,[0.5 0.5 1],'EdgeColor','none');
    P = plot(ithx,ithy);

    %%% Create the legend
    lgd=cell(nITH,1);
    for k = 1:length(P)
        l(k)=scatter(itht(:,k),ithL(:,k),[],colmat(index(k),:),'filled','MarkerEdgeColor','k');
        set(P(k),'Color',colmat(index(k),:));
        lgd{k}=strcat(num2str(isoT(k)),' °C');
    end
    legend(l,lgd,'Location','SouthWest');

    %%% Add text to figure
    text(3e5,1.05,'B','fontweight','bold');
    text(2,1.05,ITH_MODEL,'fontweight','bold');                             % isothermal decay model

    %%% Save plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(1,2,1);
    text((round(max(t)+1000)-4000),0.10,filename,'fontweight','bold');      % sample name
    TypeMeasurement(i) = records(i).typeMeasurement;                        % type of measured signal (ESR, OSL, etc.)
    TypeSignal(i) = records(i).typeSignal;                                  % type of signal measured (Al-center, Ti-center, etc.)
    text((round(max(t)+1000)-4000),0.05,[cell2mat(TypeMeasurement(i)) ' ' cell2mat(TypeSignal(i))],'fontweight','bold'); % signal measured

    print('-dpng',['./Figures/' filename '_' char(TypeMeasurement(i)) '-' char(TypeSignal(i)) '_' SAR_MODEL '_' ITH_MODEL '.png']);
    % print('-dsvg',['./Figures/' filename '_' char(TypeMeasurement(i)) '-' char(TypeSignal(i)) '_' SAR_MODEL '_' ITH_MODEL '.svg']);
    % print('-depsc',['./Figures/' filename '_' char(TypeMeasurement(i)) '-' char(TypeSignal(i)) '_' SAR_MODEL '_' ITH_MODEL '.eps']);

end


%%% Running time
tEnd = toc(tStart);
fprintf('Stage2b_PlotFit took %d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));
%%%% STAGE 2b, PlotLxTx %%%%
%%%% Plots data fits with selected models, requires input from Stage 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Georgina King, 2022 georgina.king@unil.ch %
% Melanie Bartz, 2022 melanie.bartz@unil.ch %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except filename filenamevec NITL nSARA nSARAvec SARA_fittype ITL_fittype SARA_MODEL ITL_MODEL; close all; clc;

load(['./ComputeData/' filename '_' SARA_MODEL '_' ITL_MODEL '_ESRfitpar.mat']); 
nt = length(records); %number of signals

colmat = colormap;
index = linspace(1,61,sum([ones(size(records(1).rawdata(4:3+NITL)))]));
logtick = logspace(0,9,10);

textmod = '%s: %0.2f +/- %0.2f';

for i=1:nt
	FIG=figure(i)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% Dose Response Curve (SARA) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Read in vectors and data for plotting
    labDdot = records(i).rawdata.Ddot; 
	x = records(i).plot.x;
	y = records(i).plot.y;
	d = records(i).plot.sy; 
	t = records(i).plot.t*labDdot; ok = isfinite(t); t = t(ok);
	ESR = records(i).plot.ESR; ESR=ESR(ok);
    DeGy = records(i).params.SARADeGy;
    if isnan(DeGy(1)); DeGy(1) =0; end; if DeGy(1)<0; DeGy(1)=1; end;

    ESR=reshape(ESR,length(ESR)/nSARA,nSARA);
    t = reshape(t,length(t)/nSARA,nSARA);
        
    ESR_out=NaN(size(ESR)); x_out=ESR_out; ESR_out2=ESR_out; x_out2=ESR_out;
     
        for j=1:nSARA
            ESR_y=ESR(:,j); ESR_x=t(:,j);
            normFac=max(ESR_y); ESR_y=ESR_y./normFac;                                         
            [my yIX] = (max(ESR_y)); ESR_yout=ESR_y(1:yIX); ESR_xout=t(1:yIX);
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

    %Plot dose response data
    subplot(1,2,1); 
	
    xlabel('Dose [Gy]'); ylabel('ESR Signal Intensity [a.u.]'); axis square; box on; hold on;
    xlim([0 max(t)+3000]); 
    if SARA_fittype == 1 || SARA_fittype == 2 ylim([0 1.1]); end
    
	fill(xx,d,[0.5 0.5 1], 'EdgeColor','blue');
    plot(x,y,'-r'); 
    
    scatter(t,ESR,[],'MarkerFaceColor','red','MarkerEdgeColor','black');
    scatter(x2,y2,[],'MarkerFaceColor','white','MarkerEdgeColor','black');
 	
    if SARA_fittype==1 || SARA_fittype == 2 ; 
    nN=records(i).params.SARAnNnat; scatter(0,nN(1),'o','MarkerFaceColor', 'y','MarkerEdgeColor','black');
    end

    text(max(t+1500),1.05,'A','fontweight','bold');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% ISOTHERMAL DECAY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Read in vectors and data for plotting
	itlx = records(i).plot.xitl;
	itly = records(i).plot.yitl;
	itldelta = records(i).plot.syitl;
	itlt = records(i).plot.titl;
	itlL = records(i).plot.ESRitl;
    Index = ones(length(itlt),1)*index; 
    
	d1 = (itly-itldelta); d1(d1<0)=0.0001;
	d2 = (itly+itldelta); d2(d2<0)=0.0001;
	itly = itly;
	itlL = itlL;

	itlxx = [itlx; flipud(itlx)]; 
	itld = [d1; flipud(d2)];

    %Plot isothermal decay data
    subplot(1,2,2,'xScale','log'); 
    
	xlabel('Time [s]');
	ylabel('ESR intensity [a.u.]'); 
	axis([1 1e6 0 1.1],'square'); set(gca,'XTick',[1e0 1e1 1e2 1e3 1e4 1e5 1e6]);
	box on; hold on;

	fill(itlxx,itld,[0.5 0.5 1], 'EdgeColor','none');
	P = plot(itlx,itly);
	scatter(itlt(:),itlL(:),[], Index(:),'filled','MarkerEdgeColor', 'k');
	for k = 1:length(P)
		set(P(k), 'Color',colmat(index(k),:));
        %%Edit temperatures for isothermal temperatures used
        legend(P,'130^oC','155^oC','180^oC','205^oC','Location','SouthWest');
    end
    
    text(3e5,1.05,'B','fontweight','bold'); 
%   print('-depsc',['./Figures/' filename '_' SARA_MODEL '_' ITL_MODEL '_' records(i).id '.eps'])
    print('-dpng',['./Figures/' filename '_' SARA_MODEL '_' ITL_MODEL '_' records(i).id '.png'])
end


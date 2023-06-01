%%%% STAGE 3b, PlotTt %%%%
%%%% Plot cooling histories - requires output of Stage 3a %%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frédéric Herman, 2015 frederic.herman@unil.ch %
% modified by Arnaud Duverger, 2015 arnaud.duverger@ens.fr %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except filename filenamevec NITL nSARA SARA_fittype ITL_fittype SARA_MODEL ITL_MODEL;

set(0,'DefaultFigureColormap',feval('jet'));
set(0,'defaultfigurecolor',[1 1 1])
set(gca,'FontSize',12);
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontName','Times New Roman');
close all; clc; clf;

% load inversion output and original n/N data
load(['./ComputeData/' filename '_' SARA_fittype '_' ITL_fittype '_Tt.mat']);
misfit = Tt.misfit;
Time = Tt.time; Temp = Tt.temp;
nNs = Tt.nNmod; maxnNs = max(nNs,[],2);
nNnat = Tt.nNnat'; sigmanNnat = Tt.snNnat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nAvM = 50;                      %resolution of pdf
MTemp = [1];                    %arbitrary number to plot the signal against
[m,n,nt] = size(maxnNs);
maxnNs = reshape(maxnNs,m,nt);

% to compute the PDF for the likelhood, the Tt paths are reinterpolated onto a grid (Av_matrix)
time = Time(1,:); temp = Temp(1,:);
time_max = max(time); time_min = min(time); dt = (time_max-time_min)/(nAvM-1);
Temp_max = max(temp); 
Temp_min = 0; 
dT = (Temp_max-Temp_min)/(nAvM-1);
vec_time = time_min:dt:time_max;
vec_Temp = Temp_min:dT:Temp_max;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Data selection %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = rand(m,1); 
prob=exp(-misfit); scale = max(prob); prob = prob/scale; test = prob>R;

idefix = find(test); idefix = idefix(end:-1:1); movea = length(idefix);
max_likelihood = max(prob(idefix)); min_likelihood = min(prob(idefix));                 % colour scheme bounds
index_color = max(1,floor(63*(prob-min_likelihood)/(max_likelihood-min_likelihood))+1); % scaling the value between 1 and 64

MTEMP = ones(movea,1)*MTemp(1:nt); INDEX = index_color(idefix)*ones(1,nt);
MaxnNs = maxnNs(idefix,:);

Av_matrix=zeros(nAvM);

% add accepted model to a matrix to compute the PDF, this is based on a rejection algorithm
for k = idefix';
	vec_temp = interp1(time,Temp(k,:),vec_time,'linear');
	Tpath = (0:49)*nAvM+round((vec_temp-Temp_min)/dT)+1;
	Av_matrix(Tpath) = Av_matrix(Tpath)+1;
end
X = cumsum(Av_matrix/movea);    % for computing CIs and median

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot figures %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

map = colormap;                 % creates a table with 64 rows and three columns for RGB
col = round(100*(max(misfit(idefix)):(min(misfit)-max(misfit(idefix)))/6:min(misfit)))/100; % create colorbar for misfit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure (1) Contrasts (n/N)mod with (n/N)meas %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = figure(1);

%%% contrast measured and modelled values %%%
if SARA_fittype==1 ||  SARA_fittype==2 
    ylabel('(n_n_a_t/N)'); ylim([0 1])
elseif SARA_fittype==3 
    ylabel('De (Gy)')
end
c1 = colorbar('yticklabel',col);
set(get(c1,'title'),'string','Misfit');hold on
scatter(MTEMP(:),MaxnNs(:),20,INDEX(:),'filled');                           % predicted values
errorbar(MTemp(1:nt),nNnat(1:nt),sigmanNnat(1:nt),'ko','MarkerSize', 15);   % observed values
pos1 = get(gca,'Position'); axis square; box on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure (2) Cooling histories %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f2 = figure(2);

xlabel('Time [Ma]'); ylabel('Temperature [^oC]');
axis([time_min time_max 0 Temp_max-50],'square');
xticks = time_min:time_max/5:time_max;                                      % define the tick-marks
xticklabels = time_max:-time_max/5:time_min;                                % define the tick-labels
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);                       % plot the Xaxis with values in reverse order
c2 = colorbar; set(get(c2,'title'),'string','PDF'); hold on
[cs,hc]=contourf(vec_time,vec_Temp,Av_matrix./movea,100);                   % plot PDF
set(hc,...
    'EdgeColor','none');
contour(vec_time,vec_Temp,X,[0.05,0.95],'k','LineWidth',2);                 % 90 CI  (CI=Confidence Interval)
contour(vec_time,vec_Temp,X,[0.20,0.80],'g','LineWidth',2);                 % 60 CI
contour(vec_time,vec_Temp,X,1,'r','LineWidth',2);                           % median model
caxis([0 0.1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figure (3) Modeled dose response for individual signals %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f3 = figure(3);
set(gcf,'renderer','Painters')
for k = 1:nt
        if nt>2; 
            subplot(2,(nt+mod(nt,2))/2,k); 
        end
        
        xlabel('Time [Ma]');
        
        if SARA_fittype==1 ||  SARA_fittype==2 
            ylabel('(n_n_a_t/N)'); ylim([0 1])
        elseif SARA_fittype==3 
            ylabel('De (Gy)')
        end
        
    xticks = time_min:time_max/5:time_max;                         % define the tick-marks
    xticklabels = time_max:-time_max/5:time_min;                   % define the tick-labels
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);   
    axis square; box on; hold on;
	P = plot(time,nNs(idefix,:,k)); 
	plot(time_max,nNnat(k),'ko','MarkerFaceColor','k');
	for curve = 1:movea
		set(P(curve),'Color',map(index_color(idefix(curve)),:));
    end
end

% % print(f1,'-depsc',['./Figures/' filename '_Spag_misfit_' SARA_MODEL '_' ITL_MODEL '.eps']);
% % print(f2,'-depsc',['./Figures/' filename '_Cool_misfit_' SARA_MODEL '_' ITL_MODEL '.eps']);
% % print(f3,'-depsc',['./Figures/' filename '_nN_misfit_' SARA_MODEL '_' ITL_MODEL '.eps']);

print(f1,'-dpng',['./Figures/' filename '_Spag_misfit_' SARA_MODEL '_' ITL_MODEL '.png']);
print(f2,'-dpng',['./Figures/' filename '_Cool_misfit_' SARA_MODEL '_' ITL_MODEL '.png']);
print(f3,'-dpng',['./Figures/' filename '_nN_misfit_' SARA_MODEL '_' ITL_MODEL '.png']);
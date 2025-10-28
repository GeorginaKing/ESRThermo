%% Fits ITL data to a GOK model %%

% Georgina King [georgina.king@unil.ch], 2015


function out = ITL_GOK(beta,t)

global isoT measESR

%%% Create matrix of times, one column by isothermal temperature
mat=nan(length(measESR),length(isoT));
mat(1:length(t))=t;

%%% Extract parameters
s=10.^beta(1); 
E=beta(2);
b=beta(3);
T=isoT; 

%%% Initialise output
out=[];

%%% Isothermal decay data
for j=1:length(mat(1,:))
    ok=isfinite(mat(:,j)); time=mat(:,j); 
    
    %%% Compute the isothermal decay curves
    out = [out; ThermaldecayGOK(time(ok),T(j),E,b,s)];  
end

out = out';
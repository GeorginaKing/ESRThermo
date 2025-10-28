%% Fits ITL data to a GAUSS model %%

% Georgina King [georgina.king@unil.ch], 2015


function out = ITL_GAUSS(beta,t)

global isoT measESR

%%% Create matrix of times, one column by isothermal temperature
mat=nan(length(measESR),length(isoT));
mat(1:length(t))=t;

%%% Extract parameters
s=10.^beta(1);
muEt=beta(2);
sigmaEt=beta(3); 
T=isoT; 

%%% Initialise output
out=[];

%%% Isothermal decay data
for j=1:length(mat(1,:))
    ok=isfinite(mat(:,j)); time=mat(:,j);

    %%% Compute the isothermal decay curves
    out = [out; ThermaldecayGAUSS(time(ok),T(j),muEt,sigmaEt,s)];
end

out = out';
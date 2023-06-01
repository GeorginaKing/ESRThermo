%%%% ESRThermo 
%%%% Loader of scripts 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Arnaud Duverger, 2016, arnaud.duverger@ens.fr %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic %records time of execution
clear all; close all; clc;
addpath('./Functions/'); addpath('./Data/'); 

% define global parameters for the other scripts
global NITL nSARA SARA_fittype  ITL_fittype 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SAMPLE PARAMETERS ==> TO CHANGE %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% List of file names %%%
filenamevec={'TTY-A05-Al_2022'};

%%% Number of measurements %%%
NITL=4;                         %number of isothermal decay temperatures
nSARAvec={3};                   %number of SARA measurements

%%% Select model %%%
SARA_fittype=[2];               %1=SSE; 2=GOK; 3=LIN 
ITL_fittype=[2];                %1=GOK; 2=Gauss;

if SARA_fittype==1; SARA_MODEL='SSE'; elseif SARA_fittype==2; SARA_MODEL='GOK'; elseif SARA_fittype==3; SARA_MODEL='LIN'; end;
if ITL_fittype==1; ITL_MODEL='GOK'; elseif ITL_fittype==2; ITL_MODEL='Gauss'; end;

for j=1:size(filenamevec,1); %%% Loops through the different filenames %%
    
    nSARA=cell2mat(nSARAvec(j,:));

    filename=cell2mat(filenamevec(j,:))     
     
    %%% Convert the excel file to a .mat format %%%
%      run Stage1_ExcelToStruct2022
  
    %%% Fit the data using the selected model %%%
%       run Stage2a_Fitparameters2022
 
    %%% Plot the results of the fits %%%
%      run Stage2b_PlotFit2022 

    %%% Invert the data following King et al. (2020, Geochronology) %%
    run Stage3a_Inversion2022

    %%% Plot the results of the inversion%%
    run Stage3b_PlotTt2022

%%% Invert the data following Biswas et al. (2018, EPSL) %%
%     run Stage4a_InversionExh

    %%% Plot the results of the inversion%%
%    run Stage4b_PlotExh 
end
toc

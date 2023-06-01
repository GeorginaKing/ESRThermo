%%%% STAGE1, ExcelToStructure %%%%
%%%% Converts .xlsx files to .MAT data for input to subsequent scripts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Benny Guralnik, 2014 benny.guralnik@gmail.com, modified by Georgina King, 2018 georgina.king@unil.ch %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic %records time of execution

clearvars -except filename filenamevec NITL nSARA nSARAvec SARA_fittype ITL_fittype SARA_MODEL ITL_MODEL; close all; clc;

ka = 1e3.*365.*24.*3600;                                        % seconds in ka
[exist,sheets] = xlsfinfo([filename '.xlsx']);                  % find number of sheets in excel file

for i=1:length(sheets);                                         %loop through the different excel sheets
	[~,~,raw] = xlsread([filename '.xlsx'],i,'A1:M32');         %read data within cell range A1:M32
	raw(~cellfun(@isfloat,raw)) = {NaN};                        %create cell for the raw data to be stored within

	records(i).id = sheets{i};                                  %take record id from the sheet name
	records(i).params.natT = [raw{2,4:5}];                      %take the sample temperature as input on the .xlsx
	records(i).params.natDdot = [raw{3,4:5}]./ka;               %take the sample dose rate as input on the .xlsx
	lastrow = length(raw);                                      %calculate the limits of the raw data
	firstrow = 5;                                               %fixed start limit of the raw data

	for k=firstrow:2:lastrow;                                   %data are in pairs in the .xlsx of time (or dose) and response
		j = (k-firstrow)/2+1;
		records(i).rawdata(j).T = raw{k,2};                     %extract measurement temperatures from .xlsx (i.e. isothermal decay temp)
		records(i).rawdata(j).Ddot = raw{k+1,2};                %extract measurement dose rate
		records(i).rawdata(j).t = [raw{k,4:end}].*1e3;          %convert measurement times from ks to s
		records(i).rawdata(j).ESR = [raw{k+1,4:end}];           %extract measured values (non normalised)
    end
    
end
save(['./ComputeData/' filename '.mat'],'records');             %create data output file for reading into other stages
toc

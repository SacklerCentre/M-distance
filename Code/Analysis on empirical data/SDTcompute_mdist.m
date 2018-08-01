% This function computes the type 1 and type 2 d-prime for the first 3
% quarters and last quarter of each participant's responses.  This provides
% the data for the main analysis in the paper.

function [r] = SDTcompute_mdist

% Clear the workspace
close all;
fclose all;
clear all; clc;

% Read the raw trial data from Excel and create pdata which is a cell array
% with each participant's data in a different cell.
[allData, ~, ~] = xlsread(['RawData.xlsx'], 'TrialData');
currentLoc  = 1;
currentPnum = allData(1,1); % First participant number

index = 1;
for i = 1:size(allData,1)
    if  currentPnum ~= allData(i,1) % Reached next participant
        pdata{index} = allData(currentLoc:i-1,:);
        index = index + 1;
        currentLoc = i;
        currentPnum = allData(i,1);
    elseif i == size(allData,1) % Reached end of data
        pdata{index} = allData(currentLoc:i,:);
        index = index + 1;
        currentLoc = i;
        currentPnum = allData(i,1);
    end
end
numPs = length(pdata); % Number of participants

% Column names for pdata
pnum    = 1;
order   = 2;
class   = 3;
conf    = 4;
accur   = 5;
gramt   = 6;
attrib  = 7;
fam     = 8;
binConf = 9;

r.pnum  = zeros(numPs,1);        % participant numbers
r.t1HR  = zeros(numPs,1);        % type I HR
r.t1FA  = zeros(numPs,1);        % type I FAR
r.t2HR  = zeros(numPs,1);        % type II HR
r.t2FA  = zeros(numPs,1);        % type II FAR

r.t1HR_1st3q  = zeros(numPs,1);  % as above for first 3 quarter of responses
r.t1FA_1st3q  = zeros(numPs,1);
r.t2HR_1st3q  = zeros(numPs,1);
r.t2FA_1st3q  = zeros(numPs,1);

r.t1HR_lastq  = zeros(numPs,1);  % as above for final quarter of responses
r.t1FA_lastq  = zeros(numPs,1);
r.t2HR_lastq  = zeros(numPs,1);
r.t2FA_lastq  = zeros(numPs,1);

r.t1Dp        = zeros(numPs,1);  % type I d'
r.t2Dp        = zeros(numPs,1);  % type II d'
r.t1Dp_1st3q  = zeros(numPs,1);  % type I d' first 3 quarters
r.t2Dp_1st3q  = zeros(numPs,1);  % type II d' first 3 quarters
r.t1Dp_lastq  = zeros(numPs,1);  % type I d' last quarter
r.t2Dp_lastq  = zeros(numPs,1);  % type II d' last quarter

for i = 1:numPs
    
    disp(i)
    
    numRs = length(pdata{i}); % establishes number of responses for this participant
    end3Q = (numRs/4) * 3;    % end of the third quarter of their responses
    st4Q  = end3Q+1;          % start of the fourth quarter of their responses
    
    % values to obtain:
    % dprime
    % c
    % metad
    % confidence criteria
    % % confident
    
    r.pnum(i)  = pdata{i}(1,pnum);
    r.t1HR(i)  = (sum(pdata{i}(:,class) == 1 & pdata{i}(:,gramt) == 1)+0.5) / (sum(pdata{i}(:,gramt) == 1)+1);
    r.t1FA(i)  = (sum(pdata{i}(:,class) == 1 & pdata{i}(:,gramt) == 0)+0.5) / (sum(pdata{i}(:,gramt) == 0)+1);
    r.t2HR(i)  = (sum(pdata{i}(:,binConf) == 1 & pdata{i}(:,accur) == 1)+0.5) / (sum(pdata{i}(:,accur) == 1)+1);
    r.t2FA(i)  = (sum(pdata{i}(:,binConf) == 1 & pdata{i}(:,accur) == 0)+0.5) / (sum(pdata{i}(:,accur) == 0)+1);
    
    r.t1HR_1st3q(i)  = (sum(pdata{i}(1:end3Q,class) == 1 & pdata{i}(1:end3Q,gramt) == 1)+0.5) / (sum(pdata{i}(1:end3Q,gramt) == 1)+1);
    r.t1FA_1st3q(i)  = (sum(pdata{i}(1:end3Q,class) == 1 & pdata{i}(1:end3Q,gramt) == 0)+0.5) / (sum(pdata{i}(1:end3Q,gramt) == 0)+1);
    r.t2HR_1st3q(i)  = (sum(pdata{i}(1:end3Q,binConf) == 1 & pdata{i}(1:end3Q,accur) == 1)+0.5) / (sum(pdata{i}(1:end3Q,accur) == 1)+1);
    r.t2FA_1st3q(i)  = (sum(pdata{i}(1:end3Q,binConf) == 1 & pdata{i}(1:end3Q,accur) == 0)+0.5) / (sum(pdata{i}(1:end3Q,accur) == 0)+1);
    
    r.t1HR_lastq(i)  = (sum(pdata{i}(st4Q:numRs,class) == 1 & pdata{i}(st4Q:numRs,gramt) == 1)+0.5) / (sum(pdata{i}(st4Q:numRs,gramt) == 1)+1);
    r.t1FA_lastq(i)  = (sum(pdata{i}(st4Q:numRs,class) == 1 & pdata{i}(st4Q:numRs,gramt) == 0)+0.5) / (sum(pdata{i}(st4Q:numRs,gramt) == 0)+1);
    r.t2HR_lastq(i)  = (sum(pdata{i}(st4Q:numRs,binConf) == 1 & pdata{i}(st4Q:numRs,accur) == 1)+0.5) / (sum(pdata{i}(st4Q:numRs,accur) == 1)+1);
    r.t2FA_lastq(i)  = (sum(pdata{i}(st4Q:numRs,binConf) == 1 & pdata{i}(st4Q:numRs,accur) == 0)+0.5) / (sum(pdata{i}(st4Q:numRs,accur) == 0)+1);
    
    % Compute the d-prime
    r.t1Dp(i)  = icdf('Normal',r.t1HR(i),0,1) - icdf('Normal',r.t1FA(i),0,1);
    r.t2Dp(i)  = icdf('Normal',r.t2HR(i),0,1) - icdf('Normal',r.t2FA(i),0,1);
    
    %% ADD CRITERION
    
    r.t1Dp_1st3q(i)  = icdf('Normal',r.t1HR_1st3q(i),0,1) - icdf('Normal',r.t1FA_1st3q(i),0,1);
    r.t2Dp_1st3q(i)  = icdf('Normal',r.t2HR_1st3q(i),0,1) - icdf('Normal',r.t2FA_1st3q(i),0,1);
    r.t1Dp_lastq(i)  = icdf('Normal',r.t1HR_lastq(i),0,1) - icdf('Normal',r.t1FA_lastq(i),0,1);
    r.t2Dp_lastq(i)  = icdf('Normal',r.t2HR_lastq(i),0,1) - icdf('Normal',r.t2FA_lastq(i),0,1);
    
    %% ADD c
    r.t1c_1st3q(i)  = -0.5*(icdf('Normal',r.t1HR_1st3q(i),0,1) - icdf('Normal',r.t1FA_1st3q(i),0,1));
    r.t1c_lastq(i)  = -0.5*(icdf('Normal',r.t1HR_lastq(i),0,1) - icdf('Normal',r.t1FA_lastq(i),0,1));
    
    
    % Remove infinite values.
    r.t1Dp(r.t1Dp == -Inf | r.t1Dp == Inf) = NaN;
    r.t2Dp(r.t2Dp == -Inf | r.t2Dp == Inf) = NaN;
    r.t1Dp_1st3q(r.t1Dp_1st3q == -Inf | r.t1Dp_1st3q == Inf) = NaN;
    r.t2Dp_1st3q(r.t2Dp_1st3q == -Inf | r.t2Dp_1st3q == Inf) = NaN;
    r.t1Dp_lastq(r.t1Dp_lastq == -Inf | r.t1Dp_lastq == Inf) = NaN;
    r.t2Dp_lastq(r.t2Dp_lastq == -Inf | r.t2Dp_lastq == Inf) = NaN;
    
    % meta-d
    
    % 1. get trials
    
    for whichQuarter = 1:2
        switch whichQuarter
            case 1; start = 1; finish = end3Q;
            case 2; start = st4Q; finish = numRs;
        end
        
        stim = pdata{i}(start:finish,gramt);
        resp = pdata{i}(start:finish,class);
        conf = pdata{i}(start:finish,binConf)+1;
        
        [nR_S1 , nR_S2] = trials2counts(stim,resp,conf,2,1);
        try
            out             = fit_meta_d_MLE(nR_S1 , nR_S2, 1);
        catch err; out.meta_da = nan; out.t2ca_rS1 = nan;  out.t2ca_rS2 = nan; out.meta_da = nan;%rethrow(err)
        end
        mdistR1       = ( out.meta_ca - out.t2ca_rS1 )/out.meta_da;
        mdistR2       = ( out.t2ca_rS2 - out.meta_ca )/out.meta_da;
        
        switch whichQuarter
            case 1; r.mdistR1_1st3q(i,1) = mdistR1; r.mdistR2_1st3q(i,1) = mdistR2; r.metad_1st3q(i,1) = out.meta_da; r.conf_1st3q(i,1) = mean(conf);
            case 2; r.mdistR1_lastq(i,1) = mdistR1; r.mdistR2_lastq(i,1) = mdistR2; r.metad_lastq(i,1) = out.meta_da; r.conf_lastq(i,1) = mean(conf);
        end
        
    end
    
end
save('Scott_et_al_data.mat')
results = [...
    r.pnum ...
    r.t1Dp_1st3q ...
    r.t1Dp_lastq ...
    r.t1c_1st3q ...
    r.t1c_lastq ...
    r.conf_1st3q(i) ...
    r.conf_lastq(i) ...
    r.metad_1st3q ...
    r.metad_lastq ...
    r.mdistR1_1st3q ...
    r.mdistR2_1st3q ...
    r.mdistR1_lastq ...
    r.mdistR2_lastq ...
    
    ];

xlswrite('ResultFile.xls',results);
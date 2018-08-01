% function Scott_et_al_analysis()
%
% This function runs through the analysis presented in our paper "Measuring
% metacognitive thresholds using signal detection theory".
%
% The data it uses is from:
% Scott, R. B., Dienes, Z., Barrett, A. B., Bor, D., & Seth, A. K. (2014).
% Blind insight: metacognitive discrimination despite chance task performance.
% Psychological science, 25(12), 2199-2208.
%
% The data and code from these authors can be found here:
% https://osf.io/ivdk4/files/.
%
% Maxine Sherman 23.07.2018
%


function Scott_et_al_analysis()

%% ========================================================================
%                               Initialise
%  ========================================================================

% clear workspace
clc

% Load analysed data
load scott_et_al_data r

% Collect the variables we need
data = [r.t1Dp_1st3q,...                                       % type 1 d' for the selection trials
        r.t1Dp_lastq,...                                           % type 1 d' for the test trials
        log(r.conf_lastq),...                                      % confidence for the test trials
        log(abs(r.mdistR1_lastq + r.mdistR2_lastq))];              % mdist for the test trials (we add so it's not response-conditional)


% exclude subjects with nan/inf anywhere
keep = [];
for i = 1:size(data,1)
    if isnan(prod(data(i,:))) || isinf(prod(data(i,:))) || ~isreal(prod(data(i,:))) % note: we need to exclude subjs with conf = 0 because we can't take their log
    else keep = [keep;i];
    end
end
data = real(data(keep,:));


% get type 1 dprime from the selection trials
dprime_selection = data(:,1);

% exclude subjects who's metad in selection was too low (ruins m-dist
% calculations)
bad_metad = abs(data(:,2)) < 0.1;

% of the remaining subjects, split according to their dprime in selection
% trials
d_0       = find( dprime_selection<=0 & ~bad_metad);
d_1       = find( dprime_selection>0  & ~bad_metad);


%% ========================================================================
%                               Analyse
%  ========================================================================

% Analysis 1: check that type 1 d' in the selection trials doesn't differ
% between the excluded and retained participants
disp('% ================================================== %')
disp('1. Type 1 d in selection trials for excluded vs. retained subjects')

X = data( find(bad_metad)   ,1 );     % excluded
Y = data( find(~bad_metad)  ,1 );     % retained

[p, t, df, vartype] = bootstrap(X,Y);
disp(['t(' num2str(df) ') = ' num2str(t) ', p = ' num2str(p)])
disp(['Levenes result: ' vartype])




% Analysis 2: log confidence for the "insight" vs "blind insight" group
disp('% ================================================== %')
disp('2. Log confidence for the "insight" vs "blind insight" group')

X = data( d_0  , 3 );     % "blind insight"
Y = data( d_1  , 3 );     % "insight"

[p, t, df, vartype] = bootstrap(X,Y);
disp(['t(' num2str(df) ') = ' num2str(t) ', p = ' num2str(p)])
disp(['Levenes result: ' vartype])



% Analysis 2: log mdist for the "insight" vs "blind insight" group
disp('% ================================================== %')
disp('3. Log mdist for the "insight" vs "blind insight" group')

X = data( d_0  , 4 );     % "blind insight"
Y = data( d_1  , 4 );     % "insight"

[p, t, df, vartype] = bootstrap(X,Y);
disp(['t(' num2str(df) ') = ' num2str(t) ', p = ' num2str(p)])
disp(['Levenes result: ' vartype])
disp('% ================================================== %')


%% ========================================================================
%                               Plot
%  ========================================================================

% Initialise figure
figure;

for iDV = 1:2
    
    % Get groups
    switch iDV
        case 1;
            d0   = data( d_0  , 3 );     % "blind insight"
            d1   = data( d_1  , 3 );     % "insight"
            name = 'log(prop. confident)';
        case 2;
            d0   = data( d_0  , 4 );     % "blind insight"
            d1   = data( d_1  , 4 );     % "insight"
            name = 'log(m-dist)';
    end
    
    % Subplot
    subplot(1,2,iDV)
    
    % Get mean & standard error
    M    = [nanmean(d0) , nanmean(d1)];
    SE   = [nanstd(d0)/sqrt(numel(d0)) , nanstd(d1)/sqrt(numel(d1))];
    
    % Plot errorbar
    errorbar([1:2],M,SE,'k','LineWidth',2);
    hold on;
    scatter(1:2,M,50,'k','filled');
    set(gca,'XTick',[1,2],'XTickLabel',{'Chance','Above Chance'});
    ylabel(name)
    xlabel('Type 1 d prime (selection trials)')
end
end


%% =========================================================================
%                       [p , T, df] = bootstrap(X,Y)
%  ========================================================================

function [p , T, df, vartype] = bootstrap(X,Y)

rng(0); % set rng so you can replicate our results precisely

% get group numbers
nX = numel(X);
nY = numel(Y);

% make X and Y column vectors
X  = reshape(X,nX,1);
Y  = reshape(Y,nY,1);

% determine empirical mean difference
empirical_mean_diff = mean(X)-mean(Y);

% Do Levene's test
L      = [ X ; Y ];
group  = [zeros(size(X));ones(size(Y))];
levene = vartestn(L,group,'TestType','LeveneAbsolute','Display','off');

switch levene < 0.05
    case 1; vartype = 'unequal';
    case 0; vartype = 'equal';
end

% t-test
[~,~,~,stats] = ttest2(X,Y,'vartype',vartype);
T             = stats.tstat;
df            = stats.df;

allData  = [X;Y];

nSamples = 10000;
p        = 0;

for i = 1:nSamples
    
    % assign groups
    idxA = randsample(1:(nX+nY),nX);
    xA   = allData(idxA);
    
    idxB = setdiff(1:(nX+nY),idxA);
    xB   = allData(idxB);
    
    % t-test
    [~,~,~,stats] = ttest2(xA,xB,'vartype',vartype);
    
    % is difference bigger?
    switch sign(T)
        case 1 % difference is positive
            if stats.tstat >= T
                p = p + 1;
            end
        case -1 % difference is negative
            if stats.tstat <= T
                p = p + 1;
            end
    end
    
end

p = p/nSamples;
end
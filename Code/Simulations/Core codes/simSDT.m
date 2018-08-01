% function [real, sim nR_S1 nR_S2]     = simSDT(p,plotme)
%
% This function simulates type 1 and 2 sdt measures from the parameter
% structure p. If plotme is set to 1 then it will plot the type 1 model
% from which responses are simulated.


function [real, sim nR_S1 nR_S2]     = simSDT(p,plotme)

%% ========================================================================
%           Put key parameters into the output structure 'Real'
%  ========================================================================

real.d       = p.mu2;
real.c       = p.c;
real.tauY    = p.tauY;
real.tauN    = p.tauN;
real.c_norm  = real.c/real.d;

%% ========================================================================
%                          Run simulation
%  ========================================================================

% Get values for X on the S1 trials
S1  = normrnd(-p.mu2/2,p.sigma1,p.ntrials,1);

% Get values for X on the S2 trials
S2 = normrnd(p.mu2/2,p.sigma2,p.ntrials*p.E,1);

% Get responses on the trials
[Hits,Misses,FAs,CRs,SDT,Conf, t2sdt] = gettrials(S1,S2,p);

% Load in the simulated SD of the two likelihood functions
sim.s1 = std(S1);
sim.s2 = std(S2);

% Shall we pad the cells? Padding hits, misses, false alarms and correct
% rejections by 0.5 allows estimation of SDT parameters when one of those
% is zero. This padding method acts like a Bayesian prior on d' = c = 0
if p.padCells;
    padding = 0.5;
else
    padding = 0;
end

% Load hit rate, false alarm rate, type 1 d' and type 1 c into sim.
sim.HR    = (Hits+padding)/(Hits + Misses + 2*padding);
sim.FAR   = (FAs+padding)/(FAs + CRs + 2*padding);
sim.d     = p.sigma2*norminv(sim.HR)-norminv(sim.FAR);
sim.c     = -0.5*(norminv(sim.HR)+norminv(sim.FAR));

% Also add in response-conditional confidence
sim.conf  = [t2sdt.conf_rS1 t2sdt.conf_rS2];

% Get parameters of the meta-d model
sim.output = computeMetaD(sim.c,p.padCells,p.sigma2,t2sdt);

% Plot [optional]
if plotme
    
    for i_plot = 1:2
        switch i_plot
            case 1; P = real; T = 'True';
                tauY = real.tauY;
                tauN = real.tauN;
                sigma = p.sigma2;
            case 2; P = sim;  T = 'Simulated';
                tauY = sim.output.taufit.tau_raw_min;
                tauN = sim.output.taufit.tau_raw_plus;
                sigma = 1;
        end
        
        subplot(2,1,i_plot)
        x = [-3:.1:3];
        gA = normpdf(x,-P.d/2,1);
        gP = normpdf(x,P.d/2,sigma);
        plot(x,gA,'r',x,gP,'b');
        hold on;
        plot(repmat(P.c,1,numel(0:0.1:1)),0:0.1:1,'k')
        plot(repmat(tauY,1,numel(0:0.1:1)),0:0.1:1,'k--')
        plot(repmat(tauN,1,numel(0:0.1:1)),0:0.1:1,'k--')
        xlabel('X'); ylabel('Probability')
        
        title([T ' distribution'])
    end
end
end

%% ========================================================================
%                              HELPER 1
%     [Hits,Misses,FAs,CRs,SDT,Conf t2sdt] = gettrials(Absent,Present,p)
%  ========================================================================
function [Hits,Misses,FAs,CRs,SDT,Conf t2sdt] = gettrials(S1,S2,p)




%% ========================================================================
%                             Initialise
%  ========================================================================

SDT{1}  = nan(size(S2));
Conf{1} = nan(size(S2));
SDT{2}  = nan(size(S1));
Conf{2} = nan(size(S1));

%% ========================================================================
%                             Start with S2 trials
%  ========================================================================
X = S2;

% Define X2 (signal for confidence - will be different from X under
% enhanced/degraded signal models
eval(['X2 = normrnd(' p.signalChange{2} ',size(X,1),size(X,2));']);

% Determine the type 1 report (R1 or R2) by comparing X to criterion
SDT{1}(X>p.c) = 1; % passes threshold => report R2 (hit)
SDT{1}(X<p.c) = 2; % doesnt pass threshold => report R1 (miss)

% Reported R2 (hits)
Conf{1}(X2 < p.c       & SDT{1} == 1)           = 0; % conf signal falls on other side of thresh => change of mind
Conf{1}(X2 > p.tauY    & SDT{1} == 1)           = 1; % conf signal surpasses c2 => conf
Conf{1}(X2 > p.c & X2 < p.tauY & SDT{1} == 1)   = 0; % conf signal doesnt pass c2 => guess

% Reported R1 (misses)
Conf{1}(X2 > p.c       & SDT{1} == 2)            = 0; % conf signal falls on other side of thresh => change of mind
Conf{1}(X2 < p.tauN    & SDT{1} == 2)            = 1; % conf signal surpasses c2 => conf
Conf{1}(X2 < p.c & X2 > p.tauN & SDT{1} == 2)    = 0; % conf signal doesnt pass c2 => guess


%% ========================================================================
%                             Now do S1 trials
%  ========================================================================
X = S1;
eval(['X2 = normrnd(' p.signalChange{1} ',size(X,1),size(X,2));']);

SDT{2}(X>p.c) = 3; % passes threshold => report S2 (false alarm)
SDT{2}(X<p.c) = 4; % doesnt pass threshold => report S1 (correct rejection)

% Reported yes (false alarms)
Conf{2}(X2 < p.c       & SDT{2} == 3)          = 0;% conf signal falls on other side of thresh => change of mind
Conf{2}(X2 > p.tauY    & SDT{2} == 3)          = 1;% conf signal surpasses c2 => conf
Conf{2}(X2 > p.c & X2 < p.tauY & SDT{2} == 3)  = 0;% conf signal doesnt pass c2 => guess

% Reported no (correct rejections)
Conf{2}(X2 > p.c       & SDT{2} == 4)          = 0; % conf signal falls on other side of thresh => change of mind
Conf{2}(X2 < p.tauN    & SDT{2} == 4)          = 1;% conf signal surpasses c2 => conf
Conf{2}(X2 < p.c & X2 > p.tauN & SDT{2} == 4)  = 0;% conf signal doesnt pass c2 => guess


%% ========================================================================
%                            Sum over response types
%  ========================================================================

% Count h, m, fa, cr
Hits    = sum(S2 > p.c) ;
Misses  = sum(S2 < p.c) ;
FAs     = sum(S1 > p.c)  ;
CRs     = sum(S1 < p.c)  ;

% Get prop. confidence
t2sdt.conf_rS1 = mean(Conf{1});
t2sdt.conf_rS2 = mean(Conf{2});

% Merge results over S1 and S2 trials
SDT  = [SDT{1};SDT{2}];
Conf = [Conf{1};Conf{2}];

% Get HRs and FARs
t2sdt.HR_rS1  = sum( SDT == 4 & Conf == 1 )/sum( SDT == 4); % type 2 HR for R1
t2sdt.HR_rS2  = sum(SDT  == 1 & Conf == 1)/sum(SDT == 1);   % type 2 HR for R2
t2sdt.FAR_rS1 = sum(SDT  == 2 & Conf == 1)/sum(SDT == 2);   % type 2 FAR for R1
t2sdt.FAR_rS2 = sum(SDT  == 3 & Conf == 1)/sum(SDT == 3);   % type 2 FAR for R2
t2sdt.hr      = sum(SDT == 1)/sum(ismember(SDT,[1,2]));     % type 1 HR
t2sdt.far     = sum(SDT == 3)/sum(ismember(SDT,[3,4]));     % type 2 HR

[t2sdt.nR_S1 , t2sdt.nR_S2] = trials2counts(ismember(SDT,[1,2]),ismember(SDT,[1,3]),Conf+1,2,1); % get the input for Maniscalco & Lau's code


end






%% ========================================================================
%                              HELPER 2
%                  fit = computeMetaD(c1,padCells, s ,t2sdt)
%  ========================================================================
function fit = computeMetaD(c1,padCells, s ,t2sdt)

if nargin < 4; s=1; end

% initialise
fit = struct('m_distance',[],'meta_d_fit',[]);

% get meta-d' fit (get both response unconditional and response
% conditional)
fit.meta_d_fit            = fit_meta_d_MLE(t2sdt.nR_S1, t2sdt.nR_S2, 1/s);
fit.rs_meta_d_fit         = fit_rs_meta_d_MLE(t2sdt.nR_S1, t2sdt.nR_S2, 1/s);

% compute metacognitive distance 
c2                = [ fit.meta_d_fit.t2ca_rS1 ; fit.meta_d_fit.t2ca_rS2 ];
metad             = fit.meta_d_fit.meta_da;
c1                = fit.meta_d_fit.meta_ca;

fit.mdist         = [ (c1 - c2(1))/abs(metad) , (c2(2) - c1)/abs(metad) ];

end


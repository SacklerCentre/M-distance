% This runs the analyses for the simulations presented in Sherman et al.
% Measuring Metacognitive Thresholds in Signal Detection Theory
%
% The simulations are run in the same order as they are presented in the
% paper.
%
% Maxine Sherman 2018



%% ========================================================================
%                               Setup
%  ========================================================================
clear all;close all;clc
addpath(genpath(cd));

rng(1);

% Which analyses should we run?
vary_type1_thresholds       = 1;
vary_type2_thresholds       = 1;
vary_type1_dprime           = 1;
vary_trialcount             = 1;
hierarchical_models         = 1;
dialChannel_models          = 1;

plotOn                      = 1;

% Set defaults
default.nTrials             = 100000;
default.dprime              = 2;
default.c                   = 0;
default.tauN                = -default.dprime/4; % distance between confidence threshold and criterion
default.tauY                = default.dprime/4;
default.sigma2              = 1;


% Create output folders
mkdir('Simulation_output');
mkdir('Simulation_output_figures');



%% -----------------------------------------------------------------------
%      SIMULATION 1: m-dist is invariant to criterion
%  -----------------------------------------------------------------------

if vary_type1_thresholds
    
    clc; disp('Running: vary_type1_thresholds');
    clear data;
    parameter_bounds  = -1.5:0.1:1.5;
    parameter_to_vary = 'criterion';
    
    data = simulate_parameter_effects( parameter_bounds , parameter_to_vary , 0 , default  );
    
    save('Simulation_output/vary_type1_thresholds.mat','data','parameter_bounds','parameter_to_vary')
    
    if plotOn
        plotSimulation(data,parameter_bounds,parameter_to_vary);
    end
end

%% -----------------------------------------------------------------------
%      SIMULATION 2: changing taus but fixing criterion
%  -----------------------------------------------------------------------

if vary_type2_thresholds
    
    clear data;
    parameter_bounds  = -1:0.1:1;
    
    for i = 1:2
        switch i
            case 1; parameter_to_vary = 'tauN';
            case 2; parameter_to_vary = 'tauY';
        end
        clc; disp(['Running: vary_type2_thresholds (' parameter_to_vary ')']);
        data = simulate_parameter_effects( parameter_bounds , parameter_to_vary , 0 , default );
        save(['Simulation_output/vary_type2_thresholds_' parameter_to_vary '.mat'],'data','parameter_bounds','parameter_to_vary')
        
        if plotOn
            plotSimulation(data,parameter_bounds,parameter_to_vary);
        end
    end
    
end


%% -----------------------------------------------------------------------
%      SIMULATION 3: changing type 1 sensitivity
%  -----------------------------------------------------------------------

if vary_type1_dprime
    
    clc; disp('Running: vary_type1_dprime');
    clear data;
    parameter_bounds  = 0.5:0.1:3;
    parameter_to_vary = 'dprime';
    
    data = simulate_parameter_effects( parameter_bounds , parameter_to_vary , 0 , default  );
    
    save('Simulation_output/vary_type1_dprime.mat','data','parameter_bounds','parameter_to_vary')
    
    if plotOn
        plotSimulation(data,parameter_bounds,parameter_to_vary);
    end
end   



%% -----------------------------------------------------------------------
%      SIMULATION 4: stability of m-dist over varying trial counts
%  -----------------------------------------------------------------------

if vary_trialcount
    clear data;
    clc; disp('Running: vary_trialcount');
    nSubj               = 15;
    parameter_bounds    = -2:0.2:2;
    parameter_to_vary   = 'criterion';
    nTrials             = 100:100:1000;
    
    
    for i_ntrials = 1:numel(nTrials)
        for i_subj = 1:nSubj
            x = simulate_parameter_effects( parameter_bounds , parameter_to_vary , 0 , default , 'nTrials' , nTrials(i_ntrials));
            results(i_ntrials,i_subj) = x;
        end
    end
    save('Simulation_output/vary_trialcount.mat','results','parameter_bounds','parameter_to_vary')
    if plotOn
        plot_trialCounts_simulation(data,parameter_bounds,parameter_to_vary);
    end
end



%% -----------------------------------------------------------------------
%      SIMULATION 5: Hierarchical models
%  -----------------------------------------------------------------------

if hierarchical_models
    clear data;
    a0                = [0.5:0.1:1.5]; % mu
    s0                = [0:0.1:2];  % sigma
    model_name        = 'enhancement'; % name is misleading & redundant. please ignore.
    parameter_bounds  = 0;
    parameter_to_vary = 'criterion';
    
    % set type 2 criteria at 1 instead of 0.5
    default2 = default; default2.tauN = -1; default2.tauY = 1;
    
    for a = 1:numel(a0)
        for s = 1:numel(s0)
            
            val = { [ 'X * ' num2str(a0(a)) ',' num2str(s0(s)) ] , [ 'X * ' num2str(a0(a)) ',' num2str(s0(s))  ] };
            
            data(a,s) = simulate_parameter_effects( parameter_bounds , parameter_to_vary , 0 , default2, 'donotcatch' , 'signalChange' , val , 'model', 'enhanced');
            
            sprintf(['DONE %d/%d ...'],[s+(numel(s0)*(a-1)) , numel(s0)*numel(a0)])
        end
    end
    param_values = {a0,s0};
    save('Simulation_output/hierarchical_models.mat','data','param_values')
    
    if plotOn
        plotHierarchicalModels(data,param_values);
    end
end



%% -----------------------------------------------------------------------
%      SIMULATION 6: m-dist & mratio under dual channel models
%  -----------------------------------------------------------------------

if dualChannel_models
    clear data;
    b0                = 0.4:0.2:2.2;%0.5:0.1:1.5;
    parameter_bounds  = [-1:0.1:1];
    parameter_to_vary = 'criterion';
    
    for b = 1:numel(b0)
        
        
        val = { [ num2str(-0.5* default.dprime) ',' num2str(b0(b)) ] , [ num2str(0.5*default.dprime) ',' num2str(b0(b))  ] };
        
        data(b) = simulate_parameter_effects( parameter_bounds , parameter_to_vary , 0 , default, 'donotcatch' , 'signalChange' , val , 'model', 'enhanced');
        
        sprintf(['DONE %d/%d ...'],[b numel(b0)])
    end
    
    
    save('Simulation_output/hierarchical_models.mat','data','b0','parameter_bounds')
    if plotOn
        plotDualChannelModels(data,parameter_bounds,{b0 , parameter_bouds});
    end
end





% Vary (a) criterion and (b) sensitivity. Determine effects on m-distance


function data = simulate_parameter_effects( param_values , param_name , plot_on , default , varargin )

%% Setup
addpath(genpath('InsideFunctions'))

%% Determine starting parameters
dprime      = default.dprime;
criterion   = default.c;
tauN        = default.tauN;
tauY        = default.tauY;
nTrials     = default.nTrials;
sigma2      = default.sigma2;

for i = 1:numel(param_values)
    
    % get value
    param_value = param_values(i);
    
    % Vary parameter of interest
    eval([ param_name  ' = '  num2str(param_value)  ';' ]);

    % Set parameters
    params = setParameters('nTrials', nTrials , 'mu2' , dprime , 'c' , criterion ,'tauY', tauY , 'tauN' , tauN , 'sigma2' , sigma2 , 'padCells',0);
    
    % Any other parameters to change?
    if ~isempty(varargin)
        if strcmpi(varargin{1},'donotcatch')
            eval(['params.' varargin{2} ' =  (varargin{3}) ;' ]);
        else
            for j = 1:2:numel(varargin)
                try
                    eval(['params.' varargin{j} ' = ' num2str(varargin{j+1}) ';' ]);
                catch
                    eval(['params.' varargin{j} ' =  (varargin{j+1}) ;' ]);
                end
            end
        end
    end
    
    % Simulate data
    [ real, sim ] = simSDT( params , 0 );
    
    % Load into data structure
    data.real( i )       = real;
    data.sim( i )        = sim;
    data.params( i )     = params;
    data.values          = param_values;
    
    
    % update researcher
    sprintf(['finished %d / %d \n'],[i numel(param_values)])

end
mkdir('output');
save(['output/' param_name '.mat'],'data');

%% Plot

if plot_on
    plot_results(data, param_name);
end





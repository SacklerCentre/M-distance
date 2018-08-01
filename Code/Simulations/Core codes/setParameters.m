% function p = setParameters(varargin)
%
% This function sets the SDT parameters for the simulation
%
%
% Maxine Sherman 23.07.2018


function p = setParameters(varargin)

if strcmpi(varargin,'default')
    varargin = cell(12*2,1);%{'varmodel',[],'model',[],'sigma2',[],'mu2',[],'theta',[],'tauY',[],'tauN',[],'c',[],'ntrials',[],'E',[],'degredation',[],'enhancement',[]};
else
    n = numel(varargin);
    n = 24-numel(varargin);
    varargin  = [varargin cell(1,n)];
end


% check for replacements
for k = 1:numel(varargin)
    
    k = find(strcmpi(varargin,'varmodel'));   if ~isempty(k);  varmodel   = varargin{k+1}; k = k + 1;   else varmodel = 'EV';       end
    k = find(strcmpi(varargin,'sigma1'));     if ~isempty(k);  sigma1     = varargin{k+1}; k = k + 1;   else sigma1   = 1;          end
    k = find(strcmpi(varargin,'sigma2'));     if ~isempty(k);  sigma2     = varargin{k+1}; k = k + 1;   else sigma2   = 1;          end
    k = find(strcmpi(varargin,'mu2'));        if ~isempty(k);  mu2        = varargin{k+1}; k = k + 1;   else mu2      = 2;          end
    k = find(strcmpi(varargin,'theta'));      if ~isempty(k);  theta      = varargin{k+1}; k = k + 1;   else theta    = 0;          end
    k = find(strcmpi(varargin,'c'));          if ~isempty(k);  c          = varargin{k+1}; k = k + 1;   else c        = 0;          end
    k = find(strcmpi(varargin,'tauY'));       if ~isempty(k);  tauY       = varargin{k+1}; k = k + 1;   else tauY     = .5;         end
    k = find(strcmpi(varargin,'tauN'));       if ~isempty(k);  tauN       = varargin{k+1}; k = k + 1;   else tauN     = -.5;        end
    k = find(strcmpi(varargin,'ntrials'));    if ~isempty(k);  ntrials    = varargin{k+1}; k = k + 1;   else ntrials  = 30000;      end
    k = find(strcmpi(varargin,'E'));          if ~isempty(k);  E          = varargin{k+1}; k = k + 1;   else E        = 1;          end
    k = find(strcmpi(varargin,'padCells'));   if ~isempty(k);  padCells   = varargin{k+1}; k = k + 1;   else padCells = 0;          end
    k = find(strcmpi(varargin,'whichMetad')); if ~isempty(k);  whichMetad = varargin{k+1}; k = k + 1;   else whichMetad = 'Lau';    end
    
    
    %% Set signal change
    %
    % signal degradation on a hierarchical model means having an input that
    % looks like this:
    %
    %  signalChange = {[normrnd(X - a0*mu2,  s0)],...
    %                 [normrnd(X +  a0*mu2 , s0) ]};
    %
    % This corresponds to degrading the signal with some gaussian noise
    % centred on the evidence presented. 
    %
    % To model a dual-channel model, you center X2 about d' with some
    % variance, like this:
    %
    %  signalChange = {[normrnd( -mu2/2 , s0)],...
    %                 [  normrnd( mu2/2  , s0]};
    
    k = find(strcmpi(varargin,'signalChange'));
    if ~isempty(k);
        X = varargin{k+1}; k = k + 1;
    else
        X = { ['X , 0'] ; ['X , 0'] };
    end
    signalChange = X;
    
end


% use da?
if sigma2~=1; mu2 = mu2/sqrt((0.5*(1+sigma2))); end
args = {'varmodel','sigma1','sigma2','mu2','theta','tauY','tauN','c','ntrials','E','padCells','signalChange'};


%% ---- LOAD IN ---- %
for i_arg = 1:numel(args);
    eval(['p.' args{i_arg} ' = ' args{i_arg} ';']);
end
p.tauN = p.tauN + p.c;
p.tauY = p.tauY + p.c;


end


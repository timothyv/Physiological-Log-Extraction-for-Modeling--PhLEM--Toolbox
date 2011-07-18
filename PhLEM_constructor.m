function PhLEM = PhLEM_constructor(varargin)

% PhLEM = PhLEM_constructor('input_name',input_value,...)
% Main PhLEM structure constructor.  Generates the PhLEM object as needed
% for the function PhLEM_process
%
% Inputs describe the name of the signal and a 1xN vector pulled from
% your local loader functions
%
% Inputs: 
%   input_name     = name of signal.
% 
%       Viable naming conventions are
%           'ECG'  -> phase of heartbeats
%           'RESP' -> Pressure waveform of breathing
%           'PPG'  -> Pulseox readout
%           'Time' -> Timestamps in seconds for recorded data (required)
%           'a_opts' -> object listing processing parameters for each data
%                        field (e.g. ecg, resp, PPO).  If none given then
%                        the default PhLEM_ana_opts is called
%
%   input_value     = 1xN array of datapoints for the signal
%
%   
%
% ASSUMPTIONS:
%  1) All inputs are the same length.
%  2) Time vector is in seconds.
%  3) All clipping to match to your imaging data volumes is done before
%     this stage.
%
% Written by T. Verstynen (January 2008).
% Liscensed under the GPL public license (version 3.0)

  
if length(varargin)
  
  if ~isstruct(varargin{1})
    if mod(length(varargin),2)
      error('Not an even matching of input_labels and input_values');
    end

    % Put the appropriate values in the right variables
    for i = 1:2:length(varargin) 
      optlabel = varargin{i};
      optvalue = varargin{i+1};
      dtmp.data = optvalue;
      
      switch optlabel
       case {'Time','time'}
        PhLEM.Time = optvalue;
       case 'a_opts';            
        a_opts = optvalue;
       case {'ecg', 'resp', 'ppghigh','ppglow','ppgmed','ppg'}            
        eval(sprintf('PhLEM.%s = dtmp;',upper(optlabel)));
       otherwise
        eval(sprintf('PhLEM.%s = dtmp;',optlabel));
      end
      
    end
  else
    
    % This lets me pass in an argument struct
    
    argStruct = varargin{1};
    
    for thisField = fields(argStruct)'
      optlabel = thisField{:};
      optvalue = getfield(argStruct,optlabel);
      dtmp.data = optvalue;
      
      switch lower(optlabel)
       case 'time'
        PhLEM.Time = optvalue;
       case 'a_opts'
        a_opts = optvalue;
       case {'ecg','resp','ppghigh','ppglow','ppgmed','ppg'}
        eval(sprintf('PhLEM.%s = dtmp;',upper(optlabel)));
       otherwise
        eval(sprintf('PhLEM.%s = dtmp;',optlabel));
      end
      
    end
    
  end
        
end


% Setup default options for analysis (edit PhLEM_ana_opts.m to changes
% these)
if ~exist('a_opts');
    a_opts = PhLEM_ana_opts([],500);
end

% Go through and process each field of PhLEM
fieldNames = fields(a_opts);
fieldNames = setdiff(fieldNames,'Time'); % Remove 'Time' field from analysis


% Put in each appropriate
for f = 1:length(fieldNames)
    dtmp = eval(['PhLEM.' fieldNames{f} ';']);
    dtmp.opts = eval(sprintf('a_opts.%s;',fieldNames{f}));       
    eval(sprintf('PhLEM.%s = dtmp;',fieldNames{f}));
end;

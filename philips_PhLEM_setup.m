function PhLEM = philips_PhLEM_setup(outfile,varargin);

% function PhLEM = philips_PhLEM_setup(outfile,options);
% 
% Example data prep function for data off of the Philips 3T. 
%
%
%DOES NOT DEAL WITH ECG CHANNEL 07/03/2009  - P Summers
%
%
% Options can be specified as name/value pairs as follows:
%   'logfile',<logfilepath>  = file with ecg, resp, pulse ox, and
%   triggering data
%
%   Additional Options that are also defined in PhLEM_defaults
%   'tr',<TR_inmilliseconds>  = time to acquire 1 volume
%   'ipat',<boolean>          = whether iPat (mSense, Grappa) was engaged
%   'ecg_loader_fcn',<handle> = pointer to ECG loader function.
%   'resp_loader_fcn',<handle>= pointer to RESP loader function.
%   'ppo_loader_fcn',<handle> = pointer to PPO loader function.
%
% If no options are specified, user will be prompted for file input.
%
%  
% This file is a loader and data prep function for the PhLEM 
% toolbox.  PhLEM creates a data structure by composed of different
% physiological measures taken while scanning.  This file will
% spit out a PhLEM structure which is used by the main function
% "make_physio_regresssors" adapted to fsl.
%
% The PhLEM structure will generate a set of subfields for each datatype
% indicated (e.g. ECG, RESP, PPO). For each datatype you have the
% following sub-fields
%  
%   data  = 1xN array of raw data signals
%   opts  = analysis and processing options from ana_opts function
%   events = 1xN array of indicator variables for event peaks, where 1 =
%            event, 0 inter-event times and N is the number of samples
%            in your physio signal
%   smoothed_data = gaussian smoothed data,with a kernel FWHM defined in
%                   the opts structure.
%   phaseK = the Kth phase expansion of the peak events.  
%
% In addition to fields for each physio signal, the PhLEM structure
% maintains an array of timestampes (TIME) in seconds.
%
% Written by T. Verstynen, J. Diedrichsen, and J. Schlerf (Sept 2008)
% Liscensed under the GPL public license (version 3.0)

warning off;

% Grab default values now
global defaults

PhLEM_defaults;

default_fields = fieldnames(defaults);
for i = 1:length(default_fields)
  eval(sprintf('%s=defaults.%s;',default_fields{i},default_fields{i}));
end
	
% Flag the output
if nargin < 1
  outfileMissing = 1;
else
  outfileMissing = 0;
end
  

if nargin > 1 & ~isempty(varargin)
  
  if ~isstruct(varargin{1})
  
    for i = 1:2:length(varargin)
      optlabel = varargin{i};
      optvalue = varargin{i+1};      
      
      if ischar(optvalue)
        eval(sprintf('%s = ''%s'';',optlabel,optvalue));
      else
        
        if islogical(optvalue), optvalue = double(optvalue); end
        if isnumeric(optvalue), optvalue = num2str(optvalue); end

        if isempty(optvalue) % if skipping a data type
          eval(sprintf('%s = [];',optlabel));
        else
          eval(sprintf('%s = %s;',optlabel, optvalue));
        end;
      
      end
      
    end
  
  else
    
    argstruct = varargin{1};
    
    for optlabel = fields(argstruct)'
    
      optvalue = getfield(argstruct,optlabel{:});
      if isnumeric(optvalue), optvalue = num2str(optvalue); end
      eval(sprintf('%s = %s;',lower(optlabel),optvalue));
      
    end
    
  end
    
end

% Flags about whether or not we have files
ECGfileSupplied  = exist('ecgfile','var');
PPOfileSupplied  = exist('ppofile','var');
RESPfileSupplied = exist('respfile','var');

% Redundant flagging in lieu of booleans
usingECG = ECGfileSupplied;
usingPPO = PPOfileSupplied;
usingRESP = RESPfileSupplied;

% Now prompt for log files based on what was not collected as input
% HeartRate information
if ~any([usingECG usingPPO usingRESP]) % no input files supplied ask for all
  if isequal(questdlg('Do you have a valid ECG file?','ECG',...
                      'Yes','No','No'),'Yes')
    ecgfile = localUIfileLoader('Select ECG File');
    usingECG = 1;
  end

  if isequal(questdlg('Do you have a valid PPO file?','PPO',...
                      'Yes','No','No'),'Yes')
    ppofile = localUIfileLoader('Select PPO File');
    usingPPO = 1;
  end

  if isequal(questdlg('Do you have a valid RESP file?','RESP',...
                      'Yes','No','No'),'Yes')
    respfile = localUIfileLoader('Select RESP File');
    usingRESP = 1;
  end

end

% TR time
if ~exist('tr','var')
  % get TR from user
  tr = inputdlg(['Please provide TR, in milliseconds ' ...
                 '(ie, 2 seconds would be 2000):']);
  tr = eval(tr{1});
end
  
% Find out how many dummy scans to take out
if ~exist('NDum','var')
  NDum = inputdlg(['How many dummy scans?']);
  NDum = eval(NDum{1});
end

opts_cell_argument = {}; % Will use this below

% Initialize PhLEM structure
if usingECG
  opts_cell_argument{end+1} = 'ECG';
  ecg_loader_hand = str2func(ecg_loader_fcn);
  [ecg ecg_time] = ecg_loader_hand(ecgfile,'TR',tr,'DOSAVE',save_data,...
    'Hz',ecg_hz);
  tmpPhLEM.ecg = ecg;
end

if usingPPO
  opts_cell_argument{end+1} = 'PPOhigh';
  opts_cell_argument{end+1} = 'PPOlow';
  ppo_loader_hand = str2func(ppo_loader_fcn);
  [ppo ppo_time] = ppo_loader_hand(ppofile,'TR',tr,'DOSAVE',save_data,...
				   'iPatOn',ipat,'Hz',ppo_hz);

  if exist('ecg_time','var') & (ecg_hz > ppo_hz)
    % Since ECG is sampled at a higher frequency, expand PPO to be the
    % same length
    ppo =  PhLEM_ts_match(ppo',length(ecg_time));
    ppo_time = ecg_time;
  end

  tmpPhLEM.ppohigh = ppo;
  tmpPhLEM.ppolow  = ppo;
end

if usingRESP
  opts_cell_argument{end+1} = 'RESP';
  resp_loader_hand = str2func(resp_loader_fcn);
  
  [resp resp_time] = resp_loader_hand(respfile,'TR',tr,'DOSAVE',save_data,...
				   'NDUM',NDum,'Hz',resp_hz);
  
  if exist('ecg_time','var')& (ecg_hz > ppo_hz)
    % Since ECG is sampled at a higher frequency, expand RESP to be the
    % same length
    resp = PhLEM_ts_match(resp',length(ecg_time));
    resp_time = ecg_time;
  end

  tmpPhLEM.resp = resp;
end

peak_hz = max([ecg_hz*usingECG resp_hz*usingRESP ppo_hz*usingPPO]);

if exist('ecg_time','var');
  tmpPhLEM.time = ecg_time;
elseif exist('ppo_time','var');
  tmpPhLEM.time = ppo_time;
else
  tmpPhLEM.time = resp_time;
end

% Get the (typical) analysis options
tmpPhLEM.a_opts = PhLEM_ana_opts(opts_cell_argument,peak_hz);

%PhLEM = PhLEM_constructor('ECG',ecg,'RESP',resp,...
%			  'PPOhigh',ppo,'PPOlow',ppo,...
%			  'Time',time,...
%			  'a_opts',a_opts);

% Construct the full PhLEM structure
PhLEM = PhLEM_constructor(tmpPhLEM);

%P = PhLEM_constructor(PrePhLEM);
PhLEM = PhLEM_process(PhLEM);
PhLEM = PhLEM_analysis(PhLEM);

if outfileMissing
  uisave('PhLEM');
  outfile = 'you specified.';
else
  save(outfile,'PhLEM');
end

if results_plot 

  if usingECG
    figure; PhLEM_resultsplot(PhLEM,'ECG');
  end

  if usingRESP 
    figure; PhLEM_resultsplot(PhLEM,'RESP');
  end

  if usingPPO 
    figure; PhLEM_resultsplot(PhLEM,'PPOhigh');
    figure; PhLEM_resultsplot(PhLEM,'PPOlow');
  end


  msgbox(['Results have been plotted.  Red lines are the identified event markers.  '...
    'If you like them, they are saved as ' ...
          outfile ['  If not, try adjusting some of the analysis options in' ...
       [' the ana_opts.m file.  And of course, hack around on' ...
        ' the source code as you see fit. -- PhLEM Team']]])
end

if usingRESP
    [regrRESP, regrRESP_name ] = make_physio_regressors('RESP', PhLEM, 'TR', tr, 'DOSAVE', save_data, 'ana_type', 'fourier');
end
if usingPPO
    [regrPPO, regrPPO_name ] = make_physio_regressors('PPOHIGH', PhLEM, 'TR', tr, 'DOSAVE', save_data, 'ana_type', 'fourier');
end


%regrs = [ regrRESP regrPPO ];
%regrNames = [ regrRESP_name regrPPO_name ]; 

[cnfdfilename,cnfdfilepath] = uiputfile('physcnfd.txt','Save file name');
File = fullfile(cnfdfilepath,cnfdfilename);

% save(File,'regrNames','regrs','-ASCII')

%cnfdfilename = inputdlg(['Please give the name of the file to save the harmonics to:']);
%cnfdfilename = eval(cnfdfilename{1});

fod = fopen(File, 'wt');
npts = length(regrRESP);
for i=1:npts
  fprintf(fod, '%6f \t ', regrRESP(i,:), regrPPO(i,:));
  fprintf(fod, '\n');
end;  
fclose(fod);


[cnfdfilehdrname,cnfdfilepath] = uiputfile('physcnfdhdr.txt','Save list of regressors: file name');
File = fullfile(cnfdfilepath,cnfdfilehdrname);
fud = fopen(File, 'wt');
fprintf(fud, '%s', regrRESP_name{:}, regrPPO_name{:});
fprintf(fud, '\n');
fclose(fud);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filepath = localUIfileLoader(prompt);

% Wrap UI calls for readability above
  
  if exist('spm_select') % use SPM5 loader if we have it
    filepath = spm_select(1,'Any',prompt);
    return
  end
  
  if exist('spm_get') % use SPM2 loader if it's here but not SPM5
    filepath = spm_get(1,'*',prompt);
    return
  end
  
  % if we don't have SPM, use matlab built-ins
  [file,path] = uigetfile('*.log', prompt);
  filepath = fullfile(path,file);
  return

function [signal, time] = philips_ECGload2(fname,varargin)

%function [ecg time] = ECG_load(fname, optname, optval);
%
% A loader for the ECG data coming off of the Philips 3T system.
%  
%   INPUT:
%   fname = name of physiological monitoring log file with respiration information
%   opts        = options on extracting the phyio signals
%
%     a) 'Hz': Sampling Rate of the signal (500Hz default)
%
%     b) 'TR': Sampling Rate of scans (2000ms default)
%
%     c) 'NDUM' : Number of DUMmy scans
%
%     d) 'DOSAVE' : Flag whether to SAVE the data vector
%     
%   OUTPUT
%   ecg   = vector of the ecg waveform
%   time   = vector of timestamps (in seconds)
%
%  
% Written by T. Verstynen (November 2007) and J Schlerf (August 2008)
%
% Rewritten to optimize loading time taking heavily from code written by V.
% Deshpande and J. Grinstead, Siemens Medical Solutions (March 2009) 
% Adapted to Philips physiology log files by P. Summers (March 2009) - not
% tested but is only a 2 character change from the ppo loader which seems to work! 
%
% Liscensed under the GPL public license (version 3.0)

% Liscensed under the GPL public license (version 3.0)

if nargin < 1 || isempty(fname);
  [filename filepath] = uigetfile('*.log','Select log File');
  fname = fullfile(filepath,filename);
end

Hz = 500;   % Sampling Frequency
TDyn = 2074.5;   % Acquisition time in seconds (for working around dummy scans)
DOSAVE = 1; % Saving Flag
NDUM = 3;

if nargin > 1 && ~isempty(varargin)
    for i = 1:2:length(varargin)
        optlabel = varargin{i};
        optvalue = varargin{i+1};
        
        if isnumeric(optvalue); optvalue = num2str(optvalue); end;
        eval(sprintf('%s = %s;',optlabel, optvalue));
        
    end;
end;


fclose('all');

fid=fopen(fname);

datastart = 0;
while (datastart == 0) 
	headline = fgetl(fid);
    if strcmp(headline(:),'#')   % single hash to start scan - or end of dummies?
        datastart = 1;
    elseif strcmp(headline(1),'#')  % hash with text for header - not sure if useful
		continue;	
    elseif feof(fid)
        return
    else           
        datastart = 0;    % prep and dummy data (?)
    end
end

data = textscan(fid,'%*d %*d %*d %d %*d %*d %*d %*d %*d %d'); %Read data

% Clip the time series to begin and end with the scan.
start = floor(NDUM*Hz*TDyn/1000);            % eliminate dummy scans at start
stop = find(data{:,2} == 20);  % System uses identifier 0020 as end of scan tag

signal_t = data{:,1};
signal = double(signal_t(start+1:stop));


% Filter the trigger markers from the data
t_on  = find(data{:,2} ~= 0);  % System uses identifier 0020 as start volume tag
t_trig = find(data{:,2} == 2);  % system uses identifier 0002 as start of TR (for segmented, multiple per dynamic)

% PROBLEM: the above don't seem to be consistently spaced, i.e. give a
% different numbers per TDynfor now hack it in....

time = (1:length(signal))./Hz;
time = double(time') ; 

% saving, then save
if DOSAVE
    [fpath fname fext] = fileparts(fname);
    save(fullfile(fpath,[fname '.mat']),'signal','time');
end;

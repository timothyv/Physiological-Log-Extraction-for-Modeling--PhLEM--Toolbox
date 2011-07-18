function [signal, time] = philips_RESPload2(fname,varargin)

%function [resp time] = RESP_load(fname, optname, optval);
%
% A loader for the Respiration data coming off of the Philips 3T system.
%  
%   INPUT:
%   fname = name of physiological monitoring log file with respiration information
%   opts        = options on extracting the ECG signal
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
%   resp   = vector of the respiration waveform
%   time   = vector of timestamps (in seconds)
%
%  
% Written by T. Verstynen (November 2007) and J Schlerf (August 2008)
%
% Rewritten to optimize loading time taking heavily from code written by V.
% Deshpande and J. Grinstead, Siemens Medical Solutions (March 2009) 
% Adapted to Philips physiology log files by P. Summers (March 2009)
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

data = textscan(fid,'%*d %*d %*d %*d %d %*d %*d %*d %*d %d'); %Read data

%data = zeros(4,500000,'int8');
%data = fscanf(fid,'%*d %*d %*d %*d %d %*d %*d %*d %*d %d'); %Read data
%ind = 1;
%while (feof(fid) == 0)
%    linein = fgetl(fid);
%    data = fscanf(fid,'%*i %*i %*i %d %d %d %*i %*i %*i %d')
%    sscanf(linein,'%*i %*i %*i %d %*d %*d %*i %*i %*i %*d')
%    sscanf(linein,'%*i %*i %*i %*i %d %*d %*i %*i %*i %*d')
%    sscanf(linein,'%*i %*i %*i %*i %*i %d %*i %*i %*i %*d')
%    sscanf(linein,'%*i %*i %*i %*i %*d %*d %*i %*i %*i %d')
%data = textscan(fid,'%*d %*d %*d %f %f %f %*d %*d %*d %f') %Read data
%until end of file.
%ind = ind + 1;
%end
%footer = textscan(fid,'%s');   %Read in remaining data (time stamps and statistics).
%Get time stamps from footer:
%for n=1:size(footer{1},1)
%    if strcmp(footer{1}(n),'LogStartMDHTime:')  %log start time
%        LogStartTime=str2num(footer{1}{n+1});
%    end
%    if strcmp(footer{1}(n),'LogStopMDHTime:')   %log stop time
%        LogStopTime=str2num(footer{1}{n+1});
%    end
%    if strcmp(footer{1}(n),'LogStartMPCUTime:') %scan start time
%        ScanStartTime=str2num(footer{1}{n+1});
%    end
%    if strcmp(footer{1}(n),'LogStopMPCUTime:')  %scan stop time
%        ScanStopTime=str2num(footer{1}{n+1});
%    end
%end
%ntpts = size(data(:)) / 4;
%for n=1:ntpts,1
%    m=4*(n-1)+1;
%   ecg(n) = data(m);
%    m=m+1;
%    resp(n) = data(m);
%    m=m+1;    
%    ppu(n) = data(m);
%    m=m+1;
%    mark(n) = data(m);
%    m=m+1;
%end

% Clip the time series to begin and end with the scan.
start = floor(NDUM*Hz*TDyn/1000);            % eliminate dummy scans at start
stop = find(data{:,2} == 20);  % System uses identifier 0200 as end of scan tag

signal_t = data{:,1};
signal = double(signal_t(start+1:stop));


% Filter the trigger markers from the data
t_on  = find(data{:,2} ~= 0);  % System uses identifier 0020 as start volume tag
t_trig = find(data{:,2} == 2);  % system uses identifier 0002 as start of TR (for segmented, multiple per dynamic)

% PROBLEM: the above don't seem to be consistently spaced, i.e. give a
% different numbers per TDynfor now hack it in....

time = (1:length(signal))./Hz;
time = double(time') ; 




%indx = setdiff(data_t(:),union(t_on,t_off)); %Note: depending on when the scan ends, the last size(t_off)~=size(t_on).
%signal = data{1}(indx);
% Reset vectors;
%signal = signal(start+1:end-stop)';
%time = time(start+1:end-stop);
%time = time(:)-time(1);

% saving, then save
if DOSAVE
    [fpath fname fext] = fileparts(fname);
    save(fullfile(fpath,[fname '.mat']),'signal','time');
end;

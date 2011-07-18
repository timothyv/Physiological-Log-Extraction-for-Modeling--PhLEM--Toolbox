function [events x] = PhLEM_event_id(TS,varargin);

% function events_marker = PhLEM_event_id(TS,opts);
% 
% Takes a 1xN physiological time-series (TS) and returns a 1xN
% array with 1 = event peaks and 0 otherwise.  
% 
% Inputs:
% 
%   TS          = 1xN array of physio signals;
%   opts        = options on extracting the ECG signal
%
%     a) 'FiltType': Flag for smoothing options
%        1) 'butter' -> Use a Butterworth bandpass filter (default)
%        2) 'gaussian' -> Smooth with a simple Gaussian kernel
%
%     b) 'Wn': smoothing parameters:
%        1) 1x2 vector of normalized frequencies if using 'butter' option
%           (default Wn = []: No Smoothing)
%        2) scalar with of Gaussian if using 'gaussian' options (in
%           samples)
%
%     c) 'delta': Minimum search distance between peaks (default 200
%                 samples);
%
%     d) 'peak_rise': minimum value that the peak needs to be preceded 
%         (to the left) by.  Defined as percent of peak filtered value.
%         (default 0.75).
%
%     e) 'transform': transformation on the TS data. Default 'none'
%        1) 'abs' -> take the absolute value of the TS signal
%        2) 'z' -> z-score normalize the TS signal
%
% Written by T. Verstynen, J. Diedrichsen, and J. Schlerf (Sept 2008)
% Liscensed under the GPL public license (version 3.0)


delta = 200;
FiltType = 'butter';
Wn = [];
peak_rise = 0.75;
transform = 'none';
process = 'fourier';

if nargin > 1 & ~isempty(varargin)
    for i = 1:2:length(varargin)
        optlabel = varargin{i};
        optvalue = varargin{i+1};
        
        if isnumeric(optvalue); 
            if length(optvalue) > 1
                eval(sprintf('%s = [%2.5f %2.5f];',...
                    optlabel,optvalue(1),optvalue(2)));
            else            
                eval(sprintf('%s = %2.2f;',optlabel, optvalue));
            end;
            
        elseif ischar(optvalue)
            eval(sprintf('%s = %s;',optlabel, ['''' optvalue '''']));
        end;

        
    end;
end;



switch transform
    case 'abs';
        % Do an absolute value transformation on the
        % data and zero mean the series
        x = abs(TS);
        x = x - mean(x);

    case {'z', 'z-score'}
        % Z-score transform the data
        x = (TS-mean(TS))./std(TS);

    otherwise
        % Just zero-mean the series
        x = [TS-mean(TS)]; %./std(TS);
end


% Determine which filter to use
switch FiltType
    case 'butter'
        if ~isempty(Wn)
            [b a] = butter(1,Wn);
            x = filter(b,a,x);
        end;

    case 'gaussian'
        x = smooth_kernel(x(:),Wn);
        x = x';

    otherwise
        error('Unknown filter type');
end;

% ********************************************
% Do phase analysis if required
% ********************************************
switch process
  case {'fourier','rvhr'}

    % $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % $$$ % This is not especially robust to amplitude fluctuations,
    % $$$ % especially if there is one particularly high peak:
    % $$$ 
    % $$$ % Set the threshold
    % $$$ peak_resp = max(x);
    % $$$ 
    % $$$ % Find the peaks and return
    % $$$ [mxtb mntb] = peakdet(x,peak_rise*peak_resp); 
    % $$$ events = zeros(size(TS));
    % $$$ 
    % $$$ peaks = mxtb(:,1);
    % $$$ dpeaks = diff(peaks);
    % $$$ kppeaks = find(dpeaks > delta);
    % $$$ newpeaks = peaks([1 kppeaks'+1]);
    % $$$ 
    % $$$ events(newpeaks) = 1;
    % $$$
    % $$$ % -- Schlerf (Jul 24, 2008)
    % $$$ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Slightly more robust:

    % Take a first pass at peak detection:
    [mxtb mntb] = peakdet(x,1e-14);

    % set the threshold based on the 20th highest rather than the highest:
    sorted_peaks = sort(mxtb(:,2));
    if length(sorted_peaks > 20) 
      peak_resp=sorted_peaks(end-20); 
    else 
      peak_resp = sorted_peaks(1); 
    end

    % And now do a second pass, more robustly filtered, to get the actual peaks:
    [mxtb mntb] = peakdet(x,peak_rise*peak_resp);

    events = zeros(size(TS));

    peaks = mxtb(:,1);
    dpeaks = diff(peaks);
    kppeaks = find(dpeaks > delta);
    newpeaks = peaks([1 kppeaks'+1]);

    events(newpeaks) = 1;

    % DEBUG CATCH:
    if length(newpeaks) < 5;
      keyboard;
    end;

    % -- Schlerf (Jul 24, 2008);
  otherwise
    events = [];
end;
  
      
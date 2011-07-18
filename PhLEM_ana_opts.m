function a_opts = PhLEM_ana_opts(fields,Hz);

% function a_opts = PhLEM_ana_opts(fields,Hz);
% 
% Default analysis options for the PhLEM system.
%
% Fields required for each data signal
%    target_signal = Signal field to extract
%    process       = processing method (default = 'fourier');
%       'fourier' -> Does the RETROICOR algorithm
%       'rvhr'    -> Does either the CR+HRF or RV+RRV model
%       'raw'     -> Downsamples the raw signal
%    rate          = search rate for looking for signal peaks
%    peak_rise     = amplitude threshold for looking for peaks
%    filter        = finter or smoothing type ('butter' or 'gaussian')
%    phase_expand  = 1xN array of fourier expansions (e.g. [1 2] for the 
%                     first and second fourier expansion).
%    Wn            = filter/smoothing kernel.  If bandpass filter, 
%                    normalized frequencies.  If gaussian smoothing, 
%                    then kernel fwhm.
%    Hz            = Frequency of the sampled data signal
%
% Written by T. Verstynen, J. Diedrichsen, and J. Schlerf (Sept 2008)
% Liscensed under the GPL public license (version 3.0)



if nargin < 2 | isempty(Hz);
    Hz = 400;
end;

if nargin < 1 | isempty(fields);
    fields = {'ECG','RESP','PPG','PPGhigh','PPGlow'};
end;

%fields = fieldnames(PhLEM);

for f = 1:length(fields)
  
  clear tmp_opts
  
  switch fields{f};
   case 'ECG'
    % Filter for the ECG data signal
    tmp_opts.target_signal = fields{f};
    tmp_opts.process = 'fourier';
    tmp_opts.rate = Hz/2; % search rate for peakdet
    tmp_opts.peak_rise = 0.75;
    tmp_opts.filter = 'butter';
    tmp_opts.phase_expand = [1 2];
    tmp_opts.Wn = [1 2]/Hz; % bandpass freq (normalized freqencies)
    tmp_opts.Hz = Hz;
    tmp_opts.transform='none';
    
   case 'RESP'
    % Filter for the Respiration waveform
    tmp_opts.target_signal = fields{f};
    tmp_opts.process = 'fourier';
    tmp_opts.rate = Hz*1.5; % search rate for peakdet
    tmp_opts.peak_rise = 0.75;
    tmp_opts.filter = 'butter';
    tmp_opts.phase_expand = [1];
    tmp_opts.Wn = [1/4 1/2]/Hz; % bandpass freq (normalized freqencies)
    tmp_opts.Hz = Hz;
    tmp_opts.transform='none';
    
    % RATHER THAN BANDPASS, CAN USE A GAUSSIAN KERNEL INSTEAD
    %tmp_opts.peak_rise = .10;   
    %tmp_opts.filter = 'gaussian';
    %tmp_opts.phase_expand = [1 2];
    %tmp_opts.Wn = Hz*0.6; % bandpass freq (normalized freqencies)
    %tmp_opts.Hz = Hz;
    %tmp_opts.transform='none';
    
   case 'PPGhigh'
    % High-pass filter on the pulse-ox signal
    tmp_opts.target_signal = upper(fields{f});
    tmp_opts.process = 'fourier';
    tmp_opts.rate = Hz/2; % search rate for peakdet
    tmp_opts.peak_rise = 0.1;
    tmp_opts.filter = 'butter';
    tmp_opts.phase_expand = [1 2];
    tmp_opts.Wn = [1 6]/Hz; % bandpass freq (normalized freqencies)
    tmp_opts.Hz = Hz;
    tmp_opts.transform='abs';
        
   case 'PPGlow'
    % Low-pass filter on the pulse-ox signal
    tmp_opts.target_signal = upper(fields{f});
    tmp_opts.process = 'fourier';
    tmp_opts.filter = 'gaussian';
    tmp_opts.rate = Hz*1.5; % search rate for peakdet
    tmp_opts.peak_rise = .5;   
    tmp_opts.phase_expand = [1 2];
    tmp_opts.Wn = Hz*0.4; % kernel FWHM
    tmp_opts.Hz = Hz;
    tmp_opts.transform='abs';
    
    % RATHER THAN GAUSSIAN SMOOTH, CAN USE A BANDPASS INSTEAD    
%     tmp_opts.process = 'fourier';
%     tmp_opts.rate = Hz/1.5; % search rate for peakdet
%     tmp_opts.peak_rise = 0.5;
%     tmp_opts.filter = 'butter';
%     tmp_opts.phase_expand = [1 2];
%     tmp_opts.Wn = [1/4 1/2]/Hz; % bandpass freq (normalized freqencies)
%     tmp_opts.Hz = Hz;
%     tmp_opts.transform='abs';

   case 'PPG'
    % Band-pass filter on the pulse-ox with non-fourier decomposition
    tmp_opts.target_signal = upper(fields{f});
    tmp_opts.process = 'raw';
    tmp_opts.filter = 'butter';
    tmp_opts.rate = Hz*1.5; % search rate for peakdet
    tmp_opts.peak_rise = .5;   
    tmp_opts.phase_expand = [];
    tmp_opts.Wn = [0.1 6] / Hz; % kernel FWHM
    tmp_opts.Hz = Hz;
    tmp_opts.transform='zscore';

  end;
    
  eval(sprintf('a_opts.%s=tmp_opts;',upper(fields{f})))
end;

function [C, names] = make_physio_regressors(param,PhLEM,varargin)

% function [C, names] = make_physio_regressors(param,PhLEM,varargin)
% 
% The PhLEM tool for making appropriate regressors for your GLM;
% 
% Inputs:
% 
%     1) param     -> Flag for which field of the PhLEM structure
%                     to look at.  Default is 'RESP'
%     2) PhLEM     -> Either the PhLEM stucture itself or a pointer
%                     to the file itself.
%     3) opts      -> Series of opt commands to adjust defaults.
%     
%     OPTS:  
%       a)  'TR'  -> TR in seconds.  Default is 2.
%       b)  'DOSAVE' -> 1 to save the output (default), 0 otherwise
%       c)  'ana_type' -> Analysis type to perform.  Right now the 
%                         options are:
%              i)  'fourier' : Use fourier expansions (default)
%              ii) 'rvhr'    : Use variation model for HR or RV
%              iii)'rate'    : Use the rate of events per TR;
%              iv) 'amp'     : Use the amplitude of events in each TR;
% Outputs
% 
%     1) C     = A NxD array of resampled values, where N = the number of 
%                TRs in your session.
%     2) names = Cell array of names for each column in C.
%     
%
% Written by T. Verstynen, J. Diedrichsen, and J. Schlerf (Sept 2008)
% Liscensed under the GPL public license (version 3.0)
              
             
if nargin < 2 | isempty(PhLEM)

  % Use the SPM5 loader if it is available otherwise, you are using an
  % earlier version of SPM, so go back to spm_get.  Otherwise default to
  % matlab UIGETFILE
  if ~exist('spm_get')
    PhLEM = spm_select(1,'mat','Select Physio File');
  else
    try
      PhLEM = spm_get(1,'*.mat','Select Physio File');
    catch;
      PhLEM = uigetfile('Select PhysioFile');
    end;
  end;
  
end;

if ischar(PhLEM)
    [fpath fname] = fileparts(PhLEM);
    load(PhLEM);
else
    [fdir] = fileparts(pwd);
    fname = sprintf('PhLEM_%d.m',param);
end;

if nargin < 1 | isempty(param)
  param = 'RESP';
  fprintf('No Parameter Identified: Assuming %s\n',param);
end;

% DEFAULTS
TR = 2;
DOSAVE = 1;
ana_type = 'fourier';
plotflag = 0;

% Check for options to reset DEFAULTS
for v=1:2:length(varargin),
    eval(sprintf('%s = varargin{%d};',varargin{v},v+1));
end

% Get the target field of the PhLEM Structure
dtmp = eval(sprintf('PhLEM.%s',param));

% Setup TR markers
tr_sampling = round(PhLEM.Time(end)/TR);
tr_timestamps = [0:(TR):PhLEM.Time(end)-(TR/2)]/mean(diff(PhLEM.Time));

% Determine what process to do based on the opts value
ana_type = dtmp.opts.process;

switch ana_type
 case 'fourier'        
  C = [];
  for exp = dtmp.opts.phase_expand
    nphs = interp1(eval(sprintf('dtmp.phase%d',...
				dtmp.opts.phase_expand(exp))),...
		   tr_timestamps,'linear','extrap');            
	  
    C = [C [sin(nphs)]'];
    C = [C [cos(nphs)]'];
    
    names{(exp-1)*2+1} = sprintf('FFT Exp %d: Sine %s',...
				 dtmp.opts.phase_expand(exp), param);
    names{(exp-1)*2+2} = sprintf('FFT Exp %d: Cosine %s',...
				 dtmp.opts.phase_expand(exp), param);
  end;
        
 case 'rate'
  C = histc(find(dtmp.events),tr_timestamps);
  names = sprintf('Rate: %s', param);
  
 case 'amp'
  amps = [dtmp.data./mean(dtmp.data)]-1;
  amps = dtmp.data(find(dtmp.events))./mean(dtmp.data) - 1;
  events = histc(find(dtmp.events),tr_timestamps);
  events(find(events)) = amps;
  C = events;
        
 case 'raw'
  norm_data = [dtmp.smoothed_data - ...
	       mean(dtmp.smoothed_data)] ./ std(dtmp.smoothed_data);
  C = interp1(norm_data,tr_timestamps,'linear','extrap')';
  names{1} = sprintf('Raw: %s', param);
  
 case 'rvhr'
  hr_fields = {'ECG','PPGHIGH'};
  rv_fields = {'RESP','PPGLOW'};
  
  % process the datafield appropriately
  if ismember(param,rv_fields)
    C = (local_make_rv(getfield(PhLEM,param),PhLEM.Time,TR));
    names = {sprintf('RVRRF: %s',param)};
  elseif ismember(param,hr_fields)
    C = (local_make_hr(getfield(PhLEM,param),PhLEM.Time,TR));
    names = {sprintf('HRCRF: %s',param)};
  end

end;
        
if plotflag

    time = PhLEM.Time;
    maxlen = 10000;

    switch ana_type
        case 'fourier'    
            n_exps = length(dtmp.opts.phase_expand);
            cpos = []; xtick = []; xlabel = [];
            for e = 1:n_exps
                subpos = [1 2]+([4 4]*(e-1));
                subplot(n_exps,4,subpos);                
                x = eval(sprintf('dtmp.phase%d',dtmp.opts.phase_expand(e)));
                if length(x) > maxlen
                    x = x(1:maxlen);
                end;
                
                plot(time(1:length(x)),sin(x),'r',time(1:length(x)),cos(x),'b');
                title(sprintf('Fourier Expansion: %d',dtmp.opts.phase_expand(e)));
                ylabel('Time');
                
                cpos = [cpos [3 4]+(e-1)*[4 4]];
                xlabel = [xlabel {sprintf('sin%d',e), sprintf('cos%d',e)}];
                
            end;

            subplot(n_exps,4,cpos);
            imagesc(C);
            title('Regressor Values');
            set(gca,'XTick',[1:size(C,2)],'XTickLabel',xlabel);
            
            scr = get(0,'ScreenSize');
            pos = [0.2*scr(3) 0.15*scr(4) 0.5*scr(3) 0.7*scr(4)];
            set(gcf,'Position',pos);
            
    end;
end;
            
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute RV regressor
function Col = local_make_rv(D,Time,TR)

num_tr = round(Time(end)/TR);
tr_vec = floor(Time/TR)+1;

% compute RV vector according to Chang et al:
for i=1:num_tr
  rv(1,i) = std(D.smoothed_data(ismember(tr_vec,[i-1,i,i+1])));
end

% normalize using z-score:
rv = zscore(rv);

% ...convolve with RRF (from Reference 2):
t = 0:TR:28;
RRF = 0.6*t.^2.1.*exp(-t/1.6)-0.0023*t.^3.54.*exp(-t/4.25);
Col = conv(rv,RRF);

% ...trim excess:
Col = Col(1:num_tr)';

% ...and we're done.
return
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute HR regressor
function Col = local_make_hr(D,Time,TR)

num_tr = round(Time(end)/TR);
beats = Time(D.events==1);
tr_vec = floor(beats/TR)+1;

% compute HR vector according to Chang et al:
for i=1:num_tr
  hr(1,i) = 60/mean(diff(beats(ismember(tr_vec,[i-1,i,i+1]))));
end

% normalize using z-score:
hr = zscore(hr); %hr-mean(hr);

% ...convolve with CRF (from Reference 1):
t = 0:TR:28;
CRF = 0.6*t.^2.7.*exp(-t/1.6) - 16/sqrt(2*pi*9)*exp(-0.5*(t-12).^2/9);

% And in case I need the CRF's time derivative:
CRF_TD = (1.62*t.^1/7 - (3/8)*t.^2.7).*exp(-t/1.6) + ...
	 (t./9 + 8/3).*exp(-(t-12).^2/18);

Col = conv(hr,CRF);

Col2 = conv(hr,CRF_TD);
Col2 = Col2(1:num_tr)';

% ...trim excess:
Col = Col(1:num_tr)';

% ...and we're done.
return

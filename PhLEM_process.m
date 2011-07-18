function PhLEM = PhLEM_process(PhLEM,varargin);

% 
% function PhLEM = PhLEM_process(PhLEM);
% 
% Uses the parameters layed out in the .opts of each relevant
% field in the PhLEM structure to smooth the signal and identify
% peaks in the data.  
% 
%
% Inputs:
%     1) PhLEM  -> An existing PhLEM structure after it has been 
%                  processed with PhLEM_constructor.
% 
% Outputs:
%     1) Modified PhLEM structure
%     
% Written by T. Verstynen (January 2008).
% Liscensed under the GPL public license (version 3.0)


% Get new analysis options (used later)
%transform = 'none'

for v=1:2:length(varargin),
    eval(sprintf('%s = varargin{%d};',varargin{v},v+1));
end


% Go through and process each field of PhLEM
fields = fieldnames(PhLEM);
fields = setdiff(fields,'Time'); % Remove 'Time' field from analysis


% Perform analysis defined in a_opts
for f = 1:length(fields)
    dtmp = eval(['PhLEM.' fields{f} ';']);   
    
    [dtmp.events dtmp.smoothed_data] = PhLEM_event_id(dtmp.data,'FiltType',...
        dtmp.opts.filter,...
        'Wn',dtmp.opts.Wn,...
        'delta',dtmp.opts.rate,...
        'peak_rise',dtmp.opts.peak_rise,...
        'transform',dtmp.opts.transform,...
        'process',dtmp.opts.process);

    eval(sprintf('PhLEM.%s = dtmp;',fields{f}));
    
end;
    

    


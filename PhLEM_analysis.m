function PhLEM = PhLEM_analysis(PhLEM,varargin);

% 
% function PhLEM = PhLEM_analysis(PhLEM);
% 
% Takes the event points from the PhLEM structure and estimates
% the unwarped phase between them.  
%
% Inputs:
%     1) PhLEM  -> An existing PhLEM structure after it has been 
%                  processed with PhLEM_process.
% 
% Outputs:
%     1) Modified PhLEM structure
%     
% Written by T. Verstynen (January 2008).


% Go through and process each field of PhLEM
fields = fieldnames(PhLEM);
fields = setdiff(fields,'Time'); % Remove 'Time' field from analysis

% Perform analysis defined in a_opts
for f = 1:length(fields)
    dtmp = eval(['PhLEM.' fields{f}]);
    
    switch dtmp.opts.process
        case 'fourier'
            % Do N many phase expansions
            if ~isfield(dtmp.opts,'phase_expand') | isempty(dtmp.opts.phase_expand)
                dtmp.opts.phase_expand = 1;
            end;
            
            for p = dtmp.opts.phase_expand
                tmp_phase = PhLEM_phase_expand(dtmp.events,PhLEM.Time);                
                eval(sprintf('dtmp.phase%d = %d*tmp_phase;',p,dtmp.opts.phase_expand(p)));
            end;
                        
    end;
        
    eval(sprintf('PhLEM.%s = dtmp;',fields{f}));
    
end;
    

    


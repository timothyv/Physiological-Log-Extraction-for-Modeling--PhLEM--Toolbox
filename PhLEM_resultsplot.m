function PhLEM = PhLEM_resultsplot(PhLEM,field,varargin);
 
% function PhLEM = PhLEM_resultsplot(PhLEM,field,varargin);
% 
% General plotting function for a field in the PhLEM structure
% after it has been processed with PhLEM_process.
% 
% Written by T. Verstynen, J. Diedrichsen, and J. Schlerf (Sept 2008)
% Liscensed under the GPL public license (version 3.0)



if nargin < 2 | isempty(field)
    field = 'ECG';
end;

dtmp = eval(['PhLEM.' upper(field)]);
time = PhLEM.Time;

subplot(2,1,1);
plot(time,dtmp.data,'k');
xlabel('Time (sec)','FontSize',12);
ylabel(sprintf('Field: %s',field),'FontSize',12);
title('Unfiltered Signal','FontSize',14);

subplot(2,1,2);
plot(time,dtmp.smoothed_data,'k');

switch dtmp.opts.process
  case 'raw'
    % If not using events for phase analysis, then skip this step.

  otherwise
    drawlines(time(find(dtmp.events)),'r');
end;

xlabel('Time (sec)','FontSize',12);
ylabel(sprintf('Smoothed Data: %s',field),'FontSize',12);


signalRate = length(find(dtmp.events)) / [time(end)] * 60;
title(sprintf('Event Rate = %2.2f per minute',signalRate),'FontSize',14);

Scrn = get(0,'ScreenSize');

Pos = round([0.10*Scrn(3) 0.35*Scrn(4) 0.80*Scrn(3) 0.50*(Scrn(4))]);
set(gcf,'Position',Pos);



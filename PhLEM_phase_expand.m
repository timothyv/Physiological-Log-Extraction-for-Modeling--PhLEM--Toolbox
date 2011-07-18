function phase = PhLEM_phase_expand(event_series,time);

% function phase = PhLEM_phase_expand(event_series,time);
% 
% Takes a series of events defined as 1's or 0's, and an index vector
% for real time (in seconds), and returns the unwarped phase between 
% the events. Works by assuming that the distance between consecutive
% events is 2pi.
% 
% Inputs:
%   event_series  = 1xN array of 1's (events) or 0's (non-events);
%   time          = 1xN array of timestamps in real time (seconds);
%   
% Written by T. Verstynen (November 2007)
%   
% Liscensed under the GPL public license (version 3.0)

% Time Stamps of Events
events = find(event_series);

% Phase Interpolation
uwp_phase = [0:length(events)-1];
uwp_phase = 2*pi.*uwp_phase;    
phase = interp1(time(events),uwp_phase,time,'spline','extrap');
    

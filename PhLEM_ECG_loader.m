function [signal time] = PhLEM_ECG_loader(fname,varargin);

% This is just a shadow function for you're own loader functions.  Most of the
% PhLEM scripts call this to load ECG data, but edit the line below as follows
% to use you're own loader routine.  Be sure to follow the FNAME and VARARGIN
% values in the example.
% 
% [sigma time] = YOUR_LOADER_FUNCTION(fname,varargin);
%
% Written by T. Verstynen (January 2008).

[sigma time] = nic_ECGload(fname,varargin);
function outTS = PhLEM_ts_match(TS,norm_len);

% function outTS = PhLEM_ts_match(TS,norm_len);
%
% Function for expanding a time series to a different length using
%   the INTERP1 function.  This is generally used for bringing
%   the respiration time series into the time of the ECG signal.
%  
%   INPUT:
%   TS       = vector time of time series to be adjusted
%   norm_len = the new length to normalize the vector TS to
%  
%   OUTPUT:
%   outTS    = the new 1 x NORM_LEN vector
%  
% Written by T. Verstynen (August 2007).
% Liscensed under the GPL public license (version 3.0)


time = 1:norm_len;
norm_time = 1:(norm_len/length(TS)):norm_len;
outTS = interp1(norm_time,TS,time,'linear','extrap');
	
function [signal, time] = exvolt_data_loader(input_file,varargin);

% function [signal,time] = exvolt_data_loader(input_file,optlabel,optvalue);
%
% A loader function for output from the Exvolt Script.
%
% INPUT: 
% input_file = pointer to output file from exvolt log
%   opts        = options on extracting the ECG signal
%
%     a) 'Hz': Sampling Rate of the signal (50Hz default)
%
% OUTPUT
% data       = vector of data signal
%
% Written by T. Verstynen (November 2010);
%  Updated: April, 2011
%
% Liscensed under the GPL public license (version 3.0)

Hz = 50;   % Sampling Frequency

if nargin > 1 & ~isempty(varargin)
	for i = 1:2:length(varargin)
	optlabel = varargin{i};
	optvalue = varargin{i+1};

	if isnumeric(optvalue); optvalue = num2str(optvalue); end;
		eval(sprintf('%s = %s;',optlabel, optvalue));

	end;
end;

% Just use DLMREAD
signal = dlmread(input_file,'\t',9,0);

% Estimate real time
time = (1:length(signal))./Hz;

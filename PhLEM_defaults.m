function defaults = PhLEM_defaults(PhLEM);

% Default global variables for the Phlem process.  Creates function
% pointers for each loader type.
%
% Written by T. Verstynen, J. Diedrichsen, and J. Schlerf (Sept 2008)
% Liscensed under the GPL public license (version 3.0)

global defaults

defaults.ecg_loader_fcn = 'siemens_ECGload';  % handle to the ECG loader
defaults.resp_loader_fcn = 'siemens_RESPload';% handle to the RESP loader
defaults.ppg_loader_fcn = 'siemens_PPGload';% handle to the PPG loader
defaults.ana_opts_fcn = 'PhLEM_ana_opts';
defaults.save_data = 0;  % Save data? 1 = yes, 0 = no
defaults.tr = 2000; % TR length of study
defaults.ipat = 1; % Was iPat on? (needed for loader functions)
defaults.ecg_hz  = 400; % Sampling rate of ECG signal
defaults.ppg_hz  = 50;  % Sampling rate of pulse-ox signal
defaults.resp_hz = 50;  % Sampling rate of respiration signal
defaults.results_plot = 0; % Plot results during setup?
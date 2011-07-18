function defaults = PhLEM_defaults(PhLEM);

% Default global variables for the Phlem process.  Creates function
% pointers for each loader type.
%
% Written by T. Verstynen, J. Diedrichsen, and J. Schlerf (Sept 2008)
% Liscensed under the GPL public license (version 3.0)

global defaults

defaults.ecg_loader_fcn = {'siemens_ECGload'};  % handle to the ECG loader
defaults.resp_loader_fcn = {'philips_RESPload'};% handle to the ECG loader
defaults.ppo_loader_fcn = {'philips_PPOload'};% handle to the ECG loader
defaults.save_data = 1;  % Save data? 1 = yes, 0 = no
defaults.tr = 2074; % TR length of study
defaults.ipat = 0; % Was iPat on? (needed for loader functions)
defaults.ecg_hz  = 500; % Sampling rate of ECG signal
defaults.ppo_hz  = 500;  % Sampling rate of pulse-ox signal
defaults.resp_hz = 500;  % Sampling rate of respiration signal
defaults.results_plot = 1; % Plot results during setup?
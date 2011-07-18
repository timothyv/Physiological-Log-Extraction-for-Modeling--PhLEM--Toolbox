function PhLEM(varargin)

% PhLEM (Physiological Log Extraction for Modeling) Toolkit
%  
% This file is an example loader and data prep function for the PhLEM 
% toolbox.  PhLEM creates a data structure by composed of different
% physiological measures taken while scanning.  This file shows how to
% spit out a PhLEM structure which is used by the main function
% "make_physio_regresssors" used in SPM.  
%
% The PhLEM structure can have one of several physio signal fields:
%
% ECG     : Electrocardiogram recordings
% RESP    : Pneumatic belt recordings
% PPG     : Unfiltered photoplethysmograph recordings
% PPGHIGH : High pass filtered photoplethysmograph data
% PPGLOW  : Low pass filtered photoplethysmograph data 
%
% The top level also contains a Nx1 array of time stamps called Time.
% 
% Within each signal field, the processed PhLEM structure should contain
% multiple subfields
%
%   data   : Nx1 array of raw data
%   opts   : analysis optins fields
%   events : Nx1 binary events log
%   smoothed_data : filtered and smoothed data
%   phase* : Nx1 array of "phase time" for the fourier series.  Can be
%           expanded to multiple phases dependin on info in the "opts" field
%
%
% The main PhLEM functions are
%
%   *_ana_opts  : contains all the analysis options and flags for each
%                       signal type
%
%   PhLEM_setup : creates the PhLEM object with the relevant
%                       fields. Type "help PhLEM_setup" for more
%                       information.
%   
%   make_physio_regressors : the main PhLEM function which does the expansion 
%                            on the Fourier series and puts the phases into 
% 			     the same space as your imaging data. Type
% 			     "help make_physio_regressors" for more info.
%
% See the demo folder for more information on how to setup the PhLEM
% structure.
%
% Note that if you are NOT using a Siemen's logging system, you will need
% to write your own loade routines to get the initial ECG and Respiration
% data.  Check this website (http://keck.ucsf.edu/~timothyv/) for loader
% functions to use for other systems.
%
% Written by: T. Verstynen & J. Diedrichsen (August 2007)
% Modified by: J. Schlerf (July 2008)
% Modified by: T. Verstynen (January 2011)

help PhLEM

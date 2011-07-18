This demo is designed to walk you through generating different 
physiological models using the data included in the demo.zip file

**************************************************************************
I. Introduction

We'll start with the basic concepts.  There are two files you should need to
modify to make your output work.  The first is the analysis options file 
(the *_ana_opts.m file).  This contains all the relevant analysis and decomposition
information you'll need to make the noise regressors for each signal.  The second
is the PhLEM_defaults.m file.  This contains more general information like
data loaders, scan parameters, etc.  Become comfortable with both of these
files before you start doing something serious.

A. The ana_opts file.  You can modify the parameters in the *_ana_opts.m
file to determine which model type to generate.  By default only the
RETROICOR model is produced.  To change this to another model type, modify
the 'process' flag to either 'rvhr' for the RV+RRV or HR+CRF model types 
or 'raw' for a simple downsampled model.

B. If you want to modify any other general paremters such as the data
loader or the sampling rates of the data files, change these options in 
the PhLEM_defaults.m file.

**************************************************************************
II.  Demo 

Now we'll try a demo of some model types.

1. First try processing the ECG, with the PPG signal being seperated 
based on high and low frequency patterns.  This will use the default parameters 
in the normal PhLEM_ana_opts.m file
 
>> PhLEM = PhLEM_setup('test.mat','ecgfile','PhysioLog_20080331T125957.ecg');

2. Next let's include the PPG data, which is filtered in several ways by default

>> PhLEM = PhLEM_setup('test.mat','ecgfile','PhysioLog_20080331T125957.ecg','ppgfile','PhysioLog_20080331T125957.puls')

3. Okay, now let's play with inserting a new analysis options file.  For this
   example we'll use the optiosn file that does some of the analysis used
   in Verstynen and Deshpande 2011.  This applies the RETROICOR to the 
   ECG and RESP data, the HR+CRF analysis on the PPGHIGH data, the RV+RRF
   analysis on the PPGLOW data, and a simple downsampling on the regular PPG 
   data.

>> PhLEM = PhLEM_setup('test.mat','ecgfile','PhysioLog_20080331T125957.ecg','ppgfile','PhysioLog_20080331T125957.puls','respfile','PhysioLog_20080331T125957.resp','ana_opts_fcn','verstynen_despande_ana_opts')

4.  Now let's generate a few regressors assuming a 2-second TR.  You can
    change the TR and other parameters in the PhLEM_defaults.m file

>> [ppghigh, ppghighlabel] = make_physio_regressors('PPGHIGH',PhLEM)

>> [ecg_retroicor, ecglabel] = make_physio_regressors('ECG',PhLEM)

>> [ppgraw, ppglabel] = make_physio_regressors('PPG',PhLEM)


5.  Now let's say you want to change how you model a signal source.  You can do
    that on the fly by changing the 'process' flag in that signal's 'opts' field.  For
    example, let's now model the ECG using the HR+CRF model. 

>> PhLEM.ECG.opts.process = 'rvhr';
>> [ecg_hrcrf, ecglabel] = make_physio_regressors('ECG',PhLEM);

6.  Finally if you want to save these out as text files to load into a model 
    simply use dlmwrite

>> dlmwrite('ECG_RETROICOR.txt',ecg_retroicor);

7. If you're not working with Siemens data or you're using Exvolt to sync
    your data to your scan, you can change the loader options in the PhLEM_defaults.m
    file.  For example, if you are using Exvolt to preprocess your data, 
    then change the following fields in the defaults file to:

defaults.ecg_loader_fcn = 'exvolt_data_loader';  % handle to the ECG loader
defaults.resp_loader_fcn = 'exvolt_data_loader';% handle to the RESP loader
defaults.ppg_loader_fcn = 'exvolt_data_loader';% handle to the PPG loader

    If you are not using a 2sec TR, then change the TR filed to whatever
    your TR is (in milliseconds).


This ends the simple demo.  Play around to get a sense of the full functional
utility of this script.
%% odas_version_info 
% Return the version of the current ODAS library.
%%
% <latex>\index{Functions!odas\_version\_info}</latex>
%
%%% Syntax
%   version = odas_version_info( )
%
% * []
% * [version] Numeric value of the version.
%
%%% Description
% Function called by other functions to determin the current version of
% the ODAS Libarary.  Used when checking for outdated data files.

% Version History
%
% * 2015-06-03 (WID) initial version 4.0
% * 2016-08-02 (WID) updated to version 4.1.  Added support for some JAC
%                    sensors.  Minor bug fixes.  Notable inclusion of
%                    show_ch fuction.
% * 2017-01-27 (WID) updated to version 4.2.  
%                    addition of AEM1-G analog and digital types
% * 2018-01-05 (JMM) updated to version 4.3.
%                    Bug fixes, including despiking of piezo-accelerometers
%                    This version will be the last major update to the ODAS
%                    library. Minor improvements and bug fixes will be
%                    made, but major changes will be released as part of
%                    the Zissou software package.
% * 2018-04-14 (JMM) updated to version 4.3.02
%                    Fixed some small plotting bugs and annoyances (they
%                    arose during training at Lake Tahoe).
% * 2018-04-25 (JMM) updated to version 4.3.03
%                    Fixed some small plotting bugs and annoyances (they
%                    arose during training at UNAL). Most significant
%                    change is that a plot of pressure will be generated
%                    even if no profiles are detected. 
% * 2018-05-18 (JMM) updated to version 4.3.04
%                    Fixed call to salinity_JAC to include sampling rate
%                    and mean speed. Fixed comments so manual could be
%                    updated. 
% * 2018-08-28 (JMM) updated to version 4.3.05
%                    Added a warning message in read_odas if there is no
%                    channel section in the setup file for an id number 
%                    that is in the address matrix.
% * 2018-09-27 (JMM) updated to version 4.3.06 
%                    Small fix for accelerometer plot for seagliders
% * 2018-09-27 (JMM) updated to version 4.3.07 
%                    Removed two lines that overwrote default speed and 
%                    cutoff frequency in salinity_JAC.m
% * 2019-04-19 (JMM) updated to version 4.3.08 
%                    Added support for type 'raw'. Fixed conflicting
%                    variable name in convert_odas.
% * 2019-06-03 (JMM) Updated to version 4.4
%                    - Added capability of handling horizontal profilers.
%                    - DOF equation for spectra updated to Nuttall formulation
%                    - Figure of Merit estimates now outputted from
%                       get_diss_odas.
%                    - Significant updates to cal_FP07_in_situ (better
%                    control of plots, warning message for small
%                    temperature ranges, consistency of inputs with
%                    quick_look)
%                    - Added functions to estimate the electronic noise for
%                      thermistor and shear channels
%                      from shear and thermistor channels.
%                    - show_ch: Warning messages removed and "convert" option
%                           modified to also prevent a deconvolve from being 
%                           performed.
%                    - convert_odas now has support for type "raw"
%                    - Other minor annotation corrections and documentation
%                    updates.                     
% * 2019-08-19 (JMM) updated to version 4.4.01
%                    Bug fix that was causing inclinometer channels to
%                    reverse. Will now only be implemented if instrument is
%                    a vmp or rvmp.
% * 2019-08-28 (JMM) updated to version 4.4.02
%                    quick_bench mods: Change figure sizes in to prevent 
%                    clipping of the legends when exporting. Also added
%                    support for rsijac_t, rsijac_c and DO channels 
% * 2019-09-24 (JMM) updated to version 4.4.03
%                    Bug fix to remove raw pressure signal from kinematics
%                    plot.
% * 2019-10-08 (JMM) updated to version 4.4.04
%                    Bug fixes:
%						- Code crashed when inclinometers and linear 
%						  accelerometers were sampled at different rates 
%						  (issue only in some custom instruments)
%                       - Fix in read_odas to read jac_t as unsigned. Was
%                       causing erroneous values in quick_bench.
%					 Also changed default plot_rawaccel flag to true in 
%				     quick_look, so that non-despiked signal is shown in 
%					 default version of the figure.
% * 2020-02-03 (JMM) updated to version 4.4.05
%                    Minor changes:
%						- Updated quick_bench.m to plot Rinko DO sensor
%						measurements in physical units
%                       - Added units to DO sensor in convert_odas.
% * 2020-04-16 (JMM) updated to version 4.4.06
%                    Minor changes:
%						- Reduced num_fft in get_diss_odas by 1. Previous
%						variable was incorrect. This will lead to
%						smaller Figures of Merit [prop. to sqrt(dof)]
%                       - Reduced last dimension if only one value for 
%                       quality control variables in get_diss_odas. Ensures
%                       consistency in size with other output variables. 
% * 2020-06-23 (JMM) updated to version 4.4.07
%                    Minor changes:
%                       - quick_bench: bug fix for seaglider, added plot of
%                       SBT/SBC channel, change axes limits automatically
%                       for ASTP spectra plot
%                       - salinity_JAC: output defaults if no input
%                       provided
%                       - show_spec: better position of figure on screen
%                       - show_ch: Minor improvements
%                       - setupstr: Allow for recDuration in addition to
%                       recSize variable.
% * 2020-10-13 (JMM) updated to version 4.4.08
%                    Minor changes:
%                     - Manual: now include all currently supported vehicle 
%                           types in Section 3.6
%                     - quick_look: Updated linewidths for vibration channel 
%                           plot so they are consistent between panels. Also
%                           when 'op_area' is set to 'tidal_ch', the 10x 
%                           string has been added to the legend to indicate 
%                           that the signal has been multiplied by a factor 
%                           of 10
%                      - get_diss_odas: Fixed size of data fast and data slow 
%                           when there is only one estimate of the 
%                           dissipation rate. This ensures output variables 
%                           are all the same size.
%                      - odas_p2mat: Convert single variables to double in 
%                           extract_vector part of odas_p2mat so that interp1 %                           
%                           function works. 
% * 2022-06-02 (JPJ) Update to version 4.5
% 	            Major Change:
% 	                - Adding correction factor to the Goodman Routine to clean Shear Spectra. 
% 	                    The routine has an bias due to the finite length of the fft segments used, this is now compensated. 
% 	                    See TN 061 for additional details. 
% 	            Minor Changes:
% 	                - The file start is now taken from the configuration string record. On RDL instrument this is a significant
% 	                    improvement in accuracy. 
%  	                - The get_profile function has had a bug fix. It now operates more robustly.  
% * 2022-07-19 (JPJ) Update to version 4.5.01
% 	            Minor Changes:
% 	                - Quick_look updated so that is works when users pass a data file that has no shear data
%       	        - Despiking is fixed based on bug found by Sam Kelly to find the adequate good data
%                   - Vehicle field in setup.cfg can be empty or spaces and the software will not crash
% 					- Making quick_look work with both EMC_Cur and EM_Cur

function version = odas_version_info( )

version = 4.5;

end


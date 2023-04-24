% ODAS Matlab Library version 4.5.1
% Created by Rockland Scientific Inc.
% Last Updated: Jul 15, 2022
%
% adis.m                       - Convert inclinometer
% atg.m                        - Adiabatic temperature gradient (Hydrographics Lib)
% cal_FP07_in_situ.m           - Calibrate thermistor
% channel_sampling_rate.m      - Calculate sampling rate
% check_bad_buffers.m          - Find bad buffers
% clean_shear_spec.m           - Compute spectra and apply Goodman routine
% conductivity_TPcorr.m        - Correct for thermal expansion of SBE4C
% convert_odas.m               - Convert to physical units
% correct_amt.m                - Adjust O2  for temperature, salinity, depth
% csd_matrix_odas.m            - Estimate the spectra of matrices
% csd_odas.m                   - Estimate the spectra of vectors
% deconvolve.m                 - Deconvolve signal plus derivative
% despike.m                    - Remove short-duration spikes
% extract_odas.m               - Extract a range of records from RSI binary file
% extract_setupstr.m           - Extract config string from RSI binary file
% fig_aspectratio.m            - Set the aspect ratio of a figure
% file_with_ext.m              - Find files with default extensions.
% fix_bad_buffers.m            - Fix bad buffers in RSI binary file
% fix_underscore.m             - Legacy function replaced with texstr.m
% fopen_odas.m                 - Open a RSI binary file
% get_diss_odas.m              - Calculate the dissipation rate
% get_latest_file.m            - Get name of latest RSI binary file
% get_profile.m                - Extract indices for a "profile"
% get_scalar_spectra_odas.m    - Calculate the spectra of scalar signals
% gradT_noise_odas.m           - Spectrum of noise of a temperature gradient
% hotelfile_Nortek_vector.m    - Generate hotel-file from Vector data
% hotelfile_Remus_mat.m        - Generate hotel-file from Remus data
% hotelfile_seaglider_netcdf.m - Generate hotel-file from Seaglider data
% hotelfile_slocum_netcdf.m    - Generate hotel-file from Slocum data
% make_gradC_odas.m            - Calculate gradient of conductivity
% make_gradT_odas.m            - Calculate gradient of temperature
% make_scientific.m            - Convert number to string (scientific notation)
% median_filter.m              - Remove fliers from ADV data
% nasmyth.m                    - Generate Nasmyth universal spectrum
% noise_shearchannel.m         - Noise spectrum of the shear channels
% noise_thermchannel.m         - Noise spectrum of the thermistor channels
% odas_lag.m                   - Add a time lag to a signal
% odas_p2mat.m                 - Converrt RSI binary file to .mat file
% odas_version_info.m          - Version info of ODAS Matlab library
% patch_odas.m                 - Fix bad buffers in RSI binary file
% patch_salinity.m             - Clean conductivity data from JAC-CT
% patch_setupstr.m             - Patch a config string into RSI binary data
% plot_freq_spec.m             - Plot frequency spectra
% plot_spec.m                  - Legacy function replaced with show_spec.m
% plot_VMP.m                   - Real-time plots from vertical profiler
% plot_wave_spec.m             - Plot wavenumber spectra
% query_odas.m                 - Report channels contained in RSI binary file
% quick_bench.m                - Generate plots from a bench test
% quick_look.m                 - Visualize contents on an RSI data file
% read_odas.m                  - Convert RSI binary file into mat file
% read_tomi.m                  - Read a range of records from a data file
% salinity.m                   - Salininty calculation (Hydrographics Lib)
% salinity_JAC.m               - Calculate salinity from JAC-CT data
% save_odas.m                  - Efficiently write vector into mat file
% segment_datafile.m           - Break data into segments
% setupstr.m                   - Read attributes of config string
% show_P.m                     - Extract pressure from RSI binary file
% show_ch.m                    - Extract channel data from RSI binary file
% show_spec.m                  - Plot frequency and wavenumber spectra
% sigma_p.m                    - Potential density (Hydrographics Lib)
% sigmat.m                     - Density at p = 0 (Hydrographics Lib)
% texstr.m                     - Convert a string into TeX form
% theta.m                      - Potential temperature (Hydrographics Lib)
% TI_correction.m              - Correct measured temperature for themal inertia
% visc00.m                     - Kinematic viscosity at S = 0
% visc35.m                     - Kinematic viscosity at S = 35
% viscosity.m                  - Kinemetic and dynamic viscosity of seawater

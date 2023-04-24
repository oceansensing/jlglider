%% get_diss_odas
% Calculate a profile of dissipation rate.
%%
% <latex>\index{Functions!get\_diss\_odas}</latex>
%
%%% Syntax
%   diss = get_diss_odas( shear, A, info )
%
% * [shear] Matrix of shear probe signals, in physical units, one shear
%        probe per column. 
% * [A] Matrix of acceleration signals, one acceleromter per column. Does
%        not have to be in physical units. Number of rows must be the same
%        in shear and A.
% * [info] Structure containing additional process-control parameters.
% * []
% * [diss] Structure of results including the profile of the rate of
%        dissipation of turbulent kinetic energy, for each probe, and
%        related information. 
%
%%% Description
%
% This function calculates a profile of the rate of dissipation of turbulent
% kinetic energy for each shear probe using the method described in RSI
% Technical Note 028. The shear data should be despiked (using the function
% despike), and  high-pass filtered at a frequency [Hz] that is about
% one-half of the inverse of the length [s] of the fft-segments. There is
% no limitation on the number of shear-probe and accelerometer signals. All
% results are returned in a structure.
%
%
%%% Input Structure
%
% * [fft_length] length of the FFT segment [samples]. ex, 1024
% * [diss_length] length of shear data over which to estimate the rate of
%       dissipation [samples]. Must be >= 2 * fft_length, 
%       recommended >= 3 * fft_length. 
% * [overlap] Distance [samples] of the overlap of successive dissipation
%       estimates. Should be 0 <= overlap <= diss_length. Suggested value: 
%       overlap = diss_length/2. overlap = 0 means no overlap.
% * [fs_fast] Sampling rate in Hz for shear and acceleration.
% * [fs_slow] Sampling rate in Hz for all slow channels in Data_slow.
% * [t_fast] Time vector for fast channels.
% * [speed] Profiling speed used to derive wavenumber spectra. Must have 
%       length equal to 1 or the number of rows in shear.
% * [T] Temperature vector used when calculating kinematic viscosity.
%       Length must match shear.
% * [P] Pressure data. Must have same length as the shear-probe data.
% * [goodman] Logical variable to determine if the goodman coherent noise
%       reduction algorithm should be used. Default is `true'.
%
%%% Optional Input Parameters
%
% * [Data_fast] Data matrix of fast channels, one per column. These will be
%       averaged over the time interval of each dissipation estimate. Must
%       have the same number of rows as the shear matrix. 
%       Default value = [ ]. 
% * [Data_slow] Same as Data_fast but assumes slow channel data. Each
%       column should have length = size(shear,1)*fs_slow/fs_fast. Default
%       value = [ ].
% * [fit_order] Order of the polynomial fit to the shear spectra, in log-log
%       space, that is used to estimate the spectral-minimum wavenumber.
%       Default = 3. 
% * [f_AA] Cut-off frequency of the anti-aliasing filter. Default = 98 Hz.
% * [fit_2_isr] Value of dissipation rate at which the estimation
%       will be derived from a fit to the inertial subrange rather than by
%       spectral integration. Default fit_2_isr = 1.5e-5 W/kg.
% * [f_limit] Maximum frequency to use when estimating the
%       rate of dissipation. Use this parameter if you have vibrational
%       contamination that is not removed by the Goddman coherent noise
%       removal function (clean_shear_spectrum). Default value is
%       f_limit=inf. 
% * [UVW] Matrix where every column is a velocity component. Must have the
%       same number of rows as the shear matrix. It will be used to
%       estimate the maximum angle of attack of the flow for each
%       dissipation estimate. Default = [ ].
%
%%% Output Structure
%
% A structure of the rate of dissipation estimated from the shear probes, along
%   with a lot of other information. The elements are:
%
%%% Fields associated with the rate of dissipation
%
% * [e] Rate of dissipation [W/kg]. One row for every shear
%      probe and one column for every dissipation estimate.
% * [K_max] A matrix of the maximum wavenumber used for each dissipation
%      estimate. Its size is [N P] where N is the number of shear probes
%      and P is the number of dissipation estimates.
% * [method] a matrix indicating the method used to make each dissipation 
%       estimate. A value of 0 indicates the variance method and a value of 
%       1 indicates the inertial subrange fitting method. Its size is [N P] 
%       where N is the number of shear probes and P is the number of 
%       dissipation estimates.
% * [dof_spec] degrees of freedom of the spectral estimates based on 
%      Nuttall (1971), i.e. dof_spec = 1.9*num_of_ffts. This is a constant
%      value based on the diss_length and fft_length.
% * [dof_e] a matrix of the degrees of freedom in each dissipation
%      estimate. Its size [N P] where N is the number of shear probes and P
%      is the number of estimates. It is the number of wavenumber points used to
%      estimate epsilon, divided by the Nyquist bandwidth, times the
%      diss_length (in samples). IMPORTANT: This is only an approximation
%      and work to determine a more appropriate estimate is ongoing.
% * [mad] A matrix of the mean absolute deviation of the base-10 logarithm
%       of the spectral values from the Nasmyth spectrum over the
%       wavenumber range used to make each dissipation estimate. 
% * [FM] A 'figure of merit' given by FM = mad*sqrt(dof_spec). This metric 
%       can be used for quality control of the dissipation estimates. Its 
%       size [N P] where N is the number of shear probes and P is the number 
%       of estimates.
%
%%% Fields related to the spectral estimates
%
% * [Nasmyth_spec] a 3-dimensional matrix of the Nasmyth spectra for each
%       shear probe and each dissipation estimate. Its size is [M N P]
%       where M is equal to the number of rows in diss.K and each
%       wavenumber element corresponds to diss.K, N is the number of shear
%       probes and P is the number of dissipation estimates.
% * [sh_clean] Wavenumber spectrum for each shear probe signal at each
%      dissipation estimate. This is a 4-dimensional matrix of size
%      [M N N P] where N is the number of shear probes and P is the number
%      of dissipation estimates. M is the number of wavenumber elements in
%      each cross-spectrum. The "diagonal" elements of the N by N submatrix are
%      the auto-spectra.
% * [sh] same as sh_clean but without cleaning by the Goodman coherent noise
%       removal algorithm.
% * [AA] the cross-spectral matrix of acceleration signals for each
%      dissipation estimate. This is a 4-dimensional matrix of size
%      [M N N P] where N is the number of acceleration signals and P is the
%      number of dissipation estimates. M is the number of wavenumber
%      elements in each cross-spectrum. The "diagonal" elements of the N by
%      N submatrix are the auto-spectra.
% * [UA] the cross-spectral matrix of shear probe and acceleration
%       signals. Its size is [M N_s N_a P] where N_s is the number of shear
%       probes and N_a is the number of acceleration signals.
% * [F] frequencies, with one column for every dissipation estimate.
%      Currently every column is identical.
% * [K] Wavenumbers [cpm, or cycles per metre] with one column for
%      every dissipation estimate. The number of rows equals the number of
%      wavenumber elements in the cross-spectra.
%
%%% Mean values at each dissipation estimate
%   
% * [speed] a column vector of the mean speed at each dissipation estimate.
% * [nu] column vector of kinematic viscosity [m^2/s] at each dissipation
%        estimate. 
% * [P] column vector of pressure at each dissipation estimate.
%       Derived directly from the input.
% * [T] column vector of temperature at each dissipation estimate.
%       Derived directly from the input.
% * [t] column vector of time at each dissipation estimate.
%       Derived directly from the input.
% * [AOA] a vector containing the maximum angle of attack of the
%       velocity field, if you passed the 3 velocity components U V W. The
%       output is empty if you did not pass the velocity vectors.
%
%%% Additional Parameters
%
% * [Data_fast] a matrix where every column is a fast vector that has
%      been averaged over the interval of each dissaption estimate. The number
%      of rows equals the length of diss.P. This matrix can be empty.
% * [Data_slow] same as diss.Data_fast except that it is for the
%      slow vectors.
% * [fs_fast] actual fast sampling rate of the data, in Hz
% * [fs_slow] actual slow sampling rate of the data, in Hz
% * [f_AA] frequency, in Hz, that is 90% of the anti-aliasing frequency. 
%       Note: Different from input value of f_AA. 
% * [f_limit] Maximum frequency used when estimating the rate of
%       dissipation.
% * [fit_order] Order of the polynomial fit to the shear spectra, in log-log
%       space, that was used to estimate the spectral-minimum wavenumber.
% * [diss_length] The length of the segment, in points, used to estimate 
%       each dissipation rate. 
% * [overlap] The overlap of each dissipation estimate, in points. 
% * [fft_length] The length of the fft-segments, in points, that were
%       ensemble-averaged into a spectrum. 
% 
%

%__________________________
% *Version History:*
% 2012-05-27 (RGL) started at Narita airport.
% 2012-06-26 (RGL) added while loop to adjust upper limit of integration if
%            spectra agree poorly with Nasmyth for wavenumber just beyond the
%            spectral peak.
% 2012-07-11 (RGL and WID) changed detection of spectral minimum. Added
%            mordern parser.
% 2012-10-30 (RGL) Temporary modification to make K_max go to at least 100 cpm
% 2012-10-31 (RGL) Added feature to force a fit to the Nasmyth spectrum in the
%            +1/3 slope range, if requested.
% 2012-11-08 (WID) Replaced call to ismatrix() with the logic equivalent - to
%            allow the function to run on earlier versions of Matlab.
% 2012-11-09 (WID) Documentation update.
% 2013-07-15 (RGL) First attempt at the new and improved version based on
%            RSI Tech. Note XXX
% 2013-09-06 (RGL) Added ability to pass an arbitrary "Data" matrix which
%       gets averaged into range bins matching the dissipation estimates.
%       size(Data,1) must be same as size(shear,1), i.e. fast sampled data,
%       only. Added capability to average a slow data matrix. We now pass fs_slow.
%       Set method threshold to 1.5e-5 W/kg.
% 2013-12-13 (RGL) changed name of fast data matrix to Data_fast.
% 2013-12-23 (RGL) updated documentation, but have not yet added DFP and
%       GF, etc.
% 2014-03-15 (RGL) added f_limit -the parameter for the upper limit of
%       frequency used to estimate the rate of dissipation of TKE.
% 2014-03-17 (RGL) changed the returned value of f_AA to be f_AA instead of 
%       f_AA/f_AA_limit, where f_AA_limit= 0.9. When f_limit < f_AA/f_AA_limit,
%       then f_AA = f_limit, otherwise the limit used is f_AA/f_AA_limit to
%       avoid the slight attenuation occuring near f_AA.;
% 2014-12-22 (RGL) added the ability to estimate the maximum angle of
%       attack (AOA) for each dissipation estimate, by (optionally) passing
%       the 3-D velocity [U V W] to this function. The default is and empty
%       set of velocity vectors and no AOA calculation. If not empty, then
%       all three components must be passed and they must have the same
%       length as the shear vectors.
% * 2015-10-28 (RGL) Documentation changes.
% * 2015-11-18 (RGL) Documentation changes.
% * 2015-11-19 (RGL) Finally added dof (degrees of freedom) amd mad (mean
%       absolute deviation) to this function. 
% * 2016-08-05 (RGL) Fixed problem with pwermitation when there is only 1
%       shear probe. Accelerometers were not correctly permuted.
% * 2016-08-15 (RGL) The above fix was not quite right.
% * 2016-08-30 (RGL) Modified definitions of dof. dof_spec is now the
%       degrees of freedom of the spectral values. It is
%       2*2*9/11*num_of_fft_segments. dof_e is the degrees of freedom of a
%       dissipation estimate. It is the number of wavenumber points used to
%       estimate epsilon, divided by the Nyquist bandwidth, times the
%       diss_length (in samples). This is not the best estimate but will
%       have to do for now. I will come up with a proper relative bandwidth
%       estimate, later. Corrected again on 2016-08-31.
% * 2016-12-12 (RGL) Added check for if the Goodman routine should be
%       called. Added by-pass in case that Goodman is not used.
% * 2017-01-11 (RGL) set K_range to no less than 3 points in case of a
%       rediculously slow speed, for which K=0 may be the only wavenumber
%       smaller than 10 cpm.  
% * 2017-01-30 (RGL) Also had to set Range (near line 518) to a length no
%       shorter than 2 in order for the trapz function not to fail. 
% * 2017-04-03 (RGL) Added the return of default parameters in case their
%        are no input arguments. 
% * 2018-05-09 (RGL) Added check length(speed)==1, and converted it to a
%        vector of size t_fast, so that the function does not bomb on a
%        scalar variable for speed.
% * 2019-04-10 (RGL) Changed DOF of spectrum to 1.9*num_fft to agree with 
%        Nuttall (1971).
% * 2019-06-20 (JMM) Added the figure of merit (FM) as an output field in 
%        the dissipation structure.
% * 2019-07-03 (JMM) At the suggestion of RGL, the minimum upper limit of
%       integration was changed to 7 cpm (from 10 cpm). If necessary, the
%       limit will be increased to ensure that there at least 3 spectral 
%       points are used in the integration.
% * 2019-07-03 (JMM). Removed 'warning' from the output structure since it
%       is not used.
% * 2020-03-19 (JMM). Corrected num_fft (was previously too high by 1). Also
%		fixed the size of quality control fields (FM, mad, etc) for the case
%		when num_of_estimates == 1. Fields of output structure should now all 
%		have correct dimensions.
% * 2020-10-11 (JMM). Fixed the size of Data_fast and Data_slow for the case
%		when num_of_estimates == 1. 
	 

function diss = get_diss_odas(shear, A, varargin)

% Default values for optional fields
default_fit_order = 3;
default_f_AA      = 98;
default_fit_2_isr = 1.5e-5; % W/kg
default_f_limit   = inf; % upper frequency limit of dissipation estimation
default_UVW       = []; % U, V, W are the 3-D velocity components.
default_Data_fast = [];
default_Data_slow = [];
default_goodman   = true;

p = inputParser;
p.CaseSensitive = true;
p.KeepUnmatched = true;

val_fft_length = @(x) isnumeric(x)  && isscalar(x) && (x >= 2);
val_positive   = @(x) isnumeric(x)  && isscalar(x) && (x >= 0);
val_matrix     = @(x) isnumeric(x)  && (size(x,1) > 1) && (size(x,2) >= 1);
val_speed      = @(x) isnumeric(x)  && isvector(x) && (x >= 0);
val_data       = @(x) (isnumeric(x) && (size(x,1) > 1) && (size(x,2) >= 1)) || isempty(x);
val_logical    = @(x) (islogical(x));

addRequired(  p, 'shear',       val_matrix);
addRequired(  p, 'A',           val_matrix);
addParamValue(p, 'fft_length',  val_fft_length);
addParamValue(p, 'diss_length', val_fft_length);
addParamValue(p, 'overlap',     val_fft_length);
addParamValue(p, 'fs_fast',     val_positive);
addParamValue(p, 'fs_slow',     val_positive);
addParamValue(p, 'speed',       val_speed);
addParamValue(p, 'T',           val_matrix);
addParamValue(p, 't_fast',      val_data);
addParamValue(p, 'Data_fast',   default_Data_fast, val_data);
addParamValue(p, 'Data_slow',   default_Data_slow, val_data);
addParamValue(p, 'P',           val_matrix);
addParamValue(p, 'fit_order',   default_fit_order, val_positive);
addParamValue(p, 'f_AA',        default_f_AA,      val_positive);
addParamValue(p, 'f_limit',     default_f_limit,   val_positive);
addParamValue(p, 'fit_2_isr',   default_fit_2_isr, val_positive);
addParamValue(p, 'UVW',         default_UVW,       val_data);
addParamValue(p, 'goodman',     default_goodman,   val_logical);


if ~nargin,
    for d = whos('default_*')',
        param = regexp(d.name, 'default_(.*)', 'tokens');
        result.(char(param{1})) = eval(d.name);
    end
    % Add the arguments from odas_p2mat.  They will not be required but can
    % be used as a reference when calling quick_look directly.
    p2mat = odas_p2mat();
    for name = fieldnames(p2mat)'
        result.(name{1}) = p2mat.(name{1});
    end
    return
end

% Parse the arguments.
parse(p, shear, A, varargin{:});

% Perform last stages of input validation.
if p.Results.diss_length < 2*p.Results.fft_length,
    error('Invalid size for diss_length - must be greater than 2 * fft_length.');
end

if (    size(p.Results.shear,1) ~= size(p.Results.A,1) || ...
        size(p.Results.shear,1) ~= size(p.Results.T,1) || ...
        size(p.Results.shear,1) ~= size(p.Results.P,1)),
    error('Same number of rows required for shear, A, T, P.');
end

if ~isscalar(p.Results.speed) && size(p.Results.shear,1) ~= size(p.Results.speed,1),
    error ('speed vector must have the same number of rows as shear');
end

if p.Results.diss_length > size(p.Results.shear,1)
    error('Diss_length cannot be longer than the length of the shear vectors');
end

if ~isempty(p.Results.UVW) ...
        && size(p.Results.shear,1)  ~= size(p.Results.UVW,1) ...
        && size(p.Results.UVW,2) ~= 3,
    error (...
        'Velocity matrix, [U V W], must have the same number of rows as shear, and have 3 columns.');
end


fft_length   = p.Results.fft_length;
diss_length  = p.Results.diss_length;
overlap      = p.Results.overlap;
fs_fast      = p.Results.fs_fast;
fs_slow      = p.Results.fs_slow;
speed        = p.Results.speed;
T            = p.Results.T;
Data_fast    = p.Results.Data_fast;
Data_slow    = p.Results.Data_slow;
P            = p.Results.P;
fit_order    = p.Results.fit_order;
f_AA         = p.Results.f_AA;
f_limit      = p.Results.f_limit;
t_fast       = p.Results.t_fast;
e_inertial_sr= p.Results.fit_2_isr;
UVW          = p.Results.UVW;
goodman      = p.Results.goodman;

% For backwards compatibility...  "t" was renamed "t_fast"
if isfield(p.Unmatched, 't'), t_fast = p.Unmatched.t; end

estimate_AOA = false;
if ~isempty(UVW), estimate_AOA = true; end

f_AA_limit = 0.9;
f_AA = f_AA_limit*f_AA; % constrain limit to 90% of F_AA
if f_limit < f_AA, f_AA = f_limit; end

if length(speed) == 1
    speed = speed * ones(size(t_fast));
end

% end of input argument checking.
select = (1:diss_length)';

% Calculate the number of dissipation estimates that can be made.
if overlap >= diss_length, overlap = diss_length/2; end
if overlap < 0, overlap = 0; end
num_of_estimates = 1 + floor((size(shear,1) - diss_length) / (diss_length-overlap));
F_length = 1 + floor(fft_length/2); % size of frequency and wavenumber vectors
num_of_shear = size(shear,2);

% This is a screwy fix to take care of the case when num_of_estimates == 1.
% By default, Matlab removes the last dimension if it is 1.
reduce_before_exit = false;
if num_of_estimates == 1
    reduce_before_exit = true;
    num_of_estimates = 2;
end

% pre-allocate variables
diss.e            = zeros(num_of_shear, num_of_estimates);
diss.K_max        = diss.e;
diss.method       = diss.e;
diss.dof_spec     = 0;
diss.dof_e        = diss.e;
diss.mad          = diss.e;
diss.FM           = diss.e;


% first dimension is frequency index

diss.Nasmyth_spec = zeros(F_length, num_of_shear, num_of_estimates);

% fourth dimension is estimate index
diss.sh_clean     = complex(zeros(F_length, num_of_shear, num_of_shear, num_of_estimates));
diss.sh           = diss.sh_clean;
diss.AA           = complex(zeros(F_length,   size(A,2),  size(A,2),  num_of_estimates));
diss.UA           = complex(zeros(F_length, num_of_shear, size(A,2),  num_of_estimates));
diss.F            = zeros(F_length, num_of_estimates);
diss.K            = diss.F;
diss.speed        = zeros(num_of_estimates,1);
diss.nu           = zeros(num_of_estimates,1);
diss.P            = zeros(num_of_estimates,1);
diss.T            = zeros(num_of_estimates,1);
diss.t            = zeros(num_of_estimates,1);
if estimate_AOA
    diss.AOA      = zeros(num_of_estimates,1);
else
    diss.AOA = [];
end
diss.Data_fast=[];if ~isempty(Data_fast), diss.Data_fast = zeros(size(Data_fast,2), num_of_estimates);end
diss.Data_slow=[];if ~isempty(Data_slow), diss.Data_slow = zeros(size(Data_slow,2), num_of_estimates);end


a     = 1.0774e9; % From Lueck's model for e/e_10
x_isr = 0.01; % The non-dimensional wavenumber of the upper limit of the
%               inertial subrange.
x_isr = 2*x_isr; % test pushing this upward, also in sub-function
x_95  = 0.1205; % non-dimensional wavenumber for 95% variance

Nasmyth_spectrum  = zeros(F_length, num_of_shear);

index = 1;
if isempty(Data_fast), mean_Data_fast = []; end
if isempty(Data_slow), mean_Data_slow = []; end

num_of_ffts = 2 * floor(length(select) / fft_length) - 1;
%dof_spec = 2*2*(9/11)*num_of_ffts;
dof_spec = 1.9*num_of_ffts; % Nuttall (1971)

%
%Main loop
while select(end) <= size(shear,1)
    select_slow = select(1: round(fs_fast/fs_slow): end);
    select_slow = 1 + round((select_slow - 1) * fs_slow/ fs_fast);

    if goodman
        [P_sh_clean, AA, P_sh, UA, F] = ...
        clean_shear_spec(A(select,:), shear(select,:), fft_length, fs_fast);
    else
        push = [2 3 1]; %The permutation
        method = 'parabolic';
        [AA, F] = ...
            csd_matrix_odas(...
            A(select,:),     [], fft_length, fs_fast,[], fft_length/2, method);
        [P_sh, F] = ...
            csd_matrix_odas(...
            shear(select,:), [], fft_length, fs_fast,[], fft_length/2, method);
        P_sh       = permute(P_sh, push);
        AA         = permute(AA,   push);
        
        % Remove extra dimensions
        P_sh       = squeeze(P_sh);
        P_sh_clean = P_sh; % artifical to make the cod rum
        AA         = squeeze(AA);
        UA         = AA; % artificial to make the code run
    end
    
    % Note, P_sh_clean and P_sh are [M M N] matrices where M is the number of
    % columns in shear (the number of shear probe signals) and N is the length of
    % the frequency vector F (usually fft_length/2 + 1).
    % If there is only a single shear probe, these matrices are [N 1] in size,
    % and I use the "squeeze" function to remove the redundant dimensions.
    %
    % We now move the wavenumber index into the first dimension, for memory
    % efficiency
    if size(P_sh,2) > 1
        P_sh       = permute(P_sh,       [3 1 2]);
        P_sh_clean = permute(P_sh_clean, [3 1 2]);
    end
    if size(AA,2) > 1
        AA         = permute(AA,         [3 1 2]);
    end
    if size(P_sh,2) > 1 && size(AA,2) > 1
        UA         = permute(UA,         [3 1 2]);
    else
        UA = UA';
    end
        P_sh       = squeeze(P_sh); % remove extra dimension if shear is a vector
        P_sh_clean = squeeze(P_sh_clean);
        AA         = squeeze(AA);

        % Convert frequency spectra to wavenumber spectra
    W = mean(abs(speed(select)));

    K = F/W;
    K_AA = f_AA/W; % anti-aliasing wavenumber
    junk = repmat(K,[1 size(P_sh,2),size(P_sh,2)]);
    junk = squeeze(junk); % remove extra dimensions if shear is only a vector

    correction = ones(size(junk));
    correction_index = find(junk <= 150);
    correction(correction_index) = 1 + (junk(correction_index) / 48).^2;
    % wavenumber correction stops beyond 150 cpm.
    % Wavenumber correction after Macoun & Lueck
    P_sh_clean = P_sh_clean * W; 
    P_sh       = P_sh       * W;
    P_sh_clean = P_sh_clean .* correction;
    P_sh       =       P_sh .* correction;

    % We now integrate the shear spectrum to 10 cpm and use our e /e_10
    % model to predict the actual dissipation rate.
    e      = zeros(1,num_of_shear); % pre-assign, one column per shear probe
    K_max  = e;
    method = e;
    dof_e  = ones(size(e)) * diss_length; % The maximum possible degrees 
    %                               of freediom of dissipation estimates.
    mad    = e;
    
    mean_T = mean(T(select));
    mean_P = mean(P(select));
    mean_t = mean(t_fast(select));
    nu     = visc35(mean_T);
    if estimate_AOA
        AOA    = 180*atan(...
            max(sqrt(UVW(select,2).^2 + UVW(select,3).^2) ./ abs(UVW(select,1))))/pi;
    end
    if ~isempty(Data_fast), mean_Data_fast = mean(Data_fast(select,:));end
    if ~isempty(Data_slow), mean_Data_slow = mean(Data_slow(select_slow,:));end

    for column_index = 1:num_of_shear
        if num_of_shear == 1 % we have only a single shear probe
            shear_spectrum = P_sh_clean;
        else % the auto-spectra are on the diagnonal
            shear_spectrum = P_sh_clean(:, column_index, column_index);
        end
        shear_spectrum = squeeze(shear_spectrum);
        K_range = find(K <= 10); % all wavenumber smaller than 10 cpm
        if length(K_range) < 3, K_range = [1 2 3]'; end % in case of very low speed
        e_10 = 7.5*nu*trapz(K(K_range), shear_spectrum(K_range));
        e_1 = e_10*sqrt(1 + a*e_10); % This is the first estimate of the
        % rate of dissipation, using Lueck's model.
        if e_1 < e_inertial_sr % Use the variance method
            if length(find(K.*sqrt(sqrt(nu^3 ./ e_1)) <= x_isr)) >= 20
                % refine the dissipation estimate by fitting to the
                % inertial subrange.
                % disp(['Refining at ' num2str(index)])
               [e_2, junk] = ...
                    inertial_subrange(K, shear_spectrum, e_1,nu, min([150, K_AA]));
                % dissipation from inertial subrange method
            else
                e_2 = e_1; % Use the estimate that we already have.
            end
            K_95 = x_95*sqrt(sqrt(e_2 ./ nu^3));
            valid_shear = find( K <= fit_in_range(min([K_AA, K_95]), [0 150]) );
            Index_limit = length(valid_shear);
            y = log10(shear_spectrum(2:Index_limit));
            x = log10(             K(2:Index_limit));
            % Note K(1) is always zero and has problems with logarithms. Now we fit a
            % polyninomial to the spectrum

            % Keep the fit_order reasonable. Any value from 3 to 8 could work
            fit_order = fit_in_range(fit_order, [3 8]);

            if Index_limit > fit_order + 2; % We have enough points for polyfit
                % p   - polynomial of a polynomial fit to spectrum
                % pd1 - first derivative of the polynomial
                % pr1 - roots of first derivative of the polynomial
                if sum(isnan(y))>0%JMM - to prevent getting errors if there are nans - might be better to fix clean shear spec code
                    disp([num2str(sum(isnan(y))),' NaNs in spectrum'])
                    if sum(isnan(y)) == 1 % Replace NaN with average
                        ind = find(isnan(shear_spectrum));
                        shear_spectrum(ind) = 0.5*(shear_spectrum(ind-1)+shear_spectrum(ind+1));
                        
                        y = log10(shear_spectrum(2:Index_limit));
                    else
                        continue %simply skip -- otherwise I get stuck in an infinite loop
                    end
                end %JMM (end)
                p = polyfit(x, y, fit_order);
                pd1 = polyder(p);
                pr1 = sort( roots(pd1) );
                pr1 = pr1(~imag(pr1)); % remove complex roots
                % Filter roots so that only the minima above 10cpm remain.
                if ~isempty(pr1)
                    pr1 = pr1( polyval( polyder(pd1), pr1 ) > 0 ); % minima only
                    pr1 = pr1( pr1 >= log10(10) ); % spectral min must be at 10cpm or higher
                    % Fit root within a given range.
                    if isempty(pr1), pr1 = log10(K_95);else pr1=pr1(1);end
                else
                    pr1 = log10(K_95); % if there is no minimum, use K_95
                end
            else %Not enough points for polyfit
                pr1 = log10(K_95); % if there is no minimum, use K_95
            end

            K_limit = fit_in_range( min([pr1, log10(K_95), log10(K_AA)]), ...
                [log10(7) log10(150)], log10(150)); % set to 7cpm, RGL, 2019-02-28
            Range = find( K <= 10^(K_limit)); %Integration range for shear probes.
            if K(Range(end))<7 % Integration limit is greater than 7cpm
                Range = [Range; Range(end)+1];
            end
            if length(Range) < 3 % At least 3 spectral points to integrate over
                Range = [1:3]; 
            end
            e_3 = 7.5*nu*trapz(K(Range), shear_spectrum(Range));

            % Next, the non-dimensional limit of integration
            x_limit = K(Range(end))*sqrt(sqrt(nu^3 ./ e_3));
            x_limit = x_limit^(4/3); % A more convinient form
            % Next the variance resolved according to Lueck's model
            variance_resolved = tanh(48*x_limit) - 2.9*x_limit*exp(-22.3*x_limit);
            e_new = e_3 ./ variance_resolved; % Adjust upward to correct the variance
            done = 0; % while loop control
            while ~done
                x_limit = K(Range(end))*sqrt(sqrt((nu^3 ./ e_new)));
                x_limit = x_limit^(4/3); % A more convinient form
                % Next the variance resolved according to Lueck's model
                variance_resolved = tanh(48*x_limit) - 2.9*x_limit*exp(-22.3*x_limit);
                e_old = e_new;
                e_new = e_3 ./ variance_resolved; % Adjust upward to correct the variance
                if e_new / e_old < 1.02
                    done = 1;
                    e_3 = e_new;
                end
            end
            % Correct for missing variance at bottom end of the spectrum
            phi = nasmyth(e_3, nu, K(2:3)); % you need to give it 2 or more wavenumbers
            phi = phi(1); % the Nasmyth value for the lowest non-zero wavenumber
            e_4 = e_3 + 0.25*7.5*nu*K(2)*phi; % add the missing variance
            if e_4 / e_3 > 1.1 % re-loop to find dissipation rate
                e_new = e_4 ./ variance_resolved; % Adjust upward to correct the variance
                done = 0; % while loop control
                while ~done
                    x_limit = K(Range(end))*sqrt(sqrt(nu^3 / e_new));
                    x_limit = x_limit^(4/3); % A more convinient form
                    % Next the variance resolved according to Lueck's model
                    variance_resolved = tanh(48*x_limit) - 2.9*x_limit*exp(-22.3*x_limit);
                    e_old = e_new;
                    e_new = e_4 / variance_resolved; % Adjust upward to correct the variance
                    if e_new / e_old < 1.02
                        done = 1;
                        e_4 = e_new;
                    end
                end
            end
            K_max(column_index) = K(Range(end));

        else % e >= e_inertial_sr, so fit to the inertial subrange
            %disp(['Using inertial subrange for P = ' num2str(mean(P(select)))]);
            K_limit = min([K_AA, 150]);
            [e_4, K_max(column_index), Range] = ...
                inertial_subrange(K, shear_spectrum, e_1, nu, K_limit);
            method(column_index) = 1;% identifies the method
        end
        
        e(column_index)   = e_4; % save it
        Nasmyth_spectrum(:,column_index) = nasmyth(e_4,nu,K);
%        dof(column_index) = dof(column_index) * length(Range - 1);
        dof_e(column_index) = dof_e(column_index)*length(Range)/length(K);
        This_Nasmyth        = Nasmyth_spectrum(Range(2:end), column_index);
        This_spectrum       =   shear_spectrum(Range(2:end));
        mad(column_index)   = mean(abs(log10(This_spectrum ./ This_Nasmyth)));
    end

    diss.e           (:,index)     = e;
    diss.dof_e       (:,index)     = dof_e;
    diss.mad         (:,index)     = mad;
    diss.FM          (:,index)     = mad .* sqrt(dof_spec);
    diss.Nasmyth_spec(:,:,index)   = Nasmyth_spectrum;
    diss.K_max       (:,index)     = K_max;
    diss.method      (:,index)     = method;
    diss.sh_clean    (:,:,:,index) = P_sh_clean;
    diss.sh          (:,:,:,index) = P_sh;
    diss.AA          (:,:,:,index) = AA;
    diss.UA          (:,:,:,index) = UA;
    diss.F           (:,index)     = F;
    diss.K           (:,index)     = K;
    diss.speed       (  index)     = W;
    diss.nu          (  index)     = nu;
    diss.T           (  index)     = mean_T;
    diss.t           (  index)     = mean_t;
    diss.P           (  index)     = mean_P;


    if estimate_AOA, diss.AOA         (  index)     = AOA; end
    if ~isempty(Data_fast), diss.Data_fast(:,index) = mean_Data_fast; end
    if ~isempty(Data_slow), diss.Data_slow(:,index) = mean_Data_slow; end
    index = index + 1;

    select = select + diss_length - overlap;
end
% return the parameters used to estimate the dissipation rate.
if reduce_before_exit
% Remove the last dimension because there is only one estimate
    diss.e            = diss.e             (:,1);
    diss.Nasmyth_spec = diss.Nasmyth_spec(:,:,1);
    diss.K_max        = diss.K_max         (:,1);
    diss.method       = diss.method        (:,1);
    diss.sh_clean     = diss.sh_clean  (:,:,:,1);
    diss.sh           = diss.sh        (:,:,:,1);
    diss.AA           = diss.AA        (:,:,:,1);
    diss.UA           = diss.UA         (:,:,:,1);
    diss.F            = diss.F             (:,1);
    diss.K            = diss.K             (:,1);
    diss.speed        = diss.speed           (1);
    diss.nu           = diss.nu              (1);
    diss.T            = diss.T               (1);
    diss.t            = diss.t               (1);
    diss.P            = diss.P               (1);
    diss.dof_e        = diss.dof_e         (:,1);
    diss.mad          = diss.mad           (:,1);
    diss.FM           = diss.FM            (:,1);
    if estimate_AOA, diss.AOA  = diss.AOA    (1); end
    if ~isempty(diss.Data_slow), diss.Data_slow= diss.Data_slow(:,1); end
    if ~isempty(diss.Data_fast), diss.Data_fast= diss.Data_fast(:,1); end
end

diss.fs_fast     = fs_fast;
diss.fs_slow     = fs_slow;
diss.f_AA        = f_AA;
diss.f_limit     = f_limit;
diss.fit_order   = fit_order;
diss.diss_length = diss_length;
diss.overlap     = overlap;
diss.fft_length  = fft_length;
diss.dof_spec    = dof_spec;


end

    function result = fit_in_range(value, range, default)
        % Fit first element of value into range.  Returns scalar result.
        % - range is of the form [min max]
        % - if value is empty, use default value

        if ~isempty(value),
            result = value(1);
        else
            result = default(1);
        end

        if result < range(1), result = range(1); end
        if result > range(2), result = range(2); end

    end

    function [e, K_max, fit_range] = inertial_subrange(K, shear_spectrum, e, nu, K_limit)
        % function to find the best fit to the inertial subrange by adjusting the
        % rate of dissipation.
        % K - is the wavenumbers, including zero
        % shear_spectrum - is the wavenumer spectrum of shear
        % e - starting value of the rate of dissipation, usually derived from the
        %       e/e_10 model of Lueck.
        % K_limit - wavenumber limits other than the inertial subrange, such as
        %       150 cpm or K_AA (the anti-aliasing wavenumber), or K_95.
        %
        % 2013-07-15 (RGL) first version.

        x_isr = 0.01; % nondimensional limit of the inertial subrange
        x_isr = 2*x_isr;% try pushing this up, a little
        %The inertial subrange for fitting is,
        fit_range = find(K <= min([x_isr*sqrt(sqrt(e / nu^3)), K_limit]));
        K_max = K(fit_range(end));
%         if K_max < 10*K(2) %less than 9 points in inertial subrange
%             K_max = K(2); % add a warning or error flaag here
%             return
%         end

        for loop_count = 1:3
            Nasmyth_values = nasmyth(e,nu,K(fit_range));
            fit_error = ...
                mean(log10(shear_spectrum(fit_range(2:end)) ./ Nasmyth_values(2:end)));
            e = e * 10^(3*fit_error/2);% scale up (or down) epsilon
            % Now do it all again, two times only, in case the change was large
        end

        % Now remove flyers
        Nasmyth_values = nasmyth(e,nu,K(fit_range));
        fit_error = ...
            log10(shear_spectrum(fit_range(2:end)) ./ Nasmyth_values(2:end)); % here it is a vector
        flyers_index = find(abs(fit_error) > 0.5);
        if ~isempty(flyers_index)
            [Y,index] = sort(fit_error(flyers_index),'descend'); % sort flyers in descending order
            bad_limit = ceil(0.2*length(fit_range));% Remove no more than 20%
            if length(Y) > bad_limit, index = index(1:bad_limit); end
            fit_range(index) = [];
            shear_spectrum(index) = [];
        end
        K_max = K(fit_range(end)); % in case it is reduced by flyers

        % Re-fit to the inertial subrange using reduced spectrum
        for loop_count = 1:2
            Nasmyth_values = nasmyth(e,nu,K(fit_range));
            fit_error = ...
                mean(log10(shear_spectrum(fit_range(2:end)) ./ Nasmyth_values(2:end)));
            e = e * 10^(3*fit_error/2);% scale up epsilon
            % Now do it all again, one time only, in case the change was large
        end
    end

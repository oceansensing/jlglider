%% deconvolve
% Deconvolve signal plus derivative to produce high resolution data
%%
% <latex>\index{Functions!deconvolve}</latex>
%
%%% Syntax
%   xHires = deconvolve( name, X, X_dX, fs, diff_gain, ver )
%
% * [name] Channel name that contains the 'diff_gain' parameter.  Can be
%       left empty if 'diff_gain' is explicitly provided.
% * [X] Low-resolution signal.  Data without pre-emphasis is not required
%       but will improve the initial accuracy of the output if available.
%       This is most notable for high gain channels such as pressure.
% * [X_dX] Pre-emphasized signals.
% * [fs] Sampling rate of the data in Hz.
% * [diff_gain] Differentiator gain. Provide either a numeric value, the
%       configuration string, or the configuration object generated from
%       parsing the configuration string.
% * [ver] (Depreciated - not required) ODAS header version. Included for
%       backward compatibility.
% * []
% * [xHires] Deconvolved, high-resolution signal.
%
%%% Description
% Deconvolve a vector of pre-emphasized data (temperature, conductivity, or
% pressure) to yield high-resolution data. The pre-emphasized signal
% (x+gain*dx/dt) is low-pass filtered using appropriate initial conditions
% after the method described in Mudge and Lueck, 1994.
%
% For pressure, you must pass both 'X' and 'X_dX' to this function. Both
% vectors are needed to make a good estimate of the initial conditions for
% filtering. Initial conditions are very important for pressure because the
% gain is usually $\SI{\sim20}{\s}$, and errors caused by a poor choice of initial
% conditions can last for several hundred seconds! In addition, the
% high-resolution pressure is linearly adjusted to match the low-resolution
% pressure so that the factory calibration coefficients can later be used to
% convert this signal into physical units.
%
% The gains for all signal types are provided in the calibration report of
% your instrument.
%
%%% Examples
%
%    >> C1_hres = deconvolve( 'C1_dC1', [], C1_dC1, fs_fast, setupfilestr )
%
% Deconvolve the micro-conductivity channel using the diff_gain parameter
% from the 'C1_dC1' channel section found within the supplied configuration
% string.  Because there is no channel without pre-emphasis, the argument is
% left empty.
%
%    >> T1_hres = deconvolve( '', T1, T1_dT1, fs_fast, 1.034 )
%
% Deconvole the thermistor data using both channels with and without
% pre-emphasis.  Using both channels improves the initial deconvolution
% accuracy and compensates for slight variations in the calibration
% coefficients between the two channels. Note the explicitly provided value
% for diff_gain. One typically provides the configuration string and
% channel name as shown in the previous example.
%
% Mudge, T.D. and R.G. Lueck, 1994: Digital signal processing to enhance
% oceanographic observations, _J. Atmos. Oceanogr. Techn._, 11, 825-836.
%
% @image @images/pressure_deconvolution1 @Deconvolution example. @The green
% curve is from the normal pressure channel, and the blue curve is derived from
% the pre-emphasized pressure channel. This data is from a profiler that has impacted
% the soft bottom of a lake. Both signals are shown with full bandwidth (0 - 32 Hz)
% without any smoothing. The full-scale range of the pressure transducer is 500 dBar.
%
% @image @images/pressure_deconvolution2 @Deconvolution example two. @Same
% as previous figure but with zoom-in on only the high-resolution pressure. Again,
% full bandwidth without any smoothing.
%
% @image @images/pressure_deconvolution3 @The rate of change of pressure derived
% from the normal pressure signal and the high-resolution pressure signals using
% the gradient function. @Full bandwidth signals of the rate of change of
% pressure. The normal pressure signal (green) produces a fall-rate that is
% useless without extensive smoothing because it is not even
% monotonic. The rate of change of the high-resolution pressure (blue) is smooth, always
% positive, and, therefore, the high-resolution pressure itself, can be used directly for
% plotting other variables as a function pressure (or depth). The
% high-resolution rate of change of pressure has been multiplied by 10 for visual
% clarity. The fall-rate is about $\SI{0.17}{\m\per\s}$.

% *Version History:*
%
% * 2005-02-08 (IG) based on OTL/RGL routines pressure_all, thermistor_all,
%        pressure_vel
% * 2010-01-14 (AWS) odas v6 related changes
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2010-02-28 (RGL) reverted back to simple initial conditions for the case of
%        pressure. The previous method was not right and can give very erroneous
%        results for the case of a pressure transducer on a near-surface mooring
%        with significant waves.
% * 2012-04-11 (WID) calls to inifile_with_instring changed to setupstr
% * 2012-11-05 (WID) documentation update
% * 2013-05-16 (WID) Reworked the "filtic" calls based on suggestions from
%        Rolf.
% * 2015-07-20 WID Modified algorithm - now uses channel without
%       pre-emphasis when available to improve guess of initial conditions.
% * 2015-07-20 WID Removed requirements for XMP and data file versions.
% * 2015-07-23 WID Provided some documentation for when a single input argument
%       is used. Changed the calling of "interp1_if_required".
% * 2015-10-28 (RGL) Documentation corrections.
% * 2015-11-18 (RGL) Documentation corrections.


function X_hires = deconvolve(data_type,X,X_dX,f_s,arg5,ver)

% When read_odas is used to read an ODAS data file, a simplified version of
% deconvolve can be used.  This simplified version loads the required parameters
% from the caller's workspace.  It can do this because read_odas generates
% variables with predictable names.
%
% This feature simplifes performing a deconvolve from the Matlab command prompt.
% If referencing deconvolve from another function it is advised that all input
% parameters are explicitly provided.

if nargin == 1
    setupfilestr = evalin( 'caller', 'setupfilestr' );

    % Allow for channel names to be given either with or witout pre-emphasis.
    % For example, "T1" and "T1_dT1" will both work.
    loc = strfind( data_type, '_d' );
    if isempty(loc)
        X_name    = data_type;
        X_dX_name = [data_type '_d' data_type];
    else
        X_name    = data_type(1:loc(1)-1);
        X_dX_name = data_type;
    end

    data_type = X_dX_name;

    % Not all sensors have data without preemphasis
    try, X = evalin( 'caller', X_name ); catch; X = []; end

    try
        X_dX = evalin( 'caller', X_dX_name );
        fast = evalin( 'caller', 'fs_fast' );
        f_s  = channel_sampling_rate( X_dX_name, setupfilestr, fast );
        arg5 = setupfilestr;

        X_hires = deconvolve(data_type,X,X_dX,f_s,arg5);
        return
    catch
        error(['Unable to call deconvole the easy way - you will have ' ...
               'to explicitly set the required parameters.']);
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Parameters %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5                 % Default differentiator gains
    error('Invalid number of input arguments.');
end

if isnumeric(arg5)
    diff_gain = arg5;
elseif ischar(data_type)
    %we are dealing with the setup file str being passed for arg5
    diff_gain_key = 'diff_gain';
    tmp = setupstr(arg5, data_type, diff_gain_key);
    if isempty(tmp)
        error('Can not find the following key: %s', diff_gain_key);
    end
    diff_gain = str2double(tmp);
else
    error(['Must provide either the differentiator gain or the channel' ...
           'name / configuration string in which it can be found.']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Prepare for filtering: coefficients, initial conditions %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Only perform an interpolation if the vectors are a different size.
X = interp1_if_required( X, X_dX, f_s );

f_c = 1/(2*pi*diff_gain);   % Cut-off frequency for filtering

[b,a] = butter(1,f_c/(f_s/2));       % Low-pass filtering coeffs.
if length(X) > 1
    p = polyfit(X,X_dX,1);
    if p(1)<-0.5, X_dX=-X_dX; end  % A few instruments have X_dX inverted with respect to X
    z = filtic(b,a,X(1),X_dX(1));
else
    % For T/C, initial conditions based on first portion of the record,
    % because we usually do not have the vector without pre-emphasis.
    timeV = (0:1/f_s:2*diff_gain)';   % Time vector for calucating initial filter value
    p = polyfit( timeV, X_dX(1:length(timeV)), 1 );
    previousOutput = p(2) - diff_gain * p(1);
%    t_ave = 2*diff_gain; % Base avering time scale on the differentiator gain.
%    z = filtic(b,a,mean(X_dX(1:round(t_ave*f_s))),X_dX(1));
    z = filtic(b, a, previousOutput, X_dX(1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Deconvolve to obtain the high-resolution data %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X_hires = filter(b, a, X_dX, z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% For the pressure, regress against the low-resolution vector to %%%%%
%%%%%% remove the small offset in the derivative circuit.             %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Algorithm for signals that have both with and without pre-emphasis data.
% 1) Interpolate the non pre-emphasis signal so we can perform a polyfit.
% 2) Deconvolve using the first "X" data point for the initial condition.
% 3) Polyfit the deconvolved signal to the non pre-emphasis signal.
% 4) Calculate new initial conditions after applying an offset to X based
%    from the polyfit applied in step 3.
% 5) Perform a new deconvolve with updated initial conditions.
% 6) Apply the polynomial from step 3 to the new data.  A new polyfit is
%    not required - they will be almost identicle.

if length(X) > 1
    warning('off');
    p = polyfit(X_hires, X, 1);
    p2 = [2-p(1), -p(2)];

    initialOutput = polyval(p2, X(1));
    z = filtic(b, a, initialOutput, X_dX(1));
    X_hires = filter(b, a, X_dX, z);
    X_hires = polyval(p, X_hires);
    warning('on');
end

end


% Perform an interpolation if vector "X" has length > 1 but is not equal
% to the length of "X_dX".  This is to ensure that X and X_dX can be polyfit
% against one another.
function newX = interp1_if_required( X, X_dX, fs )
    newX = X;
    if length(X) == length(X_dX) || length(X) <= 1
        return
    end

    f_slow = fs * length(X) / length(X_dX);
    t_slow = (0:length(X)-1)' / f_slow;
    t_fast = (0:length(X_dX)-1)' / fs;
    newX = interp1(t_slow, X, t_fast, 'pchip','extrap');
end

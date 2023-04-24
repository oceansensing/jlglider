%% convert_odas
% Convert data from raw counts to physical units
%%
% <latex>\index{Functions!convert\_odas}</latex>
%
%%% Syntax
% %  [phys, units] = convert_odas(X, name, empty, config, ver )
%
% * [X] Data vector to convert.
% * [name] Name of the channel containing the calibration coefficients or
%       cell array of names if coefficients exist in multiple channels.
% * [empty] Parameter ignored unless exactly 3 parameters are used in the
%       function call. When three parameters are used, it is assumed
%       that this parameter contains what would have been provided within
%       the 'config' parameter.
% * [config] Source of the calibration coefficients. This is typically the
%       variable 'setupfilestr' but can also be the name of a file
%       containing the configuration string.  Alternatively, this can be a
%       structure containing the configuration string within the field
%       'setupfilestr'.
% * [ver] Parameter ignored, can optionally be removed.  Included for
%       backwards compatibility.
% * []
% * [X_phys] X converted into physical units.
% * [units] String representation of the converted variable units.
%
%%% Description
% Convert data from raw counts into physical units.  Data acquisition does
% not modify the raw data coming from the sensors.  All data can
% be converted using this function.
%
% Information required to perform the conversion from counts to physical
% units is contained within the configuration string. Users must provide
% convert_odas with the configuration string and the channel section name.
% The channel section name can be either a string or a cell array of
% strings. When calibration coefficients exist in multiple channels, the
% cell array should be used. If duplicate parameters exist, preference is
% given to the parameter that is found first when scanning the cell array
% from front to back. The convert_odas function will then extract the
% calibration coefficients and perform the conversion.
%
% Different types of channels use different conversion algorithms. The
% ``type" parameter within the configuration-string determines which
% algorithm is used. There are functions within convert_odas.m for each
% supported type. To see exactly how conversion is performed for a specific
% type, refer to the appropriate embedded function.
%
% Users can add support for new types or change the behaviour of existing
% types by suppling their own conversion functions. See the example
% included within this function for more details.
%
% For detailed information on the conversion from raw data [counts] to
% physical units, please see Section $\ref{sec:channel_types}$.
%
%%% Examples
%
%    >> T2      = convert_odas( T2,      'T2', '', setupfilestr )
%    >> T2_fast = convert_odas( T2_fast, 'T2', '', setupfilestr )
%
% Convert data using calibration coefficients from section 'T2' found
% within the configuration file 'setupfilestr'. The 'T2_fast' data is the
% result of deconvolving the T2_dT2 channel. The correct
% calibration coefficients for this channel are found within the 'T2'
% channel section of 'setupfilestr'.
%
%    >> P_slow = convert_odas( P_slow, 'P', setupfilestr )
%
% Convert the deconvolved pressure data, P_slow, using the calibration
% coefficients found within the 'P' channel section of the configuration
% string, 'setupfilestr'. 
%

% *Version History:*
%
% * 2005-02-15 (IG) Based on RGL/OTL routines including plot_tomi, sb_t, sb_c,
%               accel_all, and pressure_all.
% * 2005-02-21 (IG) now case-insensitive
% * 2005-05-01 (IG) corrected a small typo in the coefficient loading section
% * 2005-07-12 (IG) Provided an option to find the start of coefficient lines
%               that will work with pre-R14 versions of MATLAB.  Made commas the
%               only acceptable delimiters
% * 2006-04-17 (RGL) Added 3-axis magnetometer section.
% * 2007-06-16 (RGL) Added recognition for Alec Electronics EM Current meter,
%               Altimeter and Fluorometers
% * 2007-11-05 (RGL)  Added SBE43F oxygen sensor. Only gives the frequency and
%               not the proper disolved oxygen concentration.
% * 2010-01-14 (AWS) odas v6 related changes
% * 2010-09-24 (AWS) added adis function, inclinometer x, y, t conversion
% * 2010-10-04 (AWS) added documentation tags for matlab publishing
% * 2010-12-03 (AWS) added conversion for thermistor and micro conductivity
% * 2011-03-24 (AWS) changed conversion for voltage type to include sensor Gain
% * 2011-11-19 (RGL) Corrected voltage tBizHub-BWype so that it divides by the gain
%               instead multiplying by the gain.
% * 2012-03-04 (RGL) big changes to make this whole damn process more
%               rational. You can now use a sensor (or channel) name to
%               find the type and convert it to physical units. It is no
%               longer necessary to know the name of the section in which
%               the channel and its coefficients are identified. Brother!!
% * 2012-03-21 (RGL) changed thermistor conversion to avoid taking log of
%               negative numbers if thermistor is broken.
% * 2012-03-23 (RGL) changed the way micro-conductivity is calculated.
% * 2012-03-28 (RGL) Expanded the XMP section so that it can habdle all
%               channels.
% * 2012-04-11 (WID) changed inifile_with_instring calls to setupstr
% * 2012-04-30 (WID) fix for incl?? channel conversions with pre-v6 data files
% * 2012-05-07 (WID) update to documentation
% * 2012-11-02 (WID) update to documentation
% * 2013-02-26 (WID) support for beta_1 and beta_2 coefficients on therm channels
% * 2013-04-09 (WID) support any number of coefficients in a poly channel
% * 2013-04-23 (WID) support for arbitrary types.  User provided functions can be
%                    used to convert different types.  All v6 and higher type
%                    conversion moved to external functions found at the end of
%                    this file.
% * 2013-05-23 (WID) support for unsigned data types.  Data defaults to
%                    type signed.  XMP support removed - no longer required
% * 2013-06-17 (WID) Added support for type xmp_volt - redirected to type poly.
% * 2013-07-04 (WID) removed support for sign parameter - moved to read_odas
% * 2013-07-04 (WID) some JAC_C and JAC_T fixes.
% * 2013-07-09 (WID) added support for "jac_t" type.  ODAS-RT requires the
%                    type even though it is processed as a poly.
% * 2014-02-06 (RGL) added support for MicroSquid thermometer and
%                    conductivity sensors. Types are t_ms and c_ms,
%                    respectively.
% * 2014-05-13 (WID) added support for adc_zero in voltage type
% * 2014-10-22 (RGL) added support for MicroPod thermistor and shear probes
% * 2015-03-30 (RGL) corrected problems with bitshifting in odas_jac_c_internal
% * 2015-04-28 (WID) added support for rotation and Vector types
% * 2015-05-07 (RGL) added support for AMT_O2 sensor.
% * 2015-05-08 (RGL) corrected entry for Sea-Bird 43F oxygen sensor.
% * 2015-07-23 (WID) stripped out legacy code.  No longer require all input
%                    arguments - works with just 3.
% * 2015-07-26 (WID) added support for using a structure for the config
%                    input. Modified documentation.  Added more comments.
% * 2015-10-27 (RGL) Minor documentation corrections.
% * 2015-11-12 (RGL) Added type piezo to the conversion.
% * 2015-06-24 (RGL) Added type jac_emc.
% * 2017-01-23 (WID) Removed jac_emc and replaced with two types - aem1g_a
%                    for analog output and aem1g_d for digital output. 
%                    Old type still supported but a depreciated warning is 
%                    generated.
% * 2017-01-25 (WID) Added warning for incorrect A,B coefficients in aem1g.
% * 2017-03-18 (WID) Support for rsijac_c and rsijac_t types.
% * 2017-03-23 (WID) Corrected rsijac_c and rsijac_t types - note they will
%                    change in the future.
% * 2017-04-26 (WID) Added support for RINKO FT instruments, model name
%                    ARO-FT, with types aroft_o2 and aroft_t.
% * 2017-04-25 (WID) Temperary fix for buggy RINKO RS232
% * 2017-05-17 (WID) Final corrections for RSIJAC-CT sensors.
% * 2017-05-18 (WID) Final corrections for RSIJAC-CT + default values for
%                    ADIS inclinometer thermistor.
% * 2018-04-23 (WID) Added support for type "raw"
% * 2019-02-10 (WID) Renamed the variable "V" to "Value__RSI" to prevent
%                    confusion of channel name and local variable.
% * 2020-01-30 (JMM) Added units for the Rinko DO oxygen channel (aroft_o2)
% * 2021-05-16 (WID) Bug fix for jac_C calculation


function [X_phys, units] = convert_odas(X, name, empty, config, ver )

% At least 3 parameters are required.  If only three are used, assume the
% empty parameter contains what should be in the setupfilestr parameter.
if nargin < 3
    error('Invalid number of input arguments.');
elseif nargin == 3
    config = empty;
end

% Turn name into a cell array if it is a string.  This is to facilitate
% extracting coefficients from multiple channel names - such as T and T_dT.
if ischar(name)
    name = {name};
end

% Allow the configuration string to be passed in within a structure.  This
% allows for future flexibility without changing the calling parameters
% thereby preserving backwards compatability.  Currently not being used.
if isstruct(config) && isfield(config, 'cfgobj')
    config = config.cfgobj;
elseif isstruct(config) && isfield(config, 'setupfilestr')
    config = config.setupfilestr;
end

% The configuration string can be within a file.
if ischar(config) && exist(config, 'file')
    config = fileread( config );
end

% Read / parse the configuraiton string.
cfg = setupstr( config );

% Copy the parameters from the setup file into a "coef" structure to allow
% for easier access.
for chname = name
    [S,K,Value__RSI] = setupstr(cfg, chname, '', '');
    for i = 1:length(K)
        try Value__RSI{i} = eval(Value__RSI{i}); catch, end
        try
            if ~isfield('params', K{i})
                params.(K{i}) = Value__RSI{i};
            end
        catch
        end
    end
end


% coefficient names correspond to the parameter names used in the setup file
odas_type = setupstr(cfg, params.name, 'type');

extName = ['odas_' odas_type{1}];
intName = ['odas_' odas_type{1} '_internal'];

% Look for the function required for conversion into physical units. There
% should be an internal function declared below, but first look to see if
% the user has provided one of their own.
if exist(extName, 'file')
    typef = str2func(extName);
elseif exist(intName, 'file')
    typef = str2func(intName);
else
    error(['Unknown sensor type: ' odas_type{1}]);
end

[X_phys, units] = typef( X, params );

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Example Custom Conversion Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [physical, units] = odas_radToDeg( input, params )
%     units = '[ deg ]';
%     physical = input * (180 / pi) * params.scale;
% end
%
% This function allows for a custom type, in this case "radToDeg".  The
% function will be called automatically when attempting to convert a
% channel of this type.  The user is only required to name the function
% accoringly and be sure it is within the Matlab search path.
%
% The following channel segment from a configuration file could be used for
% the previous example.
%
%   [channel]
%   id = 99
%   name = myAngleChannel
%   type = radToDeg
%   scale = 1.0;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [physical, units] = odas_jac_c_internal( input, params )
    units = '[ s ]';
    i = floor(input / 2^16); % get the MS byte-pair
    v = rem(input, 2^16);
    physical = i ./ v;
    polyVals = [params.c, params.b, params.a];
    physical = polyval(polyVals, physical);
end

function [physical, units] = odas_jac_t_internal( input, params )
    units = '[ C ]';
    input(input < 0) = input(input < 0) + 2^16;  % not required - fixed in read_odas
    polyVals = [params.f, params.e, params.d, params.c, params.b, params.a];
    physical = polyval(polyVals, input);
end

function [physical, units] = odas_rsijac_c_internal( input, params )
    units = '[ mS/cm ]';
    try df     = params.division_factor; catch, df     = 2^28;      end
    try gain   = params.gain;            catch, gain   = 10.88/100; end
    try offset = params.y_offset;        catch, offset = 2.86e-7;   end
    try a      = params.a;               catch, a      = 0;         end
    try b      = params.b;               catch, b      = 1;         end

    % Input values are in Q3.28 fixed-point format with implied decimal
    % point. Convert to floating point.
    counts = input / df;

    % Apply polynomial correction for electronics calibration.
    counts = polyval([b a], counts);

    % Conductivity formula
    physical = (counts * gain - offset) / params.cell_constant;

    % Return in units of ms/cm
    physical = physical * 10;
end

function [physical, units] = odas_rsijac_t_internal( input, params )
    units = '[ C ]';
    t0 = params.t_0;
    b1 = params.beta_1;
    b2 = params.beta_2;
    b3 = params.beta_3;
    try df = params.division_factor; catch, df = 2^28; end
    try a  = params.a;               catch, a  = 0;    end
    try b  = params.b;               catch, b  = 1;    end

    % Input values are in Q3.28 fixed-point format with implied decimal
    % point. Convert to floating point.
    counts = input / df;
    
    % Apply polynomial correction for electronics calibration.
    counts = polyval([b a], counts);

    kelvin = 1 ./ (1./t0 + ...
                   1./b1*power(log(counts),1) + ...
                   1./b2*power(log(counts),2) + ...
                   1./b3*power(log(counts),3) );
            
    physical = kelvin - 273.15;
end

function [physical, units] = odas_xmp_therm_internal( input, params )
    units = '[ ^{\circ}C ]';
    ADC_FS = 4.096;
    ADC_Bits = 16;
    Gain = 1; % bridge gain
    E_B = 4.096; % Bridge excitation voltage
    a = 0; % offset
    b = 1; % slope
    Z = ((input - a)/b) * (ADC_FS / 2^ADC_Bits)*2/(Gain*E_B);
    Z(Z >  0.6) =  0.6; % prevents taking log of negative numbers
    Z(Z < -0.6) = -0.6;
    X_phys = (1-Z) ./ (1+Z); %this is the resistance ratio
    X_phys = 1/params.coef0 + (1/params.coef1)*log(X_phys);
    physical = (1 ./X_phys) - 273.15; %temperature in deg C
end

function [physical, units] = odas_xmp_shear_internal( input, params )
    units = '[ s^{-1} ]';
    ADC_FS = 4.096;
    ADC_Bits = 16;
    physical = (ADC_FS / 2^ADC_Bits) * input;
    physical = physical ./(2*sqrt(2)*params.sens*params.diff_gain);
end

function [physical, units] = odas_xmp_pres_internal( input, params )
    units = '[ dBar ]';
    a0 = -11.719*params.coef0 / params.coef1;
    a1 = 0.12208 / params.coef1;
    physical = polyval([a1 a0],input);
end

function [physical, units] = odas_xmp_pitch_internal( input, params )
    units = '[ ^{\circ} ]';
    % coef0 is the maximum reading
    % coef1 is the minimum reading
    a0 = (params.coef0 + params.coef1)/2; % offset in counts
    a1 = (params.coef0 - params.coef1)/2; % sensitivity to gravity
    physical = (180/pi)*asin((input-a0)/a1);
end

function [physical, units] = odas_xmp_sn_internal( input, params )
    units = '';
    physical = input;
end

function [physical, units] = odas_xmp_volt_internal( input, params )
    units = '[ V ]';
    physical = odas_poly_internal( input, params );
end


function [physical, units] = odas_sbt_internal( input, params )
    units = '[deg C]';
    f = params.coef6*params.coef5./double(input);
    c = log(params.coef4./f);
    physical = 1./polyval([params.coef3 params.coef2 params.coef1 ...
                           params.coef0],c) - 273.15;
end


function [physical, units] = odas_sbc_internal( input, params )
%%%% Note that this version does not correct for thermal
%%%% expansion/compressibility; this must be done separately.
    units = '[mS/cm]';
    f = params.coef6*params.coef5./double(input)/1000;
    physical = polyval([params.coef4 params.coef3 params.coef2 ...
                        params.coef1 params.coef0],f);
end


function [physical, units] = odas_o2_43f_internal( input, params )
    units = '[ ? ]';
    f = params.f_clk*params.M ./ double(input); % in Hz
    physical = params.Soc * (f - params.Foffset); % relative to saturation [mL/L]
    % Still needs temperature and pressure correction
end


function [physical, units] = odas_poly_internal( input, params )
    units = ' ';
    polyVals = [];
    for x = 0:9
        tmpName = sprintf('coef%d', x);
        if isfield( params, tmpName )
            polyVals = [params.(tmpName) polyVals];
        else
            break;
        end
    end
    physical = polyval(polyVals,input);
end


function [physical, units] = odas_gnd_internal( input, params )
    units = '[ counts ]';
    physical = input;
end


function [physical, units] = odas_accel_internal( input, params )
    units = '[ m sec^{-2} ]';
    g = 9.81;               % Gravitational constant

    try adc_zero = params.adc_zero; catch, adc_zero = 0; end
    try adc_fs   = params.adc_fs;   catch, adc_fs   = 1; end
    try adc_bits = params.adc_bits; catch, adc_bits = 0; end
    try sig_zero = params.sig_zero; catch, sig_zero = 0; end
    
    input = input*adc_fs/2^adc_bits + adc_zero - sig_zero;
    physical = g*(input-params.coef0)/params.coef1;
end

function [physical, units] = odas_piezo_internal( input, params )
    units = '[ counts ]';

    try a_0 = params.a_0; catch, a_0 = 0; end

    physical = input - a_0;
end

function [physical, units] = odas_magn_internal( input, params )
    units = '[ micro-Tesla ]';
    physical = (input-params.coef0)/params.coef1;
end


function [physical, units] = odas_Alec_EMC_internal( input, params )
    units = '[ m / s ]';
    physical = polyval([params.coef1 params.coef0],input);
end


function [physical, units] = odas_altimeter_internal( input, params )
    units = '[ m ]';
    physical = polyval([params.coef1 params.coef0],input);
end


function [physical, units] = odas_inclxy_internal( input, params )
    units = '[ deg ]';
    [raw, old_flag, error_flag] = adis(input);
    physical = polyval([params.coef1 params.coef0],raw);
end


function [physical, units] = odas_inclt_internal( input, params )
    units = '[ deg C ]';
    % Set default coefficients for ADIS temperature that will never change.
    try a = params.coef0; catch, a = 624;   end
    try b = params.coef1; catch, b = -0.47; end
    [raw, old_flag, error_flag] = adis(input);
    physical = polyval([b a],raw);
end


function [physical, units] = odas_altimeter_error_internal( input, params )
    units = ' ';
    physical = polyval([params.coef1 params.coef0],input);
end


function [physical, units] = odas_raw_internal( input, params )
    units = '[ counts ]';
    physical = input;
end


function [physical, units] = odas_voltage_internal( input, params )
    units = '[ V ]';
    try zero = params.adc_zero; catch, zero = 0; end
    try gain = params.g;        catch, gain = 1; end
    physical = (zero + input * params.adc_fs/2^params.adc_bits) / gain;
end


function [physical, units] = odas_therm_internal( input, params )
    units = '[ deg C ]';
    Z = ((input - params.a)/params.b) * (params.adc_fs/2^params.adc_bits)*2/(params.g*params.e_b);
    Z(Z >  0.6) =  0.6; % prevents taking log of negative numbers
    Z(Z < -0.6) = -0.6;
    physical = (1-Z) ./ (1+Z); %this is the resistance ratio
    Log_R = log(physical);
    if isfield(params, 'beta')
        physical = 1/params.t_0 + (1/params.beta)*Log_R;
    elseif isfield(params, 'beta_1')
        physical = 1/params.t_0 + (1/params.beta_1)*Log_R;
    else
        warning('No beta or beta_1 parameter for this thermistor')
    end
    if isfield(params, 'beta_2')
        physical = physical + (1/params.beta_2)*Log_R.^2;
        if isfield(params, 'beta_3')
            physical = physical + (1/params.beta_3)*Log_R.^3;
        end
    end
    physical = 1 ./physical - 273.15; %temperature in deg C
end

function [physical, units] = odas_t_ms_internal( input, params )
    % For MicroSquid Thermometer
    units = '[ deg C ]';
    try zero = params.adc_zero; catch, zero = 0; end
    Z = input * (params.adc_fs/2^params.adc_bits) + zero; % Turn input into a voltage
    Z = ((Z - params.a)/params.b) *2 / (params.g*params.e_b);
    Z(Z >  0.6) =  0.6; % prevents taking log of negative numbers
    Z(Z < -0.6) = -0.6;
    physical = (1-Z) ./ (1+Z); %this is the resistance ratio
    Log_R = log(physical);
    if isfield(params, 'beta')
        physical = 1/params.t_0 + (1/params.beta)*Log_R;
    elseif isfield(params, 'beta_1')
        physical = 1/params.t_0 + (1/params.beta_1)*Log_R;
    else
        warning('No beta or beta_1 parameter for this thermistor')
    end
    if isfield(params, 'beta_2')
        physical = physical + (1/params.beta_2)*Log_R.^2;
        if isfield(params, 'beta_3')
            physical = physical + (1/params.beta_3)*Log_R.^3;
        end
    end
    physical = 1 ./physical - 273.15; %temperature in deg C
end


function [physical, units] = odas_ucond_internal( input, params )
    units = '[ mS / cm ]';
    try adc_zero = params.adc_zero; catch, adc_zero = 0; end
    try adc_fs   = params.adc_fs;   catch, adc_fs   = 1; end
    try adc_bits = params.adc_bits; catch, adc_bits = 0; end
    physical = (adc_fs / 2^adc_bits) * input + adc_zero; %turn into voltage at uC board
    physical = (physical - params.a) ./ params.b; %it is now in units of conductance [S]
    physical = physical ./ params.k; %it is now in units of conductivity [S/m]
    physical = 10*physical; %and now in [mS/cm]
end


function [physical, units] = odas_c_ms_internal( input, params )
    units = '[ mS / cm ]';
    physical = (params.adc_fs / 2^params.adc_bits) * input + params.adc_zero; %turn into voltage at MS-C board
    physical = (physical - params.a) ./ params.b; %it is now in units of conductance [S]
    physical = physical ./ params.k; %it is now in units of conductivity [S/m]
    physical = 10*physical; %and now in [mS/cm]
end


function [physical, units] = odas_shear_internal( input, params )
    units = '[ s^{-1} ]';
    try adc_zero = params.adc_zero; catch, adc_zero = 0; end
    try sig_zero = params.sig_zero; catch, sig_zero = 0; end
    physical = (params.adc_fs / 2^params.adc_bits) * input + (adc_zero - sig_zero);
    physical = physical ./(2*sqrt(2)*params.diff_gain*params.sens);
end


function [physical, units] = odas_vector_internal( input, params )
    units = '[ m / s ]';
    % First convert signal to a voltage
    try bias = params.bias; catch, bias = 0; end
    V = params.adc_zero + input * (params.adc_fs / 2^params.adc_bits);
    % Next convert voltage to velocity, in m / s.
    % Sens is in m/s per volt. Offset, the zero-velocity voltage, is in volts.
    V = (V - params.offset) * params.sens;
    physical = V - bias;
end


function [physical, units] = odas_jac_emc_internal( input, params )
    units = '[ m / s ]';
    % First convert signal to a voltage
    % zero = Voltage of digitizer virtual ground wrt the input signal
    %        ground.
    try bias = params.bias; catch, bias = 0; end
    V = params.adc_zero + input * (params.adc_fs / 2^params.adc_bits);
    % Next convert voltage to velocity, in m / s.
    % Remember to divide coefficents by 100 for cm/s to m/s.
    a = params.a / 100;
    b = params.b / 100;
    V = a + b * V;
    physical = V - bias;
    warning(['The type "jac_emc" is depreciated.  Please change your ' ...
             'configuration file to reference type "aem1g_a" in place ' ...
             'of "jac_emc".']);
end


function [physical, units] = odas_aem1g_a_internal( input, params )
    units = '[ m / s ]';
    % First convert signal to a voltage
    % zero = Voltage of digitizer virtual ground wrt the input signal
    %        ground.
    try bias = params.bias;     catch, bias = 0;                 end
    try zero = params.adc_zero; catch, zero = params.adc_fs / 2; end
    V = zero + input * (params.adc_fs / 2^params.adc_bits);    
    % Next convert voltage to velocity, in m / s.
    % Remember to divide coefficents by 100 for cm/s to m/s.
    a = params.a / 100;
    b = params.b / 100;
    V = a + b * V;
    physical = V - bias;
    if b < 1
        warning(['Provided A, B coefficients are not correct.  Each ' ...
                 'instrument provides two sets of calibration results ' ...
                 '- one for analog output and one for digital output. ' ...
                 'The analog values should be used for this channel type.']);
    end
end


function [physical, units] = odas_aem1g_d_internal( input, params )
    units = '[ m / s ]';
    % Data originates as 16 bit unsigned values but, by default, is read as
    % 16 bit signed values then converted to doubles.  Now we convert any
    % negative values back into unsigned values.
    input(input < 0) = input(input < 0) + 2^16;
    % Next convert voltage to velocity, in m / s.
    % Remember to divide coefficents by 100 for cm/s to m/s.
    a = params.a / 100;
    b = params.b / 100;
    physical = a + b * input;
    if b > 1
        warning(['Provided A, B coefficients are not correct.  Each ' ...
                 'instrument provides two sets of calibration results ' ...
                 '- one for analog output and one for digital output. ' ...
                 'The digital values should be used for this channel type.']);
    end
end


function [physical, units] = odas_aroft_o2_internal( input, params )
    % Type used for RINKO FT fast optical DO sensors with RS232 output.
    units = '[ \mumol L^{-1} ]';
    % Data originates as 16 bit unsigned values but, by default, is read as
    % 16 bit signed values then converted to doubles.  Now we convert any
    % negative values back into unsigned values.
    input(input < 0) = input(input < 0) + 2^16;
    
    % Temperary fix of a buggy RS232 receiver.  Filter obvious glitches.
    for i=2:length(input)
        if input(i) <= 50 || input(i) >= 2^16 - 50
            input(i) = input(i-1);
        end
    end
    
    physical = input/100;
end

function [physical, units] = odas_aroft_t_internal( input, params )
    % Type used for RINKO FT fast optical DO sensors with RS232 output.
    units = '[ deg C ]';
    % Data originates as 16 bit unsigned values but, by default, is read as
    % 16 bit signed values then converted to doubles.  Now we convert any
    % negative values back into unsigned values.
    input(input < 0) = input(input < 0) + 2^16;

    % Temperary fix of a buggy RS232 receiver.  Filter obvious glitches.
    for i=2:length(input)
        if input(i) <= 50 || input(i) >= 2^16 - 50
            input(i) = input(i-1);
        end
    end
    
    physical = input/1000 - 5;
end


function [physical, units] = odas_amt_o2_internal( input, params )
    % For MicroSquid AMT Dissolved oxygen sensor
    units = '[  ]';

    % process the oxygen signal
    % First convert data into volts
    % default is for a signal that is already in volts.
    try adc_zero = params.adc_zero; catch, adc_zero = 0; end
    try adc_fs   = params.adc_fs;   catch, adc_fs   = 1; end
    try adc_bits = params.adc_bits; catch, adc_bits = 0; end
    try sig_zero = params.sig_zero; catch, sig_zero = 0; end

    physical = input*adc_fs/2^adc_bits + adc_zero - sig_zero; % now in volts

    % Now get the user calibration voltage readings for 0% and 100% oxygen
    try U0   = params.u0;  catch,  U0   = 0.125; end % nominal no O2 reading
    try U100 = params.u100; catch, U100 = 2.5  ; end % typical 100% O2 reading

    % read air pressure (mbar) during calibration
    try PL = params.pl; catch, PL = 1013; end % nominal air pressure

    % read temperature (celsius) during calibration.
    try T_cal   = params.t_cal; catch, T_cal = 22; end % nominal lab temperature

    % Get the water vapour pressure at time of calibration
    Pw_cal = get_wvP(T_cal + 273.15);% Water vapor pressure during calibration

    % calculate the temperature adjustment factor for the calibration
    % temperature
    b    = zeros(1,4);
    b(1) = params.b0;
    b(2) = params.b1;
    b(3) = params.b2;
    b(4) = params.b3;
    b = fliplr(b);

    ET_cal = polyval(b, T_cal);% temperature correction for this sensor

    X_O2 = 0.2095; % The mole fraction of O2 in air

    % calcuate sensitivity coefficient
    a20 = (X_O2 / ET_cal) * (PL - Pw_cal) / (U100 - U0);% sensitivity coefficient

    physical = physical - U0;    % subtract the zero-oxygen reading from all data
    physical = a20 .* physical;  % multiply by sensitivity coefficient
    physical = physical ./ X_O2; % divide by molar fraction of O2 in air.
    % This does not make much sense but is the AMT procedure.

    % The physical value calculated here still needs correction. The correction
    % factors (.*) and divisors (./) are listed below.
    %   (1) In situ temperature, .* ET.
    %   (2) In situ water vapour pressure, ./ (PN - Pw), where PN=1013
    %       and Pw is the water vapour pressure at the in situ temperature.
    %   (3) In situ saturation concentration, .* C_S, where
    %       C_S = 1.4289*C_T .* C_S.
    %       C_T depends only on temperature, and
    %       C_S depends on salinity and temperature.
    %   (4) In situ depth, .* C_D = exp(-0.3775*Depth ./ T_K), where depth is
    %       in meteres and T_K is in kelvin.

end

function value = get_wvP(T)
% 
% value = get_wvP(T)
% function to calculate the water vapour pressure, in mbar. Temperature, T,
% in % kelvin.
% In support of AMT disolved O2 sensor.
%
% 2015-05-07, RGL, orginal version.

P = [-2.16961e5 -3.8407e3 11.8571];
value = polyval(P, 1 ./ T);
value = 1e3*exp(value);
end

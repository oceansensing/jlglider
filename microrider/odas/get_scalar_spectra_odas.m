%% get_scalar_spectra_odas
% Calculate spectra of scalar signals for an entire profile.
%%
% <latex>\index{Functions!get\_scalar\_spectra\_odas}</latex>
%
%%% Syntax
%   sp = get_scalar_spectra_odas( scalar_vectors, ref_vector, P, t, speed, ... )
%
% * [scalar_vectors] Matrix in which each column is a scalar vector.
% * [ref_vector] Matrix of column vectors of reference signals that are
%       contaminating the scalar signals. Usually this is a measure of the
%       magnetic field induced by a near-by electromagnetic current meter.
%       If present AND if the optional locical value of "goodman" is true,
%       then this function will call clean_shear_spec.m to remove the
%       contamination of the scalar signals. Can be empty, in which case
%       there is no decontamination regardless of the value of "goodman". 
% * [P] Column vector of pressure that will be used to calculate the average
%       pressure of the data used for each spectrum. Its length must match
%       the number of rows in scalar_vectors.
% * [t] Column vector of time that will be used to calculate the average
%       time of the data used for each spectrum. Its length must match
%       the number of rows in scalar_vectors.
% * [speed] Scalar or vector of profiling speed that is used to derive
%       wavenumber spectra. If it is a vector, it must have a length equal to
%       the number of rows in scalar_vectors.
% * [...] Structure, or parameter / value pair, of optional parameters.
% * []
% * [sp] Structure containing resulting spectra and anciliary information.
%       Exact contents of the structure are discussed within the
%       description section. 
%
%%% Description
% Computes the spectra of scalar signals such as the gradients of
% temperature and micro-conductivity. It is usually used in conjunction
% with $\texttt{get\_diss\_odas}$, which returns the spectra of shear
% signals. 
%
%
%%% Optional Input Parameters
%
% * [diff_gain] Vector of differentiator gains used to pre-emphasize the
%      signals. It is used to make a small correction to the spectra due to
%      the frequency characteristics of the deconvolve function. Must be
%      empty or its length must equal the number of columns in
%      scalar_vectors, if gradient_method = 'first_difference'. 
%      Default = [ ].
% * [gadient_method] Method used to create the scalar gradient vectors. 
%      'first_difference' corrects for the deconvolution and the first
%      difference operator used to estimate the scalar gradient.
%      'high_pass' applies no correction. Default = 'first_difference'
% * [fft_length] FFT segment length [samples] used to estimate the
%      auto-spectrum of scalar_vectors. Default fft_length = 512.
% * [spec_length] Length of data used for each spectrum. It should be
%      at least 2 times fft_length. The recommended value is larger than 2 times
%      fft_length, for statistical reliability. Default = 3*fft_length.
% * [overlap] Number of samples by which the auto-spectral estimates
%      overlap. Should be between 0 and spec_length. The recommended value
%      is overlap = spec_length/2, or the value used for the function
%      get_diss_odas. 
% * [fs] Sampling rate in Hz. Default value = 512.
% * [f_AA] Cut-off frequency of the anti-aliasing filter. Default = 98 Hz.
% * [goodman] Logical variable to determine if the function
%      clean_shear_spec should be called to remove coherent signal
%      contamination of the scalar signals using the signals in ref_vector.
%
%%% Resulting Output Structure
%
% * [K] Wavenumbers [cpm, or cycles per metre]. Contains one column for
%      every spectral estimate.
% * [F] Frequencies with one column for every spectral estimate. Every
%      column is identical. 
% * [scalar_spec] Wavenumber spectrum for each scalar vector.
%      This is a 3 dimensional matrix of size [M N L]. Every column is a
%      spectrum -- one for each signal. Thus, N equals the number of
%      signals. Each layer (the third index, L) corresponds to a 
%      spectral estimate. The number of layers, L, is determined by the
%      values of spec_length, and overlap, and the length of the signals.
%      The number of rows, M, is the number of wavenumber values (from 0
%      to the Nyquist wavenumber). Usually, M = 1 + fft_length/2. M equals
%      the number of rows in K. 
% * [speed] Column vector of the mean speed for each spectral estimate.
%      Derived directly from the input.
% * [P] Column vector of the mean pressure for each spectral estimate.
%      Derived directly from the input.

% *Version History:*
% 2013-12-13 (RGL) Original version, based on get_diss_odas.
% 2013-12-18 (RGL) Its ready for prime-time. William, please add the feature to
%       make this function return the default varargin if it is called
%       without any input variables.
% 2014-04-11 (RGL) Forces scalar_spec to be 3-D, even if there is only one
%       probes. Otherwise there dimensional problems if there is only one probe.
% 2014-07-25, RGL, changed "addParameter" for diff_gain to "addRequired",
%       so that this functions works with oer versions of Matlab.
% * 2015-10-29 (RGL) Documentation updates.
% * 2015-11-18 RGL Documentation updates.
% * 2016-12-13 RGL Added coherent noise removal to support cleaning the
%       spectra with the function clean_shear_spec, also known as the
%       Goodman coherent noise removal algorithm. 
% * 2018-05-09 RGL Minor annotation corrections. 

function sp = get_scalar_spectra_odas(scalar_vectors, ref_vector, P, t, speed, varargin)

% Default values for optional fields
default_fft_length  = 512;
default_spec_length = inf;
default_overlap     = inf;
default_fs          = 512;
default_gm          = 'first_difference';
default_f_AA        = 98;
default_diff_gain   = [];
default_goodman     = false;

z = inputParser;
z.CaseSensitive = true;
z.KeepUnmatched = true;

val_fft_length = @(x) isnumeric(x) && isscalar(x) && (x >= 2);
val_positive   = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
val_matrix     = @(x) isnumeric(x) && (size(x,1) > 1) && (size(x,2) >= 1);
val_speed      = @(x) isnumeric(x) && (isvector(x) || isscalar(x)) && ~any(x < 0);
val_vector     = @(x) isnumeric(x) && isvector(x);
val_string     = @(x) ischar(x);
val_diff_gain  = @(x) (isvector(x) && ~any(x < 0)) || isempty(x);
val_logical    = @(x) (islogical(x));

addRequired(  z, 'scalar_vectors', val_matrix);
addRequired(  z, 'P',              val_vector);
addRequired(  z, 't',              val_vector);
addRequired(  z, 'speed',          val_speed);
%addRequired(  z, 'diff_gain',      val_diff_gain);
addParamValue(z, 'diff_gain',      default_diff_gain,  val_diff_gain);
addParamValue(z, 'fft_length',     default_fft_length,  val_fft_length);
addParamValue(z, 'spec_length',    default_spec_length, val_fft_length);
addParamValue(z, 'overlap',        default_overlap,     val_positive);
addParamValue(z, 'fs',             default_fs,          val_positive);
addParamValue(z, 'gradient_method',default_gm,          val_string);
addParamValue(z, 'f_AA',           default_f_AA,        val_positive);
addParamValue(z, 'goodman',        default_goodman,     val_logical);

% Parse the arguments.
parse(z, scalar_vectors, P, t, speed, varargin{:});

% Perform last stages of input validation.
L = size(z.Results.scalar_vectors,1); % length of the scalar vectors
N = size(z.Results.scalar_vectors,2); % number of the scalar vectors

if z.Results.spec_length < 2*z.Results.fft_length
    error('Invalid size for spec_length - must be greater than 2 * fft_length.');
end

if (L ~= size(z.Results.t,1) || L ~= size(z.Results.P,1))
    error('Same number of rows required for scalar_vectors, t, and P.');
end

if ~isscalar(z.Results.speed) && L ~= size(z.Results.speed,1)
    error ('speed vector must have the same number of rows as scalar_vectors');
end

if z.Results.spec_length > L
    error('spec_length cannot be longer than the length of the scalar vectors');
end

if length(z.Results.diff_gain) < N ...
        && strcmpi(z.Results.gradient_method,'first_difference')
    error(['You must specify a differentiator gain ' ...
        'for each column of the scalar vector matrix, ' ...
        'when method = first_difference.']);
end
if ...
        ~(strcmpi(z.Results.gradient_method, 'first_difference') || ...
          strcmpi(z.Results.gradient_method, 'high_pass'))
    error('The only gradient methods currently supported are first_difference and high_pass');
end

scalar_vectors = z.Results.scalar_vectors;
diff_gain      = z.Results.diff_gain;
fft_length     = z.Results.fft_length;
spec_length    = z.Results.spec_length;
overlap        = z.Results.overlap;
fs             = z.Results.fs;
speed          = z.Results.speed;
P              = z.Results.P;
t              = z.Results.t;
f_AA           = z.Results.f_AA;
gradient_method = z.Results.gradient_method;
goodman         = z.Results.goodman;

if spec_length == inf, spec_length = 3*fft_length;      end
if overlap     == inf, overlap     = round(spec_length / 2);end

% Make sure that all vectors are column vectors
if isrow(t), t = t'; end
if isrow(P), P = P'; end
if isscalar(speed), speed = speed*ones(size(t)); end
if isrow(speed), speed = speed'; end
if size(scalar_vectors,1) ~= size(ref_vector,1)
    goodman = false;
end

% end of input argument checking.
select = (1:spec_length)';

% Calculate the number of spectral estimates that can be made.
if overlap >= spec_length, overlap = spec_length/2; end
if overlap < 0, overlap = 0; end
num_of_estimates = 1 + floor((size(scalar_vectors,1) - spec_length) / (spec_length-overlap));
F_length = 1 + floor(fft_length/2); % size of frequency and wavenumber vectors
num_of_vectors = size(scalar_vectors,2); % the number of vectors in the matrix

% pre-allocate matrices
% First dimension is frequency index
sp.scalar_spec  = zeros(F_length,  num_of_vectors, num_of_estimates);
sp.F            = zeros(F_length,                  num_of_estimates);
sp.K            = sp.F;
sp.speed        = zeros(num_of_estimates,1);
sp.P            = sp.speed;
sp.t            = sp.speed;

index = 1;

%%
%Main loop
while select(end) <= size(scalar_vectors,1)
    W = mean(speed(select)); % the mean speed for theis segment
    sp.speed(index)  = W;
    sp.P    (index)  = mean(P(select));
    sp.t    (index)  = mean(t(select));

    if goodman
        [P_scalar_clean, AA, P_scalar, UA, F] = clean_shear_spec(...
            ref_vector(select,:), scalar_vectors(select,:), fft_length, fs);
        if size(P_scalar_clean,2) > 1
            P_scalar_clean = permute(P_scalar_clean, [3 1 2]);
        end
        P_scalar_clean = squeeze(P_scalar_clean);
        for k = 1:num_of_vectors
            sp.scalar_spec(:,k,index) = P_scalar_clean(:, k, k) * W;
            % we return only the auto-spectra of the cleaned scalar spectra
        end
    else
        for k = 1:num_of_vectors
            [junk, F] = ...
                csd_odas(scalar_vectors(select,k),...
                [], fft_length, fs, [], fft_length/2,'linear'); % frequency spectrum
            sp.scalar_spec(:,k,index) = junk * W; % It is now a wavenumber spectrum.
        end
    end
    
    
    F = F(:);
    sp.F(:,index) = F;   % The frequency - never changes between segments.
    sp.K(:,index) = F / W; % The wavenumberas - changes between segments.

    if index == 1 && strcmpi(gradient_method, 'first_difference') % Only calculate once
        correction = ones(size(F)); % pre-fill
        correction(2:end) = ((pi*F(2:end)) ./ (fs*sin(pi*F(2:end)/fs))).^2; % Divide by zero problem
        % This is the correction for using the first difference operator to
        % estimate the gradient of these scalar vectors.
        bl_correction = ones(size(F,1),num_of_vectors); %Pre-assign
%        n = find((F > 0) & (F < 0.9*fs/2)); % index to frequency smaller than 90% of Nyquist.
        n = (1:length(F)-1)'; % all but the last value.

        for k = 1:num_of_vectors
            [b,a] = butter(1, 1/(2*pi*diff_gain(k)*fs/2)); % The LP-filter that was applied
            junk = (abs(freqz(b, a, F(n), fs))).^2; % The mag-squared of the applied LP-filter.
            H = 1 + (2*pi*F(n)*diff_gain(k)).^2; H = 1 ./ H; % What should have been applied
            bl_correction(n,k) = H ./ junk; % The bilinear transformation correction.
            bl_correction(end,k) = bl_correction(end-1,k); % is ininite at f_N
        end
    end

    if strcmpi(gradient_method, 'first_difference')
        for k = 1:num_of_vectors % apply correction to each column of spectra
        sp.scalar_spec(:,k,index) = ...
            sp.scalar_spec(:,k,index) .* correction .* bl_correction(:,k);
        end
    end

    % Note, scalar_spec is a [L N M] matrices where L is the number of frequency
    % estimates (0 to f_N), N is the number of probes (or signals), and M
    % is the number of depth (or time) estimates. The matrix is 3-D even if
    % there is only one probe.

    index = index + 1;

    select = select + spec_length - overlap;
end

sp.scalar_spec = sp.scalar_spec;
% return the parameters used to estimate the scalar spectra.
sp.fs              = fs;
sp.f_AA            = f_AA;
sp.gradient_method = gradient_method;
sp.spec_length     = spec_length;
sp.overlap         = overlap;
sp.fft_length      = fft_length;
sp.diff_gain       = diff_gain;

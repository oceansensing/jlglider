%% csd_odas
% Estimate the cross- and auto-spectrum of one or two vectors.
%%
% <latex>\index{Functions!csd\_odas}</latex>
%
%%% Syntax
%   [Cxy, F, Cxx, Cyy] = csd_odas(x, y, nFFT, rate, window, overlap, msg)
%
% * [x] Vector over which to estimate the cross-spectrum.
% * [y] Vector over which to estimate the cross-spectrum. If empty or
%       equal to x, only calculate the autospectra.  Length must match 'x'.
% * [nFFT] Length of x, y segments over which to calculate the ffts.
% * [rate] Sampling rate of the vectors.
% * [window] Window of length nFFT to place over each segment before calling 
%       the fft. If empty, window defaults to the cosine window and is 
%       normalized to have a mean-square of 1. The window is used as given and 
%       it is up to the user to make sure that the window preserves the variance
%       of the signals.
% * [overlap] Number of points to overlap between the successive segments used 
%       to calculate the fft. A good value is nFFT/2 which corresponds to 50%.
% * [msg] Message string indicating how each segment is to be detrended before 
%       calling the fft. Options are 'constant', 'linear', 'parabolic', 'cubic',
%       and 'none'. If empty, omitted, or 'none' then no detrending or mean 
%       removal is performed on the segments.
% * []
% * [Cxy] Complex cross-spectrum of x and y if x and y are different. Otherwise, 
%       the auto-spectrum of x.
% * [F] Frequency vector corresponding to the cross- and/or auto-spectrum. 
%       Ranges from 0 to the Nyquist frequency.
% * [Cxx] Auto-spectrum of x if x is different from y.
% * [Cyy] Auto-spectrum of y if x is different from y and y is not empty.
%
%%% Description
% Use this
% function only if the vectors are complex. Otherwise, use the function
% $\texttt{csd\_matrix\_odas}$. This function will become inactive and point
% to $\texttt{csd\_matrix\_odas}$ once it can handle complex signals.
%
% Estimate the cross- and auto-spectrum of one or two vectors. For best
% results nFFT should be several times shorter than the lengths of
% $\texttt{x}$ and $\texttt{y}$ so the ensemble average of the segments of
% length nFFT has some statistical significance. The cross-spectrum has the
% property that its integral from 0 to the Nyquist frequency equals the
% covariance of $\texttt{x}$ and $\texttt{y}$, if they are real variables.
% The vectors $\texttt{x}$ and $\texttt{y}$ can be complex. 
%
% This function tests if $\texttt{x}$ is identical to $\texttt{y}$ and if
% $\texttt{y}$ is empty, $\texttt{y=[ ]}$. If so, it computes only the
% auto-spectrum and places it into $\texttt{Cxy}$, and it will not return
% $\texttt{Cxx}$ and $\texttt{Cyy}$, or it returns them empty. It behaves
% likewise if $\texttt{y}$ equals the empty matrix. Thus, there is no need
% for a separate function to calculate auto-spectra. In addition, if
% $\texttt{x}$ and $\texttt{y}$ are not the same vector, then the two
% auto-spectra are optionally available and can be used to calculate the
% coherency spectrum and the transfer function of $\texttt{y}$ relative to
% $\texttt{x}$. This function $\texttt{csd\_odas}$ is very similar to the
% Matlab function $\texttt{csd}$ that is no longer supported by MathWorks.
% We do not like the new approach used by MathWorks for spectral estimation
% because their new methods do not support the linear detrending of each
% segment before calling the fft. It only allows the detrending of the
% entire vectors $\texttt{x}$ and $\texttt{y}$, which can leave undesirable
% low-frequency artifacts in the cross- and auto-spectrum. In addition, you
% do not need the Signal Processing Toolbox to use this function.  The
% equation for the squared-coherency spectrum is:
%
% $$G_{xy}(f)=\frac{|C_{xy}(f)|^2}{C_{xx}(f)C_{yy}(f)}$$
%
% where $G_{xy}$ is the squared coherency-spectrum, $C_{xy}$ is the complex 
% cross-spectrum, and $C_{xx}$ and $C_{yy}$ are the autospectra.  The complex 
% transfer function, $H_{xy}$, of $y$ relative to $x$, is:
%
% $$H_{xy}(f)=\frac{C_{xy}(f)}{C_{xx}(f)}$$
%

% *Version History:*
%
% * 2009-04-01 (RGL) Original function
% * 2011-08-23 (RGL) Modified to accept the case of y = [ ] for implicit
%   auto-spectrum. Forced Cxx and Cyy to [], for auto-spectrum
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-07-04 (RGL) added parabolic and cubic detrending feature
% * 2015-11-18 (RGL) Divided zero and Nyquist frequency components by 2.
%      Update the documentation.

function [Cxy, F, Cxx, Cyy] = csd_odas(x, y, n_fft, rate, Window, over_lap, msg)

% Here we check the input arguments
if (nargin < 6)
    error('Input arguments must be 6 or 7')
end
if (nargin == 6) ; %force detrend to 'none'
    msg = 'none';
end
if (nargin == 7); % then check the detrend message
    if ~(...
            strcmp(msg(1:3),'con') || strcmp(msg(1:3),'lin') ...
         || strcmp(msg(1:3),'par') || strcmp(msg(1:3),'cub') ...
         || strcmp(msg(1:3),'non'))
        error('detrending string must one of "none", or "constant", or "linear"')
    else
        if (strcmp(msg(1:3), 'con')); msg = 'constant';  order = 0; end;
        if (strcmp(msg(1:3), 'lin')); msg = 'linear'  ;  order = 1; end;
        if (strcmp(msg(1:3), 'non')); msg = 'none'    ;             end;
        if (strcmp(msg(1:3), 'par')); msg = 'parabolic'; order = 2; end;
        if (strcmp(msg(1:3), 'cub')); msg = 'cubic'    ; order = 3; end;
    end
end
% Check x and y
auto = 0;
if isempty(y), y=x; auto = 1; end

if (isvector(x) && isvector(y))
    if (size(x,2) ~= 1), x = x'; end % force vectors to be column vectors
    if (size(y,2) ~= 1), y = y'; end % force vectors to be column vectors
else
    error ('Only vectors are acceptable, no matrices, please')
end
if (size(x,1) ~= size(y,1)), error ('Input vectors must be of equal size'), end


if ~auto,
    if x==y, 
        auto = 1; 
    end % Only calculate the auto-spectrum
end

% Check n_fft
if ~isscalar(n_fft), ...
        error('n_fft must be a scalar and positive whole number'), end
if (n_fft < 0), error ('n_fft must be a positve whole number'), end
n_fft = floor(n_fft); % force it to be a whole number

if (size(x,1) < 2*n_fft), ...
        error('The length of the input vectors should be more than twice n_fft'), end

% Check rate
if ~isscalar(rate), ...
        error('The sampling rate must be a scalar and positive number'), end
if (rate < 0), error ('The sampling rate must be a positive number'), end

% Check Window
if (isempty(Window)) % Generate cosine window
    Window = 1 + cos(pi*(-1 + 2*(0:n_fft-1)'/n_fft));
    Window = Window ./ (sqrt(mean(Window.^2))); ...
        % Force it to have a mean-square of 1.
end
if isvector(Window)
    if (size(Window, 2) ~= 1), Window = Window'; end % force into column vector
else
    error('Window must be a vector of length n_fft')
end
if (size(Window,1) ~= n_fft), error('Window length must equal n_fft'), end

% Check over_lap
if ~isscalar(over_lap), ...
        error('over_lap must be a scalar and positive whole number'), end
if ((over_lap < 0) || (over_lap > (n_fft-1))), ...
        error ('over_lap must be a positve whole number smaller than n_fft'), end
over_lap = floor (over_lap); % force it to be a whole number

% OK, enough already. Get on with the work!

num_of_segments = floor((size(x,1) - over_lap) / (n_fft - over_lap));
% This is the total number of segments that will be processed to produce a
% spectrum

increment = n_fft - over_lap; % the shift from one segment to the next
select = 1:n_fft; % the range of indices that are selected for a segment

Cxx = zeros(n_fft,1); Cyy = Cxx; Cxy = Cxx + 1j*Cyy; % pre-assign vectors
ramp = (0:length(select)-1)';

if auto
    for i=1:num_of_segments
        if strcmp(msg,'none')
            xw = Window.*x(select);
        else
            p = polyfit(ramp,x(select),order);
            xx = x(select) - polyval(p,ramp);
            xw = Window.*xx;
        end
        select = select + increment;
        X = fft(xw,n_fft);
        Cxy = Cxy + abs(X).^2;
    end
else
    for i=1:num_of_segments
        if strcmp(msg,'none')
            xw = Window.*x(select);
            yw = Window.*y(select);
        else
            p = polyfit(ramp,x(select),order);
            xx = x(select) - polyval(p,ramp);
            xw = Window.*xx;
            p = polyfit(ramp,y(select),order);
            yy = y(select) - polyval(p,ramp);
            yw = Window.*yy;
        end
        select = select + increment;
        X = fft(xw,n_fft);
        Y = fft(yw,n_fft);
        XY = Y.*conj(X);
        Cxy = Cxy + XY;
        if (nargout > 2)
            Cxx = Cxx + abs(X).^2;
            Cyy = Cyy + abs(Y).^2;
        end
    end
end

% Select the frequency range.
info_x = whos('x');
info_y = whos('y');

if info_x.complex || info_y.complex
    % Either x or y is complex.
    thems_real = 0;
    range = 1:n_fft;
else
    % Both x and y are real.
    thems_real = 1;

    range = 1:ceil((n_fft+1)/2);   % includes DC AND Nyquist frequency
    
    Cxy = Cxy(range);
    Cxx = Cxx(range);
    Cyy = Cyy(range);
end

F = (range - 1)'*rate/n_fft;

% Now normalize the spectral variance so that:
% (1)the integral from 0 to the F_N equals the variance of the signal, and
% (2) the integral from -F_N to F_N - one point equals the variance for
% complex x and y.
%
% And, we note that the DC and Nyquist frequency (F_N)
% components appear only once.

Cxy = Cxy / num_of_segments; % find ensemble mean
if (nargout > 2)
    Cxx = Cxx / num_of_segments;
    Cyy = Cyy / num_of_segments;
end

if (thems_real)
    Cxy = Cxy / (n_fft*rate/2);
     Cxy(1)   = Cxy(1) / 2; % Uncomment if you plan to use the cumsum
     Cxy(end) = Cxy(end) / 2; % function to integrate the spectrum.
    if (nargout > 2)
        Cxx = Cxx / (n_fft*rate/2);
        Cyy = Cyy / (n_fft*rate/2);
         Cxx(1)   = Cxx(1) / 2;   
         Cyy(end) = Cyy(end) / 2; % function to integrate the spectrum.
    end
else % They are complex
    Cxy = Cxy / (n_fft*rate);
    if (nargout > 2)
        Cxx = Cxx / (n_fft*rate);
        Cyy = Cyy / (n_fft*rate);
    end
    
    n = find (F >= rate/2);
    F(n) = F(n) - rate;
    F = fftshift(F);
    Cxy = fftshift(Cxy);
    if (nargout > 2)
        Cxx = fftshift(Cxx);
        Cyy = fftshift(Cyy);
    end       
end





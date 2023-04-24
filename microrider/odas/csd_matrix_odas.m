%% csd_matrix_odas
% Estimate the cross- and auto-spectrum of one or two real matrices
%%
% <latex>\index{Functions!csd\_matrix\_odas}</latex>
%
%%% Syntax
%   [Cxy, F, Cxx, Cyy] = csd_matrix_odas(x, y, nFFT, rate, window, overlap, msg)
%
% * [x] Matrix over which to estimate the cross-spectrum.
% * [y] Matrix over which to estimate the cross-spectrum. If empty,
%       only calculate the autospectra of x.  Length must match 'x'.
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
% * [Cxy] Complex cross-spectrum matrix of x and y if x and y are different. Otherwise, 
%       the auto-spectrum matrix of x.
% * [F] Frequency vector corresponding to the cross- and/or auto-spectrum. 
%       Ranges from 0 to the Nyquist frequency.
% * [Cxx] Auto-spectrum matrix of x if x is different from y.
% * [Cyy] Auto-spectrum matrix of y if x is different from y and y is not empty.
%
%%% Description
% Estimate the cross- and auto-spectrum of one or two matrices, where every
% column is a real vector, and the number of rows is the length of the
% vectors. For best results nFFT should be several times shorter than the
% lengths of $\texttt{x}$ and $\texttt{y}$ so that the ensemble average of
% the segments of length nFFT have some statistical significance. The
% cross-spectrum has the property that its integral from 0 to the Nyquist
% frequency equals the covariance of $\texttt{x}$ and $\texttt{y}$, if they
% are real variables. Use the function $\texttt{csd\_odas}$ if either
% $\texttt{x}$ or $\texttt{y}$, or both, are complex vectors. 
%
% This function tests if $\texttt{y}$ is empty, $\texttt{y=[ ]}$. If 
% so, it computes only the auto-spectrum and places it into $\texttt{Cxy}$,
% and it will not return $\texttt{Cxx}$ and $\texttt{Cyy}$, or it returns
% them empty. Thus, there is no need for a separate function to calculate
% auto-spectra. In addition, if $\texttt{y}$ is not an empty matrix, then the two
% auto-spectra are optionally available and can be used to calculate the
% coherency spectrum and the transfer function of $\texttt{y}$ relative to
% $\texttt{x}$. This function $\texttt{csd\_matrix\_odas}$ is similar to the
% RSI function $\texttt{csd\_odas}$, when $\texttt{x}$ and $\texttt{y}$ are
% columns vectors. 
% 
% This function was motivated by the need to make the Goodman
% coeherent-noise reduction function, $\texttt{clean\_shear\_spec}$, more
% efficient. If you use $\texttt{csd\_odas}$ to build the auto- and
% cross-spectral matrices, then there are many redundant calls to the
% fft-function, and this redundancy increases with the number of columns in
% the input matracies. The function $\texttt{csd\_matrix\_odas}$ does not
% make redundant fft calls.

% *Version History:*
%
% * 2014-04-10 (RGL) Original function
% * 2015-10-28 (RGL) Documentation corrections.
% * 2015-11-18 (RGL) Documentation corrections.

function [Cxy, F, Cxx, Cyy] = csd_matrix_odas(x, y, n_fft, rate, Window, over_lap, msg)

% Here we check the input arguments
if (nargin < 6)
    error('Input arguments must be 6 or 7')
end
if (nargin == 6) ; %force detrend to 'none'
    msg = 'none';
end
if (nargin == 7); % then check the detrend message
    if ~(...
            strcmpi(msg(1:3),'con') || strcmpi(msg(1:3),'lin') ...
         || strcmpi(msg(1:3),'par') || strcmpi(msg(1:3),'cub') ...
         || strcmpi(msg(1:3),'non'))
        error('detrending string must one of "none", "constant", "linear", "parabolic", or "cubic" ')
    else
        if (strcmpi(msg(1:3), 'con')); msg = 'constant';  order = 0; end;
        if (strcmpi(msg(1:3), 'lin')); msg = 'linear'  ;  order = 1; end;
        if (strcmpi(msg(1:3), 'non')); msg = 'none'    ;             end;
        if (strcmpi(msg(1:3), 'par')); msg = 'parabolic'; order = 2; end;
        if (strcmpi(msg(1:3), 'cub')); msg = 'cubic'    ; order = 3; end;
    end
end
% Check x and y
auto = false; % auto = true means only calculate the spectral matrix for x
if isempty(y), y = x; auto = true; end
if isscalar(x) || (~isempty(y) && isscalar(y))
    error('x and y must be vectors or matrices')
end

if isrow(x), x = x'; end % force vectors to be column vectors
if ~isempty(y) && isrow(y), y = y'; end
if (size(x,1) ~= size(y,1))
    error ('number of rows of input matrices must be equal')
end

% Check n_fft
if ~isscalar(n_fft) || n_fft<=2, ...
        error('n_fft must be a positive scalar > 2 .'), end
n_fft = floor(n_fft); % force it to be a whole number

if (size(x,1) < 2*n_fft), ...
        error('The length of the input vectors should be twice n_fft or larger.'), end

% Check rate
if ~isscalar(rate) || rate<0
        error('The sampling rate must be a positive scalar.')
end

% Check Window
if isempty(Window) % Generate cosine window
    Window = 1 + cos(pi*(-1 + 2*(0:n_fft-1)'/n_fft));
    Window = Window ./ (sqrt(mean(Window.^2))); % Force it to have a mean-square of 1.
end
if isrow(Window), Window = Window'; end % force into column vector
if ~isvector(Window) || (size(Window,1) ~= n_fft)
    error('Window must be empty or a vector of length n_fft')
end

% Check over_lap
if ~isscalar(over_lap)
        error('over_lap must be a scalar and positive whole number')
end
if over_lap < 0
    over_lap = 0;
elseif over_lap > n_fft*0.9 % That is too much overlap, you silly.
    over_lap = round(n_fft/2);
end

over_lap = floor (over_lap); % force it to be a whole number

% OK, enough already. Get on with the work!

num_of_segments = floor((size(x,1) - over_lap) / (n_fft - over_lap));
% This is the total number of segments that will be processed to produce
% spectra.

increment = n_fft - over_lap; % the shift from one segment to the next
select = 1:n_fft; % the range of indices that are selected for a segment

NN = 1 + floor(n_fft/2); % The number of non-negative frequency indices

columns_in_x = size(x,2);
columns_in_y = max(size(y,2),1);

% pre-assign vectors
Cxx = complex(zeros(NN,columns_in_x,columns_in_x));
Cyy = complex(zeros(NN,columns_in_y,columns_in_y));
Cxy = complex(zeros(NN,columns_in_x,columns_in_y));
ramp = (0:length(select)-1)';

if auto
    for i=1:num_of_segments
        xx = x(select,:); %extract the selected range 
        if strcmp(msg,'none')
            for index = 1:columns_in_x
                xx(:,index) = Window.*xx(:,index);% apply window to each column
            end
        else
            for index = 1:columns_in_x
                p = polyfit(ramp,xx(:,index),order);% detrend each column
                xx(:,index) = xx(:,index) - polyval(p,ramp);
                xx(:,index) = Window.*xx(:,index);% apply window
            end
        end
        junk  = fft(xx);
        fft_x = junk(1:NN,:); % extract non-negative frequency values
        fft_x_conj = conj(fft_x);
        for page = 1:columns_in_x % construct 3-D x spectral matrix
            for index = 1:columns_in_x
                Cxy(:,index,page) = ...
                    Cxy(:,index,page) + fft_x(:,page) .* fft_x_conj(:,index);                
            end
        end
        select = select + increment;
    end
else
    for i=1:num_of_segments
        xx = x(select,:); % extract selected range
        yy = y(select,:);
        if strcmp(msg,'none')
            for index = 1:columns_in_x
                xx(:,index) = Window.*xx(:,index);% apply window to each column
            end
            for index = 1:columns_in_y
                yy(:,index) = Window.*yy(:,index);% apply window to each column
            end
        else
            for index = 1:columns_in_x
                p = polyfit(ramp,xx(:,index),order);% detrend each column
                xx(:,index) = xx(:,index) - polyval(p,ramp);
                xx(:,index) = Window.*xx(:,index);% apply window
            end
            for index = 1:columns_in_y
                p = polyfit(ramp,yy(:,index),order);% detrend each column
                yy(:,index) = yy(:,index) - polyval(p,ramp);
                yy(:,index) = Window.*yy(:,index);% apply window
            end
        end
        junk  = fft(xx);
        fft_x = junk(1:NN,:);
        fft_x_conj = conj(fft_x);
        junk  = fft(yy);
        fft_y = junk(1:NN,:);
        fft_y_conj = conj(fft_y);

        for page = 1:columns_in_x % construct 3-D spectral matrix for x
            for index = 1:columns_in_x
                Cxx(:,index,page) = ...
                    Cxx(:,index,page) + fft_x(:,page) .* fft_x_conj(:,index);                
            end
        end
        for page = 1:columns_in_y % construct 3-D spectral matrix for y
            for index = 1:columns_in_y
                Cyy(:,index,page) = ...
                    Cyy(:,index,page) + fft_y(:,page) .* fft_y_conj(:,index);                
            end
        end
        for page = 1:columns_in_y % construct 3-D spectral matrix for xy
            for index = 1:columns_in_x
                Cxy(:,index,page) = ...
                    Cxy(:,index,page) + fft_x_conj(:,index) .* fft_y(:,page);                
            end
        end

        select = select + increment;
    end
end

F = (0 : NN - 1)'*rate/n_fft;

% Now normalize the spectral variance so that:
% (1) the integral from 0 to the F_N equals the variance of the signal, 
%   and, we note that the DC and Nyquist frequency (F_N)
%   components appear only once in an fft output.

Cxy = Cxy / num_of_segments; % find ensemble mean
if (nargout > 2)
    Cxx = Cxx / num_of_segments;
    Cyy = Cyy / num_of_segments;
end

Cxy = Cxy / (n_fft*rate/2);
Cxy(1,:,:)   = Cxy(1,:,:)   / 2; 
Cxy(end,:,:) = Cxy(end,:,:) / 2; 
if (nargout > 2)
    Cxx = Cxx / (n_fft*rate/2);
    Cyy = Cyy / (n_fft*rate/2);
    Cxx(1,:,:)   = Cxx(1,:,:)   / 2;
    Cyy(end,:,:) = Cyy(end,:,:) / 2; 
end
if (nargout == 1),     F = []; Cxx = []; Cyy = []; end
if (nargout == 2),             Cxx = []; Cyy = []; end
if (nargout == 3),                       Cyy = []; end




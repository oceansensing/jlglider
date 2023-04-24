%% clean_shear_spec
% Remove acceleration contamination from all shear channels.
%%
% <latex>\index{Functions!clean\_shear\_spec}</latex>
%
%%% Syntax
%   [cleanUU, AA, UU, UA, F] = clean_shear_spec( A, U, nFFT, rate )
%
% * [A] Matrix of acceleration signals usually represented by [Ax Ay Az] 
%       where Ax Ay Az are column vectors of acceleration components.
% * [U] Matrix of shear probe signals with each column representing a 
%       single shear probe.
% * [nFFT] Length of the fft used for the calculation of auto- and cross-spectra.
% * [rate] Sampling rate of the data in Hz
% * []
% * [cleanUU] Matrix of shear probe cross-spectra after the coherent 
%             acceleration signals have been removed. The diagonal has the 
%             auto spectra of the clean shear.
% * [AA] Matrix of acceleration cross-spectra. The diagonal has the auto-spectra.
% * [UU] Matrix of shear probe cross-spectra without the removal of coherent
%        acceleration signals. The diagonal has the auto spectra of the original 
%        shear signal.
% * [UA] Matrix of cross-spectra of shear and acceleration signals.
% * [F] Column vector of frequency for the auto- and cross-spectra.
%
%%% Description
% Remove components within a shear probe signal that are coherent with 
% signals reported by the accelerometers. It is very effective at removing
% vibrational contamination from shear probe signals.
% 
% It works only in the spectral domain. 
%
% Developed by Louis Goodman (University of Massachusetts) and adapted by RSI.
%
% For best results and statistical significance, the length of the fft (nFFT) 
% should be several times shorter than the length of the vectors, [size(U,1)], 
% passed to this function. That is, several ffts have to be used to form an 
% ensemble-averaged auto- and cross-spectrum. The absolute minimum is
% L >= nFFT * (0.5*size(A,2) + 1), or L >= 2*nFFT,
% whichever is greater, where L = size(U,1), the length of the shear data.
% 
% @image @images/clean_shear_spectra.pdf @Coherent noise removal applied to shear
% probe spectrum. @Original spectrum (thin green line) and after noise removal 
% (thick green line). Data collected with Remus-600 in Buzzards Bay MA, courtesy 
% of Tim Boyd , Scottish Association for Marine Science. Only shear probe 2 was 
% installed, so, the spectrum for probe 1 falls off the deep end of the plot. 
%

% Version History
%
% * circa 2003 (Lou Goodman)
% * 2004-05-03 (RGL) tweaks
% * 2006-05-05 (RGL) "squeezed" some of the output vectors
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-05-15 (WID) updated documentation for publishing
% * 2012-07-11 (ILG) replaced csd_rolf and psd_rolf function calls with calls
%                    to csd_odas.m; for loops take advantage of additional 
%                    outputs available from csd_odas.m
% * 2014-04-10 (RGL) Modified to use csd_matrix_odas function for the
%   spectral calculations. Redundant fft calls are eliminated. Speed
%   improvement is about a factor of 2.
% * 2015-11-18 RGL, Documentation changes.
% * 2017-09-07 RGL, Applied new minimum length of shear data relative to
%   nFFT. See documentation.
% * 2022-04-13 RGL, Added code to correct for the bias of the goodman
%   cleaning based on Technical Note 61. 

function  [clean_UU, AA, UU, UA, F] = clean_shear_spec(A, U, n_fft, rate)

% do some checking
if isscalar(A) || isscalar(U)
    error('Acceleration and shear matrices must contain vectors')
end
if isrow(A), A = A'; end % Force into column vector
if isrow(U), U = U'; end 
if (size(A,1)~=size(U,1))
    error('Acceleration and shear matrices must have the same number of rows')
end
if ~isscalar(n_fft) || n_fft<2
    error('n_fft must be larger than 2.')
end

if (size(A,1) < 2*n_fft) || (size(A,1) < n_fft * (size(A,2)/2 + 1))
    error(['Vector length must be the larger of ' ...
        '(2*nFFT), OR ' ...
        '(nFFT * (size(A,2)/2 + 1)), ' ...
        'or longer.'])
end
if ~isscalar(rate) || rate <= 0
    error('Sampling rate must be positive scalar')
end
% end of checking

%method = 'parabolic';
method = 'linear';

% if new = true, we use the new and very efficient csd_matrix_odas function instead 
% of the csd_odas function. The improvement is about a factor of 2.
new = true; 

push = [2 3 1]; %The permutation 
if new
    [Cxy,F,Cxx,Cyy] = ...
        csd_matrix_odas(U, A, n_fft, rate,[], n_fft/2, method);
    UU = permute(Cxx, push);
    UA = permute(Cxy, push);
    AA = permute(Cyy, push);
    
else
    % Pre-allocate shear matrices
    AA = zeros(size(A,2),size(A,2),n_fft/2+1);
    UU = zeros(size(U,2),size(U,2),n_fft/2+1);
    UA = zeros(size(U,2),size(A,2),n_fft/2+1);
    
    
    % Compute auto- and cross-spectra
    for k = 1:size(U,2)
        % Shear-probe and acceleration cross-spectra, auto-spectra
        for m = 1:size(A,2)
            [UA(k,m,:), F, UU(k,k,:), AA(m,m,:)] = csd_odas(U(:,k),...
                A(:,m),n_fft,rate,[],n_fft/2,method);
        end
        % Shear probe cross-spectra
        for n = k+1:size(U,2)
            UU(k,n,:) = csd_odas(U(:,k),U(:,n),n_fft,rate,[],n_fft/2,method);
            UU(n,k,:) = conj(UU(k,n,:));
        end
    end
    % Acceleration cross-spectra
    for m = 1:size(A,2)
        for n = m+1:size(A,2)
            AA(m,n,:) = csd_odas(A(:,m), A(:,n), n_fft,rate,[],n_fft/2,method);
            AA(n,m,:) = conj(AA(m,n,:));
        end
    end
end
% Clean the shear spectra
clean_UU = complex(zeros(size(UU)));
for ii = 1:length(F)
   clean_UU (:,:,ii) = UU(:,:,ii) - (UA(:,:,ii)/AA(:,:,ii)) *conj(UA(:,:,ii)).';
end

% Remove extra dimensions
clean_UU = squeeze(clean_UU);
clean_UU = real(clean_UU); % probably not necessary
UU       = squeeze(UU);
AA       = squeeze(AA); 
UA       = squeeze(UA);

% Now remove the bias due to using a finite number of fft-segments.
fft_segments = floor( 2*size(U,1) / n_fft) -1;
vibration_signals = size(A,2);
R = 1 - 1.02*vibration_signals/fft_segments;
R = 1/R;

clean_UU = R*clean_UU;



   

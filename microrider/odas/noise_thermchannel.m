%% noise_thermchannel
% Compute the output noise spectrum of the FP07 thermistor channels.
%
% <latex>\index{Functions!noise\_thermchannel}</latex>
%
%%% Syntax
% [result, params] = noise_thermchannel(f, noise_info, ...)
%
% * [f] Column vector of frequencies [Hz] at which to evaluate the shear
%       circuit noise.
% * [noise_info] Structure (optional) containing configuration parameters. 
%       A template is generated and returned when noise_shearchannel is called 
%       with no input parameters. The parameters are described below.
% * [...] Optional configuration parameters to supplement, or override, 
%       those values included within noise_info. Inputs are accepted as
%       string/value pairs.
% * []
% * [result] Column vector of noise spectrum evaluated at frequencies f, in
%       units of counts^2 per Hz.
% * [params] An optional structure containing the parameters used to compute
%       the noise spectrum.
%
%%% Description
% Computes the output noise spectrum of the FP07 thermistor channels, in terms of
% raw counts from the analog-to-digital converter, including the bridge
% excitation, the differentiator, the anti-aliasing filter and the sampling
% noise. See RSI Technical Note 040 for more information.
%
% If $\texttt{noise\_thermchannel}$ is called without input arguments, it 
% returns the default parameters used by $\texttt{noise\_thermchannel}$. 
% You can then customize this structure to your particular processing 
% requirements. For example,
%
%
%   >> noise_info = noise_thermchannel
% 
%   noise_info = 
%            Bits: 16
%             E_n: 4.0000e-09
%            E_n2: 8.0000e-09
%              FS: 4.0960
%             G_D: 0.9400
%             K_B: 1.3820e-23
%             R_0: 3000
%             T_K: 295
%            f_AA: 110
%              fc: 18.7000
%            fc_2: 42
%              fs: 512
%            gain: 6
%       gamma_RSI: 3
%     make_figure: 0
%
% The configuration parameters can be provided as either a structure with
% specified fields of string/value pairs. The parameters can be separated
% into two groups: (1) general parameters that are independent of the shear
% probe circuitry and, (2) circuit-specific parameters that are determined 
% by the components used in the circuit design. 
%
%%% General parameters:
%
% * [T_K] temperature in Kelvin, default = 295.
% * [K_B] Boltzman's constant, 1.382e-23.
% * [VFS] full-scale voltage in volts of analog to digital converter
%         (sampler), default = 4.096.
% * [Bits] number of bits in the ADC, default = 16.
% * [gamma_RSI] factor by which the RSI sampler noise is higher than ideal,
%        default = 3. 
% * [fs] Sampling rate in Hz, default = 512.
% * [make_figure] logical value to toggle figures, default = false.
% 
%%% Circuit-specific parameters:
%
% * [R_0] nominal resistance, in Ohms, of thermistor at 17C, default = 3000.
% * [gain] gain of thermistor circuit, default = 6.
% * [G_D]  gain of differentiator, in seconds, default = 0.94.
% * [f_AA] half-power frequency, in Hz, of each stage of AA-filter, default = 110.
% * [E_n]  high-frequency input voltage noise of the first stage amplifier,
%       in V/sqrt(Hz), default = 4e-9.
% * [fc]  first-stage flicker noise knee frequency in Hz, default = 18.7.
% * [E_n2] high-frequency input voltage noise of the second stage amplifier,
%       in V/sqrt(Hz), default = 8e-9.
% * [fc_2]  second-stage flicker noise knee frequency in Hz, default = 42.
% 

% *Version History:*
% 2017-10-31 - Justine McMillan
% 2018-12-20 - RGL, Some clean up and testing an alternate noise model for
%       the LTC1678 (first stage). 
% 2019-04-08 - RGL, Revised to make the new noise formula intuative and
%   clear. Modified to make it a formal RSI function. 

function [result, params] = noise_thermchannel(f,varargin)

%-----------------------------------------------------------------
% ----- Default parameters ---------------------------------------
%-----------------------------------------------------------------

% set defaults (General)
default_T_K       = 295;% temperature in Kelvin
default_K_B       = 1.382e-23;% Boltzman constant
default_FS        = 4.096; % of ADC
default_Bits      = 16;
default_gamma_RSI = 3; % RSI comes  within a factor 3 of the ideal sampler-noise variance
default_fs        = 512; % Sampling rate

% Thermistor components
default_R_0   = 3000;   % Thermistor resistance
default_gain  = 6;      % thermistor channel gain
default_G_D   = 0.94;   % typical pre-emphasis gain factor
default_f_AA  = 110;    % cut-off frequency of each 4-pole Butterworth filter. 2 filters
default_E_n   = 4e-9;   % first-stage amplifier input voltage noise at high frequency, units of V/sqrt(Hz)
default_fc    = 18.7;   % first-stage flicker noise knee frequency
default_E_n2  = 8e-9;   % second-stage amplifier input voltage noise at high frequency, units of V/sqrt(Hz).
default_fc_2  = 42;     % second-stage flicker noise knee frequency
default_make_figure  = false; 

% Get defaults if there are no input parameters
if ~nargin 
    for d = whos('default_*')'
        param = regexp(d.name, 'default_(.*)', 'tokens');
        result.(char(param{1})) = eval(d.name);
    end
    return
end

%-----------------------------------------------------------------
% ----- Parse Inputs ---------------------------------------------
%-----------------------------------------------------------------
p = inputParser;
p.KeepUnmatched = true;
p.CaseSensitive = true;

val_numeric     = @(x) isnumeric(x) && isscalar(x);
val_logical     = @(x) islogical(x);
%val_string      = @(x) ischar(x);
  
addParameter(p, 'T_K',       default_T_K,        val_numeric);
addParameter(p, 'K_B',       default_K_B,        val_numeric);      
addParameter(p, 'FS',        default_FS,         val_numeric);   
addParameter(p, 'Bits',      default_Bits,       val_numeric);     
addParameter(p, 'gamma_RSI', default_gamma_RSI,  val_numeric); 
addParameter(p, 'fs',        default_fs,         val_numeric);        
addParameter(p, 'R_0',       default_R_0,        val_numeric);     
addParameter(p, 'gain',      default_gain,       val_numeric);
addParameter(p, 'G_D',       default_G_D,        val_numeric);     
addParameter(p, 'f_AA',      default_f_AA,       val_numeric); 
addParameter(p, 'E_n',       default_E_n,        val_numeric);     
addParameter(p, 'fc',        default_fc,         val_numeric); 
addParameter(p, 'E_n2',      default_E_n2,       val_numeric);     
addParameter(p, 'fc_2',      default_fc_2,       val_numeric); 
addParameter(p, 'make_figure', default_make_figure, val_logical);


% Save the input for record keeping
% ql_info_in = join_arguments(varargin); 

% Parse the arguments.
parse(p, varargin{:});

% Save into params struct
names = fieldnames(p.Results);
for name = names'
    eval(['params.' char(name) ' = p.Results.' char(name) ';']);
    eval([char(name) ' = p.Results.' char(name) ';']);
end

clear default_*



%-----------------------------------------------------------------
% ----- Compute Noise ---------------------------------------------
%-----------------------------------------------------------------
% f = f';
if size(f,2) > 1, f = f'; end % make it a column vector

delta_s  = FS / pow2(Bits); % step size of sampler
fN       = fs/ 2; % Nyquist frequency

% ---- Stage 1 ----
V1 = 2*E_n^2 * sqrt(1 + (f/fc).^2) ./ (f / fc); % An alternate but very standard model.
%disp('Using Rolf model noise for input voltage noise')

phi_R   = 4 * K_B * R_0 * T_K; % Thermal or Johnson noise
Noise_1 = gain^2 * (V1 + phi_R); % The noise output of the first stage.

% ---- Stage 2 ----
G_2 = 1 + (2 * pi * G_D * f).^2; % Pre-emphasis gain
V2  = 2*E_n2^2 .* sqrt(1 + (f/fc_2).^2) ./ (f / fc_2); % voltage noise of 
        % the LTC6240. The factor of 2 is because it is a differential
        % op-amp.
Noise_2 = G_2 .* (Noise_1 + V2); % Noise output of stage 2,


% ---- Anti-aliasing stage ----
G_AA = (1 + (f/f_AA).^8 ).^2; G_AA = 1 ./ G_AA;
Noise_3 = Noise_2 .* G_AA; % Noise output of AA-filter stage

% ---- Sampling stage ----
Noise_4 = Noise_3 + gamma_RSI * (delta_s^2) ./ (12 * fN); % Final output 
        % noise in equivalent voltage variance of the sampled data.

% ---- Total Noise -----
Noise_4_counts = Noise_4 / delta_s^2; % Now in terms of counts of the ADC.


% ---- Resistor Noise -----
noise_R = phi_R * gain^2 * G_2 / delta_s^2; % Output noise due to thermistor, alone.

result = Noise_4_counts;


%-----------------------------------------------------------------
% ----- Make figures, if requested -------------------------------
%-----------------------------------------------------------------
if make_figure
    % Stage 1
    figure(1), clf, set(gcf, 'color', 'w')
    y_lim = [1e-15 1e-12];
    x_lim = [1e-1 1e3];
    h=loglog(f, [gain^2*V1  gain^2*phi_R*ones(size(f)) Noise_1]);
    set(gca, 'ylim', y_lim)
    set(gca, 'xlim', x_lim)
    legend('\phi_{\itAmp}', '\phi_{\itR_T}', 'Sum', 'location','northeast')
    xlabel('$ f\ \mathrm{[\ Hz\ ]}$', 'interpreter', 'latex')
    ylabel('$\mathrm{[\ V^2 \, Hz^{-1}\ ]}$', 'interpreter', 'latex')
    title_string = {...
        '\rmThermistor channel output noise -- stage 1',...
        ['\itR_T\rm = ' num2str(R_0, 5) '\Omega, Gain = ' num2str(gain)]}; 
    title(title_string)
    set(h,'linewidth',2)
    set(gca, 'fontsize', 40)
    
    figure(2),clf, set(gcf, 'color', 'w')
    y_lim = [1e-15 1e-8];
    h=loglog(f, [Noise_1 Noise_2 Noise_3 Noise_4]);
    set(gca, 'ylim', y_lim)
    set(gca, 'xlim', x_lim)
    set(h,'linewidth',2)
    xlabel('$ f\ \mathrm{[\ Hz\ ]}$', 'interpreter', 'latex')
    ylabel('$\mathrm{[\ V^2 \, Hz^{-1}\ ]}$', 'interpreter', 'latex')
    legend('stage 1','stage 2','AA-filter','sampler', 'location', 'northwest')
    title_string = {...
        '\rmThermistor channel output noise',...
        ['\itR_T\rm = ' num2str(R_0, 5) '\Omega, Gain = ' num2str(gain) ', \itG_D\rm = ' num2str(G_D) 's'], ...
        ['\itV_{FS}\rm = ' num2str(FS) 'V, \itB\rm = ' num2str(Bits) ', \itf_s\rm = ' num2str(fs) 'Hz, \gamma = ' num2str(gamma_RSI)]}; 
    title(title_string)
    set(h,'linewidth',1.5)
    set(h(2),'linewidth',7)
    set(h(3),'linewidth',3)
    set(gca, 'fontsize', 40)

    % Total noise
    figure(3),clf, set(gcf, 'color', 'w')
    position = [0.13   0.1100    0.75    0.7150];
    subplot('position', position);
    % in units of Volts^2/Hz
    line(f, Noise_4, 'color', 'k', 'linewidth',3)
    ax1 = gca;
    set(ax1, 'xcolor', 'k');
    set(ax1, 'xscale', 'log');
    set(ax1, 'yscale', 'log');
    ax1_pos = get(ax1, 'position');
    y_lim = [1e-12 1e-8];
    set(ax1, 'xlim', x_lim, 'ylim', y_lim)
    xlabel('$ f\  \mathrm{[\ Hz\ ]}$', ...
        'interpreter','latex', 'color', 'k')
    ylabel('$ \Psi_{N_V} (f)\  \mathrm{[\ V^2\, Hz^{-1}\ ]}$', ...
        'interpreter','latex', 'color', 'k')
    legend('$ \Psi_{N_V} $', 'interpreter', 'latex', 'location', 'northwest')

%     title_string = {...
%         '\rmOutput noise of sampled thermistor data',...
%         ['\itR_T\rm = ' num2str(R_T, 5) '\Omega']}; 
    title(title_string)
    set(gca, 'fontsize', 30)

    % in units of counts^2/Hz
    ax2 = axes(...
        'position', ax1_pos, ...
        'xaxislocation', 'bottom', ...
        'yaxislocation','right', ...
        'color', 'none');
    line(f, Noise_4_counts, 'Parent', ax2, 'color', 'r', 'linewidth',3)
    hold all
    line(f, noise_R,'Parent', ax2, 'color', 'r', 'linewidth',3, 'linestyle','--')
    legend('$ \Psi_{N_S} $', '$ \Psi_{R_T} $', 'interpreter', 'latex')

    set(ax2, 'ycolor', 'r');
    set(ax2, 'xscale', 'log');
    set(ax2, 'yscale', 'log');
    y_lim = [1e-4 1e0];
    set(ax2, 'xlim', x_lim, 'ylim', y_lim)
    set(ax2, 'xlim', x_lim)
    ylabel('$ \Psi_{N_S} (f)\  \mathrm{[\ counts^2\, Hz^{-1}\ ]}$', ...
        'interpreter','latex', 'color', 'r')

    set(gca, 'fontsize', 30)


end



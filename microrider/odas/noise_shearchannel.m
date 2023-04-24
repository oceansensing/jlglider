%% noise_shearchannel
% Compute the output noise spectrum of the shear channels.
%
%%
% <latex>\index{Functions!noise\_shearchannel}</latex>
%
%%% Syntax
%   [result,params] = noise_shearchannel(f,noise_info,...)
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
%       units of counts^2/Hz.
% * [params] An optional structure containing the parameters used to compute
%       the noise spectrum.
%
%%% Description
% This function can be used to estimate the spectrum of the electronic
% noise of the shear probe circuit that includes a charge-transfer 
% amplifier, differentiator, AA-filter, and data sampler. See RSI
% Technical Note 042 for more details. 
%
% If $\texttt{noise\_shearchannel}$ is called without input arguments, it 
% returns the default parameters used by $\texttt{noise\_shearchannel}$. 
% You can then customize this structure to your particular processing 
% requirements. For example,
%
%    >> noise_info = noise_shearchannel
%    
%    noise_info = 
% 
%             Bits: 16
%               C1: 1.5000e-09
%               C2: 9.4000e-07
%               C3: 4.7000e-10
%               CP: 0
%              E_1: 9.0000e-09
%              I_1: 5.6000e-16
%              K_B: 1.3820e-23
%               R1: 1.0000e+09
%               R2: 499
%               R3: 1000000
%              T_K: 295
%              VFS: 4.0960
%             f_AA: 110
%               fc: 50
%               fs: 512
%        gamma_RSI: 2.5000
%     make_figures: 0
%
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
%        default = 2.5. 
% * [fs] Sampling rate in Hz, default = 512.
% * [make_figures] logical value to toggle figures, default = false.
% 
%%% Circuit-specific parameters:
%
% * [R1] first stage resistance in Ohms, default = 1e9.      
% * [C1] first stage capacitance in F, default = 1.5e-9. 
% * [E_1] first-stage amplifier input voltage noise at high frequency, 
%       in V/sqrt(Hz), default = 9e-9.
% * [fc] first-stage flicker noise knee frequency in Hz, default = 50.
% * [I_1] first-stage amplifier input current noise in A/sqrt(Hz), 
%       default = 0.56e-15.
% * [R2] second stage resistance in Ohms, default = 499.
% * [C2] second stage capacitance in F, default = 0.94e-6.
% * [R3] third stage resistance in Ohms, default = 1e6.
% * [C3] third stage capacitance in F, default = 470e-12.
% * [CP] probe capacitance in F, default = 0 (no probe).
%         Note, set to 1e-9 when probe installed. 
% * [f_AA] cut-off frequency in Hz of each antialiasing filter, default = 110. 

% *Version History:*
% 2017-10-31 (JMM) - Original version
% 2019-01-22 (JMM) - Cleaned up to try to mimic noise_thermchannel
% 2019-01-23 (RGL) - Corrected default parameters, some cleaning
% 2019-06-05 (RGL) - Added probe capacitance to the noise model. Changed
%                   scale on one figure, slightly, cleanup the titles,
%                   fonts, etc. Changed names of parameter to match the
%                   Application Note.
% 2019-06-06 (RGL) - Corrected a blunder in the first-stage noise gain.
% 2019-06-26 (JMM) - Changed the name of the function to add it to the odas
%                    library and updated documentation.


function [result,params] = noise_shearchannel(f,varargin)

%-----------------------------------------------------------------
% ----- Default parameters ---------------------------------------
%-----------------------------------------------------------------

% set defaults (General)
default_T_K       = 295;  % temperature in Kelvin
default_K_B       = 1.382e-23;% Boltzman constant
default_VFS       = 4.096; % Full-scale range of ADC
default_Bits      = 16; % number of bits in the ADC
default_gamma_RSI = 2.5;   % RSI comes  within a factor 2.5 of the ideal sampler-noise variance
default_fs        = 512; % Sampling rate

% shear channel components
default_R1      = 1e9;      % the 1 G-ohm feedback capacitor (LTC6240 Charge-transfer amplifier)
default_C1      = 1.5e-9;   % The 1.5 nF charge-transfer capacitor (LTC6240 Charge-transfer amplifier)
default_R2      = 499;      % differentiator resistor (P049R02a spec sheet)
default_C2      = 0.94e-6;  % differentiator capacitor (P049R02a spec sheet, 2*0.47 uF)
default_R3      = 1e6;      % differentiator resistor
default_C3      = 470e-12;  % differentiator capacitor 
default_CP      = 0;        % probe capacitance. ==0 no probe, or =1e-9 when probe installed. 
default_f_AA    = 110;      % cut-off frequency of each 4-pole Butterworth filter. 2 filters in series.

% Op Amp parameters
default_E_1 = 9e-9;     % first-stage amplifier input voltage noise at high frequency, units of V/sqrt(Hz).
default_fc  = 50;       % first-stage flicker noise knee frequency
default_I_1 = 0.56e-15; % first-stage amplifier input current noise -- frequency independent up to 300 Hz

% Other parameters
default_make_figures  = false; 

% Return the default values, if there are no input parameters.
if ~nargin 
    for d = whos('default_*')'
        param = regexp(d.name, 'default_(.*)', 'tokens');
        result.(char(param{1})) = eval(d.name);
    end
    return
end

% parse defaults and inputs
p = inputParser;
p.KeepUnmatched = true;
p.CaseSensitive = true;

val_logical     = @(x) islogical(x);
val_positive    = @(x) isnumeric(x) && isscalar(x) && x>=0;
  
addParameter(p, 'T_K',       default_T_K,        val_positive);
addParameter(p, 'K_B',       default_K_B,        val_positive);      
addParameter(p, 'VFS',       default_VFS,        val_positive);   
addParameter(p, 'Bits',      default_Bits,       val_positive);     
addParameter(p, 'gamma_RSI', default_gamma_RSI,  val_positive); 
addParameter(p, 'fs',        default_fs,         val_positive); 
addParameter(p, 'R1',        default_R1,         val_positive);     
addParameter(p, 'C1',        default_C1,         val_positive);
addParameter(p, 'R2',        default_R2,         val_positive);  
addParameter(p, 'C2',        default_C2,         val_positive);  
addParameter(p, 'R3',        default_R3,         val_positive);  
addParameter(p, 'C3',        default_C3,         val_positive);  
addParameter(p, 'CP',        default_CP,         val_positive);  
addParameter(p, 'f_AA',      default_f_AA,       val_positive); 
addParameter(p, 'E_1',       default_E_1,        val_positive);     
addParameter(p, 'fc',        default_fc,         val_positive); 
addParameter(p, 'I_1',       default_I_1,        val_positive); 
addParameter(p, 'make_figures',  default_make_figures,  val_logical); 

%-----------------------------------------------------------------
% ----- Parse Inputs ---------------------------------------------
%-----------------------------------------------------------------
parse(p, varargin{:});

% Save into info struct
names = fieldnames(p.Results);
for name = names'
    eval(['params.' char(name) ' = p.Results.' char(name) ';']);
    eval([char(name) ' = p.Results.' char(name) ';']);
end

clear default_*

%-----------------------------------------------------------------
% ----- Compute Noise ---------------------------------------------
%-----------------------------------------------------------------
if size(f,2) > 1, f = f'; end % make it a column vector

omega = 2*pi *f;

delta_s  = VFS / pow2(Bits); % step size of sampler
fN       = fs/ 2; % Nyquist frequency

% ---- Stage 1 ----
V_V1 = E_1^2 * (fc ./ f) .* sqrt(1 + (f/fc).^2); % voltage noise -- model based on spec sheet 
V_I1 = (I_1^2) .* R1^2 ./ (1 + (omega*R1*C1).^2); % current noise -- based on spec sheet frequency independent up to 300 Hz
V_R1 = 4*K_B*T_K* R1 ./ (1 + (omega*R1*C1).^2); % Resistor thermal (Johnson) noise.

G_1 = (1 + (omega*R1*(CP+C1)).^2) ./ (1 + (omega*R1*C1).^2  ); % Noise gain of first stage

Noise_1 = G_1 .* (V_V1 + V_I1) + V_R1;

% ---- Stage 2 ----
G_2 = (omega*R3*C2).^2 ./ (1 + (omega*R2*C2).^2 ) ./ (1 + (omega*R3*C3).^2); % Squared gain 

Noise_2 = (Noise_1 + V_V1) .* G_2; 

% ---- Anti-aliasing stage ----
G_AA = 1 + (f / f_AA).^8; % inverse of one AA-filter
G_AA = 1 ./ G_AA.^2; % cascade of 2 AA-filters

Noise_3 = Noise_2 .* G_AA; % Noise output of AA-filter stage


% ---- Sampling stage ----
Noise_4 = Noise_3 + gamma_RSI * (delta_s^2) ./ (12 * fN); % Final output units of V^2 Hz^{-1}

% ---- Total Noise -----
Noise_4_counts = Noise_4 / delta_s^2; % Now in terms of counts of the ADC.

result = Noise_4_counts;

%-----------------------------------------------------------------
% ----- Make figures, if requested -------------------------------
%-----------------------------------------------------------------
if make_figures
    title_string_1 = {...
        '\rmCharge-Transfer Amp Output Noise, LTC6240',...
        ['\itR_{\rm1}\rm = ' make_scientific(R1,2) '\Omega, ' ...
        '\itC_{\rm1}\rm = ' make_scientific(C1,2) 'F, ' ...
        '\itC_p\rm = ' make_scientific(CP,2) 'F']};
    if CP == 0
        title_string_1 = {...
            '\rmCharge-Transfer Amp Output Noise, LTC6240',...
            ['\itR_{\rm1}\rm = ' make_scientific(R1,2) '\Omega, ' ...
            '\itC_{\rm1}\rm = ' make_scientific(C1,2) 'F, ' ...
            '\itC_p\rm = 0F']};
    end
    
    title_string_2 = {...
        '\rmShear-probe circuit noise, LTC6240',...
        ['\itR_{\rm1}\rm = ' make_scientific(R1,2) '\Omega, ' ...
        '\itC_{\rm1}\rm = ' make_scientific(C1,2) 'F, ' ...
        '\itC_{\rm2}\rm = ' make_scientific(C2,2) 'F, ' ...
        '\itC_p\rm = ' make_scientific(CP,2) 'F, ' ...
        '\itf\rm_s = ' num2str(fs,3) 's^{-1}']};
    if CP == 0
        title_string_2 = {...
            '\rmShear-probe circuit noise, LTC6240',...
            ['\itR_{\rm1}\rm = ' make_scientific(R1,2) '\Omega, ' ...
            '\itC_{\rm1}\rm = ' make_scientific(C1,2) 'F, '  ...
            '\itC_{\rm2}\rm = ' make_scientific(C2,2) 'F, ' ...
            '\itC_p\rm = 0F, \itf\rm_s = ' num2str(fs,3) 's^{-1}']};
    end
    
    % Stage 1
    figure(1), clf
    set(gcf, 'color', 'w')
    y_lim = [1e-18 1e-10];
    h=loglog(f, [V_V1  V_I1 V_R1 Noise_1]);
    set(gca, 'ylim', y_lim)
    set(gca, 'fontname', 'times')
    legend(...
        '$\phi_{V_1}$', ...
        '$\phi_{I_1}$', ...
        '$\phi_{R_1}$', ...
        '$\phi_1$', ...
        'interpreter', 'latex', 'location','northeast')
    xlabel('$ f\ \mathrm{[\ Hz\ ]}$', 'interpreter', 'latex')
    ylabel('$\mathrm{[\ V^2 \, Hz^{-1}\ ]}$', 'interpreter', 'latex')
    title(title_string_1, 'fontsize', 20)
    set(h,'linewidth',2)
    
    % Noise at the output of different stages
    figure(2), clf
    set(gcf, 'color', 'w')
    y_lim = [1e-13 1e-9];
    x_lim = [1e-1 max(f)];
    h=loglog(f, [Noise_1 Noise_2 Noise_3 Noise_4]);
    set(gca, 'ylim', y_lim)
    set(gca, 'fontname', 'times')
    legend('$\phi_1$','$\phi_2$','$\phi_3$','$\phi_N$', 'interpreter', 'latex')
    xlabel('$ f\ \mathrm{[\ Hz\ ]}$', 'interpreter', 'latex')
    ylabel('$\mathrm{[\ V^2 \, Hz^{-1}\ ]}$', 'interpreter', 'latex')
    title(title_string_2, 'fontsize', 20)
    set(h,'linewidth',2)
    
    % Total noise in voltage and counts
    figure(3),clf
    set(gcf, 'color', 'w')
    position = [0.13   0.1400    0.75    0.7150];
    subplot('position', position);
    % in units of Volts^2/Hz
    line(f, Noise_4, 'color', 'k', 'linewidth',2)
    ax1 = gca;
    set(gca, 'fontname', 'times')
    set(ax1, 'xcolor', 'k');
    set(ax1, 'xscale', 'log');
    set(ax1, 'yscale', 'log');
    ax1_pos = get(ax1, 'position');
    y_lim = [1e-12 2e-10];
    set(ax1, 'xlim', x_lim, 'ylim', y_lim)
    xlabel('$ f\  \mathrm{[\ Hz\ ]}$', ...
        'interpreter','latex', 'color', 'k')
    ylabel('$ \phi_N (f)\  \mathrm{[\ V^2\, Hz^{-1}\ ]}$', ...
        'interpreter','latex', 'color', 'k')
    
    title(title_string_2, 'fontsize', 20)
    box on
    
    % in units of counts^2/Hz
    ax2 = axes(...
        'position', ax1_pos, ...
        'xaxislocation', 'bottom', ...
        'yaxislocation','right', ...
        'color', 'none');
    line(f, Noise_4_counts, 'Parent', ax2, 'color', 'r', 'linewidth',2)
    hold all
    % line(f, Noise_thermal, 'Parent', ax2, 'color', 'b', 'linewidth',2)
    
    set(gca, 'fontname', 'times')
    set(ax2, 'ycolor', 'r');
    set(ax2, 'xscale', 'log');
    set(ax2, 'yscale', 'log');
    set(ax2,'xgrid','off','ygrid','off')
    y_lim = [1e-4 2e-2];
    set(ax2, 'xlim', x_lim, 'ylim', y_lim)
    set(ax2, 'xlim', x_lim)
    ylabel('$ \Psi_N (f)\  \mathrm{[\ counts^2\, Hz^{-1}\ ]}$', ...
        'interpreter','latex', 'color', 'r')
end


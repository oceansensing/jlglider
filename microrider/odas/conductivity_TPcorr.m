%% conductivity_TPcorr
% Apply a correction for the thermal expansion and pressure contraction of
% a Sea-Bird SBE4C conductivity cell.
%%
% <latex>\index{Functions!conductivity\_TPcorr}</latex>
%
%%% Syntax
%   [sbcCorr] = conductivity_TPcorr( sbc, sbt, pres, CPcor, CTcor )
%
% * [sbc] Conductivity from Sea-Bird SBE4C conductivity sensor, in units
%         of mS/cm.
% * [sbt] Temperature from Sea-Bird SBE3F thermometer in units of Celsius.
% * [pres] Pressure in units of dBar.
% * [CPcor] Correction factor for pressure. default: -9.57e-8
% * [CTcor] Correction factor for temperature. default: 3.25e-6
% * []
% * [sbcCorr] Conductivity with correction applied, in units of mS/cm.
%
% Inputs CPcor and CTcor are an optional pair.  Both must be either
% included or excluded.
%
%%% Description
% Apply a correction for the thermal expansion and pressure contraction
% of the Sea-Bird SBE4C conductivity cell.
%
% The algorithm used:
%
% $$sbcCorr = \frac{sbc}{1 + CTcor \cdotp sbt + CPcor \cdotp pres}$$
%

% Version History
%
% * 2005-02-14 IG from OTL sb_c.m routine.
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-05-14 (WID) documentation update
% * 2012-11-02 (WID) documentation update
% * 2012-10-27 (RGL) Documentatioon changes.

function sbc_corr = conductivity_TPcorr(sbc, sbt, pres, CPcor, CTcor)

% Default coefficients for compressibility and thermal expanion (from
% Sea-Bird calibration sheet, May 2004)
if nargin<5
    CPcor = -9.5700e-8;
    CTcor = 3.2500e-6;
end

% Correct the conductivity data.  Note that this is typically a fairly
% small correction.
sbc_corr = sbc./(1+CTcor*sbt+CPcor*pres);

end

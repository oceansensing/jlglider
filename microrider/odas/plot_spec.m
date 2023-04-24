%% plot_spec
% Plot frequency and wavenumber spectra of shear, scalar gradients and
% acceleration.
%%
% <latex>\index{Depreciated!plot\_spec}</latex>
%
%%% Description
% This is a legacy function that has been renamed to show_spec. Call
% show_spec instead of plot_spec in future usage.

% 2015-10-30 RGL Replaced with show_spec

function plot_spec(diss, titleString)

warning('The plot_spec function is calling show_spec. In future, please call show_spec directly')

show_spec (diss, titleString)

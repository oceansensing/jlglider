%% fix_underscore
% Legacy function replaced with texstr
%%
% <latex>\index{Depreciated!fix\_underscore}</latex>
%
%%% Description
% Depreciated function.  Please use $\texttt{texstr}$ in place of this
% function.

% Version History
%
% * 2000-07-01 (RGL)
% * 2011-09-01 (AWS) added documentation tags for matlab publishing
% * 2012-05-14 (WID) documentation update
% * 2012-11-05 (WID) documentation update

function string_out=fix_underscore(string_in)
    warning('The function fix_underscore is depreciated.  Please use texstr.');
    string_out = texstr(string_in);
end


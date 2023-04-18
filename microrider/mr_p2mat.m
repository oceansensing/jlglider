% This matlab script runs odas_p2mat and converts all P files into MAT files
% gong@vims.edu 2023-04-18

datadirJM = '/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20221021-norse-janmayen-complete/mr1000g/';
datadirLBE = '/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20221102-norse-lofoten-complete/mr1000g/';

datadir = datadirJM
datafiles = dir([datadir '*.p']);

for ii = 1:length(datafiles)
    odas_p2mat([datadir datafiles(ii).name]);
end

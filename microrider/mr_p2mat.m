% This matlab script runs odas_p2mat and converts all P files into MAT files
% gong@vims.edu 2023-04-18

datadirJM2022 = '/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20221021-norse-janmayen-complete/mr1000g/';
datadirLBE2022 = '/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20221102-norse-lofoten-complete/mr1000g/';
datadirJM2023 = '/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/sea064-20231127-norse-complete/mr/';

datadir = datadirJM2023
datafiles = dir([datadir '*.p']);

for ii = 1:length(datafiles)
    odas_p2mat([datadir datafiles(ii).name]);
end

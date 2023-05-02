% This script loads MR data from SEA064 and data from UAF glider Freya processed by Laur Ferris.
% Time in Matlab's datenum format is converted to unix/posix time prior to resaving the data files.
% gong@vims.edu  2023-05-01

load jm_P.mat
dtime = datetime(dnum,'ConvertFrom','datenum','TimeZone','America/New_York');
unixt = posixtime(dtime);
save jm_P -nocompression

load lb_P.mat
dtime = datetime(dnum,'ConvertFrom','datenum','TimeZone','America/New_York');
unixt = posixtime(dtime);
save lb_P -nocompression

load freya.mat
freya.dtime = datetime(freya.dnum, 'ConvertFrom', 'datenum', 'TimeZone', 'UTC');
freya.unixt = posixtime(freya.dtime);
freyanav.dtime = datetime(freyanav.dnum, 'ConvertFrom', 'datenum', 'TimeZone', 'UTC');
freyanav.unixt = posixtime(freyanav.dtime);
save freya.mat -nocompression

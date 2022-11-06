# loading realtime SeaExplorer data from NORSE
# gong@vims.edu 20221021

import seaexplorer_functions: load_NAV, load_PLD
#import seaexplorer_types: NAV_RT, PLD_RT

# setting src and data directory paths
srcdir = "/Users/gong/GitHub/jlglider/";
#dataroot = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/";
#dataroot = "/Users/gong/oceansensing Dropbox/Donglai Gong/Projects/NORSE/2022_fieldwork/";
dataroot = "/Users/gong/Research/sea064/";

# define dataset loading parameters
#project = "maracoos"
#deploydate = "20220311"
#suffix = "data"

projectname = "norse"
deploydate = "20221021"
suffix = "deployment"

if (@isdefined gliderSN) != true
    gliderSN = 64
end

if (@isdefined mission) != true
    mission = 38
end

# define data load location
#datadir = dataroot * glidername * "-" * deploydate * "-" * project * "-" * suffix * "/";
datadir = dataroot;
navdir = datadir * "realtime/";

if mission < 38
    scidir = datadir * "delayed/";
else
    scidir = datadir * "realtime/";
end

(sea064nav, sea064nav1d) = load_NAV(gliderSN, mission, navdir, 1);
(sea064pld, sea064pld1d) = load_PLD(gliderSN, mission, scidir, 1); # last dataflag parameter, 0 for sub individual files, 1 for sub all, >2 for raw individual files

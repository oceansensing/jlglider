# loading realtime SeaExplorer data from NORSE
# gong@vims.edu 20221021

import seaexplorer_functions: load_NAV, load_PLD
import seaexplorer_functions: missing2nan, cleanTime, cleanAD2CPtime, cleanEPS, cleanTemp, cleanSalt, cleanPress, clean9999, cleanAD2CP, cleanFLBBCDchl, cleanFLBBCDbb700, cleanFLBBCDcdom
#import seaexplorer_types: NAV_RT, PLD_RT

# setting src and data directory paths
srcdir = "/Users/gong/GitHub/jlglider/";
#dataroot = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/";
dataroot = "/Users/gong/oceansensing Dropbox/Donglai Gong/Projects/NORSE/2022_fieldwork/";
#dataroot = "/Users/gong/Research/sea064/";

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
    mission = 37
end

# define data load location
#datadir = dataroot * glidername * "-" * deploydate * "-" * project * "-" * suffix * "/";

if mission == 37 # Jan Mayen 2022
    datadir = dataroot * "sea064-20221021-norse-janmayen-complete/";
    navdir = datadir * "nav/logs/";
    scidir = datadir * "science/logs/";
    dataflag = 2;
elseif mission == 38 # Lofoten Basin 2022
    datadir = dataroot;
    navdir = datadir * "realtime/";
    scidir = datadir * "realtime/";
    dataflag = 1;
end

(sea064nav, sea064nav1d) = load_NAV(gliderSN, mission, navdir, dataflag);
(sea064pld, sea064pld1d) = load_PLD(gliderSN, mission, scidir, dataflag); # last dataflag parameter, 0 for sub individual files, 1 for sub all, >2 for raw individual files

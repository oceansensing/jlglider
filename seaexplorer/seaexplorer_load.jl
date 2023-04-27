# loading realtime SeaExplorer data from NORSE
# gong@vims.edu 20221021

using Glob, DataFrames, CSV, Dates, Missings

import seaexplorer_functions: load_NAV, load_PLD
import seaexplorer_functions: missing2nan, cleanTime, cleanAD2CPtime, cleanEPS, cleanTemp, cleanSalt, cleanPress, clean9999, cleanAD2CP, cleanFLBBCDchl, cleanFLBBCDbb700, cleanFLBBCDcdom
#import seaexplorer_types: NAV_RT, PLD_RT

# setting src and data directory paths
srcdir = "/Users/gong/GitHub/jlglider/";
dataroot = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/";
#dataroot = "/Users/gong/oceansensing Dropbox/Donglai Gong/Projects/NORSE/2022_fieldwork/";
#dataroot = "/Users/gong/Research/sea064/";

if (@isdefined gliderSN) != true
    gliderSN = 64
end

if (@isdefined mission) != true
    mission = 37
end

# define data load location
#datadir = dataroot * glidername * "-" * deploydate * "-" * project * "-" * suffix * "/";

# define dataset loading parameters
#project = "maracoos"
#deploydate = "20220311"
#suffix = "data"

if (gliderSN == 64) & (mission == 37)
    projectname = "norse"
    deploydate = "20221021"
    suffix = "deployment"
    datadir = dataroot * "sea064-20221021-norse-janmayen-complete/";
    dataflag = "all";

elseif (gliderSN == 64) & (mission == 38)
    projectname = "norse"
    deploydate = "20221102"
    suffix = "deployment"
    datadir = dataroot * "sea064-20221102-norse-lofoten-complete/";
    dataflag = "all";
end

if (dataflag == "realtime") | (dataflag == "all") 
    navdir = datadir * "glimpse/";
    scidir = datadir * "glimpse/";
else 
    navdir = datadir * "nav/logs/";
    scidir = datadir * "science/logs/";
end

(sea064nav, sea064nav1d) = load_NAV(gliderSN, mission, navdir, dataflag);
(sea064pld, sea064pld1d) = load_PLD(gliderSN, mission, scidir, dataflag); # last dataflag parameter, 0 for sub individual files, 1 for sub all, >2 for raw individual files

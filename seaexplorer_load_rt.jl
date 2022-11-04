# loading realtime SeaExplorer data from NORSE
# gong@vims.edu 20221021

import seaexplorer_functions: load_NAV_rt, load_PLD_rt
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

project = "norse"
deploydate = "20221021"
suffix = "deployment"

gliderSN = 64
mission = 38

# define data load location
#datadir = dataroot * glidername * "-" * deploydate * "-" * project * "-" * suffix * "/";
datadir = dataroot;
navdir = datadir * "realtime/";
scidir = datadir * "realtime/";

(sea064nav, sea064nav1d) = load_NAV_rt(gliderSN, mission, navdir, 1);
(sea064pld, sea064pld1d) = load_PLD_rt(gliderSN, mission, scidir, 1);

# loading realtime SeaExplorer data from NORSE
# gong@vims.edu 20221021

import seaexplorer_functions: load_NAV_rt, load_PLD_rt
import seaexplorer_types: NAV_RT, PLD_RT

# setting src and data directory paths
srcdir = "/Users/gong/GitHub/jlglider/";
#dataroot = "/Users/gong/oceansensing Dropbox/C2PO/glider/gliderData/";
dataroot = "/Users/gong/oceansensing Dropbox/Donglai Gong/Projects/NORSE/2022_fieldwork/";

# define dataset loading parameters
#project = "maracoos"
#deploydate = "20220311"
#suffix = "data"

project = "norse"
deploydate = "20221021"
suffix = "deployment"

glidername = "sea064"
mission = "37"

# define data load location
datadir = dataroot * glidername * "-" * deploydate * "-" * project * "-" * suffix * "/";
navdir = datadir * "realtime/";
scidir = datadir * "realtime/";

(sea064nav, sea064nav1d) = load_NAV_rt(glidername, mission, navdir);
(sea064pld, sea064pld1d) = load_PLD_rt(glidername, mission, scidir);

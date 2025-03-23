module Glider
    include("/Users/gong/GitHub/jlglider/seaexplorer/seaexplorerType.jl")
    include("/Users/gong/GitHub/jlglider/plotting/gliderPlotType.jl")
    import .seaexplorerType: NAV_RT, PLD_RT, SeaExplorerData, SeaExplorerCTD, LEGATO
    import .gliderPlotType: plotSetting, plotStruct
end
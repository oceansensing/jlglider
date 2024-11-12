module Glider
    include("/Users/gong/GitHub/jlglider/seaexplorer/seaexplorerType.jl")
    include("/Users/gong/GitHub/jlglider/seaexplorer/gliderPlotType.jl")
    include("/Users/gong/GitHub/jlglider/slocum/slocumType.jl")
    import .seaexplorerType: NAV_RT, PLD_RT, SeaExplorerData, SeaExplorerCTD, LEGATO
    import .slocumType: engStruct, ctdStruct, sciStruct
    import .gliderPlotType: plotSetting, plotStruct
end
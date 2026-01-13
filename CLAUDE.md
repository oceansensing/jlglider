# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

jlglider is a Julia-based toolkit for processing ocean glider data from multiple platforms (SeaExplorer, Slocum) and sensors (AD2CP current profilers, CTD, optical sensors, microstructure sensors). The codebase processes time-series oceanographic data collected during glider missions.

## Running Julia Code

This repository contains Julia scripts (`.jl` files) that need to be run with the Julia REPL:

```bash
# Start Julia REPL
julia

# Within Julia, include and run scripts
include("/full/path/to/script.jl")

# Or run directly from command line
julia script.jl
```

**Important**: All `include()` statements in this codebase use hardcoded absolute paths (e.g., `/Users/gong/GitHub/jlglider/...`). When working on different machines, these paths will need to be updated.

## Core Architecture

### Module Structure

The codebase is organized into platform-specific and sensor-specific modules:

- **`seaexplorer/`**: SeaExplorer glider data processing (primary platform)
- **`slocum/`**: Slocum glider data processing
- **`ad2cp/`**: AD2CP current profiler data processing
- **`microrider/`**: Microstructure sensor (MR1000G) data processing for turbulence
- **`acoustics/`**: Acoustic data processing and matched filtering
- **`plotting/`**: Shared visualization utilities using GMT and GLMakie
- **`common/`**: Shared utilities in the C2PO module

### Data Type System

The codebase uses mutable structs to represent glider data with array fields for time series. Key patterns:

**SeaExplorer types** (`seaexplorer/seaexplorerType.jl`):
- `NAV_RT`: Real-time navigation data with glider metadata
- `PLD_RT`: Payload data including all sensor measurements
- `SeaExplorerData`: Processed data structure with derived quantities (N2, spice, etc.)
- `SeaExplorerCTD`: CTD-focused data structure

**Slocum types** (`slocum/slocumType.jl`):
- `engStruct`: Engineering/navigation data
- `ctdStruct`: CTD sensor data with derived oceanographic quantities
- `sciStruct`: Generic science sensor data

**Common fields across platforms**:
- `glidertype`, `gliderSN`, `glidername`, `missionID`, `project`: Metadata
- `t::Array{Float64}`: Time in Unix timestamps
- `z::Array{Float64}`: Depth in meters
- `lon`, `lat::Array{Float64}`: Position
- `temp`, `salt`, `sigma0`, `spice0`, `sndspd`: CTD-derived oceanographic variables

### Configuration System

Mission parameters are defined in YAML files organized by project:
- `seaexplorer/seaexplorer_yaml_NORSE/`: NORSE project configurations
- `seaexplorer/seaexplorer_yaml_PASSENGERS/`: PASSENGERS project configurations
- `slocum/slocum_yaml_PASSENGERS/`: Slocum PASSENGERS missions
- `slocum/slocum_yaml_MARACOOS/`: MARACOOS project missions

**YAML structure**:
```yaml
gliderType: SeaExplorer
gliderName: SEA064
gliderSN: 64
missionID: 37
project: norse-janmayen
deploydate: 20221021
dataroot: /path/to/data/
suffix: complete
dataflag: delayed
```

### Processing Pipeline

1. **Configuration**: Load mission parameters from YAML files
2. **Data Loading**:
   - SeaExplorer: `seaexplorerFunc.seaexplorerYAMLload()` loads all missions from a directory
   - Slocum: `slocumLoad.slocumYAMLload()` for Slocum missions
3. **Caching**: Processed data saved as `.jld2` files for fast reloading
4. **Analysis**: Mission-specific scripts in `run_*.jl` files
5. **Visualization**: Plotting functions generate figures using GLMakie and GMT

### Key Functions

**SeaExplorer** (`seaexplorer/seaexplorerFunc.jl`):
- `seaexplorerYAMLload(dirpath)`: Load all missions from YAML directory
- `seaexplorer_load_mission(yamlfile)`: Load single mission
- `seaexplorer_process(data)`: Process and derive oceanographic quantities

**Slocum** (`slocum/slocumLoad.jl`):
- `slocumYAMLload(dirpath)`: Load Slocum missions from YAML directory
- `load_glider_ctd(...)`: Load CTD data from dbd/ebd files
- `glider_ctd_qc(data)`: Quality control filtering

**Plotting** (`plotting/gliderPlot.jl`):
- `plotGliderCTD(array, ps, pst)`: Multi-panel CTD section plots
- `plotGliderMap(array, pst, ...)`: Geographic maps with data overlay
- `plotGliderTSarray(array, ps, pst)`: Temperature-Salinity diagrams

**Utilities** (`common/C2PO.jl`):
- `missing2nan(data)`: Convert Julia `missing` values to `NaN`
- `datetime2yearday(dt)`: Convert DateTime to year-day format
- `yearday2datetime(year, yday)`: Inverse conversion
- `gc_distance(lat1, lon1, lat2, lon2)`: Great circle distance calculation

## Workflow Examples

### Processing a SeaExplorer Mission

```julia
# Load required modules
using Glider, JLD2
include("/path/to/seaexplorerFunc.jl")

# Load all missions from YAML directory
missiondir = "/path/to/seaexplorer_yaml_NORSE/"
gliderCTDarray = seaexplorerYAMLload(missiondir)

# Save processed data
jldsave("processed_data.jld2"; gliderCTDarray)

# Plot results
include("/path/to/gliderPlot.jl")
plotGliderCTD(gliderCTDarray, ps, pst)
```

### Processing a Slocum Mission

```julia
using Glider, JLD2
include("/path/to/slocumLoad.jl")

# Load from YAML directory
missiondir = "/path/to/slocum_yaml_PASSENGERS/"
gliderCTDarray = slocumYAMLload(missiondir)

# Plot
plotGliderCTD(gliderCTDarray, ps, pst)
```

## Data Conventions

### Time Handling
- **Primary format**: Unix timestamps (Float64) stored in `t` field
- **Conversion**: Use `unix2datetime()` and `datetime2unix()` from Dates
- **All times are UTC** - be careful with local time conversions
- Year-day format available via `C2PO.datetime2yearday()` and `C2PO.yearday2datetime()`

### Missing Data
- Julia `missing` values converted to `NaN` using `C2PO.missing2nan()`
- Use `NaNMath` functions for statistics with NaN values
- Use `skipmissing()` when working with arrays containing missing values

### Coordinate Systems
- **Longitude**: Decimal degrees, typically negative for Western hemisphere
- **Latitude**: Decimal degrees
- **Depth (z)**: Negative values (meters below surface)
- **Pressure (p)**: Positive values (decibars)

### Oceanographic Variables
Derived using GibbsSeaWater.jl (TEOS-10 standard):
- `salt`: Practical salinity
- `saltA`: Absolute salinity
- `sigma0`: Potential density anomaly referenced to surface
- `spice0`: Spiciness (orthogonal T-S index)
- `sndspd`: Sound speed
- `n2`: Buoyancy frequency squared (computed from density gradients)

## Project-Specific Information

### Active Projects
- **NORSE**: Arctic Ocean missions (Jan Mayen, Lofoten Basin) using SEA064
- **PASSENGERS**: Mid-Atlantic deployments using Slocum gliders
- **NESMA**: North East Shelf Monitoring Array using SEA064, SEA094
- **MARACOOS**: Mid-Atlantic Regional Association Coastal Ocean Observing System

### Typical Gliders
- **SEA064**: SeaExplorer SN 64
- **SEA094**: SeaExplorer SN 94
- Various Slocum gliders (electa, modena, etc.)

## Important Notes

### Path Dependencies
- All `include()` statements use absolute paths starting with `/Users/gong/`
- Data paths reference Dropbox: `/Users/gong/oceansensing Dropbox/C2PO/`
- These paths must be updated when working on different systems

### Module Loading
Scripts typically use:
```julia
workdir = "/full/path/to/module/directory"
if (workdir in LOAD_PATH) == false
    push!(LOAD_PATH, workdir)
end
```

### Reload Flag Pattern
Many scripts use a `reloadflag` variable:
- `reloadflag = true`: Reprocess data from raw files (slow)
- `reloadflag = false`: Load from cached `.jld2` files (fast)

### Memory Considerations
- AD2CP current profile data can be memory-intensive (4D arrays)
- Large missions may require chunked processing
- Use `pint` parameter to decimate data for plotting if needed

### File Organization
- **Scripts**: `run_*.jl` files are mission-specific entry points
- **Types**: `*Type.jl` files define data structures
- **Functions**: `*Func.jl` files contain processing logic
- **Plotting**: `*Plot.jl` files contain visualization functions
- **Testing**: `testing/` and `nonfunc/` directories contain experimental code
- **Retired**: `retired/` directories contain old versions for reference

## Plotting System

Uses `gliderPlotType` structs for configuration:

**plotSetting**: Controls plot appearance
- `pint`: Data decimation factor
- `iday`: Day interval for time axis
- `ms`: Marker size
- `pres`: Plot resolution (width, height)
- `fs`: Font size

**plotStruct**: Controls plot ranges
- Geographic bounds: `lonmin/max`, `latmin/max`
- Depth/pressure: `zmin/max`, `pmin/max`
- Variable ranges: `tempmin/max`, `saltmin/max`, etc.

Plots are saved to configured `figoutdir` paths.

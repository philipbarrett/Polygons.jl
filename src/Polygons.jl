module Polygons

using Gadfly
using Colors
import Base: *

export

    # Types
    Polygon

    # Methods
    ptsToDirs,
    dirsToPts,
    add,
    (*),
    weightedSum,
    gScan,
    crop,
    polyPlot

# TODO: Hausdorff distance

## Source files
include("common.jl")
include("simple.jl")
include("utils.jl")
include("polyPlot.jl")

end # module

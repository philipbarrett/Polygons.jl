module Polygons

using Gadfly
using Colors
using FixedSizeArrays
using CHull2d
import Base: *, .*, +, union

export
    # Types
    Polygon,

    # Methods
    acwOrder,
    ptsToDirs,
    dirsToPts,
    (*),
    (.*),
    (+),
    setSum,
    weightedSum,
    grahamScan,
    chull,
    crop,
    polyPlot,
    deeDoop,
    hausdorff

# TODO: Hausdorff distance

## Source files
include("common.jl")
include("simple.jl")
include("convex.jl")
include("utils.jl")
include("hausdorff.jl")
include("polyPlot.jl")

end # module

using BenchmarkTools
using CHull2D
using FixedSizeArrays
using PlotlyJS
using Polygons

# Create random points we want
function gen_random_points(N::Int, seed::Int=42)
    srand(seed)

    # Generate points
    x = randn(N, 2)

    # Create FixedSizeArrays.Point array
    x_pts = Point{2, Float64}[Point(x[i, 1], x[i, 2]) for i=1:N]

    return x, x_pts
end

N = 50_000
x, x_pts = gen_random_points(N, 5);

# Compute convex hull from Polygons
@time out = Polygons.grahamScan(x)
p_x = Polygon( pts = x )
@time out2 = Polygons.convexhull(p_x)

npts_polygons = size(out.pts, 1)
ep_polygons = out.pts
println("Polygons produced a convex hull with $npts_polygons points")

# Compute convex hull from CHull2D
ch = CHull2D.convexhull(x_pts)

npts_chull2d = length(ch.extremepoints)
ep_chull2d = Array(Float64, npts_chull2d, 2);
println("CHull2D produced a convex hull with $npts_chull2d points")

for i=1:npts_chull2d
    curr_pt = ch.extremepoints[i]
    ep_chull2d[i, 1] = curr_pt[1]
    ep_chull2d[i, 2] = curr_pt[2]
end

# Order them the same way
ep_polygons_s = sortrows(ep_polygons)
ep_chull2d_s = sortrows(ep_chull2d)

if npts_polygons == npts_chull2d
    maxdiff = maxabs(ep_chull2d_s - ep_polygons_s)
    println("The max difference in points is $maxdiff")
else
    println("Produced different number of points")
end

t1 = scatter(;x=ep_polygons[:, 1], y=ep_polygons[:, 2], mode="lines+markers")
t2 = scatter(;x=ep_chull2d[:, 1], y=ep_chull2d[:, 2], mode="lines+markers")

plot([t1, t2])

bPolygons = @benchmark Polygons.grahamScan(x)
bHybrid = @benchmark Polygons.grahamScan(p_x)
bCHull2D = @benchmark CHull2D.convexhull(x_pts)

println("\n\nBenchmarking with Polygons")
println(bPolygons)
println("\n\nBenchmarking with Hybrid conversion")
println(bPolygons)
println("\n\nBenchmarking with CHull2D")
println(bCHull2D)

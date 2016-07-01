#= convex.jl
Philip Barrett, pobarrett@gmail.com
01jul2016, Chicago

Tests the (two) convex hull routines to guarantee consistency
=#

Z = [ 1 -1; 1 1; -1 1; 0 -1 ; -1 -1 ]
    # Add a redundant point
a = Polygon( pts=Z ) ;
b=grahamScan(a)

c = Point{2, Float64}[Point(Z[i, 1], Z[i, 2]) for i=1:size(Z)[1]]
d = CHull2D.convexhull(c)

Y = [ 1 -1; 1 1; -1 1; 0 -1 + 1e-06 ; -1 -1 ]
e = Polygon( pts=Y ) ;
f = grahamScan(e)

g = Point{2, Float64}[Point(Y[i, 1], Y[i, 2]) for i=1:size(Y)[1]]
h = CHull2D.convexhull(g)

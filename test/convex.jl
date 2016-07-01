#= convex.jl
Philip Barrett, pobarrett@gmail.com
01jul2016, Chicago

Tests the (two) convex hull routines to guarantee consistency
=#

Z = [ 1 -1; 1 1; -1 1; 0 -1 ; -1 -1 ]
b = Polygon( pts=Z ) ;

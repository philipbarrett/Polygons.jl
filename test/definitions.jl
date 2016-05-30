#= definitions.jl
Philip Barrett, pobarrett@gmail.com
30may2016, Chicago

Tests the polygon defitions
=#

Z = [ 1 -1; 1 1; -1 1; -1 -1 ]
G = [ 1 0; 0 1; -1 0; 0 -1 ]
m = [ 1, 1, 1, 1 ]
a = Polygon( Z, G, m ) ;
b = Polygon( pts=Z ) ;
c = Polygon( dirs=G, dists=m ) ;

@test a == b
@test a == c

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

badZ = Z[ [1, 3, 4, 2], : ]
d = Polygon(pts=badZ)

for( compare in [b c d] )
    @test a.pts == compare.pts
    @test a.dirs == compare.dirs
    @test a.dists == compare.dists
end

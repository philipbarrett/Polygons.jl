#= definitions.jl
Philip Barrett, pobarrett@gmail.com
30may2016, Chicago

Tests the polygon definitions
=#

### NEED TO RE-CHECK ALL THESE TESTS

Z = [ -1 -1; 1.0 -1; 1 1; -1 1 ]::Matrix{Float64}
G = [ 0 -1; 1.0 0; 0 1; -1 0 ]
m = [ 1.0, 1, 1, 1 ]
a = Polygon( Z, G, m ) ;
b = Polygon( pts=Z ) ;
c = Polygon( dirs=G, dists=m ) ;

badZ = Z[ [1, 3, 4, 2], : ]
    # Mis-ordered Z
d = Polygon(pts=badZ)

intZ = vcat( Z, [0.0 0.5 ] )
    # Interior point added
e = Polygon(pts=intZ)
    # Works well enough

dupeZ = Z[ [1, 1, 2, 3, 4, 1, 2 ], : ]
    # Duplicate points
f = Polygon(pts=dupeZ)
    # Also good

nearZ = vcat( Z, Z[ [3,3,3,3,3],:] +
          [0 -1e-11 ; -1e-11 0 ; -1e-11 -1e-11 ; 0 -1e-11 ; -1e-11 0 ] )
g = Polygon(pts=nearZ)
    # Works nicely

for( compare in [b c d e f g] )
    @test a.pts == compare.pts
    @test a.dirs == compare.dirs
    @test a.dists == compare.dists
end

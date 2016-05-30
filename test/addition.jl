#= addition.jl
Philip Barrett, pobarrett@gmail.com
30may2016, Chicago

Tests the addition functions for the polygon class
=#

h = Polygon( pts=[ 0 1 ; 1 0 ; -1 0 ; 0 -1 ] )
exactdirs = [ a.dirs ; h.dirs ]
    # The union of the face directions.  Will give an exact result
exact = setSum( a, h, exactdirs )

@test exact.pts == [ 2 -1; 2 1; 1 2; -1 2; -2 1; -2 -1; -1 -2; 1 -2]

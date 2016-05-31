#= addition.jl
Philip Barrett, pobarrett@gmail.com
30may2016, Chicago

Tests the addition functions for the polygon class
=#

## Binary set sum
Z = [ 1 -1; 1 1; -1 1; -1 -1 ]
a = Polygon( pts=Z )
h = Polygon( pts=[ 0 1 ; 1 0 ; -1 0 ; 0 -1 ] )
exactdirs = [ a.dirs ; h.dirs; -1 -.1 ]
    # The union of the face directions.  Will give an exact result. The
    # extra search direction guarantees that the inner approx will be exact.
exact = setSum( a, h, exactdirs )
exactinner = setSum( a, h, exactdirs, false )

@test deeDoop(exact.pts) ≈ [ 2 -1; 2 1; 1 2; -1 2; -2 1; -2 -1; -1 -2; 1 -2]
@test deeDoop(exact.pts) ≈ exactinner.pts

## Array set sum
exactarray = setSum( [a, h], exactdirs )
@test exact.pts == exactarray.pts

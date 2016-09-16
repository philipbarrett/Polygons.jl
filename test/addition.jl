#= addition.jl
Philip Barrett, pobarrett@gmail.com
30may2016, Chicago

Tests the addition functions for the polygon class
=#

## Binary set sum
Z = [ 1.0 -1; 1 1; -1 1; -1 -1 ]
a = Polygon( pts=Z )
h = Polygon( pts=[ 0.0 1 ; 1 0 ; -1 0 ; 0 -1 ] )
NN = 40
exactdirs = vcat([ [ cos(2*(i-1)*pi/NN) sin(2*(i-1)*pi/NN) ] for i in 1:NN ]...)
    # The union of the face directions.  Will give an exact result.
exact = setSum( a, h, exactdirs )
exactinner = setSum( a, h, exactdirs, false )

@test deeDoop(exact.pts) ≈ [ -2 1; -2.0 -1; -1 -2; 1 -2; 2.0 -1; 2 1; 1 2; -1 2 ]
@test deeDoop(exact.pts) ≈ [ exactinner.pts[ end, :]; exactinner.pts[ 1:(end-1), :] ]
    # Arises because we have a vertical left edge.  Annoying.

## Array set sum
exactarray = setSum( [a, h], exactdirs )
@test exact.pts == exactarray.pts

## Scratch pad

using Polygons

P = [ .8 .2; .2 .8 ]
U = [ .5 -.5 ; -.5 .5 ]
bet = 0.9
vlow = ( 1 - bet ) * ( ( eye(2) - bet * P ) \ U )

N=40
dirs = zeros( N, 2 )
dirs[:,1] = [ cos(2 * pi * (i-1)/N) for i in 1:N ]
dirs[:,2] = [ sin(2 * pi * (i-1)/N) for i in 1:N ]

# Search directions
U1(s) = [ 2 -1; 3+s s ]
U2(s) = [ 2 3-s; -1 -s ]
# The period payoff functions

function T_unc( W::Array{ Polygon, 1 }, s, a, outer = true)
    # s is the value for the state and a is a vector containing the action indices
    u = vec( [ U1(s)[a] ; U2(s)[a] ] )
        # Period payoff
    s_idx = ( s == .5 ) ? 1 : 2
        # The index of s.  Required for selecting the right transition probability row
    return ( 1 - bet ) * u + bet * weightedSum( W, vec( P[ s_idx, :] ), dirs, outer )
end

winit1 = Polygon( pts = [ 2 2; -1 3.5; -3.5 1; -0.5 0.5 ; 0.5 -0.5 ] )
winit = [ winit1, winit1 ]
ff = T_unc( winit, .5, [1 1] )
tt = weightedSum( winit, vec( P[ 1, :] ), dirs )

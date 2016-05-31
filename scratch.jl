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

function T_unc( W::Array{ Polygon, 1 }, s::Float64, a::Vector, outer = true)
    # s is the value for the state and a is a vector containing the action indices
    u = vec( [ U1(s)[a[1],a[2]] ; U2(s)[a[1],a[2]] ] )
        # Period payoff
    s_idx = ( s == .5 ) ? 1 : 2
        # The index of s.  Required for selecting the right transition probability row
    return ( ( 1 - bet ) * u ) + ( bet * ( weightedSum( W, vec( P[ s_idx, :] ), dirs, outer ) ) )
end

function T_a( W::Array{ Polygon, 1 }, s::Float64, a::Vector, outer = true)
    # calcultes the incentive-compatible value set conditional on the action a
    unc = T_unc( W, s, a, outer )
        # The unconstrained set
    pddev = [ maximum( U1(s)[ :, a[2] ] ), maximum( U2(s)[ a[1], : ] ) ]
        # The deviating period payoffs
    s_idx = ( s == .5 ) ? 1 : 2
        # The index of s.  Required for selecting the right transition probability row
    dev = ( 1 - bet ) * pddev + bet * vlow[ s_idx, : ]'
        # The deviating payoff
    ic = crop( crop( unc, 1, dev[1] ), 2, dev[2] )
        # The incentive compatible set of payoffs\
    return ic
end

A_idx = [ 1 1; 1 2 ; 2 1; 2 2 ]
    # The matrix of action combinations

function T_operator( W::Array{ Polygon, 1 }, states::Vector, outer = true)
    # Forms the T operator
    out = [ union( [ T_a( W, s, vec(A_idx[j,:]))::Polygon for j in 1:size(A_idx)[1] ] )::Polygon for s in states ]
    return out
end

# polyPlot( [ winit1, T_operator( winit, [ .5 , -.5 ] ) ] )

winit1 = grahamScan( [ 2 2; -1 3.5; 3.5 -1; -0.5 0.5 ; 0.5 -0.5 ] )
winit = [ winit1, winit1 ]
aa = [ T_a( winit, -.5, vec(A_idx[j,:]))::Polygon for j in 1:size(A_idx)[1] ]
ff = union( aa )
polyPlot( [ aa, ff ])
polyPlot( aa[2] )

s = -.5
a = vec(A_idx[2,:])
nn=T_a( winit, s, a)::Polygon
unc = T_unc( winit, s, a)
unc2 = grahamScan(unc)
pddev = [ maximum( U1(s)[ :, a[2] ] ), maximum( U2(s)[ a[1], : ] ) ]
    # The deviating period payoffs
s_idx = ( s == .5 ) ? 1 : 2
    # The index of s.  Required for selecting the right transition probability row
dev = ( 1 - bet ) * pddev + bet * vlow[ s_idx, : ]'
    # The deviating payoff
ic = crop( crop( unc, 1, dev[1] ), 2, dev[2] )


polyPlot([unc, nn])

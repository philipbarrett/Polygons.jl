#= common.jl
Philip Barrett, pobarrett@gmail.com
26may2016, Chicago

Defines polygon type and methods for construction
=#

"""
    type polygon
Defines the polygon type with three entries: the vertices, the normals of the
sides, and the projection of each face onto the normal.
"""
immutable Polygon{N}
  pts::Mat{N,2,Float64}       # The points representation
  dirs::Mat{N,2,Float64}      # The directions representation
  dists::Vec{N,Float64}       # The associated distances
end



"""
    acwOrder( ptsDirs::Matrix, dists::Vector=[NaN] )
Orients either a set of points, or a a directions/distance pair anit-clockwise.
"""
function acwOrder( ptsDirs::Mat, dists::Vec=Vec{1, Float64}(NaN) )
    N = size(ptsDirs)[1]
        # Number of points/directions
    if( isnan( dists[1] ) )
      mu = ones(1,N) * ptsDirs / N
      centered = ptsDirs - ones(N) * mu
            # Relative to centre if doing points
            # Enhances numerical stability of orientation
    else
        centered = ptsDirs
    end
    cah = [ centered[i,1] / norm( [ centered[i,1], centered[i,2] ] ) for i in 1:N ]
        # The cosine
    ord = [ ( ( centered[i,2] >= 0 ) ? cah[i] : - 2 - cah[i] ) for i in 1:N ]
        # Order based on a decreasing transform of the cosine
    ord_idx = sortperm(ord, rev=true)
        # The sorting vector
    ptsDirsOut = Mat{N,2,Float64}(
      [ ptsDirs[ ord_idx[i], j ] for i in 1:N, j in 1:2 ] )
        # Reorder the opints

    if( isnan(dists[1]) )
        return ptsDirsOut
    else
        distsOut =
            Vec{N,Float64}( [ dists[ord_idx[i]] for i in 1:N ] )
        return ptsDirsOut, distsOut
    end
end

"""
    deeDoop( pts::Matrix )
De-duplicates a collection of points.  Must already be ordered anti-clockwise.
"""
function deeDoop( pts::Mat, tol::Float64=1e-10 )
  N = size(pts)[1]
      # Number of points
  neighbors = Mat{N,2,Float64}(
                [ [ pts[ i, j ]::Float64 for i in 2:N, j in 1:2 ];
                [ pts[1,1] pts[1,2] ] ] )
      # The neighboring points
  diffce = pts - neighbors
      # Difference between the two
  distce = [ norm( [ diffce[i,1], diffce[i,2] ] ) for i in 1:N ]
      # The vector of norms
  idx = distce .> tol
      # Location
  M = sum( idx )
      # Number of points to return
  idx_loc = find(idx)
      # The indices of the return points
  out = Mat{M,2,Float64}( [ pts[idx_loc[i],j] for i in 1:M, j in 1:2 ] )
      # Output
  return out
end

function deeDoop( dirs::Mat, dists::Vec, tol=1e-8 )
  N = size(dirs)[1]
      # Number of points
  neighbors = [ [ dirs[ i, j ]::Float64 for i in 2:N, j in 1:2 ];
                [ dirs[1,1] dirs[1,2] ] ]
      # The neighboring points
  diffce = dirs - neighbors
      # Difference between the two
  distce = [ norm( [ diffce[i,1], diffce[i,2] ] ) for i in 1:N ]
      # The vector of norms
  idx = distce .> tol
      # Location
  M = sum( idx )
      # Number of points to return
  idx_loc = find(idx)
      # The indices of the return points
  outdirs = Mat{M,2,Float64}( [ dirs[idx_loc[i],j] for i in 1:M, j in 1:2 ] )
  outdists = Vec{M,Float64}( [ dists[idx_loc[i]] for i in 1:M ] )
      # Output
  return outdirs, outdists
end

# function deeDoop( poly::Polygon, tol=1e-10 )
#   dirs, dists = deeDoop( poly.dirs, poly.dists )
#   return Polygon( dirs=dirs, dists=dists )
# end
#
#
"""
    polygon( ; pts=[ NaN NaN ], dirs=[ NaN NaN ], dists=[NaN]  )
Constructor for poylgon.  Removes duplicates and orders anti-clockwise.
"""
function Polygon( ; pts::Mat=Mat{1,2,Float64}(NaN),
                    dirs::Mat=Mat{1,2,Float64}(NaN),
                    dists::Vec=Vec{1,Float64}(NaN) )

    if( !isnan( pts[1] ) && isnan( dirs[1] ) )
        pts_clean = deeDoop( acwOrder( chull( pts ) ) )
        N = size(pts_clean)[1]
        pts_clean2 = Mat{N,2,Float64}(
            [ [ pts_clean[N,1] pts_clean[N,2] ] ;
              [ pts_clean[i,j] for i in 1:(N-1), j in 1:2 ] ] )
        dirs_ord, dists_ord = ptsToDirs( pts_clean2 )
    elseif( !isnan( dirs[1] ) && isnan( pts[1] ) )
        dirs_ord, dists_ord = acwOrder( dirs, dists )
    end
    dirs_final, dists_final = deeDoop( dirs_ord, dists_ord )
        # The ultimate test of dee-dooping is on the directions.
        # This is because we can have many points satisfying the same
        # directional constraint (i.e. on a straight line)
    N_final = size(dirs_final)[1]
    pts_final = dirsToPts( dirs_final, dists_final )
    return Polygon{N_final}( pts_final, dirs_final, dists_final )
end

"""
    dirsToPts( dirs::Matrix, dists::Vector )
Given a normal/distance description of a set, returns a set of points mZ which
lie at the vertices.  This assumes that the vector of normals is ordered
counter-clockwise already.
"""
function dirsToPts( dirs::Mat, dists::Vec )
  N = length(dists)
      # Number of points to be computed
  pts = zeros(N,2)
      # Initialize the output vector
  counter = 0
      # Counter: needed in cases where multiple boundary lines intersect at only
      # one point

  # Create the extended matrices (put the first case at the end)
  dirsExt = Mat{N+1,2,Float64}( [ [ dirs[i,j] for i in 1:N, j in 1:2 ];
                                  [dirs[1,1] dirs[1,2] ] ] )
  distsExt = Vec{N+1, Float64}( [ [ dists[i] for i in 1:N ] ; dists[1] ] )

  # Main loop: Over all points
  for i in 1:N
    A = Mat{2,2,Float64}( [ dirsExt[i+k,j] for k in 0:1, j in 1:2 ] )
        # The directional vectors to be intersected
    b = Vec{2,Float64}( [ distsExt[i], distsExt[i+1] ] )
        # And their distances
    sol = inv(A) * b
        # The candidate solution

    # Now compute the error: Detects multiple intersections
    err = norm( A * sol - b )
        # Inversion error
    relerr = ( norm(b) == 0 ) ? err : err / norm(b)
        # Relative error (where possible)

    # If there is no error, then assign to the output
    if( relerr < 1e-14 )
      counter += 1
      pts[counter,1] = sol[1]
      pts[counter,2] = sol[2]
    end
  end

  out = Mat{counter, 2, Float64}([ pts[counter,:] ; pts[1:(counter-1),:] ] )
      # Formulate the output
  return out
end

"""
    ptsToDirs( pts::Matrix )
Given a set of points, returns the normal-distance representation of the edges
Assumes that the points are ordered anti-clockwise already.
"""
function ptsToDirs( pts::Mat )
  N = size(pts)[1]
      # Number of points
  flip = Mat{2,2, Float64}([ 0.0 1.0; -1.0 0.0 ])
      # Matrix to convert grdient to normal
  neighbors = Mat{N,2,Float64}(
                [ [ pts[ i, j ]::Float64 for i in 2:N, j in 1:2 ];
                [ pts[1,1] pts[1,2] ] ] )
      # The neighboring points
  unscaled = ( pts - neighbors ) * flip
      # The unscaled vaules
  lens = [ norm( [unscaled[i,1], unscaled[i,2] ] ) for i in 1:N ]
      # Lengths of the unscaled directions
  dirs = Mat{N,2,Float64}( [ unscaled[i,j] / lens[i] for i in 1:N, j in 1:2] )
      # The directions
  dists = (dirs .* pts) * Vec{2,Float64}(1)
      # The distances in each direction
  return dirs, dists
end

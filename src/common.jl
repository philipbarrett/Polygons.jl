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
type Polygon
  pts::Matrix{Float64}     # The points representation
  dirs::Matrix{Float64}    # The directions representation
  dists::Vector{Float64}   # The distances associated with the vectors
end

"""
    acwOrder( ptsDirs::Matrix, dists::Vector=[NaN] )
Orients either a set of points, or a a directions/distance pair anit-clockwise.
"""
function acwOrder( ptsDirs::Matrix, dists::Vector=[NaN] )
    N = size(ptsDirs)[1]
        # Number of poits/directions
    cah = zeros( N )
        # Initialize the cosine function
    if( isnan( dists[1] ) )
        centered = ptsDirs - ones(size(ptsDirs)[1]) * mean( ptsDirs, 1 )
            # Relative to centre if doing points
            # Enhances numerical stability of orientation
    else
        centered = ptsDirs
    end
    cah = [ centered[i,1] / norm( centered[i,:] ) for i in 1:N ]
        # The cosine
    ord = [ ( ( centered[i,2] >= 0 ) ? cah[i] : - 2 - cah[i] ) for i in 1:N ]
        # Order based on a decreasing transform of the cosine
    ord_idx = sortperm(ord, rev=true)
        # The sorting vector
    ptsDirsOut::Matrix{Float64} = ptsDirs[ ord_idx, : ]
        # Reorder the opints

    if( isnan(dists[1]) )
        return ptsDirsOut
    else
        distsOut::Vector{Float64} = dists[ord_idx]
        return ptsDirsOut, distsOut
    end
end

"""
    deeDoop( pts::Matrix )
De-duplicates a collection of points.  Must already be ordered anti-clockwise.
"""
function deeDoop( pts::Matrix, tol=1e-10 )
  N = size(pts)[1]
      # Number of points
  neighbors = [ pts[ 2:end,: ]; pts[1,:] ]
      # The neighboring points
  diff = [ norm( pts[i,:] - neighbors[i,:] ) for i in 1:N ]
      # The vector of norms
  return pts[ diff .> tol, : ]
end

function deeDoop( dirs::Matrix, dists::Vector, tol=1e-10 )
  N = size(dirs)[1]
      # Number of points
  neighbors = [ dirs[ 2:end,: ]; dirs[1,:] ]
      # The neighboring distance vectors
  diff = [ norm( dirs[i,:] - neighbors[i,:] ) for i in 1:N ]
      # The vector of norms
  return dirs[ diff .> tol, : ], dists[ diff .> tol ]
end

function deeDoop( poly::Polygon, tol=1e-10 )
  return Polygon( pts=deeDoop( poly.pts ) )
end


"""
    polygon( ; pts=[ NaN NaN ], dirs=[ NaN NaN ], dists=[NaN]  )
Constructor for poylgon.  Removes duplicates and orders anti-clockwise.
When setting points via pts, must already be a convex hull.
"""
function Polygon( ; pts=[ NaN NaN ], dirs=[ NaN NaN ], dists=[NaN]  )

    if( !isnan( pts[1] ) && isnan( dirs[1] ) )
        pts = deeDoop( acwOrder( pts ) )
        pts = [ pts[ end, : ]; pts[ 1:(end-1), : ] ]
        dirs, dists = ptsToDirs( pts )
    elseif( !isnan( dirs[1] ) && isnan( pts[1] ) )
        dirs, dists = acwOrder( dirs, dists )
        dirs, dists = deeDoop( dirs, dists )
        pts = dirsToPts( dirs, dists )
    end
    return Polygon( pts, dirs, dists )
end

"""
    dirsToPts( dirs::Matrix, dists::Vector )
Given a normal/distance description of a set, returns a set of points mZ which
lie at the vertices.  This assumes that the vector of normals is ordered
counter-clockwise already.
"""
function dirsToPts( dirs::Matrix, dists::Vector )
  nPts = length(dists)
      # Number of points to be computed
  pts = zeros(Float64, 2, nPts)
      # Initialize the (transpose) output vector
  counter = 1
      # Counter: needed in cases where multiple boundary lines intersect at only
      # one point

  # Create the extended matrices (put the first case at the end)
  dirsExt = [ dirs ; dirs[1,:] ]
  distsExt = [ dists; dists[1,:] ]

  # Main loop: Over all points
  for i in 1:nPts
    A = dirsExt[i:(i+1),:]
        # The directional vectors to be intersected
    b = distsExt[i:(i+1)]
        # And their distances
    sol = A \ b
        # The candidate solution

    # Now compute the error: Detects multiple intersections
    err = norm( A * sol - b )
        # Inversion error
    relerr = ( norm(b) == 0 ) ? err : err / norm(b)
        # Relative error (where possible)

    # If there is no error, then assign to the output
    if( relerr < 1e-14 )
      pts[:,counter] = sol
      counter += 1
    end
  end
  out = transpose(pts)
  return( [ out[ nPts, : ] ; out[ 1:(nPts-1), : ] ] )
end

"""
    ptsToDirs( pts::Matrix )
Given a set of points, returns the normal-distance representation of the edges
Assumes that the points are ordered anti-clockwise already.
"""
function ptsToDirs( pts::Matrix )
  nPts = size(pts)[1]
      # Number of points
  flip = [ 0 1; -1 0 ]
      # Matrix to convert grdient to normal
  neighbors = [ pts[ 2:nPts, : ] ; pts[1,:] ]
      # The neighboring points
  unscaled = ( pts - neighbors ) * flip
      # The unscaled vaules
  dirs::Matrix{Float64} = unscaled ./ sqrt( sum( unscaled .* unscaled, 2 ) )
      # The directions
  dists::Vector{Float64} = vec(sum( pts .* dirs, 2))
      # The distances in each direction
  return dirs, dists
end

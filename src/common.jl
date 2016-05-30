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
    polygon( ; pts=[ NaN NaN ], dirs=[ NaN NaN ], dists=[NaN]  )
Constructor for poylgon.  Maitained assumptions are that dirs is ordered
clockwise and that pts is already a convex hull.
"""
function Polygon( ; pts=[ NaN NaN ], dirs=[ NaN NaN ], dists=[NaN]  )

    # TODO: Check for oreintation of the points and normals (each a separate function)

    if( !isnan( pts[1] ) && isnan( dirs[1] ) )
        dirs, dists = ptsToDirs( pts )
    elseif( !isnan( dirs[1] ) && isnan( pts[1] ) )
        pts = dirsToPts( dirs, dists )
    end
    return Polygon( pts, dirs, dists )
end

"""
    dirsToPts( dirs::Matrix, dists::Vector )
Given a normal/distance description of a set, returns a set of points mZ which
lie at the vertices.  This assumes that the vector of normals is ordered
clockwise already.
"""
# TODO: Check that this works for anti-clockwise orientation

function dirsToPts( dirs::Matrix{Float64}, dists::Vector{Float64} )
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
  return( transpose(pts) )
end

"""
    ptsToDirs( pts::Matrix )
Given a set of points, returns the normal-distance representation of the edges
Assumes that the points are ordered anti-clockwise already.
"""
# TODO: Make sure that this works for anti-clockwise orderd points

function ptsToDirs( pts::Matrix{Float64} )

  nPts = size(pts)[1]
      # Number of points
  flip = [ 0 1; -1 0 ]
      # Matrix to convert grdient to normal
  neighbors = [ pts[ 2:nPts, : ] ; pts[1,:] ]
      # The neighboring points
  unscaled = ( neighbors - pts ) * flip
      # The unscaled vaules
  dirs::Matrix{Float64} = unscaled ./ sqrt( sum( neighbors .* neighbors, 2 ) )
      # The directions
  dists::Vector{Float64} = vec(sum( pts .* dirs, 2))
      # The distances in each direction
  return dirs, dists
end

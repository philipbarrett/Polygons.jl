#= convex.jl
Philip Barrett, pobarrett@gmail.com
27may2016, Chicago

Calculates the convex hull of a set of points
=#

"""
    acw( p1, p2, p3 )
Positive if p1 -> p2 -> p3 is anti-clockwise
"""
function acw( p1, p2, p3 )
  return (p2[1] - p1[1]) * (p3[2] - p1[2]) - (p2[2] - p1[2]) * (p3[1] - p1[1])
end

"""
    grahamScan( pts::Matrix )
Computes a convex hull using the Graham Scan algortithm NB: Returns a polygon object even when a matrix is submitted
"""

function grahamScan( pts::Matrix )

  N = size(pts)[1]
      # Number of points
  lex = sortrows( pts[ :, [2, 1 ] ] )[ :, [ 2, 1 ] ]
      # Sorts the points lexicographically by y then x
  P = lex[ 1, : ]
      # The initial point
  otherpts = lex[ 2:end, : ]
      # The other points
  orderpts = [ P; otherpts ; P ]
      # Initialize the anti-clockwise ordered points
  cah = zeros( N - 1 )
      # The cosine of the angle between each point and P.  Will be in
      # [-1,1] because of choice of P

  ## Order the points anti-clockwise ##
  for( i in 1:(N-1) )
    diff = otherpts[i,:] - P
    cah[i] = diff[1] / norm(diff)
        # Cosine is adjacent over hypotenuse
  end
  orderpts[2:N,:] = otherpts[ sortperm(cah, rev=true), : ]
      # Order the points by the angle measure

  ## Create the output ##
  out = similar( pts )
  out[ 1:2, : ] = orderpts[1:2,:]
      # Initialize
  M = 2
      # Counts number of rows in convex hull
  for( i in 3:N )
    while( acw( out[ M-1, : ], out[ M, : ], orderpts[ i, : ] ) < 0 )
      if( M > 2 )
        M -= 1
      elseif( i == N )
        break
      else
        i += 1
      end
    end
    M += 1
        # Increment counter
    out[ M, : ] = orderpts[ i, : ]
        # Add to the output if we have an acw angle
  end
  return Polygon( pts=out[ 1:M, : ] )
end

function grahamScan( poly::Polygon )
    return grahamScan( poly.pts )
end

"""
    convexhull( poly::Polygon )
Computes a convex hull using the Graham Scan algortithm from the CHull2D package
"""
function chull( poly::Polygon )
  N = size(poly.pts)[1]
      # Number of points
  oldpts = Point{2, Float64}[Point(poly.pts[i, 1], poly.pts[i, 2])
                            for i=1:N]
      # Convert poly.pts to an array of points
  ch = CHull2D.convexhull(oldpts)
      # The convex hull object
  N_newpts = length(ch.extremepoints)
  newpts = Array(Float64, N_newpts, 2)
      # Create the container for the new points
  for i=1:N_newpts
      curr_pt = ch.extremepoints[i]
      newpts[i, 1] = curr_pt[1]
      newpts[i, 2] = curr_pt[2]
  end

  return Polygon( pts = newpts )
end

#= utils.jl
Philip Barrett, pobarrett@gmail.com
27may2016, Chicago

Defines utilities for Polygon objects
=#

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

function deeDoop( poly::Polygon, tol=1e-10 )
  return deeDoop( poly.pts )
end

"""
    crop( poly::Polygon, dim::Int, dist, upper=true )
Returns a Polygon cropped in either the x or y dimension, with upper = true
retaining the part of the Polygon above the chop, and if false the part below.
The Polygon should already be oreinted counterclockwise
"""
## TODO: DEBUG ME!
## Strategy here is all wrong.  Better: Identify all the points that are to the
## left (right) of the slice.  Then throw away them (and their associated
## normals).  Finally, reintroduce the slice boundary into the dirs/dist
## framework
function crop( poly::Polygon, dim::Int, dist, upper=true )
  N = size(poly.pts)[1]

  # Compute the direction of the chop
  if( dim == 1 && upper )
    dir = [ -1, 0 ]
    dist = - dist
  elseif( dim == 1 && !upper )
    dir = [ 1, 0 ]
  elseif( dim==2 && upper )
    dir = [ 0, -1 ]
    dist = - dist
  else
    dir = [ 0, 1 ]
  end
  c = dir[1]
  s = dir[2]
      # Sign and cosine of the dir vector

  # The chop doesn't bind
  l = minimum( poly.pts[:,dim] )
  u = maximum( poly.pts[:,dim] )
  if( sum( dist * dir ) < l || sum( dist * dir ) > u )
    return(poly)
  end
  if( sum( dist * dir ) < u || sum( dist * dir ) > l )
    return(Polygon())
  end

  # Find the points which are to be chopped
  # init = 1, term = 1
  #     # Initial and terminal chopping points
  # for( i in 1:N )
  #   if(  )
  # end

#
#
#
#
#
#
#
#
#
#
#
#
#   # Now loop over the direction vectors
#   i = 1
#       # Counter
#   c1 = poly.dirs[1,1] / norm( poly.dirs[1,:] )
#   s1 = poly.dirs[1,2] / norm( poly.dirs[1,:] )
#       # Initiate the next point's sign and cos
#   while( i<(N-1) )
#     c0 = c1
#     s0 = c1
#         # Old is now new. How sad, how true.
#     c1 = poly.dirs[i+1,1] / norm( poly.dirs[i+1,:] )
#     s1 = poly.dirs[i+1,2] / norm( poly.dirs[i+1,:] )
#         # Cosine and sine of the next direction vector
#     if( min( s1, s0 ) > 0  )
#       if( c0 >= c > c1 )
#         break
#       end
#     elseif( max( s1, s0 ) < 0 )
#       if( c0 <= c < c1 )
#         break
#       end
#     elseif( s0 > 0 > s1 )
#       if( c <= min( c0, c1 ) )
#         break
#       end
#     else
#       if( c >= max( c0, c1 ))
#         break
#       end
#     end
#     i += 1
#   end
#       # So dir fits in between the ith and (i+1)th direction vectors
#   newdirs = [ poly.dirs[1:i, :] ; dir' ; poly.dirs[ (i+1):end, : ] ]
#   newdists = [ poly.dists[1:i] ; dist ; poly.dists[ (i+1):end ] ]
#       # Create the new directions and distances
#
#   ## Special cases for exact dir in poly.dirs ##
#   if( poly.dirs[i,:] == dir)
#     newdirs = poly.dirs
#     newdists = [ poly.dists[1:(i-1)] ; dist ; poly.dists[ (i+1):end ] ]
#   end
#   if( poly.dirs[i+1,:] == dir)
#     newdirs = poly.dirs
#     newdists = [ poly.dists[1:i] ; dist ; poly.dists[ (i+2):end ] ]
#   end
#
# println( "i=", i )
# println( "newdirs\n", newdirs )
# println( "newdists\n", newdists )
#
#   return( Polygon( dirs = newdirs, dists = newdists ) )

end

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

  # The crop doesn't bind
  l = minimum( poly.pts[:,dim] )
  u = maximum( poly.pts[:,dim] )
  if( upper )
      if( dist < l )
          return poly
      end
      if( dist > u )
          return Polygon()
      end
  else
      if( dist < l )
          return Polygon()
      end
      if( dist > u )
          return poly
      end
  end

  # Compute the direction of the crop
  if( dim == 1 && upper )
    dir = [ -1.0, 0.0 ]
    dist = - dist
  elseif( dim == 1 && !upper )
    dir = [ 1.0, 0.0 ]
  elseif( dim==2 && upper )
    dir = [ 0.0, -1.0 ]
    dist = - dist
  else
    dir = [ 0.0, 1.0 ]
  end
  c = dir[1]
  s = dir[2]
      # Sign and cosine of the dir vector

  # Find the points which are to be cropped
  init = 0      # Last point to be retained
  term = 0      # Last point to be cropped
      # Initial and terminal cropping points
  for( i in 1:N )
      if( ( poly.pts[i,:] * dir )[1] < dist )
          if( 0 < init < i-1 && term > 0 )
              break
          end
          init = ( init == 0 ) ? term + 1 : init + 1
      else
          if( 0 < term < i-1 && init > 0 )
              break
          end
          term = ( term == 0 ) ? init + 1 : term + 1
            # Increment the terminal counter if non-zero

      end
  end

  if( init < term )
    newdirs = [ poly.dirs[ 1:init, : ] ;
                dir' ; poly.dirs[ term:end, : ] ]
    newdists = [ poly.dists[ 1:init ]; dist;
                poly.dists[ term:end ] ]
  else
    newdirs = [ dir' ; poly.dirs[ term:init, : ] ]
    newdists = [ dist; poly.dists[ term:init ] ]
  end

  return Polygon( dirs = newdirs, dists = newdists )

end

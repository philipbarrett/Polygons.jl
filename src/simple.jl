#= simple.jl
Philip Barrett, pobarrett@gmail.com
26may2016, Chicago

Defines basic polygon operations:
    * Vector Addition
    * Setsum (of both polygon pairs and arrays)
    *
=#

"""
    add( poly::polygon, pt::Vector )
Adds a vector to all points in a polygon
"""
function add( poly::polygon, u::Vector )

  pts = similar( poly.pts )
  pts[:,1] = poly.pts[:,1] + u[1]
  pts[:,2] = poly.pts[:,2] + u[2]
      # Shift the points
  return polygon( pts = pts )
end

"""
    add( poly1::polygon, poly2::polygon, dirs, outer=true )
Addition functionality.  Provides either an inner or out approximation
"""
function add( poly1::polygon, poly2::polygon, dirs, outer=true )
  dists1 = maximum( poly1.pts * dirs', 1 )
  dists2 = maximum( poly2.pts * dirs', 1 )
      # The distances in each direction
  dists = vec( dists1 + dists2 )

  if( outer )
    return polygon( dirs=dirs, dists=dists )
  end
  return polygon( pts=dirsToPts( dirs, dists ) )
end

"""
    add( poly1::Array{polygon,1}, poly2::Array{polygon,1}, dirs, outer=true )
Adds an array of polygons
"""
function add( polys::Array{polygon,1}, dirs, outer=true )

  N = length(polys)
      # Number of polygons to add
  dists = zeros( size(dirs)[1] )
      # Initate the distances
  for( i in 1:N )
    dists += maximum( polys[i].pts * dirs', 1 )'
        # Sum the distances
  end

  if( outer )
    return polygon( dirs=dirs, dists=vec( dists ) )
  end
  return polygon( pts=dirsToPts( dirs, vec( dists ) ) )
end

"""
    times( poly1::Array{polygon,1}, poly2::Array{polygon,1}, dirs, outer=true )
Scalar multiplication
"""
function (*)( k::Number, poly::polygon )
  return polygon( k * poly.pts, poly.dirs, k * poly.dists )
end

"""
    weightedSum( polys::Array{polygon,1}, wts::Vector dirs, outer=true )
Computes a weighted sum of polygons
"""
function wtdSum( polys::Array{polygon,1}, wts::Vector, dirs, outer=true )
  N = length( polys )
      # Number of polygons
  polys2 = polys
      # Initiate scaled polygons with input array
  for( i in 1:N )
    polys2[i] = wts[i] * polys[i]
  end
  return( add( polys2, dirs, outer ) )
end

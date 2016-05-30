#= simple.jl
Philip Barrett, pobarrett@gmail.com
26may2016, Chicago

Defines basic Polygon operations:
    * Vector Addition
    * SetsetSum (of both Polygon pairs and arrays)
    * Scalar multiplication
    * Weighted setSum
=#

"""
    add( poly::Polygon, pt::Vector )
Adds a vector to all points in a Polygon
"""
function (+)( poly::Polygon, u::Vector )

  pts = similar( poly.pts )
  pts[:,1] = poly.pts[:,1] + u[1]
  pts[:,2] = poly.pts[:,2] + u[2]
      # Shift the points
  return Polygon( pts = pts )
end

"""
    setSum( poly1::Polygon, poly2::Polygon, dirs, outer=true )
Addition functionality.  Provides either an inner or out approximation
"""
function setSum( poly1::Polygon, poly2::Polygon, dirs, outer=true )
  dists1 = maximum( poly1.pts * dirs', 1 )
  dists2 = maximum( poly2.pts * dirs', 1 )
      # The distances in each direction
  dists = vec( dists1 + dists2 )

  if( outer )
    return Polygon( dirs=dirs, dists=dists )
  end
  return Polygon( pts=dirsToPts( dirs, dists ) )
end

"""
    setSum( poly1::Array{Polygon,1}, poly2::Array{Polygon,1}, dirs, outer=true )
Adds an array of Polygons
"""
function setSum( polys::Array{Polygon,1}, dirs, outer=true )

  N = length(polys)
      # Number of Polygons to add
  dists = zeros( size(dirs)[1] )
      # Initate the distances
  for( i in 1:N )
    dists += maximum( polys[i].pts * dirs', 1 )'
        # Sum the distances
  end

  if( outer )
    return Polygon( dirs=dirs, dists=vec( dists ) )
  end
  return Polygon( pts=dirsToPts( dirs, vec( dists ) ) )
end

"""
    times( poly1::Array{Polygon,1}, poly2::Array{Polygon,1}, dirs, outer=true )
Scalar multiplication
"""
function (*)( k::Number, poly::Polygon )
  return Polygon( k * poly.pts, poly.dirs, k * poly.dists )
end

"""
    weightedSum( polys::Array{Polygon,1}, wts::Vector dirs, outer=true )
Computes a weighted setSum of Polygons
"""
function wtdSum( polys::Array{Polygon,1}, wts::Vector, dirs, outer=true )
  N = length( polys )
      # Number of Polygons
  polys2 = polys
      # Initiate scaled Polygons with input array
  for( i in 1:N )
    polys2[i] = wts[i] * polys[i]
  end
  return( setSum( polys2, dirs, outer ) )
end

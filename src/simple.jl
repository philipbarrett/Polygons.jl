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
  return Polygon( pts = poly.pts + ones(size(poly.pts)[1]) * u' )
end

function (+)( u::Vector, poly::Polygon )
  return poly + u
end

"""
    setSum( poly1::Polygon, poly2::Polygon, dirs, outer=true )
Addition functionality.  Provides either an inner or out approximation
"""
# TODO: Improve the inner approximation when the serach direction is exactly
# a face normal
function setSum( poly1::Polygon, poly2::Polygon, dirs::Matrix, outer=true )
  dists1 = poly1.pts * dirs'
  dists2 = poly2.pts * dirs'
      # The distances in each direction
  maxdists1 = maximum( dists1, 1 )
  maxdists2 = maximum( dists2, 1 )
      # The maximum distances
  maxdists = vec( maxdists1 + maxdists2 )

  if( outer )
    return Polygon( dirs=dirs, dists=maxdists )
  end

  Z1 = zeros( dirs )
  Z2 = zeros( dirs )
        # Initiate output points

  for( i in 1:size(dirs)[1])
      Z1[i,:] = poly1.pts[ dists1[:,i] .== maxdists1[i], : ][1,:]
      Z2[i,:] = poly2.pts[ dists2[:,i] .== maxdists2[i], : ][1,:]
  end
  return Polygon( pts = Z1 + Z2 )
end

"""
    setSum( poly1::Array{Polygon,1}, poly2::Array{Polygon,1}, dirs, outer=true )
Adds an array of Polygons
"""
function setSum( polys::Array{Polygon,1}, dirs::Matrix, outer=true )

  N = length(polys)
      # Number of Polygons to add
  dists = zeros( size(dirs)[1] )
      # Initate the distances
  innerpts = zeros(dirs)
      # Initialize the points for the inner approximation
  for( i in 1:N )
    thisdists = polys[i].pts * dirs'
    maxdists = maximum( thisdists, 1 )'
    dists += maximum( polys[i].pts * dirs', 1 )'
        # Sum the distances
    if( !outer )
      for( j in 1:size(dirs)[1])
        innerpts[j,:] += polys[i].pts[ thisdists[:,j] .== maxdists[j], : ][1,:]
      end
    end
  end

  if( outer )
    return Polygon( dirs=dirs, dists=vec( dists ) )
  end
  return Polygon( pts=innerpts )
end

"""
    times( poly1::Array{Polygon,1}, poly2::Array{Polygon,1}, dirs, outer=true )
Scalar multiplication
"""
function (*)( k::Number, poly::Polygon )
  return Polygon( pts = k * poly.pts )
end

"""
    weightedSum( polys::Array{Polygon,1}, wts::Vector dirs, outer=true )
Computes a weighted setSum of Polygons
"""
function weightedSum( polys::Array{Polygon,1}, wts::Vector, dirs::Matrix, outer=true )
  N = length( polys )
      # Number of Polygons
  polys2 = similar(polys)
      # Initiate scaled Polygons with input array
  for( i in 1:N )
    polys2[i] = wts[i] * polys[i]
  end
  return( setSum( polys2, dirs, outer ) )
end
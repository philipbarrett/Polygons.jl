#= hausdorff.jl
Philip Barrett, pobarrett@gmail.com
30may2016, Chicago

Computes the hausdorff distance between polygons (and points)
=#

"""
    hausdorff

Computes the hausdorff distance between two polygons, or two arrays of polygons
"""
function hausdorff( poly1::Polygon, poly2::Polygon )
  return hausdorffCore( poly1.pts, poly2 )
end

function hausdorff( poly1::Array{Polygon,1}, poly2::Array{Polygon,1} )
  N = length(poly1)
  length(poly2) == N || error("Polygon arrays must have same diemnsion")
      # Error handling
  return [ hausdorffCore( poly1[i].pts, poly2[i] ) for i in 1:N ]
end

function hausdorff( poly1::Array{Polygon,2}, poly2::Array{Polygon,2} )
  N, M = size(poly1)
  size(poly2) == (N, M) || error("Polygon arrays must have same diemnsion")
      # Error handling
  return [ hausdorffCore( poly1[i,j].pts, poly2[i,j] ) for i in 1:N, j in 1:M ]
end
"""
    function hausdorffCore( newpts::Matrix, poly::Polygon  )

Core part of any hausdorff distance routine.  Measures the distance between a matrix of points newpts and a polygon poly
"""
function hausdorffCore( newpts::Matrix, poly::Polygon )

  pts, dirs, dists = poly.pts, poly.dirs, poly.dists
      # Unpack the polygon object
  N = size(newpts)[1]
      # number of points
  dline, dvert, newdist, hausdorff = 0.0, 0.0, 0.0, 0.0
      # initialize the distance measures
  for i in 1:N
    dline = setsNearestLine( vec(newpts[i,:]), dirs, dists )
        # Distance to nearest line
    dvert = setsNearestVert( vec(newpts[i,:]), pts )
        # Distance to nearest vertex
    newdist = min( dline, dvert )
        # Minimum distance
    hausdorff  = max( hausdorff, newdist )
        # Update hausdorff dist via minmax
  end
  return hausdorff
end

function setsNearestLine( pt::Vector, dirs::Matrix, dists::Vector )

  N::Int = length(dists)
      # Number of sides
  dff = dists - dirs * pt
  absdiff = abs( dff )
      # Difference between the projection directions and the boundaries (and absolute value)
  nearpts = ones( N ) * pt' + dirs .* ( dff * ones( 1, 2 ) )
      # Recover the closest points on the normal surfaces
  idx = setsInsidePtsIdx( nearpts, dirs, dists )
      # Test whether these points are inside the Polygon
  if !any(idx)
    return Inf
        # Outsize return if none of the nearest points are in the polygon
  end
  distinside = [ idx[i] ? absdiff[i] : Inf for i in 1:N ]
      # Compute the minimum distances to the points on the lines that are
      # inside the polygon
  return minimum( distinside )
      # Return the minimum
end

function setsInsidePtsIdx( cands::Matrix, dirs::Matrix,
                            dists::Vector )
  N::Int = size(cands)[1]
      # Number of candidate points
  D = dirs * cands'  - dists * ones( 1, N )
      # Negative if inside the set
  return all( D .<= 0, 1 )
end

function setsNearestVert( cand::Vector, pts::Matrix )
  N = size(pts)[1]
      # Nubmer of points
  return minimum( [ norm( cand' - pts[i,:] ) for i in 1:N ] )
end

# double sets_nearest_vert_cpp_eigen( VectorXd pt, MatrixXd mB ){
# // C++ translation of sets.nearest.vert
#
#   int iB = mB.rows() ;
#         // Number of rows in B
#   VectorXd vDists( iB ) ;
#         // Initialize the distance vector
#   for( int iII=0; iII<iB; iII++ ){
#     vDists[ iII ] = ( pt.transpose() - mB.row( iII ) ).norm() ;
#           // Compute the distance to pt
#   }
#   return vDists.minCoeff() ;
# }

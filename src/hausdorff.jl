#= hausdorff.jl
Philip Barrett, pobarrett@gmail.com
30may2016, Chicago

Computes the hausdorff distance between polygons (and points)
=#

"""
    function hausdorffCore( newpts::Matrix, poly::Polygon  )

Core part of any hausdorff distance routine.  Measures the distance between a matrix of points newpts and a polygon poly
"""
function hausdorffCore( newpts::Matrix, poly::Polygon )

  pts, dirs, dists = poly.pts, poly.dirs, poly.dists
      # Unpack the polygon object


end



# //[[Rcpp::export('sets.Hausdorff.core')]]
# double sets_Hausdorrf_core_cpp( NumericMatrix mA, NumericMatrix mZ, NumericMatrix mG, NumericVector vC ){
# // Computes the Hausdorf distance between points in mA and the set defined by mZ-mG-vC
#
#   // Eigen Conversion
#   const eMatd emA( as<eMatd>( mA ) ) ;
#   const eMatd emZ( as<eMatd>( mZ ) ) ;
#   const eMatd emG( as<eMatd>( mG ) ) ;
#   const eVecd evC( as<eVecd>( vC ) ) ;
#
#   int iRows = emA.rows() ;
#         // The number of points in mA
#   double dLine = 0 ;
#   double dVert = 0 ;
#   double dNewDist = 0;
#   double dHausdorff = 0 ;
#         // Initializing the distance variables
#
#   // Loop over points in mA
#   for( int iII=0; iII<iRows; iII++ ){
#     dLine = sets_nearest_line_cpp_eigen( emA.row( iII ), emG, evC ) ;
#           // Distance to nearest line
#     dVert = sets_nearest_vert_cpp_eigen( emA.row( iII ), emZ ) ;
#           // Distance to nearest vertex
#     dNewDist = ( dLine < dVert ? dLine : dVert ) ;
#           // The minimum distance
#     dHausdorff = ( dNewDist > dHausdorff ? dNewDist : dHausdorff ) ;
#           // Update Hausdorff with a minmax
#   }
#
#   return dHausdorff ;
# }

function setsNearestLine( pt::Vector, dirs::Matrix, dists::Vector )

  N::Int = length(dists)
      # Number of sides
  diff::Vector = dists - dirs * pt
  absdiff::Vector = abs( diff )
      # Difference between the projection directions and the boundaries (and absolute value)
  nearPts::Matrix = ones( N ) * pt' + dists .* ( diff * ones( 1, 2 ) )
      # Recover the closest points on the normal surfaces
  idx = setsInsidePtsIdx( nearPts, dirs, dists )
      # Test whether these points are inside the Polygon
  if !any(idx)
    return Inf
        # Outsize return if none of the nearest points are in the polygon
  end
  #### CONTINUE FROM HERE ###
end

# double sets_nearest_line_cpp_eigen( VectorXd pt, MatrixXd mG, VectorXd vC ){
# // Eigen-format c++ rewrite of sets.nearest.line
# // Computes distance of a point to the lines satisfying g.z = c, then picks
# // the nearest which is inside B.  Relies on mG all being unit length
#   int iNN = vC.size() ;
#         // Number of lines
#   const VectorXd vDiff = vC - mG * pt ;
#         // The difference of projections of p onto the mG vectors from vC
#   const VectorXd vDiffAbs( vDiff.cwiseAbs() ) ;
#         // Absolute value of the distance
#   const MatrixXd mNearPts = VectorXd::Ones( iNN ) * pt.transpose() +
#                               mG.cwiseProduct( vDiff * VectorXd::Ones( 2 ).transpose() ) ;
#         // Compute the matrix of nearest points to p
#   LogicalVector boIdx = sets_inside_pts_idx_cpp_eigen( mNearPts, mG, vC ) ;
#         // Vector of logical indices for mNearPts inside mG, vC set
#   if ( is_false( any( boIdx ) ) ) return( INFINITY ) ;
#         // Return outsized distance if none of mNearPts is inside the set
#   VectorXd vDistInside( iNN ) ;
#         // Initialize the vector distances of mNearPts entries inside the set
#   for( int iII=0; iII<iNN; iII++ )
#     vDistInside[ iII ] = ( boIdx[ iII ] == 1 ? vDiffAbs[ iII ] : INFINITY ) ;
#         // Compute the vector distances of mNearPts entries inside the set
#   return vDistInside.minCoeff() ;
# }

function setsInsidePtsIdx( cands::Matrix, dirs::Matrix,
                            dists::Vector )
  N::Int = size(cands)[1]
      # Number of candidate points
  D = dirs * cands'  - dists * ones( 1, N )
      # Negative if inside the set
  return all( D .<= 0, 1 )
end


# LogicalVector sets_inside_pts_idx_cpp_eigen( MatrixXd mPts, MatrixXd mG, VectorXd vC ){
# // Returns the sets of indices of points inside the set described by mG, vC
#   int iRows = mPts.rows() ;
#         // Number of points
#   LogicalVector boInside( iRows ) ;
#         // Instantiate boolean output
#   for( int iII; iII<iRows; iII++ ){
#     boInside[ iII ] = sets_inside_pt_cpp_eigen( mPts.row( iII ), mG, vC ) ;
#   }
#   return( boInside ) ;
# }

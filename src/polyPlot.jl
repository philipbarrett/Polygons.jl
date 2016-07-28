#= polyPlot.jl
Philip Barrett, pobarrett@gmail.com
27may2016, Chicago

Provides plotting functionality for the Polygon class
=#

function polyPlot( pts::Matrix )
  P = [ pts ; pts[1,:] ]
      # Wrap the matrix of points
  if ( size(pts)[1] == 3 ) & (pts[1,1]==pts[2,1]==pts[3,1]) & (pts[1,2]==pts[2,2]==pts[3,2])
    plot( x=P[:,1], y=P[:,2], Geom.point )
  end
  plot( x=P[:,1], y=P[:,2], Geom.path )
end

function polyPlot( poly::Polygon )
  polyPlot( poly.pts )
end

function polyPlot( polys::Array{Polygon,1} )
## TODO: ACCOUNT FOR SINGLETON SETS

  color_vec = ["magenta" "red" "blue" "black" "green" "cyan" "orange"  ]
  # Theme(default_color=)?
  plot( [ layer( x=[ polys[i].pts[:,1] ; polys[i].pts[1,1] ],
                 y=[ polys[i].pts[:,2] ; polys[i].pts[1,2] ],
             Geom.path,
             Theme(default_color=color(parse(Colorant, color_vec[i%7+1]))) )
             for i in 1:length(polys), line_width=8]...)
end

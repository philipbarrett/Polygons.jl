#= polyPlot.jl
Philip Barrett, pobarrett@gmail.com
27may2016, Chicago

Provides plotting functionality for the Polygon class
=#

function polyPlot( pts::Matrix )
  P = [ pts ; pts[1,:] ]
      # Wrap the matrix of points
  polyPlot( x=P[:,1], y=P[:,2], Geom.path )
end

function polyPlot( poly::Polygon )
  polyPlot( poly.pts )
end

function polyPlot( polys::Array{Polygon,1} )
  polyPlot( [ layer( x=[ polys[i].pts[:,1] ; polys[i].pts[1,1] ],
                 y=[ polys[i].pts[:,2] ; polys[i].pts[1,2] ],
             Geom.path,
             Theme(default_color=distinguishable_colors(length(polys))[i]))
             for i in 1:length(polys), line_width=4]...)
end

## Scratch pad

Z = [ 1 -1; 1 1; -1 1; -1 -1 ]
G = [ 1 0; 0 1; -1 0; 0 -1 ]
m = [ 1, 1, 1, 1 ]
a = Polygon( Z, G, m ) ;
b = Polygon( pts=Z ) ;
c = Polygon( dirs=G, dists=m ) ;

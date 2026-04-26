integer :: ix,iy
real :: height2sigma,height

height2sigma(height,ix,iy) = top*(height-topo(ix,iy))/(top-topo(ix,iy))

#! /usr/bin/env python
#
def disk_grid_regular ( n, r, c, ng ):

#*****************************************************************************80
#
## DISK_GRID_REGULAR computes grid points inside a disk.
#
#  Discussion:
#
#    The grid is defined by specifying the radius and center of the disk,
#    and the number of subintervals N into which the horizontal radius
#    should be divided.  Thus, a value of N = 2 will result in 5 points
#    along that horizontal line.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    07 April 2015
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, integer N, the number of subintervals.
#
#    Input, real R, the radius of the disk.
#
#    Input, real C(2), the coordinates of the center of the disk.
#
#    Input, integer NG, the number of grid points, as determined by
#    DISK_GRID_COUNT.
#
#    Output, real CG(2,NG), the grid points inside the disk.
#
  import numpy as np

  cg = np.zeros ( [ 2, ng ] )

  p = 0

  for j in range ( 0, n + 1 ):

    i = 0
    x = c[0]
    y = c[1] + r * float ( 2 * j ) / float ( 2 * n + 1 )
    cg[0,p] = x
    cg[1,p] = y
    p = p + 1

    if ( 0 < j ):

      cg[0,p] = x
      cg[1,p] = 2.0 * c[1] - y
      p = p + 1

    while ( True ):

      i = i + 1
      x = c[0] + r * float ( 2 * i ) / float ( 2 * n + 1 )
      if ( r * r < ( x - c[0] ) ** 2 + ( y - c[1] ) ** 2 ):
        break

      cg[0,p] = x
      cg[1,p] = y
      p = p + 1
      cg[0,p] = 2.0 * c[0] - x
      cg[1,p] = y
      p = p + 1

      if ( 0 < j ):
        cg[0,p] = x
        cg[1,p] = 2.0 * c[1] - y
        p = p + 1
        cg[0,p] = 2.0 * c[0] - x
        cg[1,p] = 2.0 * c[1] - y
        p = p + 1

  return cg

def disk_grid_regular_test ( ):

#*****************************************************************************80
#
#% DISK_GRID_REGULAR_TEST tests DISK_GRID_REGULAR.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    07 April 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  from disk_grid_regular_count import disk_grid_regular_count
  from disk_grid_display import disk_grid_display
#
#  Plot the grid, and save the plot in a file.
#
  n = 50   # Do 50 subdivisions
  r = 1300. # Our radius
  c = np.array([0.0, 0.0])
  ng = disk_grid_regular_count(n, r, c) # Figures out the total number of points given the subdivisions (n)
  cg = disk_grid_regular(n, r, c, ng) # returns the points for the given number of subdivisions (n)
  print cg.shape[1]
  filename = 'disk_grid_regular.png'
  disk_grid_display ( n, r, c, ng, cg, filename )
  return
  
if ( __name__ == '__main__' ):
  disk_grid_regular_test ( )

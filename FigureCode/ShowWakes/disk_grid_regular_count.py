#! /usr/bin/env python
#
def disk_grid_regular_count ( n, r, c ):

#*****************************************************************************80
#
#% DISK_GRID_REGULAR_COUNT counts the grid points inside a disk.
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
#    Output, integer NG, the number of grid points inside the disk.
#
  ng = 0

  for j in range ( 0, n + 1 ):

    i = 0
    x = c[0]
    y = c[1] + r * float ( 2 * j ) / float ( 2 * n + 1 )
    ng = ng + 1
    if ( 0 < j ):
      ng = ng + 1

    while ( True ):

      i = i + 1
      x = c[0] + r * float ( 2 * i ) / float ( 2 * n + 1 )

      if ( r * r < ( x - c[0] ) ** 2 + ( y - c[1] ) ** 2 ):
        break

      ng = ng + 1
      ng = ng + 1

      if ( 0 < j ):
        ng = ng + 1
        ng = ng + 1

  return ng
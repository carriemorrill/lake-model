!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                    TRIDIAG_SOLVE
!     compute the solution of a tridiagonal linear system
!
!     based on a streamlined version of subroutine trdi in the model of 
!     Schneider and Thompson (JGR, 1981).  Note: this subroutine 
!     executes satisfactorily if the input matrix is diagonally
!     dominant and non-singular. The diagonal elements are used to
!     pivot, and no tests are made to determine singularity. If a 
!     singular or numerically singular matrix is used as input, a 
!     divide zero or floating point overflow will result.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine tridiag_solve (ne, a, b, c, d, x)
  
      implicit none

      integer ne   ! the number of unknowns = lake + sed layers [count]
      real*8 a(ne) ! subdiagonals of tridiagonal matrix [unitless]
      real*8 b(ne) ! main diagonals of tridiagonal matrix [unitless]
      real*8 c(ne) ! super-diagonals of tridiagonal matrix [unitless]
      real*8 d(ne) ! right-hand side of linear system [degrees C]
      real*8 x(ne) ! solutions of the linear system [degrees C]

      real*8 alpha(ne) ! working array
      real*8 gamma(ne) ! working array

      integer i    ! counter for looping through unknowns
      integer ib   ! counter for looping through unknowns
 
      alpha(1) = 1. / b(1)
      gamma(1) = c(1) * alpha(1)

      do i = 2,ne-1
        alpha(i) = 1. / (b(i) - a(i) * gamma(i-1))
        gamma(i) = c(i) * alpha(i)
      enddo
 
      x(1) = d(1) * alpha(1)

      do i = 2,ne-1
        x(i) = (d(i) - a(i) * x(i-1)) * alpha(i)
      enddo
  
      x(ne) = (d(ne) - a(ne) * x(ne-1)) / (b(ne) - a(ne) * gamma(ne-1))
                   
      do i = 1,ne-1
         ib = ne - i
         x(ib) = x(ib) - gamma(ib) * x(ib+1)
      enddo

      return
      end 

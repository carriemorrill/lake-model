!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                 FILE_OPEN
!     open input and output files
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine file_open

!     input files
      open(unit=15,file='met-input.txt',status='old')

!     output files
      open(unit=50,file='profile.dat',status='unknown')
      open(unit=51,file='surface.dat',status='unknown')

      return
      end

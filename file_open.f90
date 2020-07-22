!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                 FILE_OPEN
!     open input and output files
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine file_open

      implicit none
      include 'lake.inc'
      
      character*11 varname(17)          ! variable names for output headers
      character*7  varformat(17)        ! variable formats for output headers
      data varname /"    tlake","   fice","    hice","   hsnow","    evap",     &
                    "   lakelev","  discharge","  mix","   swup","    lwup", &
                    "   swnet","   lwnet","    sens","  latent","   d18O", &
                    "     d2H","  salty"/
      data varformat /"(a9)","(a7)","(a8)","(a8)","(a8)","(a10)","(a11)",  &
                      "(a5)","(a7)","(a8)","(a8)","(a8)","(a8)","(a8)",   &
                      "(a7)","(a8)","(a7)"/

      integer j                         ! counter for printing

!     input files
      open(unit = 15,file = datafile,status = 'old')

!     output files
      open(unit = 50,file = 'surface.txt',status = 'unknown')
      write(50,"(a12)",advance="no") "YEAR MON DAY"
      do j = 1,nout
         write(50,varformat(output(j)),advance="no") varname(output(j))
      enddo
      write(50,*)                  ! advance to next line

      if (tlake_profile) then
         open(unit = 51,file = 'profile-laketemp.txt',status = 'unknown')
         write(51,"(a12)",advance="no") "YEAR MON DAY"
         if (tlake_print(2).eq.-999) then
           write(51,"(1000(3x,i2,2x))") (j,j=tlake_print(1),max_dep,tlake_print(3))
         else
           write(51,"(1000(3x,i2,2x))") (j,j=tlake_print(1),tlake_print(2),tlake_print(3))
         endif
      endif

      if (tsed_profile) then
         open(unit = 52,file = 'profile-sedtemp.txt',status = 'unknown')
         if (tsed_print(2).eq.-999) tsed_print(2) = nzsed
         write(52,"(a12,1000(3x,i2,2x))") "YEAR MON DAY",(j,j=tsed_print(1),tsed_print(2),tsed_print(3))
      endif

      if (d18O_profile) then
         open(unit = 53,file = 'profile-d18O.txt',status = 'unknown')
         write(53,"(a12)",advance="no") "YEAR MON DAY"
         if (d18O_print(2).eq.-999) then
           write(53,"(1000(3x,i2,2x))") (j,j=d18O_print(1),max_dep,d18O_print(3))
         else
           write(53,"(1000(3x,i2,2x))") (j,j=d18O_print(1),d18O_print(2),d18O_print(3))
         endif
      endif

      if (d2H_profile) then
         open(unit = 54,file = 'profile-d2H.txt',status = 'unknown')
         write(54,"(a12)",advance="no") "YEAR MON DAY"
         if (d2H_print(2).eq.-999) then
           write(54,"(1000(4x,i2,2x))") (j,j=d2H_print(1),max_dep,d2H_print(3))
         else
           write(54,"(1000(4x,i2,2x))") (j,j=d2H_print(1),d2H_print(2),d2H_print(3))
         endif
      endif

      if (salty_profile) then
         open(unit = 55,file = 'profile-salinity.txt',status = 'unknown')
         write(55,"(a12)",advance="no") "YEAR MON DAY"
         if (salty_print(2).eq.-999) then
           write(55,"(1000(3x,i2,2x))") (j,j=salty_print(1),max_dep,salty_print(3))
         else
           write(55,"(1000(3x,i2,2x))") (j,j=salty_print(1),salty_print(2),salty_print(3))
         endif
      endif

      return
      end

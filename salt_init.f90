!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                   SALT_INIT
!     calculate freezing point, in Gill 1982, from Millero 1978
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine salt_init (ps, tfreezeC, s)

      implicit none
  
      real ps     ! surface air pressure [Pa]	  
      real tfreezeC  ! lake surface temperature [degrees C]
      real s      ! lake surface salinity [ppt]

      real psd    ! surface air pressure [decibar]

      psd = ps / 10000.             
      tfreezeC = -0.0575 * s + 1.710523e-3 * s**(3. / 2.)               &
         - 2.154996e-4 * s** 2. - 7.53e-4 * psd
     
      return
      end 

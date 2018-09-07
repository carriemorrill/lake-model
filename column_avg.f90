!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                 COLUMN_AVG
!     average water and ice columns for fractional cover
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine column_avg (nzlk, t, t_i, tsed, tsed_i, salty, salty_i,&
                             fracprv)

      implicit none
      include 'lake.inc'

      real t(max_dep)   ! input: water temp in open water col [degrees C]
                        ! output: water temp average of two columns
      real t_i(max_dep) ! water temp in ice covered column [degrees C]
      real tsed(nzsed)  ! input: sediment temp in open water col [degrees C]
                        ! output: sediment temp average of two cols
      real tsed_i(nzsed)   ! sediment temp in ice covered col [degrees C]
      real salty(max_dep)  ! input: salinity of open water col
                           ! output: salinity average of two cols
      real salty_i(max_dep)! salinity of ice covered col [ppt]
      real fracprv      ! ice cover fraction [0-1]
  
      real dnstyw   ! water density in open water column [kg/m3]
      real dnstyi   ! water density in ice covered column [kg/m2]
      real cpw      ! water spec. heat cap. in open water col [J/kgK]
      real cpi      ! water spec. heat cap. in ice cover col [J/kgK]
      real z        ! thickness of lake water layer [m]	  
      integer j     ! counter for looping through lake and sed layers

      do j = 1,nzlk ! average temp and salinity over water layers
        call density(t(j), salty(j), dnstyw)
        call density(t_i(j), salty_i(j), dnstyi)
        call specheat(t(j), salty(j), cpw)
        call specheat(t_i(j), salty_i(j), cpi)
        z = dz
        if (j.eq.1) z = surf
        t(j) = ((1.-fracprv) * t(j) * z * (1.e3 + dnstyw) * cpw +       &
          fracprv * t_i(j) * z * (1.e3 + dnstyi) * cpi )/               &
          ((z * (1.e3 + dnstyw) * cpw + z * (1.e3 + dnstyi) * cpi) *    &
          0.5)
        salty(j) = (1. - fracprv) * salty(j) + fracprv * salty_i(j)    
      enddo

      do j = 1,nzsed
        tsed(j) = (1. - fracprv) * tsed(j) + fracprv * tsed_i(j)
      enddo

      return
      end 

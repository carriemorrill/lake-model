!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                 COLUMN_AVG
!     average water and ice columns for fractional cover
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine column_avg (nzlk, t, t_i, tsed, tsed_i, tracer, tracer_i,&
                             fracprv)

      implicit none
      include 'lake.inc'

      real t(max_dep)   ! input: water temp in open water col [degrees C]
                        ! output: water temp average of two columns
      real t_i(max_dep) ! water temp in ice covered column [degrees C]
      real tsed(nzsed)  ! input: sediment temp in open water col [degrees C]
                        ! output: sediment temp average of two cols
      real tsed_i(nzsed)   ! sediment temp in ice covered col [degrees C]
      real tracer(3,max_dep)  ! input: tracers of open water col
                           ! output: tracers average of two cols
      real tracer_i(3,max_dep)! tracers of ice covered col [ppt]
      real fracprv      ! ice cover fraction [0-1]
  
      real dnstyw   ! water density in open water column [kg/m3]
      real dnstyi   ! water density in ice covered column [kg/m2]
      real cpw      ! water spec. heat cap. in open water col [J/kgK]
      real cpi      ! water spec. heat cap. in ice cover col [J/kgK]
      real z        ! thickness of lake water layer [m]	  
      integer j     ! counter for looping through lake and sed layers
      integer i           ! counter for looping through tracers

      do j = 1,nzlk ! average temp and salinity over water layers
        call density(t(j), tracer(1,j), dnstyw)
        call density(t_i(j), tracer_i(1,j), dnstyi)
        call specheat(t(j), tracer(1,j), cpw)
        call specheat(t_i(j), tracer_i(1,j), cpi)
        z = dz
        if (j.eq.1) z = surf
        t(j) = ((1.-fracprv) * t(j) * z * (1.e3 + dnstyw) * cpw +       &
          fracprv * t_i(j) * z * (1.e3 + dnstyi) * cpi )/               &
          ((z * (1.e3 + dnstyw) * cpw + z * (1.e3 + dnstyi) * cpi) * 0.5) 
        do i = 1,3
           tracer(i,j) = (1. - fracprv) * tracer(i,j) + fracprv * tracer_i(i,j)  
        enddo
      enddo

      do j = 1,nzsed
        tsed(j) = (1. - fracprv) * tsed(j) + fracprv * tsed_i(j)
      enddo

      return
      end 

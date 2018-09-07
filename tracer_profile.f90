!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                      TRACER_PROFILE
!     compute diffusive mixing of salinity profile using Crank
!     Nicolson solution of 1-D diffusion equation
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine tracer_profile (de, nzlk, salty)
   
      implicit none
      include 'lake.inc'

      real de(nzlk)      ! eddy diffusivity [m2/s]
      real salty(max_dep)! lake layer salinity [ppt]
 
      real dzall(nzlk)   ! thickness of each lake and sediment layer [m]
      real ztop(nzlk)    ! depth of each lake and sediment layer at top [m]	  
      real*8 tall(nzlk)  ! starting temperature of layers [degrees C]
      real dei(nzlk)     ! thermal conductivity at layer interface [W/mK] 
      real dztop         ! distance from one layer node to the layer node above [m]
      real dzbot         ! distance from one layer node to the layer node below [m]
      real*8 a(nzlk)     ! "a" vector for tridiagonal matrix [unitless]
      real*8 b(nzlk)     ! "b" vector for tridiagonal matrix [unitless]
      real*8 c(nzlk)     ! "c" vector for tridiagonal matrix [unitless]
      real*8 d(nzlk)     ! "d" vector for tridiagonal matrix [degrees C]
      real*8 tnew(nzlk)  ! new layer temperature from Crank-Nicolson [degrees C]
      integer k          ! counter for looping through layers

      do k = 1, nzlk
        if (k.eq.1) then
          dzall(k) = surf
        else
          dzall(k) = dz
        endif
      enddo

      ztop(1) = 0.
      do k = 2, nzlk
        ztop(k) = ztop(k-1) + dzall(k-1)
      enddo

      do k = 1,nzlk 
        tall(k) = salty(k)
      enddo

! ************** calculate layer interface thermal conductivities ******

      do k=1,nzlk
         if (k.lt.nzlk)                                                 &
            dei(k) = de(k) * de(k+1) * (dzall(k) + dzall(k+1)) /        &
                       (de(k) * dzall(k+1) + de(k+1) * dzall(k))
         if (k.eq.nzlk)                                                 &
            dei(k) = 0.
      enddo

! ******** calculate arrays for tridiagonal matrix *********************

      k = 1   ! first lake layer, diffusion only to bottom
        dzbot = (dzall(k) + dzall(k+1)) / 2.
        c(k) = -0.5 * dei(k) / dzbot * dt / dzall(k)
        a(k) = 0.
        b(k) = 1. - c(k)  
        d(k) = tall(k) + 0.5 * dei(k) * dt / dzall(k) *                 &
                (tall(k+1) - tall(k)) / dzbot
     
      do k = 2,nzlk-1  ! middle layers, diffusion to top and bot
        dzbot = (dzall(k) + dzall(k+1)) / 2.
        dztop = (dzall(k-1) + dzall(k)) / 2.
        c(k) = -0.5 * dei(k) / dzbot * dt / dzall(k)  
        a(k) = -0.5 * dei(k-1) / dztop * dt / dzall(k)  
        b(k) = 1. - c(k) - a(k)
        d(k) = tall(k) + 0.5 * dei(k) * dt / dzall(k) *                 &
                (tall(k+1) - tall(k)) / dzbot - 0.5 * dei(k-1) *        &
                dt / dzall(k) * (tall(k) - tall(k-1)) / dztop
      enddo

      k = nzlk ! bottom sed layer, diffusion only to top
        dztop = (dzall(k-1) + dzall(k)) / 2.
        c(k) = 0.
        a(k) = -0.5 * dei(k-1) / dztop * dt / dzall(k)
        b(k) = 1. - a(k)
        d(k) = tall(k) - 0.5 * dei(k-1) * dt / dzall(k) *               &
               (tall(k) - tall(k-1)) / dztop

! ******* solve matrix and reset temp and density arrays ***************

       call tridiag_solve (nzlk, a, b, c, d, tnew)
  
       do k = 1, nzlk  ! reset temps and densities
          salty(k) = tnew(k)
       enddo

       return
       end 

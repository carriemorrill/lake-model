!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                      TEMP_PROFILE_NOSED
!     compute diffusive mixing of temperature profile using Crank
!     Nicolson solution of 1-D heat equation. This version does not
!     model lake sediments.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine temp_profile_nosed (iwater, qbot, qw, t, sw, lnet, qe, &
                                     qh, de, nzlk, salty)
   
      implicit none
      include 'lake.inc'

      integer iwater     ! open water or ice flag	  
      real qbot          ! solar radiation flux at bottom of ice [W/m2]
      real qw            ! heat flux to ice from water [W/m2]
      real t(nzlk)       ! lake layer temperature [degrees C]  
      real sw            ! downward shortwave radiation at lake surface [W/m2]
      real lnet          ! net downward longwave radiation at lake surface [W/m2]
      real qe            ! latent heat flux to lake [W/m2]
      real qh            ! sensible heat flux to lake [W/m2]
      real de(nzlk)      ! eddy diffusivity [m2/s]
      real salty(max_dep)! lake layer salinity [ppt]
 
      real dzall(nzlk)   ! thickness of each lake layer [m]
      real ztop(nzlk)    ! depth of each lake layer at top [m]	  
      real*8 tall(nzlk)  ! starting temperature of layers [degrees C]
      real dnsty(nzlk)   ! water density anomaly from 1000 [kg/m3] 
      real cpz(nzlk)     ! specific heat of lake layer [J/kgK]  
      real tk(nzlk)      ! thermal conductivity of water at layer node [W/mK]
      real factx(nzlk)   ! conversion factor from W/m2 to K [sm2K/J]  
      real*8 h(nzlk)     ! heating terms for each layer [W/m2]
      real tki(nzlk)     ! thermal conductivity at layer interface [W/mK] 
      real dztop         ! distance from one layer node to the layer node above [m]
      real dzbot         ! distance from one layer node to the layer node below [m]
      real*8 a(nzlk)     ! "a" vector for tridiagonal matrix [unitless]
      real*8 b(nzlk)     ! "b" vector for tridiagonal matrix [unitless]
      real*8 c(nzlk)     ! "c" vector for tridiagonal matrix [unitless]
      real*8 d(nzlk)     ! "d" vector for tridiagonal matrix [degrees C]
      real*8 tnew(nzlk)  ! new layer temperature from Crank-Nicolson [degrees C]
      integer k          ! counter for looping through layers

      do k = 1, nzlk
        if (k.eq.1) dzall(k) = surf
        if (k.gt.1) dzall(k) = dz
      enddo

      ztop(1) = 0.
      do k = 2, nzlk
        ztop(k) = ztop(k-1) + dzall(k-1)
      enddo

      do k = 1,nzlk 
        tall(k) = t(k)
        call density (t(k), salty(k), dnsty(k))
        call specheat (t(k), salty(k), cpz(k))
        tk(k) = de(k) * (1.e3 + dnsty(k)) * cpz(k)
        factx(k) = dt / ((1.e3 + dnsty(k)) * cpz(k) * dzall(k))
      enddo

! *************** calculate heating terms ******************************

      if (iwater.eq.1) then  ! open water
         do k=1,nzlk
            if (k.eq.1) &   ! topmost lake layer
              h(k) = (sw * beta) + (1. - beta) * sw * (1. - exp(-eta *  &
                ztop(k+1))) + (lnet + qe + qh)
            if (k.gt.1.and.k.lt.nzlk) &  ! middle lake layers          
              h(k) = (1. - beta) * sw * (exp(-eta * ztop(k))            &
                - exp(-eta * ztop(k+1))) 
            if (k.eq.nzlk) &  ! middle lake layers          
              h(k) = (1. - beta) * sw * exp(-eta * ztop(k))          
         enddo
      else   ! lake is ice covered
         do k=1,nzlk
            if (k.eq.1)                                                 &
              h(k) = (qbot * beta) + (1. - beta) * qbot * (1. - exp     &
                (-eta * ztop(k+1))) - qw
            if (k.gt.1.and.k.lt.nzlk)                                   &
              h(k) = (1. - beta) * qbot * (exp(-eta * ztop(k))          &
                   - exp(-eta * ztop(k+1))) 
            if (k.eq.nzlk) &  ! middle lake layers          
              h(k) = (1. - beta) * sw * exp(-eta * ztop(k)) 
         enddo
      endif

! ************** calculate layer interface thermal conductivities ******

      do k=1,nzlk
         if (k.lt.nzlk)                                                 &
            tki(k) = tk(k) * tk(k+1) * (dzall(k) + dzall(k+1)) /        &
                       (tk(k) * dzall(k+1) + tk(k+1) * dzall(k))
         if (k.eq.nzlk)                                                 &
            tki(k) = 0.
      enddo

! ******** calculate arrays for tridiagonal matrix *********************

      k = 1   ! first lake layer, diffusion only to bottom
        dzbot = (dzall(k) + dzall(k+1)) / 2.
        c(k) = -0.5 * tki(k) / dzbot * factx(k)
        a(k) = 0.
        b(k) = 1. - c(k)  
        d(k) = tall(k) + factx(k) *                                     &
              (h(k) + 0.5 * tki(k) * (tall(k+1) - tall(k)) / dzbot)
     
      do k = 2,nzlk-1  ! middle layers, diffusion to top and bot
        dzbot = (dzall(k) + dzall(k+1)) / 2.
        dztop = (dzall(k-1) + dzall(k)) / 2.
        c(k) = -0.5 * tki(k) / dzbot * factx(k)  
        a(k) = -0.5 * tki(k-1) / dztop * factx(k)  
        b(k) = 1. - c(k) - a(k)
        d(k) = tall(k) + factx(k) *                                     &
              (h(k) + 0.5 * tki(k) * (tall(k+1) - tall(k)) / dzbot      &
               - 0.5 * tki(k-1) * (tall(k) - tall(k-1)) / dztop)
      enddo

      k = nzlk     ! bottom layer, diffusion only to top
        dztop = (dzall(k-1) + dzall(k)) / 2.
        c(k) = 0.
        a(k) = -0.5 * tki(k-1) / dztop * factx(k)
        b(k) = 1. - a(k)
        d(k) = tall(k) + factx(k) *                                     &
              (h(k) - 0.5 * tki(k-1) * (tall(k) - tall(k-1)) / dztop)

! ******* solve matrix and reset temp and density arrays ***************

       call tridiag_solve (nzlk, a, b, c, d, tnew)
  
       do k = 1, nzlk  ! reset temps and densities
          t(k) = tnew(k)
       enddo

       return
       end 

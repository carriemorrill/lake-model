!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                      TEMP_PROFILE
!     compute diffusive mixing of temperature profile using Crank
!     Nicolson solution of 1-D heat equation. This version includes
!     both lake and sediment layers
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine temp_profile (iwater, qbot, qw, t, sw, lnet, qe, qh,   &
                               de, nzlk, salty, tsed)
   
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
      real tsed(nzsed)   ! sediment layer temperature [degrees C]
 
      real zsed(nzsed)   ! sediment node depths [m]
      real dzall(nzlk+nzsed+1)  ! thickness of each lake and sediment layer [m]
      real ztop(nzlk+nzsed+1)   ! depth of each lake and sediment layer at top [m]	  
      real*8 tall(nzlk+nzsed+1) ! starting temperature of layers [degrees C]
      real dnsty(nzlk+1) ! water density anomaly from 1000 [kg/m3] 
      real cpz(nzlk+1)   ! specific heat of lake layer [J/kgK]  
      real tk(nzlk+1)    ! thermal conductivity of water at layer node [W/mK]
      real factx(nzlk+nzsed+1)  ! conversion factor from W/m2 to K [sm2K/J]  
      real*8 h(nzlk+nzsed+1)    ! heating terms for each layer [W/m2]
      real tki(nzlk+nzsed+1)    ! thermal conductivity at layer interface [W/mK] 
      real dztop         ! distance from one layer node to the layer node above [m]
      real dzbot         ! distance from one layer node to the layer node below [m]
      real*8 a(nzlk+nzsed+1)    ! "a" vector for tridiagonal matrix [unitless]
      real*8 b(nzlk+nzsed+1)    ! "b" vector for tridiagonal matrix [unitless]
      real*8 c(nzlk+nzsed+1)    ! "c" vector for tridiagonal matrix [unitless]
      real*8 d(nzlk+nzsed+1)    ! "d" vector for tridiagonal matrix [degrees C]
      real*8 tnew(nzlk+nzsed+1) ! new layer temperature from Crank-Nicolson [degrees C]
      integer k          ! counter for looping through layers

! ********** set up combined lake and sediment arrays ******************
! ****** splitting lowest lake layer into two to improve solution ******

      do k = 1, nzsed  
        zsed(k) = 0.025 * (exp(0.5 * (k - 0.5)) - 1.)
      enddo

      do k = 1, nzlk+nzsed
        if (k.eq.1) dzall(k) = surf
        if (k.gt.1.and.k.le.nzlk-1) dzall(k) = dz
        if (k.eq.nzlk) dzall(k) = dz - (0.5 * (zsed(1) + zsed(2)))
        if (k.eq.nzlk+1) dzall(k) = 0.5 * (zsed(1) + zsed(2))
        if (k.eq.nzlk+2) dzall(k) = 0.5 * (zsed(1) + zsed(2))
        if (k.gt.nzlk+2) &
          dzall(k) = 0.5 * (zsed(k-nzlk) - zsed(k-nzlk-2))
      enddo
      dzall(nzlk+nzsed+1) = zsed(nzsed) - zsed(nzsed-1)

      ztop(1) = 0.
      do k = 2, nzlk+nzsed+1
        ztop(k) = ztop(k-1) + dzall(k-1)
      enddo

      do k = 1,nzlk 
        tall(k) = t(k)
        call density (t(k), salty(k), dnsty(k))
        call specheat (t(k), salty(k), cpz(k))
        tk(k) = de(k) * (1.e3 + dnsty(k)) * cpz(k)
        factx(k) = dt / ((1.e3 + dnsty(k)) * cpz(k) * dzall(k))
      enddo

      tall(nzlk+1) = t(nzlk)
      dnsty(nzlk+1) = dnsty(nzlk)
      cpz(nzlk+1) = cpz(nzlk) 
      tk(nzlk+1) = tk(nzlk)
      factx(nzlk+1) = factx(nzlk)

      do k=1,nzsed
        tall(nzlk+1+k) = tsed(k)
        factx(nzlk+1+k) = dt / (csed * dzall(nzlk+1+k))
      enddo

! *************** calculate heating terms ******************************

      if (iwater.eq.1) then  ! open water
         do k=1,nzlk+nzsed+1
            if (k.eq.1) &   ! topmost lake layer
              h(k) = (sw * beta) + (1. - beta) * sw * (1. - exp(-eta *  &
                ztop(k+1))) + alb_sed * (1. - beta) * sw * exp(-eta *   &
                ztop(nzlk+2)) * (exp(-eta * (ztop(nzlk+2) - ztop(k+1))) &
                -exp(-eta * (ztop(nzlk+2) - ztop(k))))                  &
                + (lnet + qe + qh)
            if (k.gt.1.and.k.le.nzlk) &  ! middle lake layers          
              h(k) = (1. - beta) * sw * (exp(-eta * ztop(k))            &
                - exp(-eta * ztop(k+1))) + alb_sed * (1. - beta) * sw * &
                exp(-eta * ztop(nzlk+2))*                               &
                (exp(-eta * (ztop(nzlk+2) - ztop(k+1)))                 &
                -exp(-eta * (ztop(nzlk+2) - ztop(k)))) 
            if (k.eq.nzlk+1) &  ! bottom lake layer
              h(k) = (1. - beta) * sw * (exp(-eta * ztop(k))            &
                - exp(-eta * ztop(k+1))) + alb_sed * (1. - beta) * sw * &
                exp(-eta * ztop(nzlk+2)) *                              &
                (1. - exp(-eta * (ztop(nzlk+2) - ztop(k)))) 
            if (k.eq.nzlk+2) &  ! topmost sed layer
              h(k) = (1. - alb_sed) * (1. - beta) * sw *                &
                     (exp(-eta * ztop(k)))  
            if (k.gt.nzlk+2)  &  ! rest of sed layers
              h(k) = 0.0
         enddo
      else   ! lake is ice covered
         do k=1,nzlk+nzsed+1
            if (k.eq.1)                                                 &
              h(k) = (qbot * beta) + (1. - beta) * qbot * (1. - exp     &
                (-eta * ztop(k+1))) + alb_sed * (1. - beta) * qbot *    &
                exp(-eta * ztop(nzlk+2)) * (exp(-eta * (ztop(nzlk+2) -  &
                ztop(k+1))) -exp(-eta * (ztop(nzlk+2) - ztop(k))))      &
                - qw 
            if (k.gt.1.and.k.le.nzlk)                                   &
              h(k) = (1. - beta) * qbot * (exp(-eta * ztop(k))          &
                   - exp(-eta * ztop(k+1))) + alb_sed * (1. - beta) *   &
                   qbot * exp(-eta * ztop(nzlk+2)) *                    &
                   (exp(-eta * (ztop(nzlk+2) - ztop(k+1)))              &
                   - exp(-eta * (ztop(nzlk+2) - ztop(k)))) 
            if (k.eq.nzlk+1)                                            &
              h(k) =  (1. - beta) * qbot * (exp(-eta * ztop(k))         &
                  - exp(-eta * ztop(k+1))) + alb_sed * (1. - beta) *    &
                  sw * exp(-eta * ztop(nzlk+2)) * (1. - exp             &
                  (-eta * (ztop(nzlk+2) - ztop(k))))  
            if (k.eq.nzlk+2)                                            &
              h(k) = (1. - alb_sed) * (1. - beta) * qbot *              &
                     (exp(-eta * ztop(k))) 
            if (k.gt.nzlk+2)                                            &
              h(k) = 0.0
         enddo
      endif

! ************** calculate layer interface thermal conductivities ******

      do k=1,nzlk+nzsed+1
         if (k.lt.nzlk)                                                 &
            tki(k) = tk(k) * tk(k+1) * (dzall(k) + dzall(k+1)) /        &
                       (tk(k) * dzall(k+1) + tk(k+1) * dzall(k))
         if (k.eq.nzlk)                                                 &
            tki(k) = tk(k) * tk(k) * (dzall(k) + dzall(k+1)) /          &
                       (tk(k) * dzall(k+1) + tk(k) * dzall(k))
         if (k.eq.nzlk+1)                                               &
            tki(k) = tk(k-1) * condsed * (dzall(k) + dzall(k+1)) /      &
                 (tk(k-1) * dzall(k+1) + condsed * dzall(k))
         if (k.gt.nzlk+1.and.k.le.nzlk+nzsed)                           &
            tki(k) = condsed * condsed * (dzall(k) +                    &
                      dzall(k+1)) / (condsed * dzall(k)                 &
                      + condsed * dzall(k+1))
         if (k.gt.nzlk+nzsed)                                           &
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
     
      do k = 2,nzlk+nzsed  ! middle layers, diffusion to top and bot
        dzbot = (dzall(k) + dzall(k+1)) / 2.
        dztop = (dzall(k-1) + dzall(k)) / 2.
        c(k) = -0.5 * tki(k) / dzbot * factx(k)  
        a(k) = -0.5 * tki(k-1) / dztop * factx(k)  
        b(k) = 1. - c(k) - a(k)
        d(k) = tall(k) + factx(k) *                                     &
              (h(k) + 0.5 * tki(k) * (tall(k+1) - tall(k)) / dzbot      &
               - 0.5 * tki(k-1) * (tall(k) - tall(k-1)) / dztop)
      enddo

      k = nzlk+nzsed+1 ! bottom sed layer, diffusion only to top
        dztop = (dzall(k-1) + dzall(k)) / 2.
        c(k) = 0.
        a(k) = -0.5 * tki(k-1) / dztop * factx(k)
        b(k) = 1. - a(k)
        d(k) = tall(k) + factx(k) *                                     &
              (h(k) - 0.5 * tki(k-1) * (tall(k) - tall(k-1)) / dztop)

! ******* solve matrix and reset temp and density arrays ***************

       call tridiag_solve (nzlk+nzsed+1, a, b, c, d, tnew)
  
       do k = 1, nzlk  ! reset temps and densities
          if (k.lt.nzlk) then
             t(k) = tnew(k)
          else  ! average two bottom lake layers back into one
             t(k) = (tnew(k) * dzall(k) + tnew(k+1) * dzall(k+1)) /   &
                      (dzall(k) + dzall(k+1))
          endif
       enddo

       do k=1,nzsed
          tsed(k) = tnew(k+nzlk+1)
       enddo

       return
       end 

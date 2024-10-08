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
      real area_top(nzlk)! area of each lake layer at top [m]
      real area_bot(nzlk)! area of each lake layer at bottom [m]
      real area_mid(nzlk)! area of each lake layer at middle [m]
      real*8 tall(nzlk)  ! starting salinity value of layers [ppt]
      real dei(nzlk)     ! tracer conductivity at layer interface [m2/s] 
      real dztop         ! distance from one layer node to the layer node above [m]
      real dzbot         ! distance from one layer node to the layer node below [m]
      real*8 a(nzlk)     ! "a" vector for tridiagonal matrix [unitless]
      real*8 b(nzlk)     ! "b" vector for tridiagonal matrix [unitless]
      real*8 c(nzlk)     ! "c" vector for tridiagonal matrix [unitless]
      real*8 d(nzlk)     ! "d" vector for tridiagonal matrix [ppt or per mil]
      real*8 tnew(nzlk)  ! new layer salinity value from Crank-Nicolson [ppt]
      integer k          ! counter for looping through layers
      integer ktop       ! index value of lake surface

      ktop = max_dep - nzlk   ! used to determine lake area from input bathymetry
      do k = 1, nzlk          ! set up lake depth and area arrays
        area_mid(k) = area(k+ktop)
        if (k.eq.1) then      ! top lake layer
          dzall(k) = surf
          ztop(k) = 0.
          area_top(k) = area(k+ktop)
          area_bot(k) = (area(k+ktop) + area(k+ktop+1))/2.
        endif
        if (k.gt.1.and.k.lt.nzlk) then
          dzall(k) = dz
          ztop(k) = ztop(k-1) + dzall(k-1)
          area_top(k) = (area(k+ktop) + area(k+ktop-1))/2.
          area_bot(k) = (area(k+ktop) + area(k+ktop+1))/2.
        endif
        if (k.eq.nzlk) then   ! bottom lake layer
          dzall(k) = dz
          ztop(k) = ztop(k-1) + dzall(k-1)
          area_top(k) = (area(k+ktop) + area(k+ktop-1))/2.
          area_bot(k) = area(k+ktop)
        endif
      enddo

      do k = 1, nzlk
        tall(k) = salty(k)
      enddo

! ************** calculate layer interface eddy conductivities ******

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
        c(k) = -0.5 * (area_bot(k)/area_mid(k)) * dei(k) /              &
                dzbot * dt / dzall(k)
        a(k) = 0.
        b(k) = 1. - c(k)  
        d(k) = tall(k) + 0.5 * (area_bot(k)/area_mid(k)) * dei(k) *     &
               dt / dzall(k) * (tall(k+1) - tall(k)) / dzbot
     
      do k = 2,nzlk-1  ! middle layers, diffusion to top and bot
        dzbot = (dzall(k) + dzall(k+1)) / 2.
        dztop = (dzall(k-1) + dzall(k)) / 2.
        c(k) = -0.5 * (area_bot(k)/area_mid(k)) * dei(k) /              &
                dzbot * dt / dzall(k)  
        a(k) = -0.5 * (area_top(k)/area_mid(k)) * dei(k-1) /            &
                dztop * dt / dzall(k)
        b(k) = 1. - c(k) - a(k)
        d(k) = tall(k) + 0.5 * (area_bot(k)/area_mid(k)) * dei(k) *     &
               dt / dzall(k) * (tall(k+1) - tall(k)) / dzbot            &
               - 0.5 * (area_top(k)/area_mid(k)) * dei(k-1) *           &
                dt / dzall(k) * (tall(k) - tall(k-1)) / dztop
      enddo

      k = nzlk ! bottom sed layer, diffusion only to top
        dztop = (dzall(k-1) + dzall(k)) / 2.
        c(k) = 0.
        a(k) = -0.5 * (area_top(k)/area_mid(k)) * dei(k-1) /            &
               dztop * dt / dzall(k)
        b(k) = 1. - a(k)
        d(k) = tall(k) - 0.5 * (area_top(k)/area_mid(k)) * dei(k-1) *   &
               dt / dzall(k) * (tall(k) - tall(k-1)) / dztop

! ******* solve matrix and reset salinity array ************************

       call tridiag_solve (nzlk, a, b, c, d, tnew)
  
       do k = 1, nzlk  ! reset salinity
          salty(k) = tnew(k)
       enddo

       return
       end 

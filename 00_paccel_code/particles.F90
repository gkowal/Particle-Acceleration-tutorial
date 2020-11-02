!!******************************************************************************
!!
!! module: particles - subroutines to prepare and advance particles
!!
!! Copyright (C) 2010 Grzegorz Kowal <grzegorz@gkowal.info>
!!
!!******************************************************************************
!!
!!  This file is part of PAccel.
!!
!!  PAccel is free software; you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation; either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  PAccel is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!
!!******************************************************************************
!!
!
module particles

  implicit none

! domain dimensions
!
  integer, dimension(3)  , save :: dm, qm

! domain bounds
!
  real   , dimension(3,2), save :: bnds

! domain size
!
  real   , dimension(3)  , save :: bsiz

! particle mass and the speed of light
!
  real(kind=8)           , save :: mrest, qom, csq, om0, fc, ln, bavg, bpar

! arrays containing the initial positions and velocities of particle
!
  real(kind=PREC), dimension(3), save :: x0, u0, p0

! array to store the dump times
!
  real(kind=8), dimension(:), save, allocatable :: tt

! global parameters
!
    real(kind=PREC), parameter :: pi2 = 6.2831853071795862319959269370884d0
!
!-------------------------------------------------------------------------------
!
  contains
!
!===============================================================================
!
! init_particle: subroutine initializes the particle position and momentum
!
!===============================================================================
!
  subroutine init_particle()

    use fields, only : get_dimensions, get_domain_bounds, bx, by, bz
    use params, only : ptype, vpar, vper, c, dens, tunit, tmulti, bunit, xc, yc, zc
    use params, only : output, tmin, tmax, ndumps
#ifdef TEST
    use params, only : bini, bshr, bamp, vamp, vrat, freq
#endif /* TEST */

    implicit none

! local variables
!
    integer      :: p, n
    real(kind=PREC) :: vp, vr, vv, va
    real(kind=PREC) :: gm, dn, mu0, om, tg, rg, mu, mp, en, ek, ba
    real(kind=PREC) :: bb, ub, uu
#ifdef ITEST
    real(kind=PREC) :: xt, yt, rt, dl, ra, rb, ec
#endif /* ITEST */

! arrays
!
    real(kind=PREC), dimension(3) :: r0
    real(kind=PREC), dimension(3) :: b, u, w

! position indices
!
    integer        , dimension(4) :: ii, jj, kk
    real(kind=8   ), dimension(4) :: cx, cy, cz
    real(kind=8   ), dimension(3) :: dr

! parameters
!
    real(kind=PREC) :: cc  = 299792457.99999998416751623153687d0   ! the speed of light [m/s]
    real(kind=PREC) :: pc  = 3.2407792896656065765177783686188d-17 ! 1 meter [pc]
    real(kind=PREC) :: sc  = 3.168876464084018437308447107767d-08  ! 1 second [yr]
!
!-------------------------------------------------------------------------------
!
! get dimain dimensions
!
    call get_dimensions(dm)

#ifdef TRICUB
    qm(:) = dm(:) + 2
#else /* TRUCUB */
    qm(:) = dm(:) + 1
#endif /* TRUCUB */
    if (dm(3) .eq. 1) qm(:) = 1

! get domain bounds
!
    call get_domain_bounds(bnds)

! calculate the domain size
!
    do p = 1, 3
      bsiz(p) = bnds(p,2) - bnds(p,1)
    end do

! compute plasma parameters
!                                                        ! c is expressed in Va
    dn   = 1.6726215850718025379202284485224d-21 * dens  ! density conversion from
                                                         ! protonmass/cm^3 to kg/m^3
    mu0  = 125.66370614359171042906382353976             ! magnetic permeability [Gs^2 m s^2 / kg]
    if (c .le. 1.0d0) then
      gm   = 1.0d0
      va   = 1.0d0 * cc
      bavg = bunit
    else
      gm   = 1.0d0 / sqrt(1.0d0 - (1.0 / c)**2)          ! Lorentz factor
      va   = gm * cc  / c                                ! Alfven speed [m/s]
      bavg = va * sqrt(mu0 * dn)                         ! magnetic field strength [Gs]
    end if
    csq  = c * c                                         ! square of the speed of light

! initialize particle parameters
!
    select case(ptype)
    case ('e')
      mrest =  0.51099890307660134070033564057667d+00    ! rest energy of electron [MeV]
      qom   = -1.75882017226579003036022186279300d+07    ! e/m [1 / Gs s]
      mp    =  9.10938188715453137087986438336060d-31    ! electron mass [kg]
    case default
      mrest =  0.93827199893682302445085952058434d+03    ! rest energy of proton   [MeV]
      qom   =  0.95788340668294185888953506946564d+04    ! e/m [1 / Gs s]
      mp    =  1.67262158507180250864766404816270d-27    ! proton mass [kg]
    end select
    vp = cc * vpar                                       ! parallel particle speed
    vr = cc * vper                                       ! perpendicular particle speed
    vv = sqrt(vpar**2 + vper**2)                         ! absolute velocity
    gm = 1.0d0 / sqrt(1.0d0 - vv * vv)
    mu = 0.5d0 * mp * vr**2 / bavg                       ! magnetic moment [kg m^2 / s^2 Gs]
    om0   = abs(qom * bavg)                              ! classical gyrofrequency
    om = om0 / gm                                        ! relativistic gyrofrequency
    tg = 1.0d0 / om                                      ! gyroperiod
    tg = pi2 * tg
    rg = vr / om                                         ! gyroradius (Larmor radius)

! print plasma parametes
!
    write( *, "('INFO      : plasma parameters:')" )
    write( *, "('INFO      : c     =',1pe15.8,' [Va]')"       ) c
    write( *, "('INFO      : Va    =',1pe15.8,' [m / s]')"    ) va
    write( *, "('INFO      : dens  =',1pe15.8,' [u / cm^3] =',1pe15.8,' [kg / m^3]')" ) dens, dn
    write( *, "('INFO      : <B>   =',1pe15.8,' [G]')"        ) bavg

! print particle parameters
!
    write( *, "('INFO      : particle parameters:')" )
    select case(ptype)
    case ('e')
      write( *, "('INFO      : trajectory for electron')" )
    case default
      write( *, "('INFO      : trajectory for proton')" )
    end select
    write( *, "('INFO      : e/m   =',1pe15.8,' [1 / G s]')" ) qom
    write( *, "('INFO      : Vpar  =',1pe15.8,' [c] =',1pe15.8,' [m / s]')" ) vpar, vp
    write( *, "('INFO      : Vper  =',1pe15.8,' [c] =',1pe15.8,' [m / s]')" ) vper, vr
    write( *, "('INFO      : |V|   =',1pe15.8,' [c] =',1pe15.8,' [m / s]')" ) vv  , vv * cc
    write( *, "('INFO      : gamma =',1pe15.8)"              ) gm
    write( *, "('INFO      : Om    =',1pe15.8,' [1 / s]')"   ) om
    write( *, "('INFO      : Tg    =',1pe15.8,' [s]')"       ) tg
    write( *, "('INFO      : Rg    =',1pe15.8,' [m] =',1pe15.8,' [pc]')" ) rg, pc * rg
    write( *, "('INFO      : mu    =',1pe15.8,' [N m / Gs]')") mu
    write( *, "('INFO      : E0    =',1pe15.8,' [MeV]')"     ) mrest

! change time unit
!
    select case(tunit)
    case('u')
      fc = 1.0d-6
    case('s')
      fc = 1.0
    case('m')
      fc = 60.0
    case('h')
      fc = 3600.0
    case('d')
      fc = 86400.0
    case('w')
      fc = 604800.0
    case('y')
      fc = 31556925.974678400903940200805664
    case default
      fc = 1.0
    end select

    fc  = tmulti * fc
    qom = qom * fc

! calculate geometry parameters
!
    ln = va * fc                                         ! the size of the box

! print geometry parameters
!
    write( *, "('INFO      : geometry parameters:')" )
    write( *, "('INFO      : T     =',1pe15.8,' [s] =',1pe15.8,' [yr]')" ) fc, sc * fc
    write( *, "('INFO      : L     =',1pe15.8,' [m] =',1pe15.8,' [pc]')" ) ln, pc * ln
    write( *, "('INFO      : Rg/L  =',1pe15.8)" ) rg / ln
    write( *, "('INFO      : Tg/T  =',1pe15.8)" ) tg / fc

    write( *, "('INFO      : code units:')" )
    write( *, "('INFO      : e/m   =',1pe25.16)" ) qom * bavg

! write parameters to info.txt
!
    open  (10, file = 'info.txt', form = 'formatted', status = 'replace')

! print plasma parametes
!
    write (10, "('INFO      : plasma parameters:')" )
    write (10, "('INFO      : c     =',1pe15.8,' [Va]')"       ) c
    write (10, "('INFO      : Va    =',1pe15.8,' [m / s]')"    ) va
    write (10, "('INFO      : dens  =',1pe15.8,' [u / cm^3] =',1pe15.8,' [kg / m^3]')" ) dens, dn
    write (10, "('INFO      : <B>   =',1pe15.8,' [G]')"        ) bavg

    write (10, "('INFO      : particle parameters:')" )
    select case(ptype)
    case ('e')
      write (10, "('INFO      : trajectory for electron')" )
    case default
      write (10, "('INFO      : trajectory for proton')" )
    end select
    write (10, "('INFO      : e/m   =',1pe15.8,' [1 / G s]')" ) qom
    write (10, "('INFO      : Om    =',1pe15.8,' [1 / s]')"   ) om
    write (10, "('INFO      : Tg    =',1pe15.8,' [s]')"       ) tg
    write (10, "('INFO      : Vpar  =',1pe15.8,' [c] =',1pe15.8,' [m / s]')" ) vpar, vp
    write (10, "('INFO      : Vper  =',1pe15.8,' [c] =',1pe15.8,' [m / s]')" ) vper, vr
    write (10, "('INFO      : |V|   =',1pe15.8,' [c] =',1pe15.8,' [m / s]')" ) vv  , vv * cc
    write (10, "('INFO      : Rg    =',1pe15.8,' [m] =',1pe15.8,' [pc]')" ) rg, pc * rg
    write (10, "('INFO      : gamma =',1pe15.8)"              ) gm
    write (10, "('INFO      : mu    =',1pe15.8,' [N m / Gs]')") mu
    write (10, "('INFO      : E0    =',1pe15.8,' [MeV]')"     ) mrest

! print geometry parameters
!
    write (10, "('INFO      : geometry parameters:')" )
    write (10, "('INFO      : T     =',1pe15.8,' [s] =',1pe15.8,' [yr]')" ) fc, sc * fc
    write (10, "('INFO      : L     =',1pe15.8,' [m] =',1pe15.8,' [pc]')" ) ln, pc * ln
    write (10, "('INFO      : Rg/L  =',1pe15.8)" ) rg / ln
    write (10, "('INFO      : Tg/T  =',1pe15.8)" ) tg / fc

    write (10, "('INFO      : code units:')" )
    write (10, "('INFO      : e/m   =',1pe15.8)" ) qom * bavg

    close (10)

! convert e/m to the units of magnetic field
!
    qom = qom * bavg

! initial position and velocity
!
    x0(:) = (/ xc, yc, zc /)

#ifdef TEST
#ifdef WTEST
    bpar = sqrt(bini**2 - bamp**2)
    b(1) = bpar
    b(2) = bamp * cos(pi2 * freq * xc)
    b(3) = bamp * sin(pi2 * freq * xc)

    u(1) = 0.0
    u(2) = 0.0
    u(3) = 1.0
#endif /* WTEST */
#ifdef ITEST
! print infor about the island problem
!
    write( *, "('PROBLEM   : motion in the contracting magnetic island:')" )

! check if the initial field is non zero
!
    if (bini .eq. 0.0d0) then
      write (*, "('ERROR     : bini must be not zero!')" )
      stop
    end if

! prepare parameters for the magnetic field topology
!
    dl = bamp / bini

    if (abs(dl) .ge. 1.0) then
      write (*, "('ERROR     : parameter bamp must be smaller than bini!')" )
      stop
    end if

    ra = 1.0d0 + dl
    rb = 1.0d0 - dl

! calculate the eccentricity
!
    ec = dsqrt(1.0d0 - (rb / ra)**2)

! print info about the problem
!
    write (*, "('INFO      : magnetic field strength          =',1pe15.8)") bini
    write (*, "('INFO      : guide field strength             =',1pe15.8)") bshr
    write (*, "('INFO      : magnetic field perturbation      =',1pe15.8)") bamp
    write (*, "('INFO      : magnetic island eccentricity     =',1pe15.8)") ec
    write (*, "('INFO      : horizontal velocity perturbation =',1pe15.8)") vamp
    write (*, "('INFO      : vertical velocity perturbation   =',1pe15.8)") vamp * vrat

! calculate the parallel and perpendicular directions with respect to the
! local magnetic field at the initial particle position
!
    xt   = xc / ra
    yt   = yc / rb

    rt   = dsqrt(xt * xt + yt * yt)

    if (rt .gt. 0.0d0) then

      b(1) =   yt / rb / rt
      b(2) = - xt / ra / rt
      b(3) = bshr

      bb = dsqrt(dot_product(b(:), b(:)))

      if (bb .gt. 0.0d0) then
        b(:) = b(:) / bb
      else
        write (*, "('ERROR     : zero magnetic field at the initial position!')" )
        stop
      end if

      u(1) = 0.0d0
      u(2) = 0.0d0
      u(3) = 1.0d0

      ub   = dot_product(u(:), b(:))

      u(1) = u(1) - b(1) * ub
      u(2) = u(2) - b(2) * ub
      u(3) = u(3) - b(3) * ub

      uu   = dsqrt(dot_product(u(:), u(:)))

      if (uu .gt. 0.0d0) then
        u(:) = u(:) / uu
      else
        write (*, "('ERROR     : cannot determine perpendicular direction!')" )
        stop
      end if
    else
      write (*, "('ERROR     : initial position cannot be located in the origin!')" )
      stop
    end if

! store the parameters in the info file
!
    open (10, file = 'info.txt', form = 'formatted', position = 'append')
    write(10, "('PROBLEM   : motion in the contracting magnetic island:')" )
    write(10, "('INFO      : magnetic field strength          =',1pe15.8)") bini
    write(10, "('INFO      : guide field strength             =',1pe15.8)") bshr
    write(10, "('INFO      : magnetic field perturbation      =',1pe15.8)") bamp
    write(10, "('INFO      : magnetic island eccentricity     =',1pe15.8)") ec
    write(10, "('INFO      : horizontal velocity perturbation =',1pe15.8)") vamp
    write(10, "('INFO      : vertical velocity perturbation   =',1pe15.8)") vamp * vrat
    close(10)
#endif /* ITEST */
#else /* TEST */
! convert position to index
!
    call pos2index(x0, r0)

! prepare coefficients for interpolation
!
    call prepare_interpolation(r0, ii, jj, kk, dr, cx, cy, cz)

! interpolate field components at the particle position
!
    b(1) = interpolate(bx, ii, jj, kk, dr, cx, cy, cz)
    b(2) = interpolate(by, ii, jj, kk, dr, cx, cy, cz)
    b(3) = interpolate(bz, ii, jj, kk, dr, cx, cy, cz)
#endif /* TEST */

! calculate the direction of the local magnetic field
!
    bb = sqrt(dot_product(b, b))
    ba = bb
    if (bb .gt. 0.0d0) then
      b(:) = b(:) / bb
    else
      write( *, "('ERROR     : ',a)" ) "B=0 at the initial position! Choose another one."
      stop
    endif

#ifndef TEST
! calculate the perpendicular unit vector
!
    if (dm(3) .eq. 1) then
      w(1) = 0.0
      w(2) = 0.0
      w(3) = 1.0
    else
      call random_number(w)
      w = w - 0.5
    end if

    bb = sqrt(dot_product(w, w))
    if (bb .gt. 0.0d0) then
      w(:) = w(:) / bb
    else
      write( *, "('ERROR     : ',a)" ) "V=0 at the initial position! Choose another one."
      stop
    end if

    u(1) = w(2) * b(3) - w(3) * b(2)
    u(2) = w(3) * b(1) - w(1) * b(3)
    u(3) = w(1) * b(2) - w(2) * b(1)

    bb = sqrt(dot_product(u, u))
    if (bb .gt. 0.0d0) then
      u(:) = u(:) / bb
    else
      write( *, "('ERROR     : ',a)" ) "V=0 at the initial position! Choose another one."
      stop
    end if
#endif /* !TEST */

! calculate the initial velocity
!
    u0(:) = (vpar * b(:) + vper * u(:)) * c

! calculate the Lorentz factor of the initial state
!
#ifdef RELAT
    gm = 1.0 / sqrt(1.0d0 - dot_product(u0, u0) / csq)
#else
    gm = 1.0
#endif

! calculate the initial particle momentuum
!
    p0(:) = gm * u0(:)

! calculate particle energy
!
#ifdef RELAT
    en = gm * mrest
    ek = en - mrest
#else
    en = 0.5 * (vpar**2 + vper**2) * csq
    ek = en
#endif

! print headers and the initial values
!
    open  (10, file = 'output.dat', form = 'formatted', status = 'replace')
    write (10, "('#',1a20,19a22)") 'Time', 'X', 'Y', 'Z', 'Vx', 'Vy', 'Vz'     &
                                 , '|V| [c]', '|Vpar| [c]', '|Vper| [c]'       &
                                 , 'gamma', 'En [MeV]', 'Ek [MeV]'             &
                                 , '<B> [Gs]', 'Omega [1/s]'                   &
                                 , 'Tg [s]', 'Rg [m]', 'Tg [T]', 'Rg [L]'      &
                                 , 'Tolerance'
    write (10, "(20(1pe22.14))") 0.0, x0(1), x0(2), x0(3), u0(1), u0(2), u0(3) &
                                    , vv, vpar, vper-vper, gm, en, ek               &
                                    , bavg * ba, om, tg, rg, tg / fc, rg / ln  &
                                    , 1.0d-16
    close (10)

! prepare dump times
!
    if (output .eq. 'l') then
      n = ndumps * (dlog10(real(tmax,kind=8)) - dlog10(real(tmin,kind=8))) + 1
      allocate(tt(n))
      do p = 1, n
        tt(p) = 10.0d0**((p - 1.0) / ndumps + dlog10(real(tmin,kind=8)))
      enddo
    endif
!
!-------------------------------------------------------------------------------
!
  end subroutine init_particle
!
!===============================================================================
!
! finit_particle: subroutine deallocates the particle variables
!
!===============================================================================
!
  subroutine finit_particle()

    implicit none
!
!-------------------------------------------------------------------------------
!
    if (allocated(tt)) deallocate(tt)

!-------------------------------------------------------------------------------
!
  end subroutine finit_particle
!
!===============================================================================
!
! integrate_trajectory_rk4: subroutine integrates particle trajectory using
!                           the 4th order RK method
!
!===============================================================================
!
  subroutine integrate_trajectory_rk4()

    use params, only : c, tmin, tmax, rho, maxtol, dtini, dtmax, ndumps,    &
                       vpar, vper

    implicit none

! local variables
!
    logical                       :: keepon = .true.
    integer                       :: n, m
    real(kind=PREC)               ::    t1, t2, t3, t4, t5
    real(kind=PREC), dimension(3) :: x, x1, x2, x3, x4, x5
    real(kind=PREC), dimension(3) :: u, u1, u2, u3, u4, u5
    real(kind=PREC), dimension(3) :: p, p1, p2, p3, p4, p5
    real(kind=PREC), dimension(3) ::    k1, k2, k3, k4, k5
    real(kind=PREC), dimension(3) ::    l1, l2, l3, l4, l5
    real(kind=PREC), dimension(3) :: a, v, b
    real(kind=PREC)               :: gm, t, dt, ds, dtn
    real(kind=PREC)               :: en, ek, ua, ba, up, ur, om, tg, rg
    real(kind=PREC)               :: tol
!
!-------------------------------------------------------------------------------
!
! initialize counters, time and timesteps
!
    n  = 0
    m  = 0
    t  = 0.0d0
    dt = dtini
    ds = qom * dt

! set the initial position, velocity, and momentum
!
    x(:) = x0(:)
    u(:) = u0(:)
    p(:) = p0(:)

! calculate the acceleration at the initial position
!
    call acceleration(t, x(:), u(:), a(:), v(:), b(:))

! separate the particle velocity into parallel and perpendicular components
!
    call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the Lorentz factor
!
    gm = lorentz_factor(p(:))

! calculate the particle energies
!
#ifdef RELAT
    en = gm * mrest
    ek = en - mrest
#else /* RELAT */
    en = 0.5d0 * ua * ua
    ek = en
#endif /* RELAT */

! print the progress
!
    write (*,"('PROGRESS  : ',a8,2x,4(a14))") 'ITER', 'TIME', 'TIMESTEP'       &
            , 'SPEED (c)', 'ENERGY (MeV)'
    write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6),a1,$)") n, t, dt, ua / c, ek    &
            , char(13)

! open the output file
!
    open  (10, file = 'output.dat', form = 'formatted', position = 'append')

!== INTEGRATION LOOP ==
!
! integrate the trajectory
!
    do while (keepon)

!! 1st step of the RK integration
!!
! integrate the position and momentum
!
      t1    = t
      x1(:) = x(:)
      p1(:) = p(:)

! calculate the Lorentz factor
!
      gm = lorentz_factor(p1(:))

! calculate the velocity
!
      u1(:) = p1(:) / gm

! calculate the acceleration for the location x1 and velocity u1
!
      call acceleration(t1, x1(:), u1(:), a(:), v(:), b(:))

! calculate the first term
!
      l1(:) = dt * u1(:)
      k1(:) = ds * a (:)

!! 2nd step of the RK integration
!!
! integrate the position and momentum
!
      t2    = t    + 0.5d0 * dt
      x2(:) = x(:) + 0.5d0 * l1(:)
      p2(:) = p(:) + 0.5d0 * k1(:)

! calculate the Lorentz factor
!
      gm = lorentz_factor(p2(:))

! calculate the velocity
!
      u2(:) = p2(:) / gm

! calculate the acceleration for the location x2 and velocity u2
!
      call acceleration(t2, x2(:), u2(:), a(:), v(:), b(:))

! calculate the second term
!
      l2(:) = dt * u2(:)
      k2(:) = ds * a (:)

!! 3rd step of the RK integration
!!
! integrate the position and momentum
!
      t3    = t    + 0.5d0 * dt
      x3(:) = x(:) + 0.5d0 * l2(:)
      p3(:) = p(:) + 0.5d0 * k2(:)

! calculate the Lorentz factor
!
      gm = lorentz_factor(p3(:))

! calculate the velocity
!
      u3(:) = p3(:) / gm

! calculate the acceleration for the location x3 and velocity u3
!
      call acceleration(t3, x3(:), u3(:), a(:), v(:), b(:))

! calculate the third term
!
      l3(:) = dt * u3(:)
      k3(:) = ds * a (:)

!! 4th step of the RK integration
!!
! integrate the position and momentum
!
      t4    = t    + dt
      x4(:) = x(:) + l3(:)
      p4(:) = p(:) + k3(:)

! calculate the Lorentz factor
!
      gm = lorentz_factor(p4(:))

! calculate the velocity
!
      u4(:) = p4(:) / gm

! calculate the acceleration for the location x4 and velocity u4
!
      call acceleration(t4, x4(:), u4(:), a(:), v(:), b(:))

! calculate the third term
!
      l4(:) = dt * u4(:)
      k4(:) = ds * a (:)

!! the final integration of the particle position and momentum
!!
      t5    = t    + dt
      x5(:) = x(:) + ( l1(:) + 2.0d0 * ( l2(:) + l3(:) ) + l4(:) ) / 6.0d0
      p5(:) = p(:) + ( k1(:) + 2.0d0 * ( k2(:) + k3(:) ) + k4(:) ) / 6.0d0

! calculate the Lorentz factor
!
      gm = lorentz_factor(p5(:))

! calculate the velocity
!
      u5(:) = p5(:) / gm

! calculate the acceleration at the updated location
!
      call acceleration(t5, x5(:), u5(:), a(:), v(:), b(:))

! estimate the error for timestep control
!
      l4(:) = l4(:) - dt * u5(:)
      k4(:) = k4(:) - ds * a (:)

      tol = sqrt(dot_product(l4(:), l4(:)) + dot_product(k4(:), k4(:))) / 6.0d0

! estimate the new timestep
!
      dtn   = dt * (rho * maxtol / tol)**0.2d0

! check if the error is below desired tolerance
!
      if (tol .gt. maxtol) then

! repeat the integration with a new timestep
!
        dt = dtn
        ds = qom * dt

      else

! update the time
!
        t   = t + dt

! check if time exceeded the maximum time
!
        if (t >= tmax) keepon = .false.

! update the new timestep
!
        dt = min(2.0d0 * dt, dtn, dtmax)
        ds = qom * dt

! update the position, velocity and momentum
!
        x(:) = x5(:)
        u(:) = u5(:)
        p(:) = p5(:)

#ifndef PERIODIC
! if the boundaries are not periodic and particle is out of the box, stop
! the integration
        if (x(1) < bnds(1,1)) keepon = .false.
        if (x(1) > bnds(1,2)) keepon = .false.
        if (x(2) < bnds(2,1)) keepon = .false.
        if (x(2) > bnds(2,2)) keepon = .false.
        if (x(3) < bnds(3,1)) keepon = .false.
        if (x(3) > bnds(3,2)) keepon = .false.
#endif /* PERIODIC */

! store the current particle state
!
        if (m .eq. ndumps) then

! separate the particle velocity into the parallel and perpendicular components
!
          call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the particle gyroperiod and gyroradius
!
          call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate the particle energies
!
#ifdef RELAT
          en = gm * mrest
          ek = en - mrest
#else /* RELAT */
          en = 0.5d0 * ua * ua
          ek = en
#endif /* RELAT */

! print the progress
!
          write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6),a1,$)") n, t, dt, ua / c  &
                  , ek, char(13)

! store the particle parameters
!
          write (10, "(20(1pe22.14))") t, x(1), x(2), x(3), u(1), u(2), u(3)   &
                                     , ua / c, up / c, ur / c - vper, gm, en, ek      &
                                     , bavg * ba, om, tg * fc, rg * ln, tg, rg &
                                     , tol
!           close (10)

! update the counters
!
          n = n + 1
          m = 0

        end if

! increase the data write counter
!
        m = m + 1

      end if

    end do

! separate the particle velocity into the parallel and perpendicular components
!
    call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the particle gyroperiod and gyroradius
!
    call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate the particle energies
!
#ifdef RELAT
    en = gm * mrest
    ek = en - mrest
#else /* RELAT */
    en = 0.5d0 * ua * ua
    ek = en
#endif /* RELAT */

! print the progress
!
    write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6))") n, t, dt, ua / c, ek

! store the particle parameters
!
    write (10, "(20(1pe22.14))") t, x(1), x(2), x(3), u(1), u(2), u(3)         &
                               , ua / c, up / c, ur / c - vper, gm, en, ek            &
                               , bavg * ba, om, tg * fc, rg * ln, tg, rg, tol

! close the output file
!
    close (10)

!-------------------------------------------------------------------------------
!
  end subroutine integrate_trajectory_rk4
!
!===============================================================================
!
! integrate_trajectory_rk4_log: subroutine integrates particle trajectory using
!                               the 4th order method with snapshots equally
!                               distributed in the logarithm of time
!
!===============================================================================
!
  subroutine integrate_trajectory_rk4_log()

    use params, only : c, tmin, tmax, rho, maxtol, dtini, dtmax, ndumps,    &
                       vpar, vper

    implicit none

! local variables
!
    logical                       :: keepon = .true.
    integer                       :: n
    real(kind=PREC)               ::    t1, t2, t3, t4, t5, tp
    real(kind=PREC), dimension(3) :: x, x1, x2, x3, x4, x5, xt
    real(kind=PREC), dimension(3) :: u, u1, u2, u3, u4, u5, ut
    real(kind=PREC), dimension(3) :: p, p1, p2, p3, p4, p5, pt
    real(kind=PREC), dimension(3) ::    k1, k2, k3, k4, k5
    real(kind=PREC), dimension(3) ::    l1, l2, l3, l4, l5
    real(kind=PREC), dimension(3) :: a, v, b
    real(kind=PREC)               :: gm, t, dt, ds, dtn
    real(kind=PREC)               :: en, ek, ua, ba, up, ur, om, tg, rg
    real(kind=PREC)               :: tol, wl, wr
!
!-------------------------------------------------------------------------------
!
! initialize counters, time and timesteps
!
    n  = 0
    t  = 0.0d0
    dt = dtini
    ds = qom * dt

! set the initial position, velocity, and momentum
!
    x(:) = x0(:)
    u(:) = u0(:)
    p(:) = p0(:)

! calculate the acceleration at the initial position
!
    call acceleration(t, x(:), u(:), a(:), v(:), b(:))

! separate the particle velocity into parallel and perpendicular components
!
    call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the Lorentz factor
!
    gm = lorentz_factor(p(:))

! calculate the particle energies
!
#ifdef RELAT
    en = gm * mrest
    ek = en - mrest
#else /* RELAT */
    en = 0.5d0 * ua * ua
    ek = en
#endif /* RELAT */

! print the progress
!
    write (*,"('PROGRESS  : ',a8,2x,4(a14))") 'ITER', 'TIME', 'TIMESTEP'       &
            , 'SPEED (c)', 'ENERGY (MeV)'
    write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6),a1,$)") n, t, dt, ua / c, ek    &
            , char(13)

!== INTEGRATION LOOP ==
!
! integrate the trajectory
!
    do while (keepon)

!! 1st step of the RK integration
!!
! integrate the position and momentum
!
      t1    = t
      x1(:) = x(:)
      p1(:) = p(:)

! calculate the Lorentz factor
!
      gm = lorentz_factor(p1(:))

! calculate the velocity
!
      u1(:) = p1(:) / gm

! calculate the acceleration for the location x1 and velocity u1
!
      call acceleration(t1, x1(:), u1(:), a(:), v(:), b(:))

! calculate the first term
!
      l1(:) = dt * u1(:)
      k1(:) = ds * a (:)

!! 2nd step of the RK integration
!!
! integrate the position and momentum
!
      t2    = t    + 0.5d0 * dt
      x2(:) = x(:) + 0.5d0 * l1(:)
      p2(:) = p(:) + 0.5d0 * k1(:)

! calculate the Lorentz factor
!
      gm = lorentz_factor(p2(:))

! calculate the velocity
!
      u2(:) = p2(:) / gm

! calculate the acceleration for the location x2 and velocity u2
!
      call acceleration(t2, x2(:), u2(:), a(:), v(:), b(:))

! calculate the second term
!
      l2(:) = dt * u2(:)
      k2(:) = ds * a (:)

!! 3rd step of the RK integration
!!
! integrate the position and momentum
!
      t3    = t    + 0.5d0 * dt
      x3(:) = x(:) + 0.5d0 * l2(:)
      p3(:) = p(:) + 0.5d0 * k2(:)

! calculate the Lorentz factor
!
      gm = lorentz_factor(p3(:))

! calculate the velocity
!
      u3(:) = p3(:) / gm

! calculate the acceleration for the location x3 and velocity u3
!
      call acceleration(t3, x3(:), u3(:), a(:), v(:), b(:))

! calculate the third term
!
      l3(:) = dt * u3(:)
      k3(:) = ds * a (:)

!! 4th step of the RK integration
!!
! integrate the position and momentum
!
      t4    = t    + 0.5d0 * dt
      x4(:) = x(:) + l3(:)
      p4(:) = p(:) + k3(:)

! calculate the Lorentz factor
!
      gm = lorentz_factor(p4(:))

! calculate the velocity
!
      u4(:) = p4(:) / gm

! calculate the acceleration for the location x4 and velocity u4
!
      call acceleration(t4, x4(:), u4(:), a(:), v(:), b(:))

! calculate the third term
!
      l4(:) = dt * u4(:)
      k4(:) = ds * a (:)

!! the final integration of the particle position and momentum
!!
      t5    = t    + dt
      x5(:) = x(:) + ( l1(:) + 2.0d0 * ( l2(:) + l3(:) ) + l4(:) ) / 6.0d0
      p5(:) = p(:) + ( k1(:) + 2.0d0 * ( k2(:) + k3(:) ) + k4(:) ) / 6.0d0

! calculate the Lorentz factor
!
      gm = lorentz_factor(p5(:))

! calculate the velocity
!
      u5(:) = p5(:) / gm

! calculate the acceleration at the updated location
!
      call acceleration(t5, x5(:), u5(:), a(:), v(:), b(:))

! estimate the error for timestep control
!
      l4(:) = l4(:) - dt * u5(:)
      k4(:) = k4(:) - ds * a (:)

      tol = sqrt(dot_product(l4(:), l4(:)) + dot_product(k4(:), k4(:))) / 6.0d0

! estimate the new timestep
!
      dtn   = dt * (rho * maxtol / tol)**0.2d0

! check if the error is below desired tolerance
!
      if (tol .gt. maxtol) then

! repeat the integration with a new timestep
!
        dt = dtn
        ds = qom * dt

      else

! update the time
!
        t   = t + dt

! check if time exceeded the maximum time
!
        if (t >= tmax) keepon = .false.

! store the intermidiate particle states
!
        do while (tt(n) .le. t .and. t .ge. tmin .and. t .lt. tmax)

! calculate the left and right weights
!
          wl = (t - tt(n)) / dt
          wr = 1.0d0 - wl

! interpolate the particle state at the proper time
!
          tp    = wl * (t - dt) + wr * t
          xt(:) = wl * x(:) + wr * x5(:)
          ut(:) = wl * u(:) + wr * u5(:)
          pt(:) = wl * p(:) + wr * p5(:)

! calculate acceleration for the location x4 and velocity v4
!
          call acceleration(tp, xt(:), ut(:), a(:), v(:), b(:))

! separate the particle velocity into the parallel and perpendicular components
!
          call separate_velocity(ut(:), b(:), ba, ua, up, ur)

! calculate the particle gyroperiod and gyroradius
!
          call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate the particle energies
!
#ifdef RELAT
          en = gm * mrest
          ek = en - mrest
#else /* RELAT */
          en = 0.5d0 * ua * ua
          ek = en
#endif /* RELAT */

! print the progress
!
          write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6),a1,$)") n, t, dt, ua / c  &
                  , ek, char(13)

! store the particle parameters
!
          open  (10, file = 'output.dat', form = 'formatted'                   &
                   , position = 'append')
          write (10, "(20(1pe22.14))") t, x(1), x(2), x(3), u(1), u(2), u(3)   &
                                     , ua / c, up / c, ur / c, gm, en, ek      &
                                     , bavg * ba, om, tg * fc, rg * ln, tg, rg &
                                     , tol
          close (10)

! update the counters
!
          n = n + 1

        end do

! update the position, velocity and momentum
!
        x(:) = x5(:)
        u(:) = u5(:)
        p(:) = p5(:)

#ifndef PERIODIC
! if the boundaries are not periodic and particle is out of the box, stop
! the integration
        if (x(1) < bnds(1,1)) keepon = .false.
        if (x(1) > bnds(1,2)) keepon = .false.
        if (x(2) < bnds(2,1)) keepon = .false.
        if (x(2) > bnds(2,2)) keepon = .false.
        if (x(3) < bnds(3,1)) keepon = .false.
        if (x(3) > bnds(3,2)) keepon = .false.
#endif /* PERIODIC */

! update the new timestep
!
        dt = min(2.0d0 * dt, dtn, dtmax, max(1.0d-16, tmax - t))
        ds = qom * dt

      end if

    end do

! separate the particle velocity into the parallel and perpendicular components
!
    call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the particle gyroperiod and gyroradius
!
    call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate the particle energies
!
#ifdef RELAT
    en = gm * mrest
    ek = en - mrest
#else /* RELAT */
    en = 0.5d0 * ua * ua
    ek = en
#endif /* RELAT */

! print the progress
!
    write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6))") n, t, dt, ua / c, ek

! store the particle parameters
!
    open  (10, file = 'output.dat', form = 'formatted', position = 'append')
    write (10, "(20(1pe22.14))") t, x(1), x(2), x(3), u(1), u(2), u(3)         &
                               , ua / c, up / c, ur / c, gm, en, ek            &
                               , bavg * ba, om, tg * fc, rg * ln, tg, rg, tol
    close (10)
!
!-------------------------------------------------------------------------------
!
  end subroutine integrate_trajectory_rk4_log
!
!===============================================================================
!
! integrate_trajectory_si4: subroutine integrates particle trajectory using
!                           the 4th order simplectic method
!
! references: "Numerical Hamiltonian Problems", J. M. Sanz-Serna & M. P. Calvo
!             Chapman & Hall, London, New York, 1994
!             "High order starting iterates for implicit Runge-Kutta methods:
!              an improvement for variable-step symplectic integrators", 2002,
!              IMA J. of Num. Ana., 22, 153
!
!===============================================================================
!
  subroutine integrate_trajectory_si4()

    use params, only : dtini, tmax, c, ndumps, vper

    implicit none

! local variables
!
    character(len=32)               :: str
    integer                         :: n, m, i, mi, ti
    real(kind=PREC), dimension(2,6) :: z, zp
    real(kind=PREC), dimension(3)   :: x , u , p , a
    real(kind=PREC), dimension(3)   :: v, b
    real(kind=PREC)                 :: gm, t, dt, dq, s, ds
    real(kind=PREC)                 :: en, ek, ua, ba, up, ur, om, tg, rg
    real(kind=PREC)                 :: tol

! local flags
!
    logical                         :: flag   = .true.
    logical                         :: keepon = .true.

! local parameters
!
    real(kind=PREC), parameter :: b1   = - dsqrt(3.0d0)                        &
                                , b2   =   dsqrt(3.0d0)
    real(kind=PREC), parameter :: c1   =   0.5d0 - dsqrt(3.0d0) / 6.0d0        &
                                , c2   =   0.5d0 + dsqrt(3.0d0) / 6.0d0
    real(kind=PREC), parameter :: b11  =   1.0d0 - 2.0d0 * dsqrt(3.0d0)        &
                                , b12  = - 6.0d0 + 4.0d0 * dsqrt(3.0d0)        &
                                , b21  = - 6.0d0 - 4.0d0 * dsqrt(3.0d0)        &
                                , b22  =   1.0d0 + 2.0d0 * dsqrt(3.0d0)
!
!-------------------------------------------------------------------------------
!
! initialize the iteration number, snapshot number, time, and time steps
!
    n  = 1
    m  = 1
    mi = 0
    ti = 0
    t  = 0.0d0
    s  = 0.0d0
    dt = dtini
    dq = qom * dt
    ds = dt * ndumps

! substitute the initial position, velocity, and momentum
!
    x(:) = x0(:)
    u(:) = u0(:)
    p(:) = p0(:)

! calculate the acceleration at the starting point
!
    call acceleration(t, x(:), u(:), a(:), v(:), b(:))

! separate particle velocity into parallel and perpendicular components
!
    call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the Lorentz factor
!
    gm = lorentz_factor(p(:))

! calculate the particle energy
!
#ifdef RELAT
    en = gm * mrest
    ek = en - mrest
#else /* RELAT */
    en = 0.5d0 * ua * ua
    ek = en
#endif /* RELAT */

! print the progress information
!
    write (*,"('PROGRESS  : ',a8,2x,4(a14))") 'ITER', 'TIME', 'TIMESTEP'       &
            , 'SPEED (c)', 'ENERGY (MeV)'
    write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6),a1,$)") n, t, dt, ua / c, ek    &
            , char(13)

! open the output file
!
    open  (10, file = 'output.dat', form = 'formatted', position = 'append')

!== INTEGRATION LOOP ==
!
! iterate until the maximum time is reached
!
    do while (keepon)

! find the initial guess for the vector Z
!
      if (flag) then

        z(1,1:3) = c1 * u(1:3)
        z(2,1:3) = c2 * u(1:3)
        z(1,4:6) = c1 * a(1:3)
        z(2,4:6) = c2 * a(1:3)

        flag = .false.
      else
        zp(:,:) = z(:,:)

        z(1,1:6) = b11 * zp(1,1:6) + b12 * zp(2,1:6)
        z(2,1:6) = b21 * zp(1,1:6) + b22 * zp(2,1:6)
      end if

! estimate the vector Z (eq. 5.3)
!
!   Z1 = dt * [ a11 * F(y + Z1) + a12 * F(y + Z2) ]
!   Z2 = dt * [ a21 * F(y + Z1) + a22 * F(y + Z2) ]
!
      call estimate_si4(x(:), p(:), z(:,:), t, dt, dq, tol, i)

! update the solution
!
!   y(n+1) = y(n) + [ b1 * Z1 + b2 * Z2 ]
!
      x(1:3) = x(1:3) + dt * (b1 * z(1,1:3) + b2 * z(2,1:3))
      p(1:3) = p(1:3) + dq * (b1 * z(1,4:6) + b2 * z(2,4:6))

#ifndef PERIODIC
! if the boundaries are not periodic and particle is out of the box, stop
! the integration
      if (x(1) < bnds(1,1)) keepon = .false.
      if (x(1) > bnds(1,2)) keepon = .false.
      if (x(2) < bnds(2,1)) keepon = .false.
      if (x(2) > bnds(2,2)) keepon = .false.
      if (x(3) < bnds(3,1)) keepon = .false.
      if (x(3) > bnds(3,2)) keepon = .false.
#endif /* PERIODIC */

! update the integration time
!
      t = s + m * dt

! check if time exceeded the maximum time
!
      if (t >= tmax) keepon = .false.

! find the maximum number of iteration in the estimator and update the counter
! of the total number of iterations
!
      mi = max(mi, i)
      ti = ti + i

! store the particle parameters at a given snapshot time
!
      if (m .eq. ndumps) then

! calculate the Lorentz factor and particle velocity
!
        gm   = lorentz_factor(p(:))
        u(:) = p(:) / gm

! calculate the acceleration at the locations x1 and x2
!
        call acceleration(t, x(:), u(:), a(:), v(:), b(:))

! separate particle velocity into parallel and perpendicular components
!
        call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the particle gyroperiod and gyroradius
!
        call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate particle energy
!
#ifdef RELAT
        en = gm * mrest
        ek = en - mrest
#else /* RELAT */
        en = 0.5d0 * ua * ua
        ek = en
#endif /* RELAT */

! update the integration time
!
        s = n * ds

! write the progress
!
        write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6),a1,$)") n, t, dt, ua / c    &
            , ek, char(13)

! write results to the output file
!
        write (10, "(20(1pe22.14),i10)") t                                     &
                                   , x(1), x(2), x(3), u(1), u(2), u(3)        &
                                   , ua / c, up / c, ur / c - vper, gm, en, ek        &
                                   , bavg * ba, om, tg * fc, rg * ln, tg, rg   &
                                   , tol, i

        n = n + 1
        m = 0

      end if

! increase data write counter
!
      m = m + 1

! end of iteration
!
    end do

! close the output file
!
    close(10)

! write the progress
!
    write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6))") n, t, dt, ua / c, ek

! write info about the estimator
!
    write(str,"(i12)") mi
    write(*,"('INFO      : maximum iterations per step = ',a)"      )          &
          trim(adjustl(str))
    write(*,"('INFO      : average iterations per step = ',1pe12.6)")          &
          real(ti, kind=8) / ((n - 1) * ndumps)

! open the info file
!
    open  (11, file = 'info.txt', form = 'formatted', position = 'append')

! write info about the estimator
!
    write(11,"('INFO      : maximum iterations per step = ',a)"      )         &
          trim(adjustl(str))
    write(11,"('INFO      : average iterations per step = ',1pe12.6)")         &
          real(ti, kind=8) / ((n - 1) * ndumps)

! close the info file
!
    close(11)
!
!-------------------------------------------------------------------------------
!
  end subroutine integrate_trajectory_si4
!
!===============================================================================
!
! integrate_trajectory_si4_log: subroutine integrates particle trajectory using
!             the 4th order simplectic method with snapshots equally
!             distributed in the logarithm of time
!
! references: "Numerical Hamiltonian Problems", J. M. Sanz-Serna & M. P. Calvo
!             Chapman & Hall, London, New York, 1994
!             "High order starting iterates for implicit Runge-Kutta methods:
!              an improvement for variable-step symplectic integrators", 2002,
!              IMA J. of Num. Ana., 22, 153
!
!===============================================================================
!
  subroutine integrate_trajectory_si4_log()

    use params, only : dtini, tmin, tmax, c

    implicit none

! local variables
!
    integer                         :: n, m, i, mi
    real(kind=PREC), dimension(2,6) :: z, zp
    real(kind=PREC), dimension(3)   :: x , u , p , a
    real(kind=PREC), dimension(3)   :: xn, pn, xt, ut, pt
    real(kind=PREC), dimension(3)   :: v, b
    real(kind=PREC)                 :: gm, t, dt, ds
    real(kind=PREC)                 :: en, ek, ua, ba, up, ur, om, tg, rg
    real(kind=PREC)                 :: tol, tp, wl, wr

! local flags
!
    logical                         :: flag   = .true.
    logical                         :: keepon = .true.

! local parameters
!
    real(kind=PREC), parameter :: b1   = - dsqrt(3.0d0)                        &
                                , b2   =   dsqrt(3.0d0)
    real(kind=PREC), parameter :: c1   =   0.5d0 - dsqrt(3.0d0) / 6.0d0        &
                                , c2   =   0.5d0 + dsqrt(3.0d0) / 6.0d0
    real(kind=PREC), parameter :: b11  =   1.0d0 - 2.0d0 * dsqrt(3.0d0)        &
                                , b12  = - 6.0d0 + 4.0d0 * dsqrt(3.0d0)        &
                                , b21  = - 6.0d0 - 4.0d0 * dsqrt(3.0d0)        &
                                , b22  =   1.0d0 + 2.0d0 * dsqrt(3.0d0)
!
!-------------------------------------------------------------------------------
!
! initialize the iteration number, snapshot number, time, and time steps
!
    n  = 1
    m  = 1
    t  = 0.0d0
    dt = dtini
    ds = qom * dt

! substitute the initial position, velocity, and momentum
!
    x(:) = x0(:)
    u(:) = u0(:)
    p(:) = p0(:)

! calculate the acceleration at the starting point
!
    call acceleration(t, x(:), u(:), a(:), v(:), b(:))

! separate particle velocity into parallel and perpendicular components
!
    call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the Lorentz factor
!
    gm = lorentz_factor(p(:))

! calculate the particle energy
!
#ifdef RELAT
    en = gm * mrest
    ek = en - mrest
#else /* RELAT */
    en = 0.5 * ua * ua
    ek = en
#endif /* RELAT */

! print the progress information
!
    write (*,"('PROGRESS  : ',a8,2x,4(a14))") 'ITER', 'TIME', 'TIMESTEP'       &
            , 'SPEED (c)', 'ENERGY (MeV)'
    write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6),a1,$)") n, t, dt, ua / c, ek    &
            , char(13)

!== INTEGRATION LOOP ==
!
! iterate until the maximum time is reached
!
    do while (keepon)

! find the initial guess for the vector Z
!
!      if (flag) then

        z(1,1:3) = c1 * u(1:3)
        z(2,1:3) = c2 * u(1:3)
        z(1,4:6) = c1 * a(1:3)
        z(2,4:6) = c2 * a(1:3)

!        flag = .false.
!      else
!        zp(:,:) = z(:,:)

!        z(1,1:6) = b11 * zp(1,1:6) + b12 * zp(2,1:6)
!        z(2,1:6) = b21 * zp(1,1:6) + b22 * zp(2,1:6)
!      end if

! estimate the vector Z (eq. 5.3)
!
!   Z1 = dt * [ a11 * F(y + Z1) + a12 * F(y + Z2) ]
!   Z2 = dt * [ a21 * F(y + Z1) + a22 * F(y + Z2) ]
!
      call estimate_si4(x(:), p(:), z(:,:), t, dt, ds, tol, i)

! update the solution
!
!   y(n+1) = y(n) + [ b1 * Z1 + b2 * Z2 ]
!
      xn(1:3) = x(1:3) + dt * (b1 * z(1,1:3) + b2 * z(2,1:3))
      pn(1:3) = p(1:3) + ds * (b1 * z(1,4:6) + b2 * z(2,4:6))

! update the integration time
!
      t = t + dt

! check if time exceeded the maximum time
!
      if (t >= tmax) keepon = .false.

! store the intermediate snapshots
!
      do while (tt(n) .le. t .and. t .ge. tmin .and. t .lt. tmax)

! calculate the left and right weights
!
        wl = (t - tt(n)) / dt
        wr = 1.0d0 - wl

! interpolate the particle state at the proper time
!
        tp    = wl * (t - dt) + wr * t
        xt(:) = wl * x(:) + wr * xn(:)
        pt(:) = wl * p(:) + wr * pn(:)

! calculate the Lorentz factor and particle velocity
!
        gm    = lorentz_factor(pt(:))
        ut(:) = pt(:) / gm

! calculate acceleration for the location x4 and velocity v4
!
        call acceleration(tp, xt(:), ut(:), a(:), v(:), b(:))

! separate particle velocity into parallel and perpendicular components
!
        call separate_velocity(ut(:), b(:), ba, ua, up, ur)

! calculate the particle gyroperiod and gyroradius
!
        call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate particle energy
!
#ifdef RELAT
        en = gm * mrest
        ek = en - mrest
#else
        en = 0.5 * ua * ua
        ek = en
#endif

! write the progress
!
        write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6),a1,$)") n, tp, dt, ua / c   &
            , ek, char(13)

! write results to the output file
!
        open  (10, file = 'output.dat', form = 'formatted', position = 'append')
        write (10, "(20(1pe22.14))") tp, x(1), x(2), x(3), u(1), u(2), u(3)    &
                                   , ua / c, up / c, ur / c, gm, en, ek        &
                                   , bavg * ba, om, tg * fc, rg * ln, tg, rg   &
                                   , tol
        close (10)

! increate the snapshot index
!
        n = n + 1

      end do

! substitute the new particle position and momentum
!
      x(:) = xn(:)
      p(:) = pn(:)

#ifndef PERIODIC
! if the boundaries are not periodic and particle is out of the box, stop
! the integration
      if (x(1) < bnds(1,1)) keepon = .false.
      if (x(1) > bnds(1,2)) keepon = .false.
      if (x(2) < bnds(2,1)) keepon = .false.
      if (x(2) > bnds(2,2)) keepon = .false.
      if (x(3) < bnds(3,1)) keepon = .false.
      if (x(3) > bnds(3,2)) keepon = .false.
#endif /* PERIODIC */

! update timestep
!
      dt = min(dt, max(1.0d-16, tmax - t))
      ds = qom * dt

! end of iteration
!
    end do

! separate particle velocity into parallel and perpendicular components
!
    call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the particle gyroperiod and gyroradius
!
    call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate particle energy
!
#ifdef RELAT
    en = gm * mrest
    ek = en - mrest
#else
    en = 0.5d0 * ua * ua
    ek = en
#endif

! print the progress info
!
    write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6))") n, t, dtini, ua / c, ek

! store the last snapshot
!
    open  (10, file = 'output.dat', form = 'formatted', position = 'append')
    write (10, "(20(1pe22.14))") t, x(1), x(2), x(3), v(1), v(2), v(3)         &
                               , ua / c, up / c, ur / c, gm, en, ek            &
                               , bavg * ba, om, tg * fc, rg * ln, tg, rg, tol
    close (10)
!
!-------------------------------------------------------------------------------
!
  end subroutine integrate_trajectory_si4_log
!
!===============================================================================
!
! estimate_si4: subroutine estimates the solution for the equation of motion
!               using a simple functional iteration (SI4 version)
!
! references: "Numerical Hamiltonian Problems", J. M. Sanz-Serna & M. P. Calvo
!             Chapman & Hall, London, New York, 1994
!
! description: This subroutines find the solution of the equation 5.3 for
!              the increment Z using the functional iteration
!
!===============================================================================
!
  subroutine estimate_si4(x, p, z, t, dt, dq, tol, it)

    use params, only : maxit, maxeps

    implicit none

! subroutine arguments
!
    real(kind=PREC), dimension(3)  , intent(in)    :: x, p
    real(kind=PREC), dimension(2,6), intent(inout) :: z
    real(kind=PREC)                , intent(in)    :: t
    real(kind=PREC)                , intent(inout) :: dt, dq, tol
    integer                        , intent(inout) :: it

! local variables
!
    real(kind=PREC), dimension(2,6) :: zn
    real(kind=PREC), dimension(6)   :: dh
    real(kind=PREC), dimension(3)   :: x1, p1, u1, a1
    real(kind=PREC), dimension(3)   :: x2, p2, u2, a2
    real(kind=PREC), dimension(3)   :: v, b
    real(kind=PREC)                 :: g1, g2, eps

! local parameter
!
    real(kind=PREC), parameter :: a11 = 1.0d0 / 4.0d0                          &
                                , a12 = 1.0d0 / 4.0d0 - dsqrt(3.0d0) / 6.0d0   &
                                , a21 = 1.0d0 / 4.0d0 + dsqrt(3.0d0) / 6.0d0   &
                                , a22 = 1.0d0 / 4.0d0
    real(kind=PREC), parameter :: e1  = - dsqrt(3.0d0)                         &
                                , e2  =   dsqrt(3.0d0)
!
!-------------------------------------------------------------------------------
!
! initiate the iteration control parameters
!
    it  = 0
    eps = 1.0d+16

! perform the simple functional iteration until the conditions are met
!
    do while (eps .gt. maxeps .and. it .lt. maxit)

! prepare the particle position and momentum for the current iteration
!
      x1(:) = x(:) + dt * z(1,1:3)
      x2(:) = x(:) + dt * z(2,1:3)
      p1(:) = p(:) + dq * z(1,4:6)
      p2(:) = p(:) + dq * z(2,4:6)

! calculate the Lorentz factors and particle velocity
!
      g1    = lorentz_factor(p1(:))
      g2    = lorentz_factor(p2(:))
      u1(:) = p1(:) / g1
      u2(:) = p2(:) / g2

! calculate the accelerations
!
      call acceleration(t, x1(1:3), u1(1:3), a1(1:3), v(1:3), b(1:3))
      call acceleration(t, x2(1:3), u2(1:3), a2(1:3), v(1:3), b(1:3))

! update the increment
!
      zn(1,1:3) = a11 * u1(1:3) + a12 * u2(1:3)
      zn(1,4:6) = a11 * a1(1:3) + a12 * a2(1:3)
      zn(2,1:3) = a21 * u1(1:3) + a22 * u2(1:3)
      zn(2,4:6) = a21 * a1(1:3) + a22 * a2(1:3)

! calculate the maximum of residuum of the increment
!
      eps = maxval(abs(zn - z))

! substitute the new solution of the increment
!
      z = zn

! increase the iteration counter
!
      it = it + 1

    end do

! estimate the integration error
!
    dh(1:3) = dt * (e1 * u1(:) + e2 * u2(:))
    dh(4:6) = dq * (e1 * a1(:) + e2 * a2(:))
    tol     = sqrt(sum(dh(:) * dh(:)))

! if the convergence was not reached write the warning about it
!
    if (it .ge. maxit) then
      open (11, file = 'info.txt', form = 'formatted', position = 'append')
      write(11,"('WARNING   : convergence not reached at t =',1pe12.5," //     &
               "' eps =',1pe12.5,' tol =',1pe12.5)") t, eps, tol
      close(11)
    end if
!
!-------------------------------------------------------------------------------
!
  end subroutine estimate_si4
!
!===============================================================================
!
! integrate_trajectory_si6: subroutine integrates particle trajectory using
!                           the 6th order simplectic method
!
! references: "Numerical Hamiltonian Problems", J. M. Sanz-Serna & M. P. Calvo
!             Chapman & Hall, London, New York, 1994
!             "High order starting iterates for implicit Runge-Kutta methods:
!              an improvement for variable-step symplectic integrators", 2002,
!              IMA J. of Num. Ana., 22, 153
!
!===============================================================================
!
  subroutine integrate_trajectory_si6()

    use params, only : dtini, tmax, c, ndumps, vper

    implicit none

! local variables
!
    logical                         :: keepon = .true.
    character(len=32)               :: str
    integer                         :: n, m, i, mi, ti
    real(kind=PREC), dimension(3,6) :: z
    real(kind=PREC), dimension(3)   :: x , u , p , a
    real(kind=PREC), dimension(3)   :: v, b
    real(kind=PREC)                 :: gm, t, dt, s, ds, dq
    real(kind=PREC)                 :: en, ek, ua, ba, up, ur, om, tg, rg
    real(kind=PREC)                 :: tol

! local parameters
!
    real(kind=PREC), parameter :: b1   =   5.0d0 / 3.0d0                       &
                                , b2   = - 4.0d0 / 3.0d0                       &
                                , b3   =   5.0d0 / 3.0d0
    real(kind=PREC), parameter :: c1   =   0.5d0 - 0.1d0 * dsqrt(15.0d0)       &
                                , c2   =   0.5d0                               &
                                , c3   =   0.5d0 + 0.1d0 * dsqrt(15.0d0)
!
!-------------------------------------------------------------------------------
!
! initialize the iteration number, snapshot number, time, and time steps
!
    n  = 1
    m  = 1
    mi = 0
    ti = 0
    t  = 0.0d0
    s  = 0.0d0
    dt = dtini
    dq = qom * dt
    ds = dt * ndumps

! substitute the initial position, velocity, and momentum
!
    x(:) = x0(:)
    u(:) = u0(:)
    p(:) = p0(:)

! calculate the acceleration at the starting point
!
    call acceleration(t, x(:), u(:), a(:), v(:), b(:))

! separate particle velocity into parallel and perpendicular components
!
    call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the Lorentz factor
!
    gm = lorentz_factor(p(:))

! calculate the particle energy
!
#ifdef RELAT
    en = gm * mrest
    ek = en - mrest
#else /* RELAT */
    en = 0.5d0 * ua * ua
    ek = en
#endif /* RELAT */

! print the progress information
!
    write (*,"('PROGRESS  : ',a8,2x,4(a14))") 'ITER', 'TIME', 'TIMESTEP'       &
            , 'SPEED (c)', 'ENERGY (MeV)'
    write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6),a1,$)") n, t, dt, ua / c, ek    &
            , char(13)

! open the output file
!
    open  (10, file = 'output.dat', form = 'formatted', position = 'append')

!== INTEGRATION LOOP ==
!
! iterate until the maximum time is reached
!
    do while (keepon)

! obtain velocity and acceleration for the initial guess of Z
!
      gm   = lorentz_factor(p(:))
      u(:) = p(:) / gm
      call acceleration(t, x(:), u(:), a(:), v(:), b(:))

! find the initial guess for the vector Z (linear estimation)
!
      z(1,1:3) = c1 * dt * u(1:3)
      z(2,1:3) = c2 * dt * u(1:3)
      z(3,1:3) = c3 * dt * u(1:3)
      z(1,4:6) = c1 * dq * a(1:3)
      z(2,4:6) = c2 * dq * a(1:3)
      z(3,4:6) = c3 * dq * a(1:3)

! estimate the vector Z (eq. 5.3)
!
!   Z1 = [ a11 * F(y + Z1) + a12 * F(y + Z2) + a13 * F(y + Z3) ]
!   Z2 = [ a21 * F(y + Z1) + a22 * F(y + Z2) + a23 * F(y + Z3) ]
!   Z3 = [ a31 * F(y + Z1) + a32 * F(y + Z2) + a33 * F(y + Z3) ]
!
      call estimate_si6(x(:), p(:), z(:,:), t, dt, dq, tol, i)

! update the solution
!
!   y(n+1) = y(n) + [ d1 * Z1 + d2 * Z2 + d3 * Z3 ]
!
      x(1:3) = x(1:3) + (b1 * z(1,1:3) + b2 * z(2,1:3) + b3 * z(3,1:3))
      p(1:3) = p(1:3) + (b1 * z(1,4:6) + b2 * z(2,4:6) + b3 * z(3,4:6))

#ifndef PERIODIC
! if the boundaries are not periodic and particle is out of the box, stop
! the integration
      if (x(1) < bnds(1,1)) keepon = .false.
      if (x(1) > bnds(1,2)) keepon = .false.
      if (x(2) < bnds(2,1)) keepon = .false.
      if (x(2) > bnds(2,2)) keepon = .false.
      if (x(3) < bnds(3,1)) keepon = .false.
      if (x(3) > bnds(3,2)) keepon = .false.
#endif /* PERIODIC */

! update the integration time
!
      t = s + m * dt

! check if time exceeded the maximum time
!
      if (t >= tmax) keepon = .false.

! find the maximum number of iteration in the estimator and update the counter
! of the total number of iterations
!
      mi = max(mi, i)
      ti = ti + i

! store the particle parameters at a given snapshot time
!
      if (m .eq. ndumps) then

! calculate the Lorentz factor and particle velocity
!
        gm   = lorentz_factor(p(:))
        u(:) = p(:) / gm

! calculate the acceleration at the locations x1 and x2
!
        call acceleration(t, x(:), u(:), a(:), v(:), b(:))

! separate particle velocity into parallel and perpendicular components
!
        call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the particle gyroperiod and gyroradius
!
        call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate particle energy
!
#ifdef RELAT
        en = gm * mrest
        ek = en - mrest
#else /* RELAT */
        en = 0.5d0 * ua * ua
        ek = en
#endif /* RELAT */

! update the integration time
!
        s = n * ds

! write the progress
!
        write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6),a1,$)") n, t, dt, ua / c    &
            , ek, char(13)

! write results to the output file
!
        write (10, "(20(1pe22.14),i10)") t                                     &
                                   , x(1), x(2), x(3), u(1), u(2), u(3)        &
                                   , ua / c, up / c, ur / c - vper, gm, en, ek        &
                                   , bavg * ba, om, tg * fc, rg * ln, tg, rg   &
                                   , tol, i

        n = n + 1
        m = 0

      end if

! increase data write counter
!
      m = m + 1

! end of iteration
!
    end do

! close the output file
!
    close(10)

! write the progress
!
    write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6))") n, t, dt, ua / c, ek

! write info about the estimator
!
    write(str,"(i12)") mi
    write(*,"('INFO      : maximum iterations per step = ',a)"      )          &
          trim(adjustl(str))
    write(*,"('INFO      : average iterations per step = ',1pe12.6)")          &
          real(ti, kind=8) / ((n - 1) * ndumps)

! open the info file
!
    open  (11, file = 'info.txt', form = 'formatted', position = 'append')

! write info about the estimator
!
    write(11,"('INFO      : maximum iterations per step = ',a)"      )         &
          trim(adjustl(str))
    write(11,"('INFO      : average iterations per step = ',1pe12.6)")         &
          real(ti, kind=8) / ((n - 1) * ndumps)

! close the info file
!
    close(11)
!
!-------------------------------------------------------------------------------
!
  end subroutine integrate_trajectory_si6
!
!===============================================================================
!
! integrate_trajectory_si6_log: subroutine integrates particle trajectory using
!             the 6th order simplectic method with snapshots equally
!             distributed in the logarithm of time
!
! references: "Numerical Hamiltonian Problems", J. M. Sanz-Serna & M. P. Calvo
!             Chapman & Hall, London, New York, 1994
!             "High order starting iterates for implicit Runge-Kutta methods:
!              an improvement for variable-step symplectic integrators", 2002,
!              IMA J. of Num. Ana., 22, 153
!
!===============================================================================
!
  subroutine integrate_trajectory_si6_log()

    use params, only : dtini, tmin, tmax, c

    implicit none

! local variables
!
    integer                         :: n, m, i, mi, ti
    real(kind=PREC), dimension(3,6) :: z, zp
    real(kind=PREC), dimension(3)   :: x , u , p , a
    real(kind=PREC), dimension(3)   :: xn, pn, xt, ut, pt
    real(kind=PREC), dimension(3)   :: v, b
    real(kind=PREC)                 :: gm, t, dt, ds
    real(kind=PREC)                 :: en, ek, ua, ba, up, ur, om, tg, rg
    real(kind=PREC)                 :: tol, tp, wl, wr

! local flags
!
    logical                         :: flag   = .true.
    logical                         :: keepon = .true.

! local parameters
!
    real(kind=PREC), parameter :: b1   =   5.0d0 / 3.0d0                       &
                                , b2   = - 4.0d0 / 3.0d0                       &
                                , b3   =   5.0d0 / 3.0d0
    real(kind=PREC), parameter :: c1   =   0.5d0 - 0.1d0 * dsqrt(15.0d0)       &
                                , c2   =   0.5d0                               &
                                , c3   =   0.5d0 + 0.1d0 * dsqrt(15.0d0)
!
!-------------------------------------------------------------------------------
!
! initialize the iteration number, snapshot number, time, and time steps
!
    n  = 1
    m  = 1
    t  = 0.0d0
    dt = dtini
    ds = qom * dt

! substitute the initial position, velocity, and momentum
!
    x(:) = x0(:)
    u(:) = u0(:)
    p(:) = p0(:)

! calculate the acceleration at the starting point
!
    call acceleration(t, x(:), u(:), a(:), v(:), b(:))

! separate particle velocity into parallel and perpendicular components
!
    call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the Lorentz factor
!
    gm = lorentz_factor(p(:))

! calculate the particle energy
!
#ifdef RELAT
    en = gm * mrest
    ek = en - mrest
#else /* RELAT */
    en = 0.5 * ua * ua
    ek = en
#endif /* RELAT */

! print the progress information
!
    write (*,"('PROGRESS  : ',a8,2x,4(a14))") 'ITER', 'TIME', 'TIMESTEP'       &
            , 'SPEED (c)', 'ENERGY (MeV)'
    write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6),a1,$)") n, t, dt, ua / c, ek    &
            , char(13)

!== INTEGRATION LOOP ==
!
! iterate until the maximum time is reached
!
    do while (keepon)

! find the initial guess for the vector Z
!
      z(1,1:3) = c1 * u(1:3)
      z(2,1:3) = c2 * u(1:3)
      z(3,1:3) = c3 * u(1:3)
      z(1,4:6) = c1 * a(1:3)
      z(2,4:6) = c2 * a(1:3)
      z(3,4:6) = c3 * a(1:3)

! estimate the vector Z (eq. 5.3)
!
!   Z1 = [ a11 * F(y + Z1) + a12 * F(y + Z2) + a13 * F(y + Z3) ]
!   Z2 = [ a21 * F(y + Z1) + a22 * F(y + Z2) + a23 * F(y + Z3) ]
!   Z3 = [ a31 * F(y + Z1) + a32 * F(y + Z2) + a33 * F(y + Z3) ]
!
      call estimate_si6(x(:), p(:), z(:,:), t, dt, ds, tol, i)

! update the solution
!
!   y(n+1) = y(n) + [ d1 * Z1 + d2 * Z2 + d3 * Z3 ]
!
      xn(1:3) = x(1:3) + dt * (b1 * z(1,1:3) + b2 * z(2,1:3) + b3 * z(3,1:3))
      pn(1:3) = p(1:3) + ds * (b1 * z(1,4:6) + b2 * z(2,4:6) + b3 * z(3,4:6))

! update the integration time
!
      t = t + dt

! check if time exceeded the maximum time
!
      if (t >= tmax) keepon = .false.

! store the intermediate snapshots
!
      do while (tt(n) .le. t .and. t .ge. tmin .and. t .lt. tmax)

! calculate the left and right weights
!
        wl = (t - tt(n)) / dt
        wr = 1.0d0 - wl

! interpolate the particle state at the proper time
!
        tp    = wl * (t - dt) + wr * t
        xt(:) = wl * x(:) + wr * xn(:)
        pt(:) = wl * p(:) + wr * pn(:)

! calculate the Lorentz factor and particle velocity
!
        gm    = lorentz_factor(pt(:))
        ut(:) = pt(:) / gm

! calculate acceleration for the location x4 and velocity v4
!
        call acceleration(tp, xt(:), ut(:), a(:), v(:), b(:))

! separate particle velocity into parallel and perpendicular components
!
        call separate_velocity(ut(:), b(:), ba, ua, up, ur)

! calculate the particle gyroperiod and gyroradius
!
        call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate particle energy
!
#ifdef RELAT
        en = gm * mrest
        ek = en - mrest
#else
        en = 0.5 * ua * ua
        ek = en
#endif

! write the progress
!
        write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6),a1,$)") n, tp, dt, ua / c   &
            , ek, char(13)

! write results to the output file
!
        open  (10, file = 'output.dat', form = 'formatted', position = 'append')
        write (10, "(20(1pe22.14))") tp, x(1), x(2), x(3), u(1), u(2), u(3)    &
                                   , ua / c, up / c, ur / c, gm, en, ek        &
                                   , bavg * ba, om, tg * fc, rg * ln, tg, rg   &
                                   , tol
        close (10)

! increate the snapshot index
!
        n = n + 1

      end do

! substitute the new particle position and momentum
!
      x(:) = xn(:)
      p(:) = pn(:)

#ifndef PERIODIC
! if the boundaries are not periodic and particle is out of the box, stop
! the integration
      if (x(1) < bnds(1,1)) keepon = .false.
      if (x(1) > bnds(1,2)) keepon = .false.
      if (x(2) < bnds(2,1)) keepon = .false.
      if (x(2) > bnds(2,2)) keepon = .false.
      if (x(3) < bnds(3,1)) keepon = .false.
      if (x(3) > bnds(3,2)) keepon = .false.
#endif /* PERIODIC */

! update timestep
!
      dt = min(dt, max(1.0d-16, tmax - t))
      ds = qom * dt

! end of iteration
!
    end do

! separate particle velocity into parallel and perpendicular components
!
    call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the particle gyroperiod and gyroradius
!
    call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate particle energy
!
#ifdef RELAT
    en = gm * mrest
    ek = en - mrest
#else
    en = 0.5d0 * ua * ua
    ek = en
#endif

! print the progress info
!
    write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6))") n, t, dtini, ua / c, ek

! store the last snapshot
!
    open  (10, file = 'output.dat', form = 'formatted', position = 'append')
    write (10, "(20(1pe22.14))") t, x(1), x(2), x(3), v(1), v(2), v(3)         &
                               , ua / c, up / c, ur / c, gm, en, ek            &
                               , bavg * ba, om, tg * fc, rg * ln, tg, rg, tol
    close (10)
!
!-------------------------------------------------------------------------------
!
  end subroutine integrate_trajectory_si6_log
!
!===============================================================================
!
! estimate_si6: subroutine estimates the solution for the equation of motion
!               using a simple functional iteration (SI6 version)
!
! references: "Numerical Hamiltonian Problems", J. M. Sanz-Serna & M. P. Calvo
!             Chapman & Hall, London, New York, 1994
!
! description: This subroutines find the solution of the equation 5.3 for
!              the increment Z using the functional iteration
!
!===============================================================================
!
  subroutine estimate_si6(x, p, z, t, dt, dq, tol, it)

    use params, only : maxit, maxeps

    implicit none

! subroutine arguments
!
    real(kind=PREC), dimension(3)  , intent(in)    :: x, p
    real(kind=PREC), dimension(3,6), intent(inout) :: z
    real(kind=PREC)                , intent(in)    :: t
    real(kind=PREC)                , intent(inout) :: dt, dq, tol
    integer                        , intent(inout) :: it

! local variables
!
    real(kind=PREC), dimension(3,6) :: zn
    real(kind=PREC), dimension(6)   :: dh
    real(kind=PREC), dimension(3)   :: x1, p1, u1, a1
    real(kind=PREC), dimension(3)   :: x2, p2, u2, a2
    real(kind=PREC), dimension(3)   :: x3, p3, u3, a3
    real(kind=PREC), dimension(3)   :: v, b
    real(kind=PREC)                 :: g1, g2, g3, eps

! local parameter
!
    real(kind=PREC), parameter :: a11 = 5.0d0 / 36.0d0                         &
                                , a12 = 2.0d0 /  9.0d0 - dsqrt(15.0d0) / 15.0d0&
                                , a13 = 5.0d0 / 36.0d0 - dsqrt(15.0d0) / 30.0d0&
                                , a21 = 5.0d0 / 36.0d0 + dsqrt(15.0d0) / 24.0d0&
                                , a22 = 2.0d0 /  9.0d0                         &
                                , a23 = 5.0d0 / 36.0d0 - dsqrt(15.0d0) / 24.0d0&
                                , a31 = 5.0d0 / 36.0d0 + dsqrt(15.0d0) / 30.0d0&
                                , a32 = 2.0d0 /  9.0d0 + dsqrt(15.0d0) / 15.0d0&
                                , a33 = 5.0d0 / 36.0d0
!     real(kind=PREC), parameter :: e1  =   1.0d1 / 3.0d0                        &
!                                 , e2  = - 2.0d1 / 3.0d0                        &
!                                 , e3  =   1.0d1 / 3.0d0
!
!-------------------------------------------------------------------------------
!
! initiate the iteration control parameters
!
    it  = 0
    eps = 1.0d+16

! perform the simple functional iteration until the conditions are met
!
    do while (eps .gt. maxeps .and. it .lt. maxit)

! prepare the particle position and momentum for the current iteration
!
      x1(:) = x(:) + z(1,1:3)
      x2(:) = x(:) + z(2,1:3)
      x3(:) = x(:) + z(3,1:3)
      p1(:) = p(:) + z(1,4:6)
      p2(:) = p(:) + z(2,4:6)
      p3(:) = p(:) + z(3,4:6)

! calculate the Lorentz factors and particle velocity
!
      g1    = lorentz_factor(p1(:))
      g2    = lorentz_factor(p2(:))
      g3    = lorentz_factor(p3(:))
      u1(:) = p1(:) / g1
      u2(:) = p2(:) / g2
      u3(:) = p3(:) / g3

! calculate the accelerations
!
      call acceleration(t, x1(1:3), u1(1:3), a1(1:3), v(1:3), b(1:3))
      call acceleration(t, x2(1:3), u2(1:3), a2(1:3), v(1:3), b(1:3))
      call acceleration(t, x3(1:3), u3(1:3), a3(1:3), v(1:3), b(1:3))

! update the increment
!
      zn(1,1:3) = dt * (a11 * u1(1:3) + a12 * u2(1:3) + a13 * u3(1:3))
      zn(2,1:3) = dt * (a21 * u1(1:3) + a22 * u2(1:3) + a23 * u3(1:3))
      zn(3,1:3) = dt * (a31 * u1(1:3) + a32 * u2(1:3) + a33 * u3(1:3))
      zn(1,4:6) = dq * (a11 * a1(1:3) + a12 * a2(1:3) + a13 * a3(1:3))
      zn(2,4:6) = dq * (a21 * a1(1:3) + a22 * a2(1:3) + a23 * a3(1:3))
      zn(3,4:6) = dq * (a31 * a1(1:3) + a32 * a2(1:3) + a33 * a3(1:3))

! calculate the maximum of residuum of the increment
!
      eps = maxval(abs(zn - z))

! substitute the new solution of the increment
!
      z = zn

! increase the iteration counter
!
      it = it + 1

    end do
!
! ! estimate the integration error
! !
!     dh(1:3) = dt * (e1 * u1(:) + e2 * u2(:) + e3 * u3(:))
!     dh(4:6) = dq * (e1 * a1(:) + e2 * a2(:) + e3 * a3(:))
!     tol     = sqrt(sum(dh(:) * dh(:)))

! if the convergence was not reached write the warning about it
!
    if (it .ge. maxit) then
      open (11, file = 'info.txt', form = 'formatted', position = 'append')
      write(11,"('WARNING   : convergence not reached at t =',1pe12.5," //     &
               "' eps =',1pe12.5,' tol =',1pe12.5)") t, eps, tol
      close(11)
    end if
!
!-------------------------------------------------------------------------------
!
  end subroutine estimate_si6
!
!===============================================================================
!
! integrate_trajectory_si6: subroutine integrates particle trajectory using
!                           the 6th order simplectic method
!
! references: "Numerical Hamiltonian Problems", J. M. Sanz-Serna & M. P. Calvo
!             Chapman & Hall, London, New York, 1994
!             "High order starting iterates for implicit Runge-Kutta methods:
!              an improvement for variable-step symplectic integrators", 2002,
!              IMA J. of Num. Ana., 22, 153
!
!===============================================================================
!
  subroutine integrate_trajectory_si8()

    use params, only : dtini, tmax, c, ndumps, vper

    implicit none

! local variables
!
    logical                         :: keepon = .true.
    character(len=32)               :: str
    integer                         :: n, m, i, mi, ti
    real(kind=PREC), dimension(4,6) :: z
    real(kind=PREC), dimension(3)   :: x , u , p , a
    real(kind=PREC), dimension(3)   :: v, b
    real(kind=PREC)                 :: gm, t, dt, s, ds, dq
    real(kind=PREC)                 :: en, ek, ua, ba, up, ur, om, tg, rg
    real(kind=PREC)                 :: tol

! local parameters, Butcher's coefficients c_i and Sans-Serna & Calvo's
! coefficients d_i
!
    real(kind=8), parameter :: c1 =  6.9431844202973712388026755553596d-02     &
                             , c2 =  3.3000947820757186759866712044838d-01     &
                             , c3 =  6.6999052179242813240133287955162d-01     &
                             , c4 =  9.3056815579702628761197324444641d-01
    real(kind=8), parameter :: d1 = -1.6407053217392567182070402516331d+00     &
                             , d2 =  1.2143939697985776653621798588684d+00     &
                             , d3 = -1.2143939697985776653621798588684d+00     &
                             , d4 =  1.6407053217392567182070402516331d+00
!
!-------------------------------------------------------------------------------
!
! initialize the iteration number, snapshot number, time, and time steps
!
    n  = 1
    m  = 1
    mi = 0
    ti = 0
    t  = 0.0d0
    s  = 0.0d0
    dt = dtini
    dq = qom * dt
    ds = dt * ndumps

! substitute the initial position, velocity, and momentum
!
    x(:) = x0(:)
    u(:) = u0(:)
    p(:) = p0(:)

! calculate the acceleration at the starting point
!
    call acceleration(t, x(:), u(:), a(:), v(:), b(:))

! separate particle velocity into parallel and perpendicular components
!
    call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the Lorentz factor
!
    gm = lorentz_factor(p(:))

! calculate the particle energy
!
#ifdef RELAT
    en = gm * mrest
    ek = en - mrest
#else /* RELAT */
    en = 0.5d0 * ua * ua
    ek = en
#endif /* RELAT */

! print the progress information
!
    write (*,"('PROGRESS  : ',a8,2x,4(a14))") 'ITER', 'TIME', 'TIMESTEP'       &
            , 'SPEED (c)', 'ENERGY (MeV)'
    write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6),a1,$)") n, t, dt, ua / c, ek    &
            , char(13)

! open the output file
!
    open  (10, file = 'output.dat', form = 'formatted', position = 'append')

!== INTEGRATION LOOP ==
!
! iterate until the maximum time is reached
!
    do while (keepon)

! obtain velocity and acceleration for the initial guess of Z
!
      gm   = lorentz_factor(p(:))
      u(:) = p(:) / gm
      call acceleration(t, x(:), u(:), a(:), v(:), b(:))

! find the initial guess for the vector Z (linear estimation)
!
      z(1,1:3) = c1 * u(1:3)
      z(2,1:3) = c2 * u(1:3)
      z(3,1:3) = c3 * u(1:3)
      z(4,1:3) = c4 * u(1:3)
      z(1,4:6) = c1 * a(1:3)
      z(2,4:6) = c2 * a(1:3)
      z(3,4:6) = c3 * a(1:3)
      z(4,4:6) = c4 * a(1:3)

! estimate the vector Z (eq. 5.3)
!
!   Z1 = [ a11 * F(y + Z1) + a12 * F(y + Z2) + a13 * F(y + Z3) ]
!   Z2 = [ a21 * F(y + Z1) + a22 * F(y + Z2) + a23 * F(y + Z3) ]
!   Z3 = [ a31 * F(y + Z1) + a32 * F(y + Z2) + a33 * F(y + Z3) ]
!
      call estimate_si8(x(:), p(:), z(:,:), t, dt, dq, tol, i)

! update the solution
!
!   y(n+1) = y(n) + [ d1 * Z1 + d2 * Z2 + d3 * Z3 ]
!
      x(1:3) = x(1:3) + dt * (d1 * z(1,1:3) + d2 * z(2,1:3) + d3 * z(3,1:3)    &
                            + d4 * z(4,1:3))
      p(1:3) = p(1:3) + dq * (d1 * z(1,4:6) + d2 * z(2,4:6) + d3 * z(3,4:6)    &
                            + d4 * z(4,4:6))

#ifndef PERIODIC
! if the boundaries are not periodic and particle is out of the box, stop
! the integration
      if (x(1) < bnds(1,1)) keepon = .false.
      if (x(1) > bnds(1,2)) keepon = .false.
      if (x(2) < bnds(2,1)) keepon = .false.
      if (x(2) > bnds(2,2)) keepon = .false.
      if (x(3) < bnds(3,1)) keepon = .false.
      if (x(3) > bnds(3,2)) keepon = .false.
#endif /* PERIODIC */

! update the integration time
!
      t = s + m * dt

! check if time exceeded the maximum time
!
      if (t >= tmax) keepon = .false.

! find the maximum number of iteration in the estimator and update the counter
! of the total number of iterations
!
      mi = max(mi, i)
      ti = ti + i

! store the particle parameters at a given snapshot time
!
      if (m .eq. ndumps) then

! calculate the Lorentz factor and particle velocity
!
        gm   = lorentz_factor(p(:))
        u(:) = p(:) / gm

! calculate the acceleration at the locations x1 and x2
!
        call acceleration(t, x(:), u(:), a(:), v(:), b(:))

! separate particle velocity into parallel and perpendicular components
!
        call separate_velocity(u(:), b(:), ba, ua, up, ur)

! calculate the particle gyroperiod and gyroradius
!
        call gyro_parameters(gm, ba, ur, om, tg, rg)

! calculate particle energy
!
#ifdef RELAT
        en = gm * mrest
        ek = en - mrest
#else /* RELAT */
        en = 0.5d0 * ua * ua
        ek = en
#endif /* RELAT */

! update the integration time
!
        s = n * ds

! write the progress
!
        write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6),a1,$)") n, t, dt, ua / c    &
            , ek, char(13)

! write results to the output file
!
        write (10, "(20(1pe22.14),i10)") t                                     &
                                   , x(1), x(2), x(3), u(1), u(2), u(3)        &
                                   , ua / c, up / c, ur / c - vper, gm, en, ek        &
                                   , bavg * ba, om, tg * fc, rg * ln, tg, rg   &
                                   , tol, i

        n = n + 1
        m = 0

      end if

! increase data write counter
!
      m = m + 1

! end of iteration
!
    end do

! close the output file
!
    close(10)

! write the progress
!
    write (*,"('PROGRESS  : ',i8,2x,4(1pe14.6))") n, t, dt, ua / c, ek

! write info about the estimator
!
    write(str,"(i12)") mi
    write(*,"('INFO      : maximum iterations per step = ',a)"      )          &
          trim(adjustl(str))
    write(*,"('INFO      : average iterations per step = ',1pe12.6)")          &
          real(ti, kind=8) / ((n - 1) * ndumps)

! open the info file
!
    open  (11, file = 'info.txt', form = 'formatted', position = 'append')

! write info about the estimator
!
    write(11,"('INFO      : maximum iterations per step = ',a)"      )         &
          trim(adjustl(str))
    write(11,"('INFO      : average iterations per step = ',1pe12.6)")         &
          real(ti, kind=8) / ((n - 1) * ndumps)

! close the info file
!
    close(11)
!
!-------------------------------------------------------------------------------
!
  end subroutine integrate_trajectory_si8
!
!===============================================================================
!
! estimate_si6: subroutine estimates the solution for the equation of motion
!               using a simple functional iteration (SI6 version)
!
! references: "Numerical Hamiltonian Problems", J. M. Sanz-Serna & M. P. Calvo
!             Chapman & Hall, London, New York, 1994
!
! description: This subroutines find the solution of the equation 5.3 for
!              the increment Z using the functional iteration
!
!===============================================================================
!
  subroutine estimate_si8(x, p, z, t, dt, dq, tol, it)

    use params, only : maxit, maxeps

    implicit none

! subroutine arguments
!
    real(kind=PREC), dimension(3)  , intent(in)    :: x, p
    real(kind=PREC), dimension(4,6), intent(inout) :: z
    real(kind=PREC)                , intent(in)    :: t
    real(kind=PREC)                , intent(inout) :: dt, dq, tol
    integer                        , intent(inout) :: it

! local variables
!
    real(kind=PREC), dimension(4,6) :: zn
    real(kind=PREC), dimension(6)   :: dh
    real(kind=PREC), dimension(3)   :: x1, p1, u1, a1
    real(kind=PREC), dimension(3)   :: x2, p2, u2, a2
    real(kind=PREC), dimension(3)   :: x3, p3, u3, a3
    real(kind=PREC), dimension(3)   :: x4, p4, u4, a4
    real(kind=PREC), dimension(3)   :: v, b
    real(kind=PREC)                 :: g1, g2, g3, g4, eps

! local parameters, Butcher's coefficients a_ij
!
    real(kind=8), parameter :: a11 =  8.6963711284363464343265987305500d-02    &
                             , a12 = -2.6604180084998793313385130476953d-02    &
                             , a13 =  1.2627462689404724515056880574618d-02    &
                             , a14 = -3.5551496857956831569109818495695d-03    &
                             , a21 =  1.8811811749986807165068554508717d-01    &
                             , a22 =  1.6303628871563653565673401269450d-01    &
                             , a23 = -2.7880428602470895224151106418997d-02    &
                             , a24 =  6.7355005945381555153986690857040d-03    &
                             , a31 =  1.6719192197418877317113330552530d-01    &
                             , a32 =  3.5395300603374396653761913180800d-01    &
                             , a33 =  1.6303628871563653565673401269450d-01    &
                             , a34 = -1.4190694931141142964153570476171d-02    &
                             , a41 =  1.7748257225452261184344295646057d-01    &
                             , a42 =  3.1344511474186834679841114481438d-01    &
                             , a43 =  3.5267675751627186462685315586596d-01    &
                             , a44 =  8.6963711284363464343265987305500d-02
    real(kind=PREC), parameter :: e1  =   1.0d1 / 3.0d0                        &
                                , e2  = - 2.0d1 / 3.0d0                        &
                                , e3  =   1.0d1 / 3.0d0
!
!-------------------------------------------------------------------------------
!
! initiate the iteration control parameters
!
    it  = 0
    eps = 1.0d+16

! perform the simple functional iteration until the conditions are met
!
    do while (eps > maxeps .and. it < maxit)

! prepare the particle position and momentum for the current iteration
!
      x1(:) = x(:) + dt * z(1,1:3)
      x2(:) = x(:) + dt * z(2,1:3)
      x3(:) = x(:) + dt * z(3,1:3)
      x4(:) = x(:) + dt * z(4,1:3)
      p1(:) = p(:) + dq * z(1,4:6)
      p2(:) = p(:) + dq * z(2,4:6)
      p3(:) = p(:) + dq * z(3,4:6)
      p4(:) = p(:) + dq * z(4,4:6)

! calculate the Lorentz factors and particle velocity
!
      g1    = lorentz_factor(p1(:))
      g2    = lorentz_factor(p2(:))
      g3    = lorentz_factor(p3(:))
      g4    = lorentz_factor(p4(:))
      u1(:) = p1(:) / g1
      u2(:) = p2(:) / g2
      u3(:) = p3(:) / g3
      u4(:) = p4(:) / g4

! calculate the accelerations
!
      call acceleration(t, x1(1:3), u1(1:3), a1(1:3), v(1:3), b(1:3))
      call acceleration(t, x2(1:3), u2(1:3), a2(1:3), v(1:3), b(1:3))
      call acceleration(t, x3(1:3), u3(1:3), a3(1:3), v(1:3), b(1:3))
      call acceleration(t, x4(1:3), u4(1:3), a4(1:3), v(1:3), b(1:3))

! update the increment
!
      zn(1,1:3) = a11 * u1(1:3) + a12 * u2(1:3) + a13 * u3(1:3) + a14 * u4(1:3)
      zn(1,4:6) = a11 * a1(1:3) + a12 * a2(1:3) + a13 * a3(1:3) + a14 * a4(1:3)
      zn(2,1:3) = a21 * u1(1:3) + a22 * u2(1:3) + a23 * u3(1:3) + a24 * u4(1:3)
      zn(2,4:6) = a21 * a1(1:3) + a22 * a2(1:3) + a23 * a3(1:3) + a24 * a4(1:3)
      zn(3,1:3) = a31 * u1(1:3) + a32 * u2(1:3) + a33 * u3(1:3) + a34 * u4(1:3)
      zn(3,4:6) = a31 * a1(1:3) + a32 * a2(1:3) + a33 * a3(1:3) + a34 * a4(1:3)
      zn(4,1:3) = a41 * u1(1:3) + a42 * u2(1:3) + a43 * u3(1:3) + a44 * u4(1:3)
      zn(4,4:6) = a41 * a1(1:3) + a42 * a2(1:3) + a43 * a3(1:3) + a44 * a4(1:3)

! calculate the maximum of residuum of the increment
!
      eps = maxval(abs(zn - z))

! substitute the new solution of the increment
!
      z = zn

! increase the iteration counter
!
      it = it + 1

    end do

! estimate the integration error
!
    dh(1:3) = dt * (e1 * u1(:) + e2 * u2(:) + e3 * u3(:))
    dh(4:6) = dq * (e1 * a1(:) + e2 * a2(:) + e3 * a3(:))
    tol     = sqrt(sum(dh(:) * dh(:)))

! if the convergence was not reached write the warning about it
!
    if (it .ge. maxit) then
      open (11, file = 'info.txt', form = 'formatted', position = 'append')
      write(11,"('WARNING   : convergence not reached at t =',1pe12.5," //     &
               "' eps =',1pe12.5,' tol =',1pe12.5)") t, eps, tol
      close(11)
    end if
!
!-------------------------------------------------------------------------------
!
  end subroutine estimate_si8
!
!===============================================================================
!
! pos2index: subroutine converts a given position to the array index
!
!===============================================================================
!
  subroutine pos2index(x, r)

    use fields, only : ng

    implicit none

! input and output arguments
!
    real(kind=PREC), dimension(3), intent(in)  :: x
    real(kind=PREC), dimension(3), intent(out) :: r

! local variables
!
    integer         :: i
    real(kind=PREC) :: t
!
!------------------------------------------------------------------------------
!
    do i = 1, DIMS
      t    = (x(i) - bnds(i,1)) / bsiz(i) + 0.5 / dm(i)
      t    = t - floor(t)
      r(i) = dm(i) * t + ng
    end do
!
!------------------------------------------------------------------------------
!
  end subroutine pos2index
!
!===============================================================================
!
! prepare_interpolation: subroutine prepares coeafficients for interpolation
!
!===============================================================================
!
  subroutine prepare_interpolation(x, ii, jj, kk, dr, cx, cy, cz)

    implicit none

! input and output arguments
!
    real(kind=PREC), dimension(3), intent(in)  :: x
    integer        , dimension(4), intent(out) :: ii, jj, kk
    real(kind=8)   , dimension(3), intent(out) :: dr
    real(kind=8)   , dimension(4), intent(out) :: cx, cy, cz
!
!------------------------------------------------------------------------------
!
#if !defined TRILIN && !defined TRICUB
!
!= nearest interpolation =
!
! calculate indices
!
    ii(1) = nint(x(1))
    jj(1) = nint(x(2))
#if DIMS == 3
    kk(1) = nint(x(3))
#endif /* DIMS == 3 */
#endif /* !TRILIN & !TRICUB */
#ifdef TRILIN
!
!= trilinear interpolation =
!
! calculate indices
!
    ii(1) = floor(x(1))
    ii(2) = ii(1) + 1
    jj(1) = floor(x(2))
    jj(2) = jj(1) + 1
#if DIMS == 3
    kk(1) = floor(x(3))
    kk(2) = kk(1) + 1
#endif /* DIMS == 3 */

! calculate intercell position
!
    dr(1) = x(1) - ii(1)
    dr(2) = x(2) - jj(1)
#if DIMS == 3
    dr(3) = x(3) - kk(1)
#endif /* DIMS == 3 */
#endif /* TRILIN */
#ifdef TRICUB
!
!= tricubic interpolation =
!
! calculate indices
!
    ii(2) = floor(x(1))
    ii(1) = ii(2) - 1
    ii(3) = ii(2) + 1
    ii(4) = ii(2) + 2
    jj(2) = floor(x(2))
    jj(1) = jj(2) - 1
    jj(3) = jj(2) + 1
    jj(4) = jj(2) + 2
#if DIMS == 3
    kk(2) = floor(x(3))
    kk(1) = kk(2) - 1
    kk(3) = kk(2) + 1
    kk(4) = kk(2) + 2
#endif /* DIMS == 3 */

! calculate intercell position
!
    dr(1) = x(1) - ii(2)
    dr(2) = x(2) - jj(2)
#if DIMS == 3
    dr(3) = x(3) - kk(2)
#endif /* DIMS == 3 */

! coefficients for dx, dy, and dz
!
    call coefficients_cubic(dr(1), cx(:))
    call coefficients_cubic(dr(2), cy(:))
#if DIMS == 3
    call coefficients_cubic(dr(3), cz(:))
#endif /* DIMS == 3 */
#endif /* TRICUB */
!
!------------------------------------------------------------------------------
!
  end subroutine prepare_interpolation
!
!===============================================================================
!
! interpolate: subroutine interpolates field value for a given position
!
!===============================================================================
!
  real(kind=8) function interpolate(f, ii, jj, kk, dr, cx, cy, cz) result(q)

    implicit none

! input and output arguments
!
    real(kind=4), dimension(:,:,:), intent(in) :: f
    integer     , dimension(4)    , intent(in) :: ii, jj, kk
    real(kind=8), dimension(3)    , intent(in) :: dr
    real(kind=8), dimension(4)    , intent(in) :: cx, cy, cz

! local variables
!
#ifdef TRILIN
    real(kind=4) :: q11, q12, q21, q22, q1, q2
#endif /* TRILIN */
#ifdef TRICUB
    real(kind=4) :: q11, q12, q13, q14, q21, q22, q23, q24                     &
                  , q31, q32, q33, q34, q41, q42, q43, q44, q1, q2, q3, q4
#endif /* TRICUB */
!
!------------------------------------------------------------------------------
!
#if !defined TRILIN && !defined TRICUB
!
!= nearest interpolation =
!
#if DIMS == 2
    q = f(ii(1),jj(1),1)
#else /* DIMS == 2 */
    q = f(ii(1),jj(1),kk(1))
#endif /* DIMS == 2 */
#endif /* !TRILIN & !TRICUB */
#ifdef TRILIN
#if DIMS == 2
!= bilinear interpolation =
!
! interpolate along the Y direction
!
    q1 = plinear(dr(2), f(ii(1),jj(1),1), f(ii(1),jj(2),1))
    q2 = plinear(dr(2), f(ii(2),jj(1),1), f(ii(2),jj(2),1))
#else /* DIMS == 2 */
!= trilinear interpolation =
!
! interpolate along the Z direction
!
    q11 = plinear(dr(3), f(ii(1),jj(1),kk(1)), f(ii(1),jj(1),kk(2)))
    q12 = plinear(dr(3), f(ii(1),jj(2),kk(1)), f(ii(1),jj(2),kk(2)))
    q21 = plinear(dr(3), f(ii(2),jj(1),kk(1)), f(ii(2),jj(1),kk(2)))
    q22 = plinear(dr(3), f(ii(2),jj(2),kk(1)), f(ii(2),jj(2),kk(2)))

! interpolate along the Y direction
!
    q1 = plinear(dr(2), q11, q12)
    q2 = plinear(dr(2), q21, q22)
#endif /* DIMS == 2 */

! interpolate the value at a given position
!
    q  = plinear(dr(1), q1 , q2 )
#endif /* TRILIN */
#ifdef TRICUB
#if DIMS == 2
!= bicubic interpolation =
!
! interpolate along the Y direction
!
    q1  = pcubic(cy, f(ii(1),jj(1),1), f(ii(1),jj(2),1)                        &
                   , f(ii(1),jj(3),1), f(ii(1),jj(4),1))
    q2  = pcubic(cy, f(ii(2),jj(1),1), f(ii(2),jj(2),1)                        &
                   , f(ii(2),jj(3),1), f(ii(2),jj(4),1))
    q3  = pcubic(cy, f(ii(3),jj(1),1), f(ii(3),jj(2),1)                        &
                   , f(ii(3),jj(3),1), f(ii(3),jj(4),1))
    q4  = pcubic(cy, f(ii(4),jj(1),1), f(ii(4),jj(2),1)                        &
                   , f(ii(4),jj(3),1), f(ii(4),jj(4),1))
#else /* DIMS == 2 */
!= tricubic interpolation =
!
! interpolate along Z direction
!
    q11 = pcubic(cz, f(ii(1),jj(1),kk(1)), f(ii(1),jj(1),kk(2))                &
                   , f(ii(1),jj(1),kk(3)), f(ii(1),jj(1),kk(4)))
    q12 = pcubic(cz, f(ii(1),jj(2),kk(1)), f(ii(1),jj(2),kk(2))                &
                   , f(ii(1),jj(2),kk(3)), f(ii(1),jj(2),kk(4)))
    q13 = pcubic(cz, f(ii(1),jj(3),kk(1)), f(ii(1),jj(3),kk(2))                &
                   , f(ii(1),jj(3),kk(3)), f(ii(1),jj(3),kk(4)))
    q14 = pcubic(cz, f(ii(1),jj(4),kk(1)), f(ii(1),jj(4),kk(2))                &
                   , f(ii(1),jj(4),kk(3)), f(ii(1),jj(4),kk(4)))

    q21 = pcubic(cz, f(ii(2),jj(1),kk(1)), f(ii(2),jj(1),kk(2))                &
                   , f(ii(2),jj(1),kk(3)), f(ii(2),jj(1),kk(4)))
    q22 = pcubic(cz, f(ii(2),jj(2),kk(1)), f(ii(2),jj(2),kk(2))                &
                   , f(ii(2),jj(2),kk(3)), f(ii(2),jj(2),kk(4)))
    q23 = pcubic(cz, f(ii(2),jj(3),kk(1)), f(ii(2),jj(3),kk(2))                &
                   , f(ii(2),jj(3),kk(3)), f(ii(2),jj(3),kk(4)))
    q24 = pcubic(cz, f(ii(2),jj(4),kk(1)), f(ii(2),jj(4),kk(2))                &
                   , f(ii(2),jj(4),kk(3)), f(ii(2),jj(4),kk(4)))

    q31 = pcubic(cz, f(ii(3),jj(1),kk(1)), f(ii(3),jj(1),kk(2))                &
                   , f(ii(3),jj(1),kk(3)), f(ii(3),jj(1),kk(4)))
    q32 = pcubic(cz, f(ii(3),jj(2),kk(1)), f(ii(3),jj(2),kk(2))                &
                   , f(ii(3),jj(2),kk(3)), f(ii(3),jj(2),kk(4)))
    q33 = pcubic(cz, f(ii(3),jj(3),kk(1)), f(ii(3),jj(3),kk(2))                &
                   , f(ii(3),jj(3),kk(3)), f(ii(3),jj(3),kk(4)))
    q34 = pcubic(cz, f(ii(3),jj(4),kk(1)), f(ii(3),jj(4),kk(2))                &
                   , f(ii(3),jj(4),kk(3)), f(ii(3),jj(4),kk(4)))

    q41 = pcubic(cz, f(ii(4),jj(1),kk(1)), f(ii(4),jj(1),kk(2))                &
                   , f(ii(4),jj(1),kk(3)), f(ii(4),jj(1),kk(4)))
    q42 = pcubic(cz, f(ii(4),jj(2),kk(1)), f(ii(4),jj(2),kk(2))                &
                   , f(ii(4),jj(2),kk(3)), f(ii(4),jj(2),kk(4)))
    q43 = pcubic(cz, f(ii(4),jj(3),kk(1)), f(ii(4),jj(3),kk(2))                &
                   , f(ii(4),jj(3),kk(3)), f(ii(4),jj(3),kk(4)))
    q44 = pcubic(cz, f(ii(4),jj(4),kk(1)), f(ii(4),jj(4),kk(2))                &
                   , f(ii(4),jj(4),kk(3)), f(ii(4),jj(4),kk(4)))

! interpolate along the Y direction
!
    q1  = pcubic(cy, q11, q12, q13, q14)
    q2  = pcubic(cy, q21, q22, q23, q24)
    q3  = pcubic(cy, q31, q32, q33, q34)
    q4  = pcubic(cy, q41, q42, q43, q44)
#endif /* DIMS == 2 */

! interpolate along the X direction
!
    q   = pcubic(cx, q1 , q2 , q3 , q4 )
#endif /* TRICUB */
!
!------------------------------------------------------------------------------
!
  end function interpolate
!
!===============================================================================
!
! plinear: subroutine performs one dimensional linear interpolation
!
!===============================================================================
!
  real(kind=8) function plinear(x, fl, fr) result(q)

    implicit none

! input and output arguments
!
    real(kind=8), intent(in)  :: x
    real(kind=4), intent(in)  :: fl, fr
!
!------------------------------------------------------------------------------
!
    q = fl + x * (fr - fl)
!
!------------------------------------------------------------------------------
!
  end function plinear
!
!===============================================================================
!
! pcubic: subroutine performs one dimensional cubic interpolation
!
!===============================================================================
!
  real(kind=8) function pcubic(c, fk, fl, fr, fq) result(q)

    implicit none

! input and output arguments
!
    real(kind=8), dimension(4), intent(in)  :: c
    real(kind=4)              , intent(in)  :: fk, fl, fr, fq

#ifdef TVD
! local parameters
!
    real(kind=4) :: dfl, dfr, ds, dl, dr
#endif /* TVD */
!
!------------------------------------------------------------------------------
!
#ifdef TVD
    ds  = (fr - fl)
    dl  = (fl - fk)
    dr  = (fq - fr)
    dfl = sign(1.0, ds) * min(abs(ds), abs(dl))
    dfr = sign(1.0, ds) * min(abs(ds), abs(dr))

    if ((dl * ds) .le. 0.0) then
      dfl = 0.0
    end if

    if ((dr * ds) .le. 0.0) then
      dfr = 0.0
    end if

    q = c(1) * fl + c(2) * fr + c(3) * dfl + c(4) * dfr
#else /* TVD */
    q = c(1) * fk + c(2) * fl + c(3) * fr + c(4) * fq
#endif /* TVD */
!
!------------------------------------------------------------------------------
!
  end function pcubic
!
!===============================================================================
!
! coefficients_cubic: subroutine prepares coeafficients for cubic interpolation
!
!===============================================================================
!
  subroutine coefficients_cubic(x, c)

    implicit none

! input and output arguments
!
    real(kind=8)              , intent(in)  :: x
    real(kind=8), dimension(4), intent(out) :: c

! local variables
!
    real(kind=8) :: x1, x2, x3, xd
!
!------------------------------------------------------------------------------
!
! prepare local variables
!
    x1 = x - 1.0
    x2 = x * x
#ifdef TVD
    x3 = x1 * x1
    xd = 2.0 * x
#else /* TVD */
#endif /* TVD */

! calculate coefficients
!
#ifdef TVD
    c(1) = x3 * (xd + 1.0)
    c(2) = x2 * (3.0 - xd)
    c(3) = x  * x3
    c(4) = x2 * x1
#else /* TVD */
    c(1) = 0.5 * x * ( ( 2.0 - x ) * x - 1.0 )
    c(2) = 0.5 * x2 * ( 3.0 * x - 5.0 ) + 1.0
    c(3) = 0.5 * x * ( ( 4.0 - 3.0 * x ) * x + 1.0 )
    c(4) = 0.5 * x1 * x2
#endif /* TVD */
!
!------------------------------------------------------------------------------
!
  end subroutine coefficients_cubic
!
!===============================================================================
!
! lorentz_factor: subroutine calculates the Lorentz factor
!
!===============================================================================
!
  real(kind=PREC) function lorentz_factor(p) result(gm)

    implicit none

! input and output arguments
!
    real(kind=PREC), dimension(3), intent(in)  :: p
!
!------------------------------------------------------------------------------
!
#ifdef RELAT
    gm = sqrt(1.0d0 + dot_product(p, p) / csq)
#else
    gm = 1.0
#endif
!
!-------------------------------------------------------------------------------
!
  end function lorentz_factor
!
!===============================================================================
!
! acceleration: subroutine calculates acceleration vector at a given location
!
!===============================================================================
!
  subroutine acceleration(t, x, v, a, u, b)

    use fields, only : ux, uy, uz, bx, by, bz
    use params, only : nghost
#ifdef TEST
    use params, only : bini, bshr, bamp, vamp, vrat, freq
#endif /* TEST */

    implicit none

! input and output arguments
!
    real(kind=PREC)              , intent(in)  :: t
    real(kind=PREC), dimension(3), intent(in)  :: x, v
    real(kind=PREC), dimension(3), intent(out) :: a, u, b

#ifdef TEST
! local variables
!
    real(kind=PREC), dimension(3) :: w
#ifdef ITEST
    real(kind=PREC)               :: dl, ra, rb, xt, yt, rt
#endif /* ITEST */
#else /* TEST */
! local variables
!
    real(kind=PREC), dimension(3) :: r, w

! position indices
!
    integer                       :: dist
    integer        , dimension(4) :: ii, jj, kk
    real(kind=8   ), dimension(4) :: cx, cy, cz
    real(kind=8   ), dimension(3) :: dr
#endif /* TEST */
!
!-------------------------------------------------------------------------------
!
#ifndef TEST
! convert position to index
!
      call pos2index(x, r)

#ifdef BNDRY
      dist = min(minval(dm(1:DIMS) - r(1:DIMS)), minval(r(1:DIMS)))
      if (dist .gt. nghost) then
#endif /* BNDRY */

! prepare coefficients for interpolation
!
        call prepare_interpolation(r, ii, jj, kk, dr, cx, cy, cz)
#endif /* !TEST */

! interpolate field components at the particle position
!
#ifdef TEST
#ifdef WTEST
        u(1) = 0.0
        u(2) = - vamp * sin(pi2 * freq * x(1))
        u(3) = 0.0

        b(1) = bpar
        b(2) = bamp * cos(pi2 * freq * x(1))
        b(3) = bamp * sin(pi2 * freq * x(1))
#endif /* WTEST */
#ifdef ITEST
! calculate the local velocity
!
        u(1) =      - vamp * x(1)
        u(2) = vrat * vamp * x(2)
        u(3) = 0.0d0

! calculate the local magnetic field
!
        dl   = bamp / bini
        ra   = 1.0d0 + dl
        rb   = 1.0d0 - dl

        xt   = x(1) / ra
        yt   = x(2) / rb

        rt   = dsqrt(xt * xt + yt * yt)

        if (rt .gt. 0.0d0) then
          b(1) =   yt / rb / rt * bini
          b(2) = - xt / ra / rt * bini
          b(3) = bshr
        else
          b(1) = 0.0d0
          b(2) = 0.0d0
          b(3) = bshr
        end if
#endif /* ITEST */
#else /* TEST */
        u(1) = interpolate(ux, ii, jj, kk, dr, cx, cy, cz)
        u(2) = interpolate(uy, ii, jj, kk, dr, cx, cy, cz)
        u(3) = interpolate(uz, ii, jj, kk, dr, cx, cy, cz)
        b(1) = interpolate(bx, ii, jj, kk, dr, cx, cy, cz)
        b(2) = interpolate(by, ii, jj, kk, dr, cx, cy, cz)
        b(3) = interpolate(bz, ii, jj, kk, dr, cx, cy, cz)
#endif /* TEST */

! subtract the fluid velocity
!
        w(:) = v(:) - u(:)

! compute the acceleration
!
        a(1) = w(2) * b(3) - w(3) * b(2)
        a(2) = w(3) * b(1) - w(1) * b(3)
        a(3) = w(1) * b(2) - w(2) * b(1)
#if !defined TEST && defined BNDRY
      else
        a(:) = 0.0
      endif
#endif /* !TEST & BNDRY */
!
!-------------------------------------------------------------------------------
!
  end subroutine acceleration
!
!===============================================================================
!
! separate_velocity: subroutine separates velocity into two components, parallel
!                    and perpendicular to the local magnetic field
!
!===============================================================================
!
  subroutine separate_velocity(v, b, ba, vu, vp, vr)

    implicit none

! input and output arguments
!
    real(kind=PREC), dimension(3), intent(in)  :: v, b
    real(kind=PREC)              , intent(out) :: ba, vu, vp, vr

! local variables
!
    real(kind=PREC), dimension(3) :: p
    real(kind=PREC)               :: pp, vv
!
!------------------------------------------------------------------------------
!
! calculate amplitude of magnetic field
!
    ba = sqrt(dot_product(b, b))

! calculate unit vector parallel to B
!
    if (ba .gt. 0.0) then
      p  = b / ba
    else
      p  = 0.0
    end if

! calculate component parallel to B
!
    pp = dot_product(v, p)**2

! calculate amplitude of velocity
!
    vv = dot_product(v, v)

! calculate amplitude of the parallel and perpendicular components of velocity
!
    vu = sqrt(vv)
    vp = sqrt(pp)
    vr = sqrt(vv - pp)
!
!-------------------------------------------------------------------------------
!
  end subroutine separate_velocity
!
!===============================================================================
!
! gyro_parameters: subroutine calculates particle gyrofrequency, gyroperiod and
!                  gyroradius
!
!===============================================================================
!
  subroutine gyro_parameters(gm, ba, vr, om, tg, rg)

    implicit none

! input and output arguments
!
    real(kind=PREC), intent(in)  :: gm
    real(kind=PREC), intent(in)  :: ba, vr
    real(kind=PREC), intent(out) :: om, tg, rg
!
!------------------------------------------------------------------------------
!
    om = om0 * ba / gm
    tg = pi2 / om / fc
    rg = vr / om / fc
!
!-------------------------------------------------------------------------------
!
  end subroutine gyro_parameters
!
end module particles

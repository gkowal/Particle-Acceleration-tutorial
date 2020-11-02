!!******************************************************************************
!!
!! module: params - subroutines to read parameters file.
!!
!! Copyright (C) 2007-2010 Grzegorz Kowal <grzegorz@gkowal.info>
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
module params

  implicit none

! input data parameters
!
  character(len = 128), save :: idir    = "./"     ! input directory
  character(len = 128), save :: odir    = "./"     ! output directory
  character(len =   4), save :: fformat = 'fits'   ! file format:
                                                   ! 'fits' - FITS
                                                   ! 'hdf5' - HDF5
  character(len =   1), save :: ftype   = 'r'      ! file type: 'r', 'p', 'f'
  integer             , save :: fnumber = 0        ! file number

! geometry parameters
!
  character(len =   1), save :: tunit   = 's'      ! time units:
                                                   ! 's' - second
                                                   ! 'm' - minute
                                                   ! 'h' - hour
                                                   ! 'd' - day
                                                   ! 'w' - week
                                                   ! 'y' - year
  real(kind=PREC)     , save :: tmulti  = 1.0      ! time unit count
  real(kind=PREC)     , save :: bunit   = 1.0      ! magnetic field unit in Gs
  integer             , save :: nghost  = 8        ! number of ghost pixels near the boundary

! initial particle state parameters
!
  character(len =   1), save :: ptype   = 'p'      ! particle type:
                                                   ! 'p' - proton
                                                   ! 'e' - electron
  real(kind=PREC)     , save :: xc      = 0.0      ! initial position
  real(kind=PREC)     , save :: yc      = 0.0
  real(kind=PREC)     , save :: zc      = 0.0
  real(kind=PREC)     , save :: vpar    = 0.0      ! initial parallel speed [in c]
  real(kind=PREC)     , save :: vper    = 0.1      ! initial perpendicular speed [in c]
  real(kind=PREC)     , save :: rho     = 0.5      ! safety coefficient

! plasma parameters
!
  real(kind=PREC)     , save :: c       = 1.0      ! the speed of light in Va
  real(kind=PREC)     , save :: dens    = 1.0      ! density [1 u/cm^3]

! integration quality parameters
!
  character(len =   4), save :: method  = 'rk4'    ! the integration method: rk4 or si4
  real(kind=PREC)     , save :: maxtol  = 1.0e-4   ! the maximi integration tolerance
  real(kind=PREC)     , save :: maxeps  = 1.0e-15  ! the maximum iteration error
  real(kind=PREC)     , save :: dtini   = 1.0e-8   ! the initial time step
  real(kind=PREC)     , save :: dtmax   = 1.0      ! maximum allowed step size
  integer             , save :: maxit   = 1000     ! the limit of iterations

! output data parameters
!
  character(len =   1), save :: output  = 'i'      ! the type of output:
                                                   ! 'i' - by the iteration number
                                                   ! 'l' - by the logarithmic time
  integer             , save :: ndumps  = 1000     ! number of steps between subsequent dumps if the output is 'i'
                                                   ! or number of dumps per time decade if the output is 'l'
  real(kind=PREC)     , save :: tmin    = 1.0e-3   ! minimum time of writing data
  real(kind=PREC)     , save :: tmax    = 1.0      ! maximum time for integration

#ifdef TEST
! test problem parameters
!
  real(kind=PREC)     , save :: bini    = 1.0      ! the mean magnetic field
  real(kind=PREC)     , save :: bshr    = 0.0      ! the guilde field
  real(kind=PREC)     , save :: bamp    = 0.0      ! the amplitude of the magnetic field fluctuations
  real(kind=PREC)     , save :: vamp    = 0.0      ! the amplitude of the velocity field fluctuations
  real(kind=PREC)     , save :: vrat    = 1.0      ! the ratio between velocity fluctuations amplitudes in different directions
  real(kind=PREC)     , save :: freq    = 1.0      ! the frequency of the field fluctuations
  real(kind=PREC)     , save :: epar    = 0.0      ! the constant electric field along the parallel direction
#endif /* TEST */
!
!-------------------------------------------------------------------------------
!
  contains
!
!===============================================================================
!
! read_config: subroutine to read config file and to fill proper parameters
!
!===============================================================================
!
  subroutine read_params
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
    character(len =   *), parameter :: configfile = './params.in' ! config file with parameters
    character(len=255) :: line, name, value
    integer            :: l, i, ios
    real               :: vv
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
    open(unit=1, file=configfile, err=100)

10  read(unit=1, fmt="(a)", end=20, err=200 ) line

    call parse_line(line, name, value)

    select case(name)
      case ('fformat')
        l = len_trim(value)
        write(fformat, "(a)" ) value(2:l-1)
      case ('ftype')
        l = len_trim(value)
        write(ftype  , "(a)" ) value(2:l-1)
      case ('idir')
        l = len_trim(value)
        write(idir   , "(a)" ) value(2:l-1)
      case ('method')
        l = len_trim(value)
        write(method , "(a)" ) value(2:l-1)
      case ('odir')
        l = len_trim(value)
        write(odir   , "(a)" ) value(2:l-1)
      case ('output')
        l = len_trim(value)
        write(output , "(a)" ) value(2:l-1)
      case ('ptype')
        l = len_trim(value)
        write(ptype  , "(a)" ) value(2:l-1)
      case ('tunit')
        l = len_trim(value)
        write(tunit  , "(a)" ) value(2:l-1)

      case ('fnumber')
        read (value  , "(i6)") fnumber
      case ('maxit')
        read (value  , "(i9)") maxit
      case ('ndumps')
        read (value  , "(i9)") ndumps
      case ('nghost')
        read (value  , "(i6)") nghost

      case ('c')
        read (value  , *     ) c
      case ('dens')
        read (value  , *     ) dens
      case ('dtini')
        read (value  , *     ) dtini
      case ('dtmax')
        read (value  , *     ) dtmax
      case ('maxeps')
        read (value  , *     ) maxeps
      case ('maxtol')
        read (value  , *     ) maxtol
      case ('rho')
        read (value  , *     ) rho
      case ('tmulti')
        read (value  , *     ) tmulti
      case ('vpar')
        read (value  , *     ) vpar
      case ('vper')
        read (value  , *     ) vper
      case ('tmin')
        read (value  , *     ) tmin
      case ('tmax')
        read (value  , *     ) tmax
      case ('xc')
        read (value  , *     ) xc
      case ('yc')
        read (value  , *     ) yc
      case ('zc')
        read (value  , *     ) zc

#ifdef TEST
      case ('bini')
        read (value  , *     ) bini
      case ('bshr')
        read (value  , *     ) bshr
      case ('bamp')
        read (value  , *     ) bamp
      case ('vamp')
        read (value  , *     ) vamp
      case ('vrat')
        read (value  , *     ) vrat
      case ('freq')
        read (value  , *     ) freq
      case ('epar')
        read (value  , *     ) epar
#endif /* TEST */
      case default
    end select

    go to 10

20  close(1)

! check input parameters
!
    vv = sqrt(vpar * vpar + vper * vper)
    if (vv .ge. 1.0) then
      write( *, "('ERROR     : ',a)" ) "absolute speed of the particle is larger than c!"
      write( *, "('ERROR     : ',a,1pe15.8)" ) "|v| = ", vv
      stop
    end if
!
    return
!
100 print *, 'Error opening file ', configfile
    stop
200 print *, 'Error reading file ', configfile
    stop
!
    return
!
  end subroutine read_params
!
!===============================================================================
!
! parse_line: subroutine to parse line into name and value of parameter
!
!===============================================================================
!
  subroutine parse_line(line, name, value)
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
    character(len=*), intent(in)  :: line
    character(len=*), intent(out) :: name, value

    integer :: l, i
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
    l = len_trim(line)

    i = index( line, '=' )

    name  = trim(adjustl(line(1:i-1)))
    value = trim(adjustl(line(i+1:l)))
!
  end subroutine parse_line

end module params

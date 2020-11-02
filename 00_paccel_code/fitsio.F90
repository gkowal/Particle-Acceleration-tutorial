!!******************************************************************************
!!
!! module: fitsio - subroutines to read data from FITS files
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
module fitsio

  implicit none

  integer, dimension(3), save :: dm
  integer              , save :: nl

!-------------------------------------------------------------------------------
!
  contains
!
!===============================================================================
!
! hdf4_init: subroutine reads attributes from a HDF file
!
!===============================================================================
!
  subroutine fits_init(nm)

    use params, only : idir

    implicit none

! input parameters
!
    character(len=4), intent(in) :: nm

! local variables
!
    character(len=255) :: fl
    logical            :: info

! FITS variables
!
    integer            :: status, iunit, bsize, bpix, naxes
!
!-------------------------------------------------------------------------------
!
! set default values
!
    dm (:) = 1

! generate filename
!
    write(fl,"(a,a4,'.fits.gz')") trim(idir), trim(nm)

! check if the file exists
!
    inquire(file = fl, exist = info)

    if (.not. info) then
      write (*,'("ERROR   : file ",a," does not exist!")') trim(fl)
      stop
    endif

! check if the file is written in FITS format
!
    status  = 0

#ifndef TEST
    call ftgiou(iunit, status)
    call ftopen(iunit, fl, 0, bsize, status)
    call ftgipr(iunit, 3, bpix, naxes, dm, status)
    call ftclos(iunit, status)
    call ftfiou(iunit, status)
#endif /* !TEST */

! calculate number of elements in array
!
    nl = product(dm)

  end subroutine fits_init
!
!===============================================================================
!
! fits_get_dims: subroutine reads dimensions of stored arrays from a FITS file
!
!===============================================================================
!
  subroutine fits_get_dims(dims)

    implicit none

! arguments
!
    integer, dimension(:), intent(inout) :: dims
!
!-------------------------------------------------------------------------------
!
    dims(:) = dm(:)

  end subroutine fits_get_dims
!
!===============================================================================
!
! fits_get_bounds: subroutine returns the limits of the box
!
!===============================================================================
!
  subroutine fits_get_bounds(xmin, xmax, ymin, ymax, zmin, zmax)

    implicit none

! arguments
!
    real, intent(out) :: xmin, xmax, ymin, ymax, zmin, zmax
!
!-------------------------------------------------------------------------------
!
    xmin = -0.5
    xmax =  0.5
    ymin = -0.5
    ymax =  0.5
    zmin = -0.5
    zmax =  0.5

  end subroutine fits_get_bounds
!
!===============================================================================
!
! fits_get_gridsize: subroutine returns the size of the cell
!
!===============================================================================
!
  subroutine fits_get_gridsize(dx, dy, dz)

    implicit none

! arguments
!
    real, intent(out) :: dx, dy, dz
!
!-------------------------------------------------------------------------------
!
    dx = 1.0 / dm(1)
    dy = 1.0 / dm(2)
    dz = 1.0 / dm(3)

  end subroutine fits_get_gridsize
!
!===============================================================================
!
! fits_get_timestep: subroutine returns the timestep
!
!===============================================================================
!
  subroutine fits_get_timestep(dt)

    implicit none

! arguments
!
    real, intent(out) :: dt
!
!-------------------------------------------------------------------------------
!
    dt = 1.0

  end subroutine fits_get_timestep
!
!===============================================================================
!
! fits_read_var: subroutine reads variables from HDF files
!
!===============================================================================
!
  subroutine fits_read_var(var, qty)

    implicit none

! arguments
!
    character(len=*)      , intent(in)    :: var
    real, dimension(:,:,:), intent(inout) :: qty

! local variables
!
    real, dimension(:,:,:), allocatable :: a
!
!-------------------------------------------------------------------------------
!
! allocate temporary arrays
!
    allocate(a(dm(1),dm(2),dm(3)))

    select case(trim(var))
    case("dlog")
      call fits_read_data("dens", qty)
      qty = alog10(qty)
    case("momx")
      call fits_read_data("dens", a)
      call fits_read_data("velx", qty)
      qty = a * qty
    case("momy")
      call fits_read_data("dens", a)
      call fits_read_data("vely", qty)
      qty = a * qty
    case("momz")
      call fits_read_data("dens", a)
      call fits_read_data("velz", qty)
      qty = a * qty
    case("dvlx")
      call fits_read_data("dens", a)
      call fits_read_data("velx", qty)
      qty = a**(1./3.) * qty
    case("dvly")
      call fits_read_data("dens", a)
      call fits_read_data("vely", qty)
      qty = a**(1./3.) * qty
    case("dvlz")
      call fits_read_data("dens", a)
      call fits_read_data("velz", qty)
      qty = a**(1./3.) * qty
    case default
      call fits_read_data(var, qty)
    end select

! deallocate temporary arrays
!
    deallocate(a)

  end subroutine fits_read_var
!
!===============================================================================
!
! fits_read_data: subroutine reads selected array from an HDF file
!
!===============================================================================
!
  subroutine fits_read_data(var, qty)

    use params, only : idir

! arguments
!
    character(len=*)      , intent(in)    :: var
    real, dimension(:,:,:), intent(inout) :: qty

! local variables
!
    character(len=255) :: fl

! FITS variables
!
    logical            :: info
    integer            :: status, iunit, bsize, bpix, naxes
!
!-------------------------------------------------------------------------------
!
! generate filename
!
    write(fl,"(a,a4,'.fits.gz')") trim(idir), trim(var)
    write( *,"(a,a4,'.fits.gz')") 'INFO      : reading from ', trim(var)

! check if the file exists
!
    inquire(file = fl, exist = info)

    if (.not. info) then
      write (*,'("ERROR   : file ",a," does not exist!")') trim(fl)
      stop
    endif

! check if the file is written in FITS format
!
    status  = 0

#ifndef TEST
    call ftgiou(iunit, status)
    call ftopen(iunit, fl, 0, bsize, status)
    call ftgpve(iunit, 0, 1, nl, 0.0, qty, info, status)
    call ftclos(iunit, status)
    call ftfiou(iunit, status)
#endif /* !TEST */

  end subroutine fits_read_data
!
!===============================================================================
!
! FITS_GET_DATA: subroutine reads data from FITS file
!
!===============================================================================
!
  subroutine fits_get_data(rfile, qty, status)

    implicit none

    character(len=*)            , intent(in)    :: rfile
    real, dimension(:,:,:)      , intent(inout) :: qty
    integer                     , intent(out)   :: status

    logical :: anyf, ext
    integer :: nelems, iunit, bsize, bpix, naxis, naxes(3)
!
!-------------------------------------------------------------------------------
!
    inquire(file = rfile, exist = ext)

    if (.not. ext) then
      write (*,'("ERROR   : file ",a," does not exist!")') trim(rfile)
      stop
    endif

    status = 0
    naxes(:) = 1

#ifndef TEST
    call ftgiou(iunit, status)
    call ftopen(iunit, rfile, 0, bsize, status)
    call ftgipr(iunit, 3, bpix, naxis, naxes, status)
    nelems = product(naxes(:))
    call ftgpve(iunit, 0, 1, nelems, 0.0, qty, anyf, status)
    call ftclos(iunit, status)
    call ftfiou(iunit, status)
#endif /* !TEST */


  end subroutine fits_get_data

!===============================================================================
!
! FITS_GET_DATA_2D: subroutine reads data from FITS file
!
!===============================================================================
!
  subroutine fits_get_data_2d(rfile, qty, status)

    implicit none

    character(len=*)            , intent(in)    :: rfile
    real, dimension(:,:)        , intent(inout) :: qty
    integer                     , intent(out)   :: status

    logical :: anyf, ext
    integer :: nelems, iunit, bsize, bpix, naxis, naxes(2)
!
!-------------------------------------------------------------------------------
!
    inquire(file = rfile, exist = ext)

    if (.not. ext) then
      write (*,'("ERROR   : file ",a," does not exist!")') trim(rfile)
      stop
    endif

    status = 0
    naxes(:) = 1

    qty(:,:) = 0.0

#ifndef TEST
    call ftgiou(iunit, status)
    call ftopen(iunit, rfile, 0, bsize, status)
    call ftgipr(iunit, 2, bpix, naxis, naxes, status)
    nelems = naxes(1) * naxes(2)
    call ftgpve(iunit, 0, 1, nelems, 0.0, qty, anyf, status)
    call ftclos(iunit, status)
    call ftfiou(iunit, status)
#endif /* !TEST */

  end subroutine fits_get_data_2d

!===============================================================================
!
! FITS_GET_DATA_1D: subroutine reads data from FITS file
!
!===============================================================================
!
  subroutine fits_get_data_1d(rfile, qty, status)

    implicit none

    character(len=*)            , intent(in)    :: rfile
    real, dimension(:)          , intent(inout) :: qty
    integer                     , intent(out)   :: status

    logical :: anyf, ext
    integer :: nelems, iunit, bsize, bpix, naxis, naxes(1)
!
!-------------------------------------------------------------------------------
!
    inquire(file = rfile, exist = ext)

    if (.not. ext) then
      write (*,'("ERROR   : file ",a," does not exist!")') trim(rfile)
      stop
    endif

    status = 0
    naxes(:) = 1

    qty(:) = 0.0

#ifndef TEST
    call ftgiou(iunit, status)
    call ftopen(iunit, rfile, 0, bsize, status)
    call ftgipr(iunit, 1, bpix, naxis, naxes, status)
    nelems = naxes(1)
    call ftgpve(iunit, 0, 1, nelems, 0.0, qty, anyf, status)
    call ftclos(iunit, status)
    call ftfiou(iunit, status)
#endif /* !TEST */

  end subroutine fits_get_data_1d

!===============================================================================
!
! FITS_GET_FOURIER_DATA: subroutine reads data from two files containing
!                        the real and imaginery parts and merges them into
!                        one complex array
!
!===============================================================================
!
  subroutine fits_get_fourier_data(refile, imfile, fft, status)

    implicit none

    character(len=*)            , intent(in)    :: refile, imfile
    complex, dimension(:,:,:)   , intent(inout) :: fft
    integer                     , intent(out)   :: status

    logical :: anyf, ext
    integer :: nelems, iunit, bsize, bpix, naxis, naxes(3)

    real, dimension(:,:,:), allocatable :: tmp
    integer                     , parameter :: FFTW_ESTIMATE = 64
!
!-------------------------------------------------------------------------------
!
    inquire(file = refile, exist = ext)

    if (.not. ext) then
      write (*,'("ERROR   : file ",a," does not exist!")') trim(refile)
      stop
    endif

    inquire(file = imfile, exist = ext)

    if (.not. ext) then
      write (*,'("ERROR   : file ",a," does not exist!")') trim(imfile)
      stop
    endif

    status = 0
    naxes(:) = 1

    fft(:,:,:) = cmplx(0.0, 0.0)

#ifndef TEST
    call ftgiou(iunit, status)
    call ftopen(iunit, refile, 0, bsize, status)
    call ftgipr(iunit, 3, bpix, naxis, naxes, status)
    nelems = product(naxes(:))
    allocate(tmp(naxes(1),naxes(2),naxes(3)))
    call ftgpve(iunit, 0, 1, nelems, 0.0, tmp, anyf, status)
    call ftclos(iunit, status)
    call ftfiou(iunit, status)

!     print *, minval(tmp), maxval(tmp)

    fft(:,:,:) = fft(:,:,:) + cmplx(tmp(:,:,:), 0.0)

!     print *, minval(real(fft)), maxval(real(fft))
!     print *, minval(imag(fft)), maxval(imag(fft))

    call ftgiou(iunit, status)
    call ftopen(iunit, imfile, 0, bsize, status)
    call ftgipr(iunit, 3, bpix, naxis, naxes, status)
    nelems = product(naxes(:))
    call ftgpve(iunit, 0, 1, nelems, 0.0, tmp, anyf, status)
    call ftclos(iunit, status)
    call ftfiou(iunit, status)
#endif /* !TEST */

!     print *, minval(tmp), maxval(tmp)

    fft(:,:,:) = fft(:,:,:) + cmplx(0.0, tmp(:,:,:))

!     print *, minval(real(fft)), maxval(real(fft))
!     print *, minval(imag(fft)), maxval(imag(fft))

    if (allocated(tmp)) deallocate(tmp)

  end subroutine fits_get_fourier_data

!===============================================================================
!
! FITS_PUT_DATA: subroutine writes data to a file
!
!===============================================================================
!
  subroutine fits_put_data(rfile, qty, status)

    implicit none

    character(len=*)            , intent(in)    :: rfile
    real, dimension(:,:,:)      , intent(in)    :: qty
    integer                     , intent(out)   :: status

    integer :: nelems, iunit, naxes(3)
!
!-------------------------------------------------------------------------------
!
    status = 0

    naxes(1) = size(qty, 1)
    naxes(2) = size(qty, 2)
    naxes(3) = size(qty, 3)

#ifndef TEST
    call ftgiou(iunit, status)
    call ftinit(iunit, rfile, 1, status)
    nelems = product(naxes(:))
    call ftphpr(iunit, .true., -32, 3, naxes, 0, 1, .true., status)
    call ftppre(iunit, 1, 1, nelems, qty, status)
    call ftclos(iunit, status)
    call ftfiou(iunit, status)
#endif /* !TEST */

  end subroutine fits_put_data

!===============================================================================
!
! FITS_PUT_DATA_1D: subroutine writes 1D data to a file
!
!===============================================================================
!
  subroutine fits_put_data_2d(rfile, qty)

    implicit none

    character(len=*)    , intent(in)    :: rfile
    real, dimension(:,:), intent(in)    :: qty

    integer :: nelems, iunit, naxes(2), status
!
!-------------------------------------------------------------------------------
!
    status = 0
    naxes(1) = size(qty, 1)
    naxes(2) = size(qty, 2)

#ifndef TEST
    call ftgiou(iunit, status)
    call ftinit(iunit, rfile, 1, status)
    nelems = naxes(1) * naxes(2)
    call ftphpr(iunit, .true., -32, 2, naxes, 0, 1, .true., status)
    call ftppre(iunit, 1, 1, nelems, qty, status)
    call ftclos(iunit, status)
    call ftfiou(iunit, status)
#endif /* !TEST */

  end subroutine fits_put_data_2d

!===============================================================================
!
! FITS_PUT_DATA_1D: subroutine writes 1D data to a file
!
!===============================================================================
!
  subroutine fits_put_data_1d(rfile, qty)

    implicit none

    character(len=*)            , intent(in)    :: rfile
    real, dimension(:)          , intent(in)    :: qty

    integer :: nelems, iunit, naxes(1), status
!
!-------------------------------------------------------------------------------
!
    status = 0
    naxes(1) = size(qty, 1)

#ifndef TEST
    call ftgiou(iunit, status)
    call ftinit(iunit, rfile, 1, status)
    nelems = naxes(1)
    call ftphpr(iunit, .true., -32, 1, naxes, 0, 1, .true., status)
    call ftppre(iunit, 1, 1, nelems, qty, status)
    call ftclos(iunit, status)
    call ftfiou(iunit, status)
#endif /* !TEST */

  end subroutine fits_put_data_1d

end module fitsio

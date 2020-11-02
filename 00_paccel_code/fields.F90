!!******************************************************************************
!!
!! module: fields - subroutines to access velocity and magnetic fields
!!
!! Copyright (C) 2009-2010 Grzegorz Kowal <grzegorz@gkowal.info>
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
module fields

  implicit none

! the dimesion of the field components
!
  integer, dimension(3), save :: dm, qm

! the domain bounds
!
  real            , save :: xmin, xmax, ymin, ymax, zmin, zmax

! the indices needed to extend arrays
!
  integer, private, save :: ib, jb, kb, ie, je, ke
  integer, private, save :: il, jl, kl, iu, ju, ku
  integer, private, save :: is, js, ks, it, jt, kt

! the allocatable arrays for storing the electric and magnetic field components
!
  real, dimension(:,:,:), allocatable, save :: ux, uy, uz, bx, by, bz

! the number of ghost layers for interpolation
!
#ifdef TRICUB
  integer, parameter :: ng = 3
#else /* TRICUB */
#ifdef TRILIN
  integer, parameter :: ng = 2
#else /* TRILIN */
  integer, parameter :: ng = 1
#endif /* TRILIN */
#endif /* TRICUB */
!
!-------------------------------------------------------------------------------
!
  contains
!
!===============================================================================
!
! init_fields: subroutine initializes the field variables
!
!===============================================================================
!
  subroutine init_fields()

    use params, only : fformat
    use fitsio, only : fits_init, fits_get_dims, fits_get_bounds, fits_read_var
#ifdef HDF5
    use hdf5io, only : hdf5_init, hdf5_get_dims, hdf5_get_bounds, hdf5_read_var
#endif /* HDF5 */

    implicit none

! the local temporary allocatable arrays
!
    real, dimension(:,:,:), allocatable :: tt
!
!-------------------------------------------------------------------------------
!
! initialize dimensions
!
    dm(:) = 1

! obtain array dimensions
!
    select case(fformat)
    case('fits')
      call fits_init('magx')
      call fits_get_dims(dm)
      call fits_get_bounds(xmin, xmax, ymin, ymax, zmin, zmax)
#ifdef HDF5
    case('hdf5')
      call hdf5_init()
      call hdf5_get_dims(dm)
      call hdf5_get_bounds(xmin, xmax, ymin, ymax, zmin, zmax)
#endif /* HDF5 */
    case default
      write( *, "('ERROR     : ',a,1x,a)" ) "unsupported data format:", fformat
      stop
    end select

! prepare extended dimensions
!
    qm(:) = dm(:) + 2 * ng
    if (dm(3) .eq. 1) qm(3) = 1

! prepare indices for array extension
!
    ib =  1 + ng
    ie = ib + dm(1) - 1
    il = ib + ng - 1
    iu = ie - ng + 1
    is = ib - 1
    it = ie + 1

    jb =  1 + ng
    je = jb + dm(2) - 1
    jl = jb + ng - 1
    ju = je - ng + 1
    js = jb - 1
    jt = je + 1

    if (dm(3) .gt. 1) then
      kb =  1 + ng
      ke = kb + dm(3) - 1
      kl = kb + ng - 1
      ku = ke - ng + 1
      ks = kb - 1
      kt = ke + 1
    else
      kb = 1
      ke = 1
      kl = 1
      ku = 1
      ks = 1
      kt = 1
    end if

! allocate space for the field components
!
    allocate(tt(dm(1),dm(2),dm(3)))
    allocate(ux(qm(1),qm(2),qm(3)))
    allocate(uy(qm(1),qm(2),qm(3)))
    allocate(uz(qm(1),qm(2),qm(3)))
    allocate(bx(qm(1),qm(2),qm(3)))
    allocate(by(qm(1),qm(2),qm(3)))
    allocate(bz(qm(1),qm(2),qm(3)))

! read the field components from the file
!
    write( *, "('INFO      : ',a)" ) "reading velocity and magnetic field"
    select case(fformat)
    case('fits')
      call fits_read_var('velx', tt)
      call expand_array(tt, ux)
      call fits_read_var('vely', tt)
      call expand_array(tt, uy)
      call fits_read_var('velz', tt)
      call expand_array(tt, uz)
      call fits_read_var('magx', tt)
      call expand_array(tt, bx)
      call fits_read_var('magy', tt)
      call expand_array(tt, by)
      call fits_read_var('magz', tt)
      call expand_array(tt, bz)
#ifdef HDF5
    case('hdf5')
      call hdf5_read_var('velx', tt)
      call expand_array(tt, ux)
      call hdf5_read_var('vely', tt)
      call expand_array(tt, uy)
      call hdf5_read_var('velz', tt)
      call expand_array(tt, uz)
      call hdf5_read_var('magx', tt)
      call expand_array(tt, bx)
      call hdf5_read_var('magy', tt)
      call expand_array(tt, by)
      call hdf5_read_var('magz', tt)
      call expand_array(tt, bz)
#endif /* HDF5 */
    end select

! deallocate the temporary local array
!
    if (allocated(tt)) deallocate(tt)
!
!-------------------------------------------------------------------------------
!
  end subroutine init_fields
!
!===============================================================================
!
! finit_fields: subroutine deallocates the field variables
!
!===============================================================================
!
  subroutine finit_fields()

    implicit none
!
!-------------------------------------------------------------------------------
!
    if (allocated(ux)) deallocate(ux)
    if (allocated(uy)) deallocate(uy)
    if (allocated(uz)) deallocate(uz)
    if (allocated(bx)) deallocate(bx)
    if (allocated(by)) deallocate(by)
    if (allocated(bz)) deallocate(bz)
!
!-------------------------------------------------------------------------------
!
  end subroutine finit_fields
!
!===============================================================================
!
! get_dimensions: subroutine returns the array dimensions
!
!===============================================================================
!
  subroutine get_dimensions(dims)

    implicit none

! output arguments
!
    integer, dimension(3), intent(out) :: dims
!
!-------------------------------------------------------------------------------
!
    dims(:) = dm(:)
!
!-------------------------------------------------------------------------------
!
  end subroutine get_dimensions
!
!===============================================================================
!
! get_domain_bounds: subroutine returns the domain bounds
!
!===============================================================================
!
  subroutine get_domain_bounds(bounds)

    implicit none

! output arguments
!
    real, dimension(3,2), intent(out) :: bounds
!
!-------------------------------------------------------------------------------
!
    bounds(1,1) = xmin
    bounds(1,2) = xmax
    bounds(2,1) = ymin
    bounds(2,2) = ymax
    bounds(3,1) = zmin
    bounds(3,2) = zmax
!
!-------------------------------------------------------------------------------
!
  end subroutine get_domain_bounds
!
!===============================================================================
!
! expand_array: subroutine expands an input array by adding the boundary layers
!
!===============================================================================
!
  subroutine expand_array(a, b)

    implicit none

! output arguments
!
    real, dimension(dm(1),dm(2),dm(3)), intent(in)  :: a
    real, dimension(qm(1),qm(2),qm(3)), intent(out) :: b
!
!-------------------------------------------------------------------------------
!
! copy the interior
!
    b(ib:ie,jb:je,kb:ke) = a(1:dm(1),1:dm(2),1:dm(3))

! copy the X boundaries
!
    b( 1:is   ,jb:je,kb:ke) = b(iu:ie,jb:je,kb:ke)
    b(it:qm(1),jb:je,kb:ke) = b(ib:il,jb:je,kb:ke)

! copy the Y boundaries
!
    b(1:qm(1), 1:js   ,kb:ke) = b(1:qm(1),ju:je,kb:ke)
    b(1:qm(1),jt:qm(2),kb:ke) = b(1:qm(1),jb:jl,kb:ke)

! copy the Z boundaries
!
    if (dm(3) .gt. 1) then
      b(1:qm(1),1:qm(2), 1:ks   ) = b(1:qm(1),1:qm(2),ku:ke)
      b(1:qm(1),1:qm(2),kt:qm(3)) = b(1:qm(1),1:qm(2),kb:kl)
    end if
!
!-------------------------------------------------------------------------------
!
  end subroutine expand_array
!
end module fields

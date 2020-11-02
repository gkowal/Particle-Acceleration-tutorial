!!******************************************************************************
!!
!! module: hdf5io - subroutines to read data from HDF files
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
module hdf5io

#ifdef HDF5
  use hdf5
#endif /* HDF5 */

  implicit none

  integer, dimension(3), save :: dm, pdm, cdm, fdm, qb, qe
  integer              , save :: nf, ng
  logical              , save :: mhd, mag
  real                 , save :: xmn, xmx, ymn, ymx, zmn, zmx, dx, dy, dz, dxi, dyi, dzi, dtc

!-------------------------------------------------------------------------------
!
  contains
!
!===============================================================================
!
! hdf5_init: subroutine reads attributes from a HDF file
!
!===============================================================================
!
  subroutine hdf5_init

#ifdef HDF5
    use params, only : idir, ftype, fnumber

    implicit none

! local variables
!
    character(len=255) :: fl
    character(len= 32) :: nm
    logical            :: info
    integer            :: err, nattr, nvars, i, tp

! HDF5 variables
!
    integer(hid_t)                 :: fid, gid, aid
    integer(hsize_t)               :: aln = 32
    integer(hsize_t), dimension(1) :: am = (/1/), cm = (/3/)
#endif /* HDF5 */
!
!-------------------------------------------------------------------------------
!
! set default values
!
    nf     = 1
    ng     = 0
    dm (:) = 1
    pdm(:) = 1
    qb (:) = 1
    qe (:) = 1
    mhd    = .false.
    mag    = .false.
    xmn    = 0.0
    xmx    = 1.0
    ymn    = 0.0
    ymx    = 1.0
    zmn    = 0.0
    zmx    = 1.0
    dx     = 1.0
    dy     = 1.0
    dz     = 1.0
    dxi    = 1.0
    dyi    = 1.0
    dzi    = 1.0
    dtc    = 1.0

#ifdef HDF5
! generate filename
!
    write(fl,"(a,a1,i6.6,'_',i5.5,'.h5')") trim(idir), ftype, fnumber, 0

! check if the file exists
!
    inquire(file = fl, exist = info)

    if (.not. info) then
      write (*,'("ERROR   : file ",a," does not exist!")') trim(fl)
      stop
    endif

! check if the file is written in HDF5 format
!
    call h5fis_hdf5_f(fl, info, err)


    if (info) then

! initialize FORTRAN interface
!
      call h5open_f(err)

! open HDF5 file
!
      call h5fopen_f(fl, H5F_ACC_RDONLY_F, fid, err)

! open group 'attributes'
!
      call h5gopen_f(fid, 'attributes', gid, err)

! get number of attributes
!
      call h5aget_num_attrs_f(gid, nattr, err)

! read attributes
!
      do i = 0, nattr - 1
        call h5aopen_idx_f(gid, i, aid, err)
        call h5aget_name_f(aid, aln, nm, err)


        select case(nm)
        case('nprocs')
          call h5aread_f(aid, H5T_NATIVE_INTEGER, nf, am, err)
        case('ng')
          if (ftype .eq. 'r') then
            call h5aread_f(aid, H5T_NATIVE_INTEGER, ng, am, err)
          end if
        case('pdims')
          call h5aread_f(aid, H5T_NATIVE_INTEGER, pdm, cm, err)
        case('dims')
          call h5aread_f(aid, H5T_NATIVE_INTEGER, dm, cm, err)
        case('xmin')
          call h5aread_f(aid, H5T_NATIVE_REAL, xmn, am, err)
        case('xmax')
          call h5aread_f(aid, H5T_NATIVE_REAL, xmx, am, err)
        case('ymin')
          call h5aread_f(aid, H5T_NATIVE_REAL, ymn, am, err)
        case('ymax')
          call h5aread_f(aid, H5T_NATIVE_REAL, ymx, am, err)
        case('zmin')
          call h5aread_f(aid, H5T_NATIVE_REAL, zmn, am, err)
        case('zmax')
          call h5aread_f(aid, H5T_NATIVE_REAL, zmx, am, err)
        case('dx')
          call h5aread_f(aid, H5T_NATIVE_REAL, dx , am, err)
        case('dy')
          call h5aread_f(aid, H5T_NATIVE_REAL, dy , am, err)
        case('dz')
          call h5aread_f(aid, H5T_NATIVE_REAL, dz , am, err)
        case('dxi')
          call h5aread_f(aid, H5T_NATIVE_REAL, dxi, am, err)
        case('dyi')
          call h5aread_f(aid, H5T_NATIVE_REAL, dyi, am, err)
        case('dzi')
          call h5aread_f(aid, H5T_NATIVE_REAL, dzi, am, err)
        case('dt')
          call h5aread_f(aid, H5T_NATIVE_REAL, dtc, am, err)
        end select
        call h5aclose_f(aid, err)
      end do

! close group 'attributes'
!
      call h5gclose_f(gid, err)

! get the number of variables
!
      call h5gn_members_f(fid, 'variables', nvars, err)

! check if file contains magnetic field or vector potential
!
      do i = 0, nvars - 1
        call h5gget_obj_info_idx_f(fid, 'variables', i, nm, tp, err)
        select case(trim(nm))
        case('potx', 'poty', 'potz')
          mhd = .true.
          mag = .false.
        case('magx', 'magy', 'magz')
          mhd = .true.
          mag = .true.
        end select
      enddo

! terminate access to the file
!
      call h5fclose_f(fid, err)

! close FORTRAN interface
!
      call h5close_f(err)

    endif

! calculate chunk dimensions and full dimension
!
    cdm = dm
    where(dm .gt. 1) cdm = dm - 2 * ng
    fdm = cdm * pdm

! calculate dxi, dyi, dzi
!
    if (dxi .eq. 1.0) dxi = fdm(1) / (xmx - xmn)
    if (dyi .eq. 1.0) dyi = fdm(2) / (ymx - ymn)
    if (dzi .eq. 1.0) dzi = fdm(3) / (zmx - zmn)

! calculate dxi, dyi, dzi
!
    if (dx  .eq. 1.0) dx  = (xmx - xmn) / fdm(1)
    if (dy  .eq. 1.0) dy  = (ymx - ymn) / fdm(2)
    if (dz  .eq. 1.0) dz  = (zmx - zmn) / fdm(3)

! calculate source indices
!
    where(dm .gt. 1) qb(:) = ng + 1
    qe(:) = qb(:) + cdm(:) - 1

#endif /* HDF5 */
  end subroutine hdf5_init
!
!===============================================================================
!
! hdf5_get_dims: subroutine reads dimensions of stored arrays from a HDF file
!
!===============================================================================
!
  subroutine hdf5_get_dims(dims)

    implicit none

! arguments
!
    integer, dimension(:), intent(inout) :: dims
!
!-------------------------------------------------------------------------------
!
    dims(:) = fdm(:)

  end subroutine hdf5_get_dims
!
!===============================================================================
!
! hdf5_get_bounds: subroutine returns the limits of the box
!
!===============================================================================
!
  subroutine hdf5_get_bounds(xmin, xmax, ymin, ymax, zmin, zmax)

    implicit none

! arguments
!
    real, intent(out) :: xmin, xmax, ymin, ymax, zmin, zmax
!
!-------------------------------------------------------------------------------
!
    xmin = xmn
    xmax = xmx
    ymin = ymn
    ymax = ymx
    zmin = zmn
    zmax = zmx

  end subroutine hdf5_get_bounds
!
!===============================================================================
!
! hdf5_get_gridsize: subroutine returns the size of the cell
!
!===============================================================================
!
  subroutine hdf5_get_gridsize(gdx, gdy, gdz)

    implicit none

! arguments
!
    real, intent(out) :: gdx, gdy, gdz
!
!-------------------------------------------------------------------------------
!
    gdx = dx
    gdy = dy
    gdz = dz

  end subroutine hdf5_get_gridsize
!
!===============================================================================
!
! hdf5_get_timestep: subroutine returns the timestep
!
!===============================================================================
!
  subroutine hdf5_get_timestep(dt)

    implicit none

! arguments
!
    real, intent(out) :: dt
!
!-------------------------------------------------------------------------------
!
    dt = dtc

  end subroutine hdf5_get_timestep
!
!===============================================================================
!
! hdf5_get_coord: subroutine reads the coordinate of subarray
!
!===============================================================================
!
  subroutine hdf5_get_coord(n, coord)

    use params, only : idir, ftype, fnumber

    implicit none

! arguments
!
    integer               , intent(in)    :: n
    integer, dimension(3) , intent(out)   :: coord

#ifdef HDF5
! local variables
!
    character(len=255) :: fl
    character(len= 32) :: at
    logical            :: info
    integer            :: err, nattr, i

! HDF5 variables
!
    integer(hid_t)                 :: fid, gid, aid
    integer(hsize_t)               :: aln = 32
    integer(hsize_t), dimension(1) :: am = (/1/)
#endif /* HDF5 */
!
!-------------------------------------------------------------------------------
!
! set default values
!
    coord(:) = 0

#ifdef HDF5
! generate filename
!
    write(fl,"(a,a1,i6.6,'_',i5.5,'.h5')") trim(idir), ftype, fnumber, n

! check if the file exists
!
    inquire(file = fl, exist = info)

    if (.not. info) then
      write (*,'("ERROR   : file ",a," does not exist!")') trim(fl)
      stop
    endif

! check if the file is written in HDF5 format
!
    call h5fis_hdf5_f(fl, info, err)

    if (info) then

! initialize FORTRAN interface
!
      call h5open_f(err)

! open HDF5 file
!
      call h5fopen_f(fl, H5F_ACC_RDONLY_F, fid, err)

! open group 'attributes'
!
      call h5gopen_f(fid, 'attributes', gid, err)

! get number of attributes
!
      call h5aget_num_attrs_f(gid, nattr, err)

! read attributes
!
      do i = 0, nattr - 1
        call h5aopen_idx_f(gid, i, aid, err)
        call h5aget_name_f(aid, aln, at, err)
        select case(trim(at))
        case('pcoords')
          am(1) = 3
          call h5aread_f(aid, H5T_NATIVE_INTEGER, coord, am, err)
        end select
        call h5aclose_f(aid, err)
      end do

! close group 'attributes'
!
      call h5gclose_f(gid, err)

! terminate access to the file
!
      call h5fclose_f(fid, err)

! close FORTRAN interface
!
      call h5close_f(err)

    endif
#endif /* HDF5 */

  end subroutine hdf5_get_coord
!
!===============================================================================
!
! hdf5_read_var: subroutine reads variables from HDF files
!
!===============================================================================
!
  subroutine hdf5_read_var(var, qty)

    use params, only : ftype, fnumber

    implicit none

! arguments
!
    character(len=*)      , intent(in)    :: var
    real, dimension(:,:,:), intent(inout) :: qty

! local variables
!
    integer                             :: n
    integer, dimension(3)               :: coord, pb, pe
    real, dimension(:,:,:), allocatable :: q
!
!-------------------------------------------------------------------------------
!
    qty(:,:,:) = 0.0

#ifdef HDF5
! allocate small subarray for variable chunks
!
    allocate(q(cdm(1),cdm(2),cdm(3)))

! iterate over all files
!
    do n = 0, nf - 1

! print info
!
    write( *,"(a,a4,a,a1,i6.6,'_',i5.5,'.h5',a1,$)") 'INFO      : reading ', trim(var), ' from ', ftype, fnumber, n, char(13)

! get coordinate of subarray
!
      call hdf5_get_coord(n, coord)

! get variable chunk from each file
!
      call hdf5_get_var(n, var, q)

! calculate destination indices
!
      pb(:) = coord(:) * cdm(:) + 1
      pe(:) =    pb(:) + cdm(:) - 1

! insert a chunk into variable array
!
      qty(pb(1):pe(1),pb(2):pe(2),pb(3):pe(3)) = q(:,:,:)

    end do

! deallocate subarray
!
    deallocate(q)

    write(*,"('')")
#endif /* HDF5 */

  end subroutine hdf5_read_var
!
!===============================================================================
!
! hdf5_get_var: subroutine reads subarray from a HDF file
!
!===============================================================================
!
  subroutine hdf5_get_var(n, var, qty)

    use params, only : ftype

! arguments
!
    integer               , intent(in)    :: n
    character(len=*)      , intent(in)    :: var
    real, dimension(:,:,:), intent(inout) :: qty

! local variables
!
    integer                             :: i, j, k
    real, dimension(:,:,:), allocatable :: a, q
!
!-------------------------------------------------------------------------------
!
! in this subroutine we read conserved or primitive variables and construct desired variable
!
    allocate(q(dm(1),dm(2),dm(3)))
    if (ftype .eq. 'r') then
      select case(trim(var))
      case("dlog")
        call hdf5_read_data(n, "dens", q)
        q = alog10(q)
      case("velx")
        allocate(a(dm(1),dm(2),dm(3)))
        call hdf5_read_data(n, "dens", a)
        call hdf5_read_data(n, "momx", q)
        q = q / a
        deallocate(a)
      case("vely")
        allocate(a(dm(1),dm(2),dm(3)))
        call hdf5_read_data(n, "dens", a)
        call hdf5_read_data(n, "momy", q)
        q = q / a
        deallocate(a)
      case("velz")
        allocate(a(dm(1),dm(2),dm(3)))
        call hdf5_read_data(n, "dens", a)
        call hdf5_read_data(n, "momz", q)
        q = q / a
        deallocate(a)
      case("dvlx")
        allocate(a(dm(1),dm(2),dm(3)))
        call hdf5_read_data(n, "dens", a)
        call hdf5_read_data(n, "momx", q)
        q = q / a**(2./3.)
        deallocate(a)
      case("dvly")
        allocate(a(dm(1),dm(2),dm(3)))
        call hdf5_read_data(n, "dens", a)
        call hdf5_read_data(n, "momy", q)
        q = q / a**(2./3.)
        deallocate(a)
      case("dvlz")
        allocate(a(dm(1),dm(2),dm(3)))
        call hdf5_read_data(n, "dens", a)
        call hdf5_read_data(n, "momz", q)
        q = q / a**(2./3.)
        deallocate(a)
      case("magx")
        if (mag) then
          call hdf5_read_data(n, "magx", q)
        else
          allocate(a(dm(1),dm(2),dm(3)))

          call hdf5_read_data(n, "potz", a)

          do j = 2, dm(2)
            q(:,j,:) =              dyi * (a(:,j,:) - a(:,j-1,:))
          end do

          if (dm(3) .gt. 1) then
            call hdf5_read_data(n, "poty", a)

            do k = 2, dm(3)
              q(:,:,k) = q(:,:,k) - dzi * (a(:,:,k) - a(:,:,k-1))
            end do
          end if

          deallocate(a)
        endif
      case("magy")
        if (mag) then
          call hdf5_read_data(n, "magy", q)
        else
          allocate(a(dm(1),dm(2),dm(3)))

          call hdf5_read_data(n, "potz", a)

          do i = 2, dm(1)
            q(i,:,:) =            - dxi * (a(i,:,:) - a(i-1,:,:))
          end do

          if (dm(3) .gt. 1) then
            call hdf5_read_data(n, "potx", a)

            do k = 2, dm(3)
              q(:,:,k) = q(:,:,k) + dzi * (a(:,:,k) - a(:,:,k-1))
            end do
          end if

          deallocate(a)
        endif
      case("magz")
        if (mag) then
          call hdf5_read_data(n, "magz", q)
        else
          allocate(a(dm(1),dm(2),dm(3)))

          call hdf5_read_data(n, "poty", a)

          do i = 2, dm(1)
            q(i,:,:) =            dxi * (a(i,:,:) - a(i-1,:,:))
          end do

          call hdf5_read_data(n, "potx", a)

          do j = 2, dm(2)
            q(:,j,:) = q(:,j,:) - dyi * (a(:,j,:) - a(:,j-1,:))
          end do

          deallocate(a)
        endif
      case default
        call hdf5_read_data(n, var, q)
      end select
    else
      select case(trim(var))
      case("dlog")
        call hdf5_read_data(n, "dens", q)
        q = alog10(q)
      case("momx")
        allocate(a(dm(1),dm(2),dm(3)))
        call hdf5_read_data(n, "dens", a)
        call hdf5_read_data(n, "velx", q)
        q = a * q
        deallocate(a)
      case("momy")
        allocate(a(dm(1),dm(2),dm(3)))
        call hdf5_read_data(n, "dens", a)
        call hdf5_read_data(n, "vely", q)
        q = a * q
        deallocate(a)
      case("momz")
        allocate(a(dm(1),dm(2),dm(3)))
        call hdf5_read_data(n, "dens", a)
        call hdf5_read_data(n, "velz", q)
        q = a * q
        deallocate(a)
      case("dvlx")
        allocate(a(dm(1),dm(2),dm(3)))
        call hdf5_read_data(n, "dens", a)
        call hdf5_read_data(n, "velx", q)
        q = a**(1./3.) * q
        deallocate(a)
      case("dvly")
        allocate(a(dm(1),dm(2),dm(3)))
        call hdf5_read_data(n, "dens", a)
        call hdf5_read_data(n, "vely", q)
        q = a**(1./3.) * q
        deallocate(a)
      case("dvlz")
        allocate(a(dm(1),dm(2),dm(3)))
        call hdf5_read_data(n, "dens", a)
        call hdf5_read_data(n, "velz", q)
        q = a**(1./3.) * q
        deallocate(a)
      case default
        call hdf5_read_data(n, var, q)
      end select
    end if

! copy variable
!
    qty(:,:,:) = q(qb(1):qe(1),qb(2):qe(2),qb(3):qe(3))

! deallocate variable sub-array
!
    deallocate(q)

  end subroutine hdf5_get_var
!
!===============================================================================
!
! hdf5_read_data: subroutine reads selected array from an HDF file
!
!===============================================================================
!
  subroutine hdf5_read_data(n, var, qty)

    use params, only : idir, ftype, fnumber

! arguments
!
    integer               , intent(in)    :: n
    character(len=*)      , intent(in)    :: var
    real, dimension(:,:,:), intent(inout) :: qty

#ifdef HDF5
! local variables
!
    character(len=255) :: fl
    character(len= 32) :: nm
    logical            :: info
    integer            :: err, nds, i

! local arrays
!
    real(kind=4), dimension(:,:,:), allocatable :: tmp4
    real(kind=8), dimension(:,:,:), allocatable :: tmp8

! HDF variables
!
    integer(hid_t)                 :: fid, gid, aid, did, tid
    integer(hsize_t), dimension(3) :: cm
!
!-------------------------------------------------------------------------------
!
! generate filename
!
    write(fl,"(a,a1,i6.6,'_',i5.5,'.h5')") trim(idir), ftype, fnumber, n

! check if the file exists
!
    inquire(file = fl, exist = info)

    if (.not. info) then
      write (*,'("ERROR   : ",a," could not be found!")') trim(fl)
      stop
    endif

! check if the file is written in HDF5 format
!
    call h5fis_hdf5_f(fl, info, err)

    if (info) then

! initialize FORTRAN interface
!
      call h5open_f(err)

! open HDF5 file
!
      call h5fopen_f(fl, H5F_ACC_RDONLY_F, fid, err)

! open group 'variables'
!
      call h5gopen_f(fid, 'variables', gid, err)

! open dataset
!
      call h5dopen_f(gid, trim(var), did, err)

! get datatype
!
      call h5dget_type_f(did, tid, err)

      cm(:) = dm(:)

      call h5tequal_f(tid, H5T_NATIVE_DOUBLE, info, err)

      if (info) then
        allocate(tmp8(dm(1),dm(2),dm(3)))
        call h5dread_f(did, H5T_NATIVE_DOUBLE, tmp8, cm(1:3), err)
        qty = real(tmp8,4)
        deallocate(tmp8)
      else
        allocate(tmp4(dm(1),dm(2),dm(3)))
        call h5dread_f(did, H5T_NATIVE_REAL, tmp4, cm(1:3), err)
        qty = real(tmp4,4)
        deallocate(tmp4)
      endif

! close dataset
!
      call h5dclose_f(did, err)

! close group 'variables'
!
      call h5gclose_f(gid, err)

! terminate access to the file
!
      call h5fclose_f(fid, err)

! close FORTRAN interface
!
      call h5close_f(err)

    end if
#endif /* HDF5 */

  end subroutine hdf5_read_data
!
!===============================================================================
!
! hdf5_get_coords: subroutine reads selected array from an HDF file
!
!===============================================================================
!
  subroutine hdf5_get_coords(x, y, z)

    implicit none

! arguments
!
    real, dimension(:), intent(inout) :: x, y, z

! local variables
!
    integer            :: i, j, k
!
!-------------------------------------------------------------------------------
!
    x = ((/(i, i = 1, fdm(1))/) - 0.5) * dx + xmn
    y = ((/(j, j = 1, fdm(2))/) - 0.5) * dy + ymn
    z = ((/(k, k = 1, fdm(3))/) - 0.5) * dz + zmn
    if (dm(3) .eq. 1) z(:) = 0.0

  end subroutine hdf5_get_coords

end module hdf5io

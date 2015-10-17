#include "cctk.h"

program CAR_test

  CCTK_INT :: nx, ny, nz
  CCTK_INT, dimension(:,:,:), allocatable :: mask
  CCTK_INT, dimension(:), allocatable :: sum_x, sum_y, sum_z, &
       sig_x, sig_y, sig_z
  CCTK_INT, dimension(3,2) :: bbox, newbbox1, newbbox2
  CCTK_INT :: min_width
  CCTK_REAL :: min_density
  CCTK_INT :: didit

  CCTK_INT :: ierr, i, j

  nx = 10
  ny = 12
  nz = 8
!!$  nx = 9
!!$  ny = 8
!!$  nz = 1

  min_width = 3
!!$  min_density = 0.8
  min_density = 0.9

  didit = -1

  bbox(:,1) = 1
  bbox(1,2) = nx
  bbox(2,2) = ny
  bbox(3,2) = nz

  newbbox1 = 0
  newbbox2 = 0

  allocate(mask(nx,ny,nz),&
       sum_x(nx),sig_x(nx),&
       sum_y(ny),sig_y(ny),&
       sum_z(nz),sig_z(nz),STAT=ierr)
  if (ierr .ne. 0) call CCTK_WARN(0, "Allocation error")

  sum_x = 0
  sum_y = 0
  sum_z = 0
  sig_x = 0
  sig_y = 0
  sig_z = 0

  mask = 0

!!$  Pruning check:
!!$  mask(3:7,4:8,1:6) = 1

!!$  Split at hole check (x)
!!$  mask(:4,:,:) = 1
!!$  mask(7:,:,:) = 1

!!$  Split at hole check (y)
!!$  mask(:,:5,:) = 1
!!$  mask(:,8:,:) = 1

!!$  Split at hole check (z)
!!$  mask(:,:,:3) = 1
!!$  mask(:,:,6:) = 1

!!$  Density check
!!$  mask = 1
!!$  mask(4:5,7:9,3:5) = 0

!!$  Signature check
!!$  This should have a density of 2/3
  mask(:4,:8,:) = 1
  mask(5:,5:,:) = 1

!!$  Signature check
!!$  This should have a density of 3/4
!!$  mask(:5,:,:) = 1
!!$  mask(6:,5:,:) = 1

!!$  Signature check
!!$  This should have a density of just under 0.8
!!$  mask(:5,:,:) = 1
!!$  mask(6:,6:,:) = 1

!!$  Test from
!!$ http://www-lmc.imag.fr/IDOPT/AGRIF/documentation/doc_AGRIF3DUser/node9.html
!!$  Requires
!!$  nx = 9
!!$  ny = 8
!!$  min_density > 0.681
!!$  This one is a bit pathological; I don't think the website description
!!$  is correct
!!$  mask(1:2,:,:) = 1
!!$  mask(3,1:2,:) = 1; mask(3,4:5,:) = 1; mask(3,7:8,:) = 1
!  mask(3,1:5,:) = 1; mask(3,7:8,:) = 1
!  mask(3,:,:)   = 1
!!$  mask(4,:,:)   = 1
!!$  mask(5,3:6,:) = 1
!  mask(6,3,:)   = 1; mask(6,5:6,:) = 1
!!$  mask(6,3:6,:) = 1
!!$  mask(7:9,3:6,:) = 1

  do j = ny, 1, -1
    do i = 1, nx
      write(*,'(i1," ")',advance="no") mask(i,j,1)
    end do
    write(*,*) " "
  end do

  call check_box(nx, ny, nz, mask, sum_x, sum_y, sum_z, &
       sig_x, sig_y, sig_z, bbox, newbbox1, newbbox2, &
       min_width, min_density, didit)
  
  write(*,*) "Done. Did it:", didit

  select case (didit)
    case(0)
      write(*,*) "Unchanged"
    case(1)
      write(*,*) "Single new bbox"
      write(*,'(a10,"(",i2,":",i2,"):(",i2,":",i2,"):(",i2,":",i2,")")') &
           "new bbox: ", newbbox1(1,:),newbbox1(2,:),newbbox1(3,:)
    case(2)
      write(*,*) "Two new bboxes"
      write(*,'(a11,"(",i2,":",i2,"):(",i2,":",i2,"):(",i2,":",i2,")")') &
           "new bbox1: ", newbbox1(1,:),newbbox1(2,:),newbbox1(3,:)
      write(*,'(a11,"(",i2,":",i2,"):(",i2,":",i2,"):(",i2,":",i2,")")') &
           "new bbox2: ", newbbox2(1,:),newbbox2(2,:),newbbox2(3,:)
    case default
      write(*,*) "Error in return; didit is",didit
  end select

end program CAR_test

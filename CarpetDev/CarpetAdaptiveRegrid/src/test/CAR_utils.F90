#include "cctk.h"

!!$ Given the mask set the 1D sums

subroutine count_points(nx, ny, nz, mask, sum_x, sum_y, sum_z)

  implicit none

  CCTK_INT :: nx, ny, nz
  CCTK_INT :: mask(nx, ny, nz)
  CCTK_INT :: sum_x(nx), sum_y(ny), sum_z(nz)
  
  CCTK_INT :: i, j, k
  
  do i = 1, nx
    
    sum_x(i) = 0
    
  end do
  
  do j = 1, ny
    
    sum_y(j) = 0
    
  end do
  
  do k = 1, nz
    
    sum_z(k) = 0
    
  end do
  
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        
        sum_x(i) = sum_x(i) + mask(i, j, k)
        sum_y(j) = sum_y(j) + mask(i, j, k)
        sum_z(k) = sum_z(k) + mask(i, j, k)
        
      end do
    end do
  end do
  
end subroutine count_points

!!$ Given a sum compute the signature

subroutine signature(nx, lsum, sig)

  implicit none
  
  CCTK_INT :: nx
  CCTK_INT :: lsum(nx), sig(nx)
  
  CCTK_INT :: i
  
  sig(1)  = 0
  sig(nx) = 0
  
  do i = 2, nx - 1
    
    sig(i) = lsum(i - 1) - 2 * lsum(i) + lsum(i + 1)
    
  end do
  
end subroutine signature

!!$ Given a sum prune the zeros.
!!$ That is, find the minimum and maximum non-zero indices

subroutine prune(nx, lsum, ilo, ihi)

  implicit none
  
  CCTK_INT :: nx
  CCTK_INT :: lsum(nx)
  CCTK_INT :: ilo, ihi
  
  CCTK_INT :: i
  
  ilo = 0
  ihi = 0
  
  i = 0
  
11 if (ilo .eq. 0) then
    i = i + 1
    if (i > nx) then
      call CCTK_WARN(0, "Error in prune; sum is all zero!")
    end if
    if (lsum(i) > 0) then
      ilo = i
    end if
    goto 11
  end if
  
  i = nx + 1
12 if (ihi .eq. 0) then
    i = i - 1
!!$    write(*,*) "prune:",nx,i,lsum(i)
    if (i < 1) then
      call CCTK_WARN(0, "Error in prune; sum is all zero!")
    end if
    if (lsum(i) > 0) then
      ihi = i
    end if
    goto 12
  end if
  
end subroutine prune

!!$ Prune box.
!!$ didit is 0 is nothing was done, otherwise 1
!!$ newbbox is always set.

subroutine prune_box(nx, ny, nz, sum_x, sum_y, sum_z, &
                     bbox, newbbox, didit)

  implicit none
  
  CCTK_INT :: nx, ny, nz
  CCTK_INT :: sum_x(nx), sum_y(ny), sum_z(nz)
  CCTK_INT :: bbox(3,3), newbbox(3,3)
  CCTK_INT :: didit
  
  CCTK_INT :: lo, hi
  
  didit = 0

  newbbox(:,3) = bbox(:,3)

!!$  x direction

  call prune(nx, sum_x, lo, hi)
  
  newbbox(1,1) = bbox(1,1) + lo - 1
  newbbox(1,2) = bbox(1,2) + hi - nx

!!$  write(*,*) "x dirn",lo,hi,nx,bbox(1,:)
  
  if ( (lo > 1) .or. (hi < nx) ) then
    didit = 1
  end if
  
!!$  y direction

  call prune(ny, sum_y, lo, hi)
  
  newbbox(2,1) = bbox(2,1) + lo - 1
  newbbox(2,2) = bbox(2,2) + hi - ny

!!$  write(*,*) "y dirn",lo,hi,ny,bbox(2,:)
  
  if ( (lo > 1) .or. (hi < ny) ) then
    didit = 1
  end if
  
!!$  z direction

  call prune(nz, sum_z, lo, hi)
  
  newbbox(3,1) = bbox(3,1) + lo - 1
  newbbox(3,2) = bbox(3,2) + hi - nz

!!$  write(*,*) "z dirn",lo,hi,nz,bbox(3,:)
  
  if ( (lo > 1) .or. (hi < nz) ) then
    didit = 1
  end if
  
end subroutine prune_box

!!$ Find holes.
!!$ That is, find the first location within the grid 
!!$ (taking into account the minimum width parameters)
!!$ that contains a zero
!!$
!!$ ihole is set to zero if there is no hole

subroutine find_hole(nx, lsum, min_width, ihole)
  
  implicit none
  
  CCTK_INT :: nx
  CCTK_INT :: lsum(nx)
  CCTK_INT :: min_width, ihole
  
  CCTK_INT :: i
  
  ihole = 0
  
  if (nx < 2 * min_width + 1) return
  
!!$  write(*,*) "find hole",nx,min_width

  i = min_width 
13 if ( (ihole .eq. 0) .and. (i < nx - min_width + 1) ) then
    i = i + 1
!!$    write(*,*) "fh",i,lsum(i)
    if (lsum(i) .eq. 0) then
      ihole = i
    end if
    goto 13
  end if

!!$  write(*,*) "Done find hole",ihole

end subroutine find_hole

!!$ Split a box at a hole if one exists
!!$ didit is zero if nothing is done and 1 otherwise
!!$ newbbox1 and 2 are always set, but newbbox2 is 
!!$ meaningless if nothing is done.

subroutine split_box_at_hole(nx, ny, nz, sum_x, sum_y, sum_z, &
                             bbox, newbbox1, newbbox2, min_width, &
                             didit)

  implicit none

  CCTK_INT :: nx, ny, nz
  CCTK_INT :: sum_x(nx), sum_y(ny), sum_z(nz)
  CCTK_INT :: bbox(3,3), newbbox1(3,3), newbbox2(3,3)
  CCTK_INT :: min_width
  CCTK_INT :: didit
  
  CCTK_INT :: ihole
  
  didit = 0
  
  newbbox1(1,1) = bbox(1,1)
  newbbox1(2,1) = bbox(2,1)
  newbbox1(3,1) = bbox(3,1)
  newbbox1(1,2) = bbox(1,2)
  newbbox1(2,2) = bbox(2,2)
  newbbox1(3,2) = bbox(3,2)
  
  newbbox2(1,1) = bbox(1,1)
  newbbox2(2,1) = bbox(2,1)
  newbbox2(3,1) = bbox(3,1)
  newbbox2(1,2) = bbox(1,2)
  newbbox2(2,2) = bbox(2,2)
  newbbox2(3,2) = bbox(3,2)
  
  newbbox1(:,3) = bbox(:,3)
  newbbox2(:,3) = bbox(:,3)

!!$  write(*,*) nx,ny,nz

!!$  Try splitting in z first

  call find_hole(nz, sum_z, min_width, ihole)

!!$  write(*,*) "Find hole in z:", ihole

  if (ihole > 0) then

    didit = 1
    
    newbbox1(3,2) = bbox(3,1) + ihole - 2
    newbbox2(3,1) = newbbox1(3,2) + 1
    
    return
    
  end if
  
!!$  Then try splitting in y next

  call find_hole(ny, sum_y, min_width, ihole)

!!$  write(*,*) "Find hole in y:", ihole

  if (ihole > 0) then
    
    didit = 1
    
    newbbox1(2,2) = bbox(2,1) + ihole - 2
    newbbox2(2,1) = newbbox1(2,2) + 1
    
    return
    
  end if

!!$  Then try splitting in x last

  call find_hole(nx, sum_x, min_width, ihole)

!!$  write(*,*) "Find hole in x:", ihole

  if (ihole > 0) then
    
    didit = 1
    
    newbbox1(1,2) = bbox(1,1) + ihole - 2
    newbbox2(1,1) = newbbox1(1,2) + 1
    
    return
    
  end if
  
!!$  We should only reach here if we did not find any holes
!!$  So set newbbox2 to dummy values

  newbbox2(1,1) = 0
  newbbox2(2,1) = 0
  newbbox2(3,1) = 0
  newbbox2(1,2) = 0
  newbbox2(2,2) = 0
  newbbox2(3,2) = 0
  
end subroutine split_box_at_hole

!!$ Find the density.
!!$ That is, find the ratio of marked points to total points.

subroutine find_density(nx, ny, nz, mask, density)

  implicit none

  CCTK_INT :: nx, ny, nz
  CCTK_INT :: mask(nx, ny, nz)
  CCTK_REAL :: density
  
  CCTK_INT :: marked, i, j, k
  
  marked = 0
  
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        
        marked = marked + mask(i, j, k)
        
      end do
    end do
  end do
  
  density = dble(marked) / dble(nx * ny * nz)
  
end subroutine find_density

!!$ Split along a direction.
!!$ That is, find the zero crossing in the signature
!!$ (taking into account the minimum width parameters)
!!$ with the largest first difference in the signature.
!!$
!!$ If the grid should not be split, isplit is zero
!!$ max_jump is passed in so that we can try all directions.

subroutine split(nx, signature, min_width, isplit, max_jump)

  implicit none
  
  CCTK_INT :: nx
  CCTK_INT :: signature(nx)
  CCTK_INT :: min_width, isplit, max_jump
  
  CCTK_INT :: i
  CCTK_INT :: jump
  
  isplit = 0
  
  if (nx < 2 * min_width + 1) return
  
  do i = min_width + 1, nx - min_width
    
!!$    write(*,*) "split:",i,signature(i),signature(i-1)

    if (signature(i) * signature(i - 1) < 0) then
      
      jump = abs( signature(i) - signature(i - 1) )
      
!!$      write(*,*) "jump at",i,"is",jump
      
      if (jump > max_jump) then
        
        isplit = i
        max_jump = jump
        
      end if
      
    end if
    
  end do
  
end subroutine split

!!$ Split the box at turning points in the signature
!!$ Arguments and general behaviour are identical to split_box_at_hole

subroutine split_box_at_sig(nx, ny, nz, sig_x, sig_y, sig_z, &
                            bbox, newbbox1, newbbox2, min_width, &
                            didit)

  implicit none

  CCTK_INT :: nx, ny, nz
  CCTK_INT :: sig_x(nx), sig_y(ny), sig_z(nz)
  CCTK_INT :: bbox(3,3), newbbox1(3,3), newbbox2(3,3)
  CCTK_INT :: min_width
  CCTK_INT :: didit
  
  CCTK_INT :: isplit, max_jump
  
  didit = 0
  max_jump = 0
  
  newbbox1(1,1) = bbox(1,1)
  newbbox1(2,1) = bbox(2,1)
  newbbox1(3,1) = bbox(3,1)
  newbbox1(1,2) = bbox(1,2)
  newbbox1(2,2) = bbox(2,2)
  newbbox1(3,2) = bbox(3,2)
  
  newbbox2(1,1) = bbox(1,1)
  newbbox2(2,1) = bbox(2,1)
  newbbox2(3,1) = bbox(3,1)
  newbbox2(1,2) = bbox(1,2)
  newbbox2(2,2) = bbox(2,2)
  newbbox2(3,2) = bbox(3,2)
  
  newbbox1(:,3) = bbox(:,3)
  newbbox2(:,3) = bbox(:,3)

!!$  Try splitting in z first

  call split(nz, sig_z, min_width, isplit, max_jump)
  
!!$  write(*,*) "Split in z at sig:",isplit
  
  if (isplit > 0) then
    
    didit = 1
    
    newbbox1(3,2) = bbox(3,1) + isplit - 2
    newbbox2(3,1) = newbbox1(3,2) + 1
    
  end if
  
!!$  Then try splitting in y next

  call split(ny, sig_y, min_width, isplit, max_jump)
  
!!$  write(*,*) "Split in y at sig:",isplit
  
  if (isplit > 0) then
    
    didit = 1
    
    newbbox1(3,2) = bbox(3,2)
    newbbox2(3,1) = bbox(3,1)
    newbbox1(2,2) = bbox(2,1) + isplit - 2
    newbbox2(2,1) = newbbox1(2,2) + 1
    
  end if
  
!!$  Then try splitting in x last

  call split(nx, sig_x, min_width, isplit, max_jump)
  
!!$  write(*,*) "Split in x at sig:",isplit
  
  if (isplit > 0) then
    
    didit = 1
    
    newbbox1(3,2) = bbox(3,2)
    newbbox2(3,1) = bbox(3,1)
    newbbox1(2,2) = bbox(2,2)
    newbbox2(2,1) = bbox(2,1)
    newbbox1(1,2) = bbox(1,1) + isplit - 2
    newbbox2(1,1) = newbbox1(1,2) + 1
    
  end if
  
  if (didit > 0) return

!!$  We should only reach here if we did not find any holes
!!$  So set newbbox2 to dummy values

  newbbox2(1,1) = 0
  newbbox2(2,1) = 0
  newbbox2(3,1) = 0
  newbbox2(1,2) = 0
  newbbox2(2,2) = 0
  newbbox2(3,2) = 0
  
end subroutine split_box_at_sig

!!$ Prune a new box
!!$ This is to be done after the grid is split

subroutine prune_new_box(nx, ny, nz, mask, bbox, newbbox)

  implicit none

  CCTK_INT :: nx, ny, nz
  CCTK_INT :: mask(nx, ny, nz)
  CCTK_INT :: bbox(3,3), newbbox(3,3)

  CCTK_INT, dimension(:,:,:), allocatable :: tmp_mask
  CCTK_INT, dimension(:), allocatable :: tmp_sum_x, tmp_sum_y, tmp_sum_z
  CCTK_INT, dimension(:), allocatable :: tmp_sig_x, tmp_sig_y, tmp_sig_z
  CCTK_INT :: tmp_nx, tmp_ny, tmp_nz
  CCTK_INT :: tmp_bbox(3,3)
  CCTK_INT :: tmp_didit
  CCTK_INT :: ierr

  tmp_nx = newbbox(1,2) - newbbox(1,1) + 1
  tmp_ny = newbbox(2,2) - newbbox(2,1) + 1
  tmp_nz = newbbox(3,2) - newbbox(3,1) + 1

  newbbox(:,3) = bbox(:,3)
  
  allocate(tmp_mask(tmp_nx,tmp_ny,tmp_nz),&
           tmp_sum_x(tmp_nx),tmp_sum_y(tmp_ny),tmp_sum_z(tmp_nz),&
           tmp_sig_x(tmp_nx),tmp_sig_y(tmp_ny),tmp_sig_z(tmp_nz),&
           STAT=ierr)
  if (ierr .ne. 0) call CCTK_WARN(0, "Allocation error!")
  
  call copy_mask(nx, ny, nz, mask, bbox, &
                 tmp_nx, tmp_ny, tmp_nz, tmp_mask, newbbox)

  call count_points(tmp_nx, tmp_ny, tmp_nz, tmp_mask, &
                    tmp_sum_x, tmp_sum_y, tmp_sum_z)

  call prune_box(tmp_nx, tmp_ny, tmp_nz, &
                 tmp_sum_x, tmp_sum_y, tmp_sum_z, &
                 newbbox, tmp_bbox, tmp_didit)

!!$  write(*,*) "Pruning split box. Did it:", tmp_didit

  if (tmp_didit > 0) then
    newbbox = tmp_bbox
  end if
  
  deallocate(tmp_mask, tmp_sum_x, tmp_sum_y, tmp_sum_z, &
             tmp_sig_x, tmp_sig_y, tmp_sig_z, &
             STAT=ierr)
  if (ierr .ne. 0) call CCTK_WARN(0, "Deallocation error!")

end subroutine prune_new_box

!!$ Check a box
!!$ This function runs through the entire clustering algorithm for a 
!!$ single box.
!!$
!!$ The return value didit has 3 states:
!!$
!!$     0: nothing happened so this box is now fine
!!$     1: the box requires pruning
!!$     2: the box requires splitting
!!$
!!$ As the new box would require memory, pass it back up to the C++ code.

subroutine check_box(nx, ny, nz, mask, sum_x, sum_y, sum_z, &
                                       sig_x, sig_y, sig_z, &
                     bbox, newbbox1, newbbox2, &
                     min_width, min_density, &
                     didit)

  implicit none

  CCTK_INT :: nx, ny, nz
  CCTK_INT :: mask(nx, ny, nz)
  CCTK_INT :: sum_x(nx), sum_y(ny), sum_z(nz)
  CCTK_INT :: sig_x(nx), sig_y(ny), sig_z(nz)
  CCTK_INT :: bbox(3,3), newbbox1(3,3), newbbox2(3,3)
  CCTK_INT :: min_width
  CCTK_REAL :: min_density
  CCTK_INT :: didit
  
  CCTK_INT :: i

  CCTK_REAL :: density

!!$  write(*,*) nx, ny, nz

!!$  First set up the sums

  call count_points(nx, ny, nz, mask, sum_x, sum_y, sum_z)

!!$  Then prune the box

  call prune_box(nx, ny, nz, sum_x, sum_y, sum_z, bbox, newbbox1, &
                 didit)

!!$  If it needs pruning go back to C++.

!!$  write(*,*) "Pruned. Did it:", didit

  if (didit > 0) return

!!$  Otherwise the box doesn't need pruning. Try finding a hole.

  call split_box_at_hole(nx, ny, nz, sum_x, sum_y, sum_z, &
                         bbox, newbbox1, newbbox2, min_width, didit)

!!$  If it needs splitting then go back to C++

!!$  write(*,*) "Split box at hole. Did it:", didit
  
  if (didit > 0) then

!!$    Prune the split boxes

!!$    Box 1

    call prune_new_box(nx, ny, nz, mask, bbox, newbbox1)

!!$    Box 2

    call prune_new_box(nx, ny, nz, mask, bbox, newbbox2)

    didit = 2
    return
  end if

!!$  Otherwise there are no holes. Check the density.

!!$  write(*,*) "Finding density"
      
  call find_density(nx, ny, nz, mask, density)
  
!!$  write(*,*) "Found density:", density
  
!!$  If the density is sufficient then this box is done.

  if (density > min_density) then
    didit = 0
    return
  end if

!!$  Otherwise we try and split the box.

!!$  Set up the signature arrays.

!!$  write(*,*) "Setting up signatures"

  call signature(nx, sum_x, sig_x)
  call signature(ny, sum_y, sig_y)
  call signature(nz, sum_z, sig_z)

!!$  write(*,*) "Sums:"
!!$
!!$  do i = ny, 1, -1
!!$    write(*,'(i3)') sum_y(i)
!!$  end do
!!$  
!!$  do i = 1, nx
!!$    write(*,'("   ",i3," ")',advance="no") sum_x(i)
!!$  end do
!!$  write(*,*)
!!$
!!$  write(*,*) "Sigs:"
!!$
!!$  do i = ny, 1, -1
!!$    write(*,'(i3)') sig_y(i)
!!$  end do
!!$  
!!$  do i = 1, nx
!!$    write(*,'("   ",i3," ")',advance="no") sig_x(i)
!!$  end do
!!$  write(*,*)

!!$  write(*,*) "Sum x:"
!!$  write(*,*) sum_x
!!$  write(*,*) "Sig x:"
!!$  write(*,*) sig_x
!!$  write(*,*) "Sum y:"
!!$  write(*,*) sum_y
!!$  write(*,*) "Sig y:"
!!$  write(*,*) sig_y
!!$  write(*,*) "Sum z:"
!!$  write(*,*) sum_z
!!$  write(*,*) "Sig z:"
!!$  write(*,*) sig_z

!!$  Try splitting by the signature

!!$  write(*,*) "Splitting at signatures"

  call split_box_at_sig(nx, ny, nz, sig_x, sig_y, sig_z, &
                        bbox, newbbox1, newbbox2, min_width, &
                        didit)

!!$  This is the last trick up our sleeve. So we just check if
!!$  it succeeded so that we can correct the value of didit.

!!$  write(*,*) "Split box at signature. Did it:", didit

  if (didit > 0) then

!!$    Prune the split boxes

!!$    Box 1

    call prune_new_box(nx, ny, nz, mask, bbox, newbbox1)

!!$    Box 2

    call prune_new_box(nx, ny, nz, mask, bbox, newbbox2)

    didit = 2
  end if

!!$  Then we are done

end subroutine check_box

!!$ Copy the mask.
!!$ Given two masks, one strictly contained within the other, copy
!!$ the intersection to the smaller mask.
!!$ 
!!$ All "s" quantities are the "source" mask - the larger one
!!$ All "d" quantities are the "destination" mask
!!$ The location of the masks in "grid space" is denoted by 
!!$ the bboxes.
!!$ ?bbox(:,1) is the lower boundary
!!$ ?bbox(:,2) is the upper boundary
!!$ ?bbox(:,3) is the stride and is irrelevant

subroutine copy_mask(snx, sny, snz, smask, sbbox, &
                     dnx, dny, dnz, dmask, dbbox)

  implicit none

  CCTK_INT :: snx, sny, snz
  CCTK_INT :: smask(snx, sny, snz)
  CCTK_INT :: sbbox(3, 2)
  CCTK_INT :: dnx, dny, dnz
  CCTK_INT :: dmask(dnx, dny, dnz)
  CCTK_INT :: dbbox(3, 2)
  
  CCTK_INT :: i, j, k, si, sj, sk, di, dj, dk
  
  if ( (dbbox(1,1) < sbbox(1,1)) .or. &
       (dbbox(2,1) < sbbox(2,1)) .or. &
       (dbbox(3,1) < sbbox(3,1)) .or. &
       (dbbox(1,2) > sbbox(1,2)) .or. &
       (dbbox(2,2) > sbbox(2,2)) .or. &
       (dbbox(3,2) > sbbox(3,2)) ) then
    call CCTK_WARN(0, &
        "The destination mask is not contained in the source mask!")
  end if

  do k = dbbox(3,1), dbbox(3,2)
    do j = dbbox(2,1), dbbox(2,2)
      do i = dbbox(1,1), dbbox(1,2)
        
        si = 1 + i - sbbox(1,1)
        sj = 1 + j - sbbox(2,1)
        sk = 1 + k - sbbox(3,1)
        
        di = 1 + i - dbbox(1,1)
        dj = 1 + j - dbbox(2,1)
        dk = 1 + k - dbbox(3,1)
        
        dmask(di, dj, dk) = smask(si, sj, sk)

!!$        write(*,*) "copying mask: loc",i,j,k
!!$        write(*,*) "copying mask:d", di, dj, dk
!!$        write(*,'("(",i2,":",i2,"):(",i2,":",i2,"):(",i2,":",i2,")")') &
!!$             dbbox
!!$        write(*,*) "copying mask:s", si, sj, sk
!!$        write(*,'("(",i2,":",i2,"):(",i2,":",i2,"):(",i2,":",i2,")")') &
!!$             sbbox
        
      end do
    end do
  end do

end subroutine copy_mask

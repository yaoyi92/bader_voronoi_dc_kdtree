!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module for analyzing the charge with a voronoi analysis
!-----------------------------------------------------------------------------------!
MODULE voronoi_mod
  USE kind_mod
  USE matrix_mod
  USE ions_mod
  USE charge_mod
  USE options_mod
  USE io_mod
  use kdtree2_module
  IMPLICIT NONE

! Public, allocatable variables
  TYPE voronoi_obj
    REAL(q2),ALLOCATABLE,DIMENSION(:) :: vorchg
    REAL(q2),ALLOCATABLE,DIMENSION(:) :: ionvol
  END TYPE

  PRIVATE
  PUBLIC :: voronoi_obj
  PUBLIC :: voronoi
  PUBLIC :: voronoi_divide_conquer
  PUBLIC :: voronoi_divide_conquer_kdtree

  CONTAINS

!-----------------------------------------------------------------------------------!
!  voronoi:  Calculate the charge density populations based upon Voronoi polyhedra.
!    In this scheme each element of charge density is associated with the atom that
!    it is closest to.
!-----------------------------------------------------------------------------------!

  SUBROUTINE voronoi(vor,ions,chg)

    TYPE(voronoi_obj) :: vor
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    REAL(q2),DIMENSION(3) :: r_lat, dr_lat, dr_car
    REAL(q2) :: dist, min_dist, vol
    INTEGER :: i, n1, n2, n3, closest, tenths_done, cr, count_max, t1, t2

    CALL SYSTEM_CLOCK(t1, cr, count_max)

    WRITE(*,'(/,2x,A)') 'CALCULATING VORONOI CHARGE DISTRIBUTION'
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'

    ALLOCATE(vor%vorchg(ions%nions))

    vor%vorchg = 0._q2
    tenths_done = 0

    ALLOCATE(vor%ionvol(ions%nions))
    vor%ionvol = 0.0

    DO n1 = 1, chg%npts(1)
      r_lat(1) = REAL(n1,q2)
      IF ((n1*10/chg%npts(1)) > tenths_done) THEN
        tenths_done = (n1*10/chg%npts(1))
        WRITE(*,'(A,$)') '**'
      END IF
      DO n2 = 1, chg%npts(2)
        r_lat(2) = REAL(n2,q2)
        DO n3 = 1, chg%npts(3)
          r_lat(3) = REAL(n3,q2)
          closest = 1
!          r_car = lat2car(chg, r_lat)
!          dr_car = r_car - ions%r_car(1,:)
          dr_lat = r_lat - ions%r_lat(1,:)
          dr_car = MATMUL(chg%lat2car, dr_lat)
          CALL dpbc_car(ions, dr_car)
          min_dist = DOT_PRODUCT(dr_car,dr_car)
          DO i=2,ions%nions
!            dr_car = r_car - ions%r_car(i,:)
            dr_lat = r_lat - ions%r_lat(i,:)
            dr_car = MATMUL(chg%lat2car, dr_lat)
            CALL dpbc_car(ions, dr_car)
            dist = DOT_PRODUCT(dr_car, dr_car)
            IF (dist < min_dist) THEN
              min_dist = dist
              closest = i
            END IF
          END DO
          vor%ionvol(closest) = vor%ionvol(closest) + 1._q2
          vor%vorchg(closest) = vor%vorchg(closest) + rho_val(chg,n1,n2,n3)
        END DO
      END DO
    END DO

    vor%vorchg(:) = vor%vorchg(:)/REAL(chg%nrho,q2)

    vol = matrix_volume(ions%lattice)
    vol = vol/chg%nrho
    DO i = 1, ions%nions
      vor%ionvol(i) = vor%ionvol(i)*vol
    END DO

    CALL SYSTEM_CLOCK(t2, cr, count_max)
    WRITE(*,'(A12,F7.2,A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'

    WRITE(*,'(/,2X,A,/)') 'VORONOI ANALYSIS RESULT'
    WRITE(*,556) '#','X','Y','Z','CHARGE','ATOMIC VOL'
    556 FORMAT(4X,1A1,9X,1A1,2(11X,1A1),8X,1A6,6X,1A10)

    WRITE(*,'(A)') '  ----------------------------------------------------------------------'
    DO i=1,ions%nions
       WRITE(*,777) i,ions%r_car(i,:),vor%vorchg(i),vor%ionvol(i)
       777 FORMAT(1I5,4F12.4,3X,1F12.4)
    END DO
    WRITE(*,'(A)') '  -----------------------------------------------------------------------'

    WRITE(*,'(2x,A,2X,1F12.5)')  '         NUMBER OF ELECTRONS: ', &
    &                                        SUM(vor%vorchg(1:ions%nions))

    
  RETURN
  END SUBROUTINE voronoi

  RECURSIVE SUBROUTINE voronoi_cell(vor,ions,chg,sub_cell)
    TYPE(voronoi_obj) :: vor
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    INTEGER,DIMENSION(6) :: sub_cell

    INTEGER :: x, y, z
    INTEGER :: ix, iy, iz
    INTEGER :: xlo, xhi, ylo, yhi, zlo, zhi
    INTEGER,DIMENSION(6) :: sub_sub_cell
    REAL(q2),DIMENSION(3) :: r_lat, dr_lat, dr_car
    REAL(q2) :: dist, min_dist, vol
    INTEGER,DIMENSION(2,2,2) :: closest
    INTEGER :: i, n1, n2, n3, tenths_done, cr, count_max, t1, t2, cell_size
    LOGICAL :: same_volume
    closest = 0
    !WRITE(*,*) sub_cell(1), sub_cell(2), sub_cell(3), sub_cell(4), sub_cell(5), sub_cell(6)
    DO ix = 1,2
      DO iy = 1,2
        DO iz = 1,2
          n1 = sub_cell(ix)
          n2 = sub_cell(iy+2)
          n3 = sub_cell(iz+4)
          r_lat(1) = REAL(n1,q2)
          r_lat(2) = REAL(n2,q2)
          r_lat(3) = REAL(n3,q2)
          min_dist = 9.9E9
          IF (ix == 2 .AND. sub_cell(2) == sub_cell(1)) THEN
            closest(ix,iy,iz) = closest(1,iy,iz)
          ELSE IF (iy == 2 .AND. sub_cell(3) == sub_cell(4)) THEN
            closest(ix,iy,iz) = closest(ix,1,iz)
          ELSE IF (iz == 2 .AND. sub_cell(5) == sub_cell(6)) THEN
            closest(ix,iy,iz) = closest(ix,iy,1)
          ELSE
            DO i=1,ions%nions
              dr_lat = r_lat - ions%r_lat(i,:)
              !WRITE(*,*) ions%r_lat(i,:)
              dr_car = MATMUL(chg%lat2car, dr_lat)
              CALL dpbc_car(ions, dr_car)
              dist = DOT_PRODUCT(dr_car,dr_car)
              !WRITE(*,*) n1, n2, n3, r_lat(1), r_lat(2), r_lat(3),i, dist(ix,iy,iz)
              IF (dist < min_dist) THEN
                min_dist = dist
                closest(ix,iy,iz) = i
                !closest(ix,iy,iz) = 1
              END IF
            END DO
          END IF
          !WRITE(*,*) closest(ix,iy,iz)
        END DO
      END DO
    END DO
    same_volume = .True.
    DO ix = 1,2
      DO iy = 1,2
        DO iz = 1,2
          IF (closest(ix,iy,iz) /= closest(1,1,1)) THEN
            same_volume = .False.
          END IF
        END DO
      END DO
    END DO
    !WRITE(*,*) same_volume

    cell_size = sub_cell(2) - sub_cell(1) + sub_cell(4) - sub_cell(3) + sub_cell(6) - sub_cell(5)
    !WRITE(*,*) "cell_size",cell_size
    IF (same_volume) THEN
    !IF (cell_size == 0) THEN

      DO n1 = sub_cell(1), sub_cell(2)
        DO n2 = sub_cell(3), sub_cell(4)
          DO n3 = sub_cell(5), sub_cell(6)
            !WRITE(*,*) "hello"
            vor%ionvol(closest(1,1,1)) = vor%ionvol(closest(1,1,1)) + 1._q2
            vor%vorchg(closest(1,1,1)) = vor%vorchg(closest(1,1,1)) + rho_val(chg,n1,n2,n3)
            !WRITE(*,*) n1,n2,n3
          END DO
        END DO
      END DO

    ELSE

    DO ix = 1,2
      IF (ix == 1) THEN
        xlo = sub_cell(1)
        xhi = sub_cell(1) + ( sub_cell(2) - sub_cell(1)) / 2
      ELSE
        xlo = sub_cell(1) + ( sub_cell(2) - sub_cell(1)) / 2 + 1
        xhi = sub_cell(2) 
      END IF
      DO iy = 1,2
        IF (iy == 1) THEN
          ylo = sub_cell(3)
          yhi = sub_cell(3) + ( sub_cell(4) - sub_cell(3)) / 2
        ELSE
          ylo = sub_cell(3) + ( sub_cell(4) - sub_cell(3)) / 2 + 1
          yhi = sub_cell(4) 
        END IF
        DO iz = 1,2
          IF (iz == 1) THEN
            zlo = sub_cell(5)
            zhi = sub_cell(5) + ( sub_cell(6) - sub_cell(5)) / 2
          ELSE
            zlo = sub_cell(5) + ( sub_cell(6) - sub_cell(5)) / 2 + 1
            zhi = sub_cell(6) 
          END IF
          sub_sub_cell(1) = xlo
          sub_sub_cell(2) = xhi
          sub_sub_cell(3) = ylo
          sub_sub_cell(4) = yhi
          sub_sub_cell(5) = zlo
          sub_sub_cell(6) = zhi
          IF (.NOT. (( ix == 2 .AND. sub_cell(1) == sub_cell(2) ) .OR. &
                     ( iy == 2 .AND. sub_cell(3) == sub_cell(4) ) .OR. &
                     ( iz == 2 .AND. sub_cell(5) == sub_cell(6) ))) THEN
            !WRITE(*,*) 'origin',sub_cell(1), sub_cell(2), sub_cell(3), sub_cell(4), sub_cell(5), sub_cell(6)
            !WRITE(*,*) xlo, xhi, ylo, yhi, zlo, zhi
            CALL voronoi_cell(vor,ions,chg,sub_sub_cell)
          END IF
        END DO
      END DO
    END DO

    END IF

  END SUBROUTINE voronoi_cell

  RECURSIVE SUBROUTINE voronoi_cell_kdtree(vor,ions,chg,sub_cell,tree)
    TYPE(voronoi_obj) :: vor
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg
    type(kdtree2), pointer :: tree

    INTEGER,DIMENSION(6) :: sub_cell

    INTEGER :: x, y, z
    INTEGER :: ix, iy, iz
    INTEGER :: xlo, xhi, ylo, yhi, zlo, zhi
    INTEGER,DIMENSION(6) :: sub_sub_cell
    REAL(q2),DIMENSION(3) :: r_lat, dr_lat, r_car
    REAL(q2) :: dist, min_dist, vol
    INTEGER,DIMENSION(2,2,2) :: closest
    INTEGER :: i, n1, n2, n3, tenths_done, cr, count_max, t1, t2, cell_size
    LOGICAL :: same_volume

    type(kdtree2_result),allocatable :: results_kdtree(:)
    allocate(results_kdtree(1))

    closest = 0
    !WRITE(*,*) sub_cell(1), sub_cell(2), sub_cell(3), sub_cell(4), sub_cell(5), sub_cell(6)
    DO ix = 1,2
      DO iy = 1,2
        DO iz = 1,2
          n1 = sub_cell(ix)
          n2 = sub_cell(iy+2)
          n3 = sub_cell(iz+4)
          r_lat(1) = REAL(n1,q2)
          r_lat(2) = REAL(n2,q2)
          r_lat(3) = REAL(n3,q2)
          min_dist = 9.9E9
          IF (ix == 2 .AND. sub_cell(2) == sub_cell(1)) THEN
            closest(ix,iy,iz) = closest(1,iy,iz)
          ELSE IF (iy == 2 .AND. sub_cell(3) == sub_cell(4)) THEN
            closest(ix,iy,iz) = closest(ix,1,iz)
          ELSE IF (iz == 2 .AND. sub_cell(5) == sub_cell(6)) THEN
            closest(ix,iy,iz) = closest(ix,iy,1)
          ELSE
            r_car = MATMUL(chg%lat2car, r_lat) * 1000.0 ! 1000.0 to reduce roundoff error
            call kdtree2_n_nearest(tp=tree,qv=r_car,nn=1,results=results_kdtree) 
            !WRITE(*,*) results_kdtree%idx
            !WRITE(*,*) ions%nions
            !WRITE(*,*) MOD(results_kdtree%idx-1, ions%nions) + 1
            closest(ix,iy,iz) = MOD(results_kdtree(1)%idx-1, ions%nions) + 1
            !i = MOD(results_kdtree%idx+1, ions%nions) + 1
            !DO i=1,ions%nions
            !  dr_lat = r_lat - ions%r_lat(i,:)
            !  !WRITE(*,*) ions%r_lat(i,:)
            !  dr_car = MATMUL(chg%lat2car, dr_lat)
            !  CALL dpbc_car(ions, dr_car)
            !  dist = DOT_PRODUCT(dr_car,dr_car)
            !  !WRITE(*,*) n1, n2, n3, r_lat(1), r_lat(2), r_lat(3),i, dist(ix,iy,iz)
            !  IF (dist < min_dist) THEN
            !    min_dist = dist
            !    closest(ix,iy,iz) = i
            !    !closest(ix,iy,iz) = 1
            !  END IF
            !END DO
            !WRITE(*,*) closest(ix,iy,iz)
          END IF
          !WRITE(*,*) closest(ix,iy,iz)
        END DO
      END DO
    END DO
    same_volume = .True.
    DO ix = 1,2
      DO iy = 1,2
        DO iz = 1,2
          IF (closest(ix,iy,iz) /= closest(1,1,1)) THEN
            same_volume = .False.
          END IF
        END DO
      END DO
    END DO
    !WRITE(*,*) same_volume

    cell_size = sub_cell(2) - sub_cell(1) + sub_cell(4) - sub_cell(3) + sub_cell(6) - sub_cell(5)
    !WRITE(*,*) "cell_size",cell_size
    IF (same_volume) THEN
    !IF (cell_size == 0) THEN

      DO n1 = sub_cell(1), sub_cell(2)
        DO n2 = sub_cell(3), sub_cell(4)
          DO n3 = sub_cell(5), sub_cell(6)
            !WRITE(*,*) "hello"
            vor%ionvol(closest(1,1,1)) = vor%ionvol(closest(1,1,1)) + 1._q2
            vor%vorchg(closest(1,1,1)) = vor%vorchg(closest(1,1,1)) + rho_val(chg,n1,n2,n3)
            !WRITE(*,*) n1,n2,n3
          END DO
        END DO
      END DO

    ELSE

    DO ix = 1,2
      IF (ix == 1) THEN
        xlo = sub_cell(1)
        xhi = sub_cell(1) + ( sub_cell(2) - sub_cell(1)) / 2
      ELSE
        xlo = sub_cell(1) + ( sub_cell(2) - sub_cell(1)) / 2 + 1
        xhi = sub_cell(2) 
      END IF
      DO iy = 1,2
        IF (iy == 1) THEN
          ylo = sub_cell(3)
          yhi = sub_cell(3) + ( sub_cell(4) - sub_cell(3)) / 2
        ELSE
          ylo = sub_cell(3) + ( sub_cell(4) - sub_cell(3)) / 2 + 1
          yhi = sub_cell(4) 
        END IF
        DO iz = 1,2
          IF (iz == 1) THEN
            zlo = sub_cell(5)
            zhi = sub_cell(5) + ( sub_cell(6) - sub_cell(5)) / 2
          ELSE
            zlo = sub_cell(5) + ( sub_cell(6) - sub_cell(5)) / 2 + 1
            zhi = sub_cell(6) 
          END IF
          sub_sub_cell(1) = xlo
          sub_sub_cell(2) = xhi
          sub_sub_cell(3) = ylo
          sub_sub_cell(4) = yhi
          sub_sub_cell(5) = zlo
          sub_sub_cell(6) = zhi
          IF (.NOT. (( ix == 2 .AND. sub_cell(1) == sub_cell(2) ) .OR. &
                     ( iy == 2 .AND. sub_cell(3) == sub_cell(4) ) .OR. &
                     ( iz == 2 .AND. sub_cell(5) == sub_cell(6) ))) THEN
            !WRITE(*,*) 'origin',sub_cell(1), sub_cell(2), sub_cell(3), sub_cell(4), sub_cell(5), sub_cell(6)
            !WRITE(*,*) xlo, xhi, ylo, yhi, zlo, zhi
            CALL voronoi_cell_kdtree(vor,ions,chg,sub_sub_cell,tree)
          END IF
        END DO
      END DO
    END DO

    END IF
    deallocate(results_kdtree)

  END SUBROUTINE voronoi_cell_kdtree

!-----------------------------------------------------------------------------------!
  SUBROUTINE voronoi_divide_conquer(vor,ions,chg)

    TYPE(voronoi_obj) :: vor
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    REAL(q2),DIMENSION(3) :: r_lat, dr_lat, dr_car
    REAL(q2) :: dist, min_dist, vol
    INTEGER :: i, n1, n2, n3, closest, tenths_done, cr, count_max, t1, t2

    INTEGER,DIMENSION(6) :: cell, sub_cell
    INTEGER :: xlo, xhi, ylo, yhi, zlo, zhi, ix, iy, iz

    CALL SYSTEM_CLOCK(t1, cr, count_max)

    WRITE(*,'(/,2x,A)') 'CALCULATING VORONOI CHARGE DISTRIBUTION'
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'

    ALLOCATE(vor%vorchg(ions%nions))

    vor%vorchg = 0._q2
    tenths_done = 0

    ALLOCATE(vor%ionvol(ions%nions))
    vor%ionvol = 0.0
    cell(1) = 1
    cell(2) = chg%npts(1)
    cell(3) = 1
    cell(4) = chg%npts(2)
    cell(5) = 1
    cell(6) = chg%npts(3)

    DO ix = 1,2
      IF (ix == 1) THEN
        xlo = cell(1)
        xhi = cell(1) + ( cell(2) - cell(1)) / 2
      ELSE
        xlo = cell(1) + ( cell(2) - cell(1)) / 2 + 1
        xhi = cell(2) 
        IF (xlo > xhi) THEN
          xlo = xhi
        END IF
      END IF
      DO iy = 1,2
        IF (iy == 1) THEN
          ylo = cell(3)
          yhi = cell(3) + ( cell(4) - cell(3)) / 2
        ELSE
          ylo = cell(3) + ( cell(4) - cell(3)) / 2 + 1
          yhi = cell(4) 
          IF (ylo > yhi) THEN
            ylo = yhi
          END IF
        END IF
        DO iz = 1,2
          IF (iz == 1) THEN
            zlo = cell(5)
            zhi = cell(5) + ( cell(6) - cell(5)) / 2
          ELSE
            zlo = cell(5) + ( cell(6) - cell(5)) / 2 + 1
            zhi = cell(6) 
            IF (zlo > zhi) THEN
              zlo = zhi
            END IF
          END IF
          WRITE(*,'(A,$)') '**'
          sub_cell(1) = xlo
          sub_cell(2) = xhi
          sub_cell(3) = ylo
          sub_cell(4) = yhi
          sub_cell(5) = zlo
          sub_cell(6) = zhi
          !WRITE(*,*) 'origin',cell(1), cell(2), cell(3), cell(4), cell(5), cell(6)
          !WRITE(*,*) xlo, xhi, ylo, yhi, zlo, zhi
          CALL voronoi_cell(vor,ions,chg,sub_cell)
        END DO
      END DO
    END DO


    vor%vorchg(:) = vor%vorchg(:)/REAL(chg%nrho,q2)

    vol = matrix_volume(ions%lattice)
    vol = vol/chg%nrho
    DO i = 1, ions%nions
      vor%ionvol(i) = vor%ionvol(i)*vol
    END DO

    CALL SYSTEM_CLOCK(t2, cr, count_max)
    WRITE(*,'(A12,F7.2,A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'


    !WRITE(*,'(/,2X,A,/)') 'VORONOI ANALYSIS RESULT'
    WRITE(*,*) 'VORONOI ANALYSIS RESULT'
    WRITE(*,556) '#','X','Y','Z','CHARGE','ATOMIC VOL'
    556 FORMAT(4X,1A1,9X,1A1,2(11X,1A1),8X,1A6,6X,1A10)

    WRITE(*,'(A)') '  ----------------------------------------------------------------------'
    DO i=1,ions%nions
       WRITE(*,777) i,ions%r_car(i,:),vor%vorchg(i),vor%ionvol(i)
       777 FORMAT(1I5,4F12.4,3X,1F12.4)
    END DO
    WRITE(*,'(A)') '  -----------------------------------------------------------------------'

    WRITE(*,'(2x,A,2X,1F12.5)')  '         NUMBER OF ELECTRONS: ', &
    &                                        SUM(vor%vorchg(1:ions%nions))

    
  RETURN
  END SUBROUTINE voronoi_divide_conquer

  SUBROUTINE voronoi_divide_conquer_kdtree(vor,ions,chg)

    TYPE(voronoi_obj) :: vor
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chg

    REAL(q2),DIMENSION(3) :: r_lat, ir_lat, ir_car
    REAL(q2) :: dist, min_dist, vol
    INTEGER :: i, n1, n2, n3, closest, tenths_done, cr, count_max, t1, t2

    INTEGER,DIMENSION(6) :: cell, sub_cell
    INTEGER :: xlo, xhi, ylo, yhi, zlo, zhi, ix, iy, iz
    type(kdtree2), pointer :: tree
    INTEGER :: ii, jj, kk ! for periodic index
    REAL(q2),DIMENSION(3,ions%nions*27) :: extended_coords
    
    extended_coords = 0.0

    CALL SYSTEM_CLOCK(t1, cr, count_max)
    !! generate kdtree
    DO ii = -1,1
      DO jj = -1,1
        DO kk = -1,1
          DO i = 1, ions%nions
            ir_lat = (/ ions%r_lat(i,1) + chg%npts(1) * ii, &
                        ions%r_lat(i,2) + chg%npts(2) * jj, &
                        ions%r_lat(i,3) + chg%npts(3) * kk /)
            ir_car = MATMUL(chg%lat2car, ir_lat) * 1000.0  ! 1000.0 to reduce roundoff error
            extended_coords(1,i+((ii+1)*9+(jj+1)*3+(kk+1))*ions%nions) = ir_car(1)
            extended_coords(2,i+((ii+1)*9+(jj+1)*3+(kk+1))*ions%nions) = ir_car(2)
            extended_coords(3,i+((ii+1)*9+(jj+1)*3+(kk+1))*ions%nions) = ir_car(3)
          END DO
        END DO
      END DO
    END DO
   
    !WRITE(*,*) ions%r_lat(:,:)
    tree => kdtree2_create(extended_coords,sort=.false.,rearrange=.false.)  ! this is how you create a tree. 
    !WRITE (*,*) extended_coords

    WRITE(*,'(/,2x,A)') 'CALCULATING VORONOI CHARGE DISTRIBUTION'
    WRITE(*,'(2x,A)')   '               0  10  25  50  75  100'
    WRITE(*,'(2x,A,$)') 'PERCENT DONE:  **'

    ALLOCATE(vor%vorchg(ions%nions))

    vor%vorchg = 0._q2
    tenths_done = 0

    ALLOCATE(vor%ionvol(ions%nions))
    vor%ionvol = 0.0
    cell(1) = 1
    cell(2) = chg%npts(1)
    cell(3) = 1
    cell(4) = chg%npts(2)
    cell(5) = 1
    cell(6) = chg%npts(3)

    DO ix = 1,2
      IF (ix == 1) THEN
        xlo = cell(1)
        xhi = cell(1) + ( cell(2) - cell(1)) / 2
      ELSE
        xlo = cell(1) + ( cell(2) - cell(1)) / 2 + 1
        xhi = cell(2) 
        IF (xlo > xhi) THEN
          xlo = xhi
        END IF
      END IF
      DO iy = 1,2
        IF (iy == 1) THEN
          ylo = cell(3)
          yhi = cell(3) + ( cell(4) - cell(3)) / 2
        ELSE
          ylo = cell(3) + ( cell(4) - cell(3)) / 2 + 1
          yhi = cell(4) 
          IF (ylo > yhi) THEN
            ylo = yhi
          END IF
        END IF
        DO iz = 1,2
          IF (iz == 1) THEN
            zlo = cell(5)
            zhi = cell(5) + ( cell(6) - cell(5)) / 2
          ELSE
            zlo = cell(5) + ( cell(6) - cell(5)) / 2 + 1
            zhi = cell(6) 
            IF (zlo > zhi) THEN
              zlo = zhi
            END IF
          END IF
          WRITE(*,'(A,$)') '**'
          sub_cell(1) = xlo
          sub_cell(2) = xhi
          sub_cell(3) = ylo
          sub_cell(4) = yhi
          sub_cell(5) = zlo
          sub_cell(6) = zhi
          !WRITE(*,*) 'origin',cell(1), cell(2), cell(3), cell(4), cell(5), cell(6)
          !WRITE(*,*) xlo, xhi, ylo, yhi, zlo, zhi
          CALL voronoi_cell_kdtree(vor,ions,chg,sub_cell,tree)
        END DO
      END DO
    END DO


    vor%vorchg(:) = vor%vorchg(:)/REAL(chg%nrho,q2)

    vol = matrix_volume(ions%lattice)
    vol = vol/chg%nrho
    DO i = 1, ions%nions
      vor%ionvol(i) = vor%ionvol(i)*vol
    END DO


    call kdtree2_destroy(tree)  

    CALL SYSTEM_CLOCK(t2, cr, count_max)
    WRITE(*,'(A12,F7.2,A8)') 'RUN TIME: ',(t2-t1)/REAL(cr,q2),' SECONDS'


    !WRITE(*,'(/,2X,A,/)') 'VORONOI ANALYSIS RESULT'
    WRITE(*,*) 'VORONOI ANALYSIS RESULT'
    WRITE(*,556) '#','X','Y','Z','CHARGE','ATOMIC VOL'
    556 FORMAT(4X,1A1,9X,1A1,2(11X,1A1),8X,1A6,6X,1A10)

    WRITE(*,'(A)') '  ----------------------------------------------------------------------'
    DO i=1,ions%nions
       WRITE(*,777) i,ions%r_car(i,:),vor%vorchg(i),vor%ionvol(i)
       777 FORMAT(1I5,4F12.4,3X,1F12.4)
    END DO
    WRITE(*,'(A)') '  -----------------------------------------------------------------------'

    WRITE(*,'(2x,A,2X,1F12.5)')  '         NUMBER OF ELECTRONS: ', &
    &                                        SUM(vor%vorchg(1:ions%nions))

    
  RETURN
  END SUBROUTINE voronoi_divide_conquer_kdtree

END MODULE voronoi_mod

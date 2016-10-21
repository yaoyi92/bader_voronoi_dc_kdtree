  MODULE critpoints_mod
    USE kind_mod
    USE matrix_mod
    USE bader_mod
    USE charge_mod 
    USE options_mod
    USE ions_mod
    USE io_mod
    IMPLICIT NONE

    PRIVATE 
    PUBLIC :: critpoint_find

    CONTAINS

!-----------------------------------------------------------------------------------!
!critpoint_find: find critical points on the edge of the Bader volumes
!NOTE: this subroutine should be called after refine_edge
!      in order to restrict the calculation to edge points
!-----------------------------------------------------------------------------------!
  SUBROUTINE critpoint_find(bdr,chg,opts)

    TYPE hessian
      REAL(q2),ALLOCATABLE,DIMENSION(:,:,:) :: rho, dx, dy, dz
      REAL(q2) :: dxdx, dydy, dzdz, dxdy, dxdz, dydz
      REAL(q2) :: r1, r2, r3, eigval1, eigval2, eigval3
      ! eigval and eigvec are eigenvalues and eigvectors of hessian matrix
    END TYPE

    TYPE(hessian) :: hes
    TYPE(bader_obj) :: bdr
    TYPE(charge_obj) :: chg
    TYPE(options_obj) :: opts

! for points, 1 and 2 are +1, -1
    INTEGER,DIMENSION(3) :: p, pt, ptt, ptx1, ptx2, pty1, pty2, ptz1, ptz2
    REAL(q2),DIMENSION(3,3) :: dM ! the deviatoric matrix
    REAL(q2),DIMENSION(3,3) :: dMSQ, dMVSS ! dM squared and dM in vss basis
    REAL(q2) :: j2, j3 ! second and third invariant of dM
    INTEGER :: n1, n2, n3, d1, d2, d3, cptnum, switch
    REAL(q2),DIMENSION(3) :: eigvec1, eigvec2, eigvec3, cartX, cartY, cartZ, tempVec
    REAL(q2),DIMENSION(3,3) :: nIdentity, identityM, devSubNIden
    ! these are vectors orthogonal to eigenvectors
    REAL(q2),DIMENSION(3) :: orthoR1, orthoR2, orthoR3, S1, S2, orT2, orT3
    REAL(q2) :: x, y, z, xx, yy, zz, xy, xz, yz, denom, nomx, nomy, nomz
    REAL(q2) :: norm, trace, traceOver3, temp, alpha
    REAL(q2) :: yita1, yita2, yita3 ! variables for degenerate eigenvalues

    WRITE(*,'(A)')  'FINDING CRITICAL POINTS'
    ALLOCATE (hes%rho(chg%npts(1),chg%npts(2),chg%npts(3)))
    ALLOCATE (hes%dx(chg%npts(1),chg%npts(2),chg%npts(3)))
    ALLOCATE (hes%dy(chg%npts(1),chg%npts(2),chg%npts(3)))
    ALLOCATE (hes%dz(chg%npts(1),chg%npts(2),chg%npts(3)))
    switch = 0
    cptnum = 0
    identityM = 0._q2
    identityM(1,1) = 1._q2
    identityM(2,2) = 1._q2
    identityM(3,3) = 1._q2
    cartX = (/1._q2,0._q2,0._q2/)
    cartY = (/0._q2,1._q2,0._q2/)
    cartZ = (/0._q2,0._q2,1._q2/)

    OPEN(97,FILE='CPF.dat',STATUS='REPLACE',ACTION='WRITE')
    DO n1 = 1,chg%npts(1)
      DO n2 = 1,chg%npts(2)
        DO n3 = 1,chg%npts(3)
            p = (/n1,n2,n3/)
            IF (rho_val(chg,p(1),p(2),p(3)) <= opts%vacval) CYCLE
!            IF (is_vol_edge(bdr,chg,p) .AND. (.NOT.is_max(chg,p))) THEN

!-----------------------------------------------------------------------------------!
! now that this subroutine can find the correct amount of edge points, lets have
! it find the hessian
!-----------------------------------------------------------------------------------!

             ! this loop finds all the neighboring points that will be needed and stores their rho
             DO d1=-2,2
               DO d2=-2,2
                 DO d3=-2,2
                   pt = p + (/d1,d2,d3/)
                   CALL pbc(pt,chg%npts)
                   ptt = p + (/d1,d2,d3/)
                   CALL pbc(ptt,chg%npts)
                   hes%rho(pt(1),pt(2),pt(3)) = rho_val(chg,pt(1),pt(2),pt(3))
                 END DO
               END DO
             END DO

! now there is a problem.  ptt(1)+1 may get out of boundary. How to do that
! boundary check? --- by creating points for every point that could possibily be
! needed and run pbc on it. Freezing is because of points out of boundary.
! However, the compiler flags seems to not catch this problem.

             DO d1=-1,1
               DO d2=-1,1
                 DO d3=-1,1
                   ptt = p + (/d1,d2,d3/)
                   CALL pbc(ptt,chg%npts)
                   ptx1 = ptt + (/1,0,0/)
                   ptx2 = ptt + (/-1,0,0/)
                   pty1 = ptt + (/0,1,0/)
                   pty2 = ptt + (/0,-1,0/)
                   ptz1 = ptt + (/0,0,1/)
                   ptz2 = ptt + (/0,0,-1/)
                   CALL pbc(ptx1,chg%npts)
                   CALL pbc(ptx2,chg%npts)
                   CALL pbc(pty1,chg%npts)
                   CALL pbc(pty2,chg%npts)
                   CALL pbc(ptz1,chg%npts)
                   CALL pbc(ptz2,chg%npts)
                   hes%dx(ptt(1),ptt(2),ptt(3)) = 0.5_q2 * &
                         (hes%rho(ptx1(1),ptt(2),ptt(3)) - &
                          hes%rho(ptx2(1),ptt(2),ptt(3)))
                   hes%dy(ptt(1),ptt(2),ptt(3)) = 0.5_q2 * &
                         (hes%rho(ptt(1),pty1(2),ptt(3)) - &
                          hes%rho(ptt(1),pty2(2),ptt(3)))
                   hes%dz(ptt(1),ptt(2),ptt(3)) = 0.5_q2 * &
                         (hes%rho(ptx1(1),ptt(2),ptz1(3)) - &
                          hes%rho(ptt(1),ptt(2),ptz2(3)))
                 END DO
               END DO
             END DO
             ptx1 = p + (/1,0,0/)
             ptx2 = p + (/-1,0,0/)
             pty1 = p + (/0,1,0/)
             pty2 = p + (/0,-1,0/)
             ptz1 = p + (/0,0,1/)
             ptz2 = p + (/0,0,-1/)
             CALL pbc(ptx1,chg%npts)
             CALL pbc(ptx2,chg%npts)
             CALL pbc(pty1,chg%npts)
             CALL pbc(pty2,chg%npts)
             CALL pbc(ptz1,chg%npts)
             CALL pbc(ptz2,chg%npts)

             hes%dx(p(1),p(2),p(3)) = 0.5_q2*(hes%rho(ptx1(1),p(2),p(3)) - &
                                              hes%rho(ptx2(1),p(2),p(3)))
             hes%dy(p(1),p(2),p(3)) = 0.5_q2*(hes%rho(p(1),pty1(2),p(3)) - &
                                              hes%rho(p(1),pty2(2),p(3)))
             hes%dz(p(1),p(2),p(3)) = 0.5_q2*(hes%rho(p(1),p(2),ptz1(3)) - &
                                              hes%rho(p(1),p(2),ptz2(3)))

             ! now calculate the curvature of rho at the point p
             hes%dxdx = 0.5_q2 * (hes%dx(ptx1(1),p(2),p(3)) - hes%dx(ptx2(1),p(2),p(3)))
             hes%dydy = 0.5_q2 * (hes%dy(p(1),pty1(2),p(3)) - hes%dy(p(1),pty2(2),p(3)))
             hes%dzdz = 0.5_q2 * (hes%dz(p(1),p(2),ptz1(3)) - hes%dz(p(1),p(2),ptz2(3)))
             hes%dxdy = 0.25_q2 * (hes%rho(ptx1(1),pty1(2),p(3)) - hes%rho(ptx2(1),pty1(2),p(3)) - &
                                   hes%rho(ptx1(1),pty2(2),p(3)) + hes%rho(ptx2(1),pty2(2),p(3)))
             hes%dxdz = 0.25_q2 * (hes%rho(ptx1(1),p(2),ptz1(3)) - hes%rho(ptx2(1),p(2),ptz1(3)) - &
                                   hes%rho(ptx1(1),p(2),ptz2(3)) + hes%rho(ptx2(1),p(2),ptz2(3)))
             hes%dydz = 0.25_q2 * (hes%rho(p(1),pty1(2),ptz1(3)) - hes%rho(p(1),pty2(2),ptz1(3)) - &
                                   hes%rho(p(1),pty1(2),ptz2(3)) + hes%rho(p(1),pty2(2),ptz2(3)))
!solutions for the vector
!x:
!(-dxdz dy dydz + dx dydz^2 + dxdz dydy dz -  dxdy dydz dz + dxdy dy dzdz - dx
!dydy dzdz)/(dxdz^2 dydy - 2 dxdy dxdz dydz + dxdx dydz^2 + dxdy^2 dzdz -  dxdx
!dydy dzdz)
!y:
!(dxdz^2 dy - dx dxdz dydz - dxdy dxdz dz + dxdx dydz dz +  dx dxdy dzdz - dxdx
!dy dzdz)/(dxdz^2 dydy - 2 dxdy dxdz dydz + dxdx dydz^2 + dxdy^2 dzdz -  dxdx
!dydy dzdz)
!z:
!(-dxdy dxdz dy + dx dxdz dydy - dx dxdy dydz + dxdx dy dydz +  dxdy^2 dz - dxdx
!dydy dz)/(dxdz^2 dydy - 2 dxdy dxdz dydz + dxdx dydz^2 + dxdy^2 dzdz -  dxdx
!dydy dzdz)
               x = hes%dx(p(1),p(2),p(3))
               y = hes%dy(p(1),p(2),p(3))
               z = hes%dz(p(1),p(2),p(3))
               xx = hes%dxdx
               yy = hes%dydy
               zz = hes%dzdz
               xy = hes%dxdy
               xz = hes%dxdz
               yz = hes%dydz
               denom = xz*xz*yy - 2.0_q2*xy*xz*yz + xx*yz*yz + xy*xy*zz - xx*yy*zz
               nomx = -xz*y*yz + x*yz*yz + xz*yy*z - xy*yz*z + xy*y*zz - x*yy*zz
               nomy = xz*xz*y - x*xz*yz - xy*xz*z + xx*yz*z + x*xy*zz - xx*y*zz
               nomz = -xy*xz*y + x*xz*yy - x*xy*yz + xx*y*yz + xy*xy*z - xx*yy*z
               hes%r1 = nomx/denom
               hes%r2 = nomy/denom
               hes%r3 = nomz/denom

               IF (ABS(hes%r1) <= 0.5_q2*bdr%stepsize) THEN
                 IF (ABS(hes%r2) <= 0.5_q2*bdr%stepsize) THEN
                   IF (ABS(hes%r3) <= 0.5_q2*bdr%stepsize) THEN
                     cptnum = cptnum + 1
                     ! eigenvalue calculation below
                     trace = hes%dxdx + hes%dydy + hes%dzdz
                     traceOver3 = (hes%dxdx + hes%dydy + hes%dzdz)/3.0_q2
                     dM(1,1) = hes%dxdx - traceOver3
                     dM(2,2) = hes%dydy - traceOver3
                     dM(3,3) = hes%dzdz - traceOver3
                     dM(1,2) = hes%dxdy
                     dM(1,3) = hes%dxdz
                     dM(2,1) = hes%dxdy
                     dM(2,3) = hes%dydz
                     dM(3,1) = hes%dxdz
                     dM(3,2) = hes%dydz
                     ! Now a loop to calculate dMSQ
                     temp = 0._q2
                     ! this loop is functional
                     dMSQ = MATMUL(dM,dM)
                     ! j2 is 1/2 tr(dM dot dM)
                     j2 = 0.5_q2*(dMSQ(1,1) + dMSQ(2,2) + dMSQ(3,3))
                     j3 = determinant(dM)
                     alpha = ACOS(j3/2._q2*(3._q2/j2)**(3._q2/2._q2))/3._q2
                     yita1 = 2._q2*SQRT(j2/3._q2)*COS(alpha)
                     ! Find the most distinct eigenvalue first
                     ! every eigenvalue are equaly distinct when alpha = pi/6
                     ! can use the code when alpha < pi/6
                     IF (alpha>pi/6._q2) THEN
                       ! Start with eigval 3 first
                       temp = hes%eigval1
                       hes%eigval1 = hes%eigval3
                       hes%eigval3 = temp
                       temp = 0
                     END IF
                     ! this is the default when yita1 is the most distinct
                     nIdentity = identityM*yita1
                     devSubNIden = dM - nIdentity
                     orthoR1 = MATMUL(devSubNIden,cartX)
                     orthoR2 = MATMUL(devSubNIden,cartY)
                     orthoR3 = MATMUL(devSubNIden,cartZ)
                     ! assume r1 is the largest, normalize all orthoR vectors
                     norm = 1._q2/SQRT(SUM(orthoR1(:)**2))
                     S1 = orthoR1*norm
                     tempVec = S1*DOT_PRODUCT(S1, orthoR2)
                     orT2 = orthoR2 - tempVec
                     tempVec = S1*DOT_PRODUCT(S1, orthoR3)
                     orT3 = orthoR3 - tempVec
                     ! assume t2 is the larger one, normalize it
                     norm = 1.0_q2/SQRT(SUM(orT2(:)**2))
                     s2 = orT2*norm
                     eigvec1 = cross_product(S1, S2)
                     ! now that I have an eigenvalue and an eigenvector,
                     ! write the deviatoric matrix using the v1,s1,s2 basis
                     ! hes%eigval1 |         0        |      0
                     ! 0           | s1 dot dM dot s1 | s1 dot dM dot s2
                     ! 0           | s2 dot dM dot s1 | s2 dot dM dot s2
                     dMVSS = 0._q2
                     dMVSS(1,1) = yita1
                     tempVec = MATMUL(dM,S1)
                     dMVSS(2,2) = DOT_PRODUCT(tempVec, S1)
                     dMVSS(2,3) = DOT_PRODUCT(tempVec, S2)
                     tempVec = MATMUL(dM,S2)
                     dMVSS(3,2) = DOT_PRODUCT(tempVec, S1)
                     dMVSS(3,3) = DOT_PRODUCT(tempVec, S2)
                     ! This is a sign function
                     temp = SIGN(1._q2,dMVSS(2,2)-dMVSS(3,3))
                     IF (dMVSS(2,2) == dMVSS(3,3)) temp = 0
                     yita2 = (dMVSS(2,2) + dMVSS(3,3))/2._q2 - 1._q2/2._q2*temp*&
                             SQRT((dMVSS(2,2) - dMVSS(3,3))**2 + 4._q2*dMVSS(2,3)*dMVSS(3,2))
                     yita3 = dMVSS(2,2) + dMVSS(3,3) - yita2
                     ! these eigenvalues are shifted by traceOver3
                     hes%eigval1 = yita1 + traceOver3
                     hes%eigval2 = yita2 + traceOver3
                     hes%eigval3 = yita3 + traceOver3

                     CALL eigenvectors(yita2, identityM, dM, s1, s2, eigvec1, eigvec2, eigvec3)

                     WRITE(97,*) '*********** A NEW ENTRY *************'
                     WRITE(97,*) 'Critical point number: ', cptnum
                     WRITE(97,*) p(1),p(2),p(3)
                     WRITE(97,*) 'Threshold: ', 0.5_q2*bdr%stepsize
                     WRITE(97,*) 'r1' , hes%r1
                     WRITE(97,*) 'r2' , hes%r2
                     WRITE(97,*) 'r3' , hes%r3
                     WRITE(97,'(3(1X,E18.11))') 0._q2, hes%rho(ptz1(1),ptz1(2),ptz1(3)), hes%rho(ptx2(1),ptx2(2),ptx2(3))
                     WRITE(97,'(3(1X,E18.11))') hes%rho(pty2(1),pty2(2),pty2(3)), hes%rho(p(1),p(2),p(3)), &
                                                hes%rho(pty1(1),pty1(2),pty1(3))
                     WRITE(97,'(3(1X,E18.11))') hes%rho(ptx1(1), ptx1(2),ptx1(3)), hes%rho(ptz2(1),ptz2(2),ptz2(3)),0._q2
                     WRITE(97,*) '---------------------------'
                     WRITE(97,*),'Hessian matrix:'
                     WRITE(97,'(3(1X,E18.11))') hes%dxdx, hes%dxdy, hes%dxdz
                     WRITE(97,'(3(1X,E18.11))') hes%dxdy, hes%dydy, hes%dydz
                     WRITE(97,'(3(1X,E18.11))') hes%dxdz, hes%dydz, hes%dzdz
                     WRITE(97,*),'Eigenvalues: '
                     WRITE(97,'(3(1X,E18.11))') ,hes%eigval1, hes%eigval2, hes%eigval3
                     WRITE(97,*) 'Eigenvectors:'
                     WRITE(97,*) eigvec1
                     WRITE(97,*) eigvec2
                     WRITE(97,*) eigvec3
                   END IF
                 END IF
               END IF

        END DO
      END DO
    END DO
    CLOSE(97)
    PRINT *, "CRITICAL POINTS FOUND: ", cptnum 

    DEALLOCATE (hes%rho)
    DEALLOCATE (hes%dx)
    DEALLOCATE (hes%dy)
    DEALLOCATE (hes%dz)
    END SUBROUTINE critpoint_find

  END MODULE


! Copyright 2006-2015
! Wenjie Tang, Andri Arnaldsson, Wenrui Chai, Samuel T. Chill, and Graeme Henkelman
!
! Bader is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! A copy of the GNU General Public License is available at
! http://www.gnu.org/licenses/

!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!
! Authors:
!   Wenjie Tang, Andri Arnaldsson, Wenrui Chai, Samuel T. Chill, and Graeme Henkelman
!
! Authors of the multipole code:
!   Sebastien Lebegue <Sebastien.Lebegue@crm2.uhp-nancy.fr>
!   Angyan Janos <Janos.Angyan@crm2.uhp-nancy.fr>
!   Emmanuel Aubert <emmanuel.aubert@crm2.uhp-nancy.fr>
!
! Contributers:
!   Johannes Voss (DTU), Erik McNellis (FHI), Matthew Dyer (Liverpool),
!   Sören Wohlthat (Sydney)
!
! Based on algorithms described in the following publications:
!
!   A fast and robust algorithm for Bader decomposition of charge density
!   G. Henkelman, A. Arnaldsson, and H. Jonsson
!   Comput. Mater. Sci. 36, 254-360 (2006).
!
!   An improved grid-based algorithm for Bader charge allocation
!   E. Sanville, S. Kenny, R. Smith, and G. Henkelman
!   J. Comput. Chem. 28, 899-908 (2007).
!
!   A grid-based Bader analysis algorithm without lattice bias
!   W. Tang, E. Sanville, and G. Henkelman
!   J. Phys.: Condens. Matter 21, 084204 (2009)
!
!   Accurate and efficient algorithm for Bader charge integration
!   M. Yu and D. Trinkle
!   J. Chem. Phys. 134, 064111 (2011).
!
!-----------------------------------------------------------------------------------!

  PROGRAM Charge

     USE options_mod
     USE ions_mod
     USE charge_mod
     USE io_mod
     USE bader_mod 
     USE voronoi_mod
     USE chgcar_mod
     USE chgcar_mod
     USE multipole_mod
     USE critpoints_mod
     USE weight_mod    
     IMPLICIT NONE

     ! Variables
     TYPE(options_obj) :: opts
     TYPE(ions_obj) :: ions
     TYPE(charge_obj) :: chgval
     TYPE(bader_obj) :: bdr
     TYPE(voronoi_obj) :: vor
    ! Write the version number
     WRITE(*,'(/,2X,A)') 'GRID BASED BADER ANALYSIS  (Version 0.95a 02/26/16)'

     ! Get the control variables
     CALL get_options(opts)

     ! Call the read routines from io_mod
     CALL read_charge(ions,chgval,opts)
     IF (opts%bader_flag) THEN
       IF (opts%bader_opt == opts%bader_weight) THEN
         CALL bader_weight_calc(bdr,ions,chgval,opts)  ! Yu-Trinkle weight method
       ELSE
         CALL bader_calc(bdr,ions,chgval,opts)         ! grid-based Bader methods
       ENDIF
       CALL bader_mindist(bdr,ions,chgval)
       CALL bader_output(bdr,ions,chgval)
       IF (opts%find_critpoints_flag) CALL critpoint_find(bdr,chgval,opts)
     END IF

     IF (opts%print_all_bader) CALL write_all_bader(bdr,opts,ions,chgval)
     IF (opts%print_all_atom) CALL write_all_atom(bdr,opts,ions,chgval)
     IF (opts%print_sel_atom) CALL write_sel_atom(bdr,opts,ions,chgval)
     IF (opts%print_sel_bader) CALL write_sel_bader(bdr,opts,ions,chgval)
     IF (opts%print_sum_atom) CALL write_sum_atom(bdr,opts,ions,chgval)
     IF (opts%print_sum_bader) CALL write_sum_bader(bdr,opts,ions,chgval)
     IF (opts%print_bader_index) CALL write_bader_index(bdr,opts,ions,chgval)
     IF (opts%print_atom_index) CALL write_atom_index(bdr,opts,ions,chgval)
     IF (opts%dipole_flag) CALL multipole_calc(bdr,ions,chgval,opts)
     !IF (opts%voronoi_flag) CALL voronoi(vor,ions,chgval)
     !IF (opts%voronoi_flag) CALL voronoi_divide_conquer(vor,ions,chgval)
     IF (opts%voronoi_flag) CALL voronoi_divide_conquer_kdtree(vor,ions,chgval)

    WRITE(*,*)
  END PROGRAM Charge


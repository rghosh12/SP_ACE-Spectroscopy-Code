!    This module contains the routines that upload the GCOGs.
!
!    It is part of the program SP_Ace, which derives stellar parameters, 
!    such as gravity, temperature, and element abundances from optical
!    stellar spectra, assuming Local Thermodynamic Equilibrium (LTE) 
!    and 1D stellar atmosphere models.
!
!    Copyright (C) 2016 Corrado Boeche
!    
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.

!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!    
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.


MODULE read_GCOG
        USE num_type
        USE space_pars, ONLY: GCOGlib
        USE read_sp_ll
        USE share, ONLY: coeff_4deg
        IMPLICIT NONE
        CHARACTER(len=180) :: file_GCOG
        
        CONTAINS
        
        SUBROUTINE load_GCOG_4deg(temp,logg,met)
        IMPLICIT NONE
        INTEGER(I4B) :: ierror
        INTEGER(I2B) :: i, j, k
        REAL(DP), INTENT(IN) :: temp, logg, met
        CHARACTER(LEN=16) :: logg_label
        CHARACTER(LEN=6) :: temp_label
        CHARACTER(LEN=5) :: met_label
        CHARACTER(LEN=10) :: label


        write(label,fmt='(I2.2)') NINT(logg*10,I1B)
        logg_label='g'// label(1:2) // '-4degpoly.dat'
        write(label,fmt='(I4.4)') NINT(temp,I2B)
        temp_label='t'// label(1:4) // '-'
        IF (met<0) THEN
        write(label,fmt='(I2.2)') NINT(-met*10,I1B)
        met_label='am'// label(1:2) // '-'
        ELSE
        write(label,fmt='(I2.2)') NINT(met*10,I1B)
        met_label='ap'// label(1:2) // '-'
        END IF

        i=INT(INDEX(GCOGlib,' ',.FALSE.),I2B)
        file_GCOG=GCOGlib(1:i-1) // met_label // temp_label // logg_label
        OPEN(unit=10,file=file_GCOG,status='OLD',action='READ',iostat=ierror)
        IF(ierror/=0) CALL error_msg(10_I1B,'I cannot open the GCOG library!')

        i=1
        DO k=1,INT(SIZE(llist_mask,1),I2B)
          IF(llist_mask(k)) THEN
            READ(unit=10,fmt=*,iostat=ierror) (coeff_4deg(j,i),j=1,70)
            i=i+1_I2B
          ELSE
            READ(unit=10,fmt='(A10)',advance='yes',iostat=ierror)
          END IF
          IF(ierror/=0) EXIT
        END DO
        CLOSE(unit=10)

        END SUBROUTINE load_GCOG_4deg
!##########################
        SUBROUTINE load_GCOG_4deg_quick
        IMPLICIT NONE
        INTEGER(I4B) :: ierror
        INTEGER(I2B) :: i, j, k

        i=INT(INDEX(GCOGlib,' ',.FALSE.),I2B)
        file_GCOG=GCOGlib(1:i-1) // 'space_6degpoly.dat'

        OPEN(unit=10,file=file_GCOG,status='OLD',action='READ',iostat=ierror)
        IF(ierror/=0) CALL error_msg(11_I1B,'I cannot open the space_6degpoly.dat file!')

        i=1
        DO k=1,INT(SIZE(llist_mask,1),I2B)
          IF(llist_mask(k)) THEN
            READ(unit=10,fmt=*,iostat=ierror) (coeff_4deg_quick(j,i),j=1,84)
            i=i+1_I2B
          ELSE
            READ(unit=10,fmt='(A10)',advance='yes',iostat=ierror)
          END IF
          IF(ierror/=0) EXIT
        END DO
        CLOSE(unit=10)

        END SUBROUTINE load_GCOG_4deg_quick

END MODULE read_GCOG

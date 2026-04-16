!    
!    This module contains a routine that reads the parameter 
!    file 'space.par'.
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

MODULE space_pars
        USE num_type
        USE share, ONLY: obs_sp_file,sigma,norm_rad,ele2write,n_ele_symb
        USE data_lib, ONLY: ELE_symb

        IMPLICIT NONE
        CHARACTER(120) :: GCOGlib
        CHARACTER(120) :: llist, llist_rej
        CHARACTER(5) :: null_val
        REAL(DP) :: fwhm
        REAL(DP), DIMENSION(10) :: w_inf,w_sup
        REAL(DP) :: sn_ratio,rv_ini
        REAL(DP), DIMENSION(2) :: TGM_force
        INTEGER(I1B) :: N_w_int
        LOGICAL :: error_est,flag_rej,sn_flag,flag_norm,&
        &flag_alpha,flag_ABD_loop,flag_salaris_MH

        CONTAINS

        SUBROUTINE read_space_pars(file_pars)
        CHARACTER(120), INTENT(IN) :: file_pars
        CHARACTER(120) :: string,array,keyword
        CHARACTER(1) :: ele_str1
        CHARACTER(2) :: ele_str2
        CHARACTER(3) :: ele_str3
        CHARACTER(4) :: ele_str4
        INTEGER(I2B) :: i
        INTEGER(I2B) :: j
        INTEGER(I4B) :: ierror
        INTEGER(I2B) :: AllocateStatus,k,n_ele
        INTEGER(I2B), DIMENSION(21), PARAMETER :: ele_allowed=&
        &INT((/6,11,12,13,14,20,21,22,23,24,25,27,28,29,30,39,40,56,57,58,60/),1)
        INTEGER(I2B), DIMENSION(101) :: ele,ele_all
        LOGICAL :: flag_ele2write
        LOGICAL,DIMENSION(21) :: mask_2write
        LOGICAL,DIMENSION(101) :: mask_any_ele
        LOGICAL,DIMENSION(4) :: flag_keys

        !some default values
        n_ele=0
        n_ele_symb = SIZE(ELE_symb)
        sn_ratio=100._dp
        rv_ini=0.
        sn_flag=.TRUE.
        TGM_force = REAL((/-99.,-99./),DP)
        error_est=.FALSE.
        flag_rej=.FALSE.
        flag_ele2write=.FALSE.
        flag_alpha=.FALSE.
        flag_ABD_loop=.TRUE.
        flag_salaris_MH=.TRUE.
        norm_rad=30._dp
        sigma=0.4_dp
        null_val='null'
        mask_2write=.FALSE.
        mask_any_ele=.FALSE.
        flag_keys=.FALSE.
        flag_norm=.TRUE.

        j=0
        DO i=1,n_ele_symb
         IF(i/=26) THEN
          j=j+1_I1B
          ele_all(j)=j
         END IF
        END DO

        OPEN(unit=10,file=TRIM(file_pars),status='OLD',action='READ',iostat=ierror)
        IF(ierror/=0) THEN
         write(*,*) 'SPACE cannot open ', file_pars, ' STOP!'
         STOP 1
        END IF

        DO

          read(unit=10,fmt='(A120)',iostat=ierror) string
          IF(ierror/=0) EXIT

          i=INT(INDEX(string,' ',.FALSE.),I2B)
          keyword= string(1:i-1)
          array = string(i:)

          IF (string(1:1) .eq. '#') THEN
            CYCLE
          ELSEIF (keyword .eq. 'obs_sp_file') THEN
            flag_keys(1)=.TRUE.
            READ (array,fmt=*,iostat=ierror) obs_sp_file
            IF(ierror/=0) CALL stop_local('obs_sp_file missing!')
          ELSEIF (keyword .eq. 'sn_ratio') THEN
           READ (array,fmt=*,iostat=ierror) sn_ratio
           IF(ierror/=0.OR.sn_ratio<1) THEN
            write(*,*) 'sn_ratio wrong! SP_Ace estimates it itself!'
           ELSE
            sn_flag=.FALSE.
           END IF
          ELSEIF (keyword .eq. 'GCOGlib') THEN
            flag_keys(2)=.TRUE.
            READ (array,fmt=*,iostat=ierror) GCOGlib 
            IF(ierror/=0) CALL stop_local('GCOG missing!')
            k=INT(INDEX(GCOGlib,' ',.FALSE.),I1B)
            llist=GCOGlib(1:k-1) // 'linelist.dat'
          ELSEIF (keyword .eq. 'llist_rej') THEN
            READ (array,fmt=*,iostat=ierror) llist_rej
            IF(ierror/=0) CALL stop_local('llist_rej missing!')
            flag_rej=.TRUE.
          ELSEIF (keyword .eq. 'fwhm') THEN
            flag_keys(3)=.TRUE.
            READ (array,fmt=*,iostat=ierror) fwhm
            IF(ierror/=0) CALL stop_local('fwhm missing!')
            sigma=fwhm/2.35_dp
          ELSEIF (keyword .eq. 'wave_lims') THEN
            flag_keys(4)=.TRUE.
            N_w_int=0
            DO j=1,50
             array=ADJUSTL(array)
             i=INT(INDEX(array,' ',.FALSE.),I2B)

             IF(LEN(array(1:i-1))>1) THEN
              N_w_int=N_w_int+1_I1B
              READ (array(1:i-1),fmt=*,iostat=ierror) w_inf(N_w_int)
              IF(ierror/=0) CALL stop_local('w_inf missing!')
              array = ADJUSTL(array(i:))
              i=INT(INDEX(array,' ',.FALSE.),I2B)
              READ (array(1:i-1),fmt=*,iostat=ierror) w_sup(N_w_int)
              IF(ierror/=0) CALL stop_local('w_sup missing!')
              array=array(i:)
             END IF
            END DO
            DO j=1,N_w_int
             IF(w_inf(j)>w_sup(j)) CALL stop_local('wrong wavelenght limits!')
            END DO
          ELSEIF (keyword .eq. 'RV_ini') THEN
           READ (array,fmt=*,iostat=ierror) rv_ini
           IF(ierror/=0.OR.abs(rv_ini)>9999) THEN
            write(*,*) 'Rv beyond the limit of +-9999 km/sec! SP_Ace starts from Rv=0!'
            rv_ini=0
           END IF
          ELSEIF (keyword .eq. 'ele2write') THEN
          flag_ele2write=.TRUE.
           n_ele=0
           DO j=2,n_ele_symb
            IF(j<10) THEN
             WRITE (ele_str1,fmt='(I1)') j
             ele_str3=' ' // ele_str1 // ' '
             i=INT(INDEX(string,ele_str3,.FALSE.),I2B)
             IF(i/=0) THEN
              n_ele=n_ele+1_I1B
              READ (string(i+1:i+2),fmt='(I1)',iostat=ierror) ele(n_ele)
              IF(ierror/=0) CALL stop_local('ele2write missing!')
             END IF
            ELSE
             WRITE (ele_str2,fmt='(I2)') j
             ele_str4=' ' // ele_str2 // ' '
             i=INT(INDEX(string,ele_str4,.FALSE.),I2B)
             IF(i/=0) THEN
              n_ele=n_ele+1_I1B
              READ (string(i+1:i+3),fmt='(I2)') ele(n_ele)
             END IF
            END IF
           END DO
          ELSEIF (keyword .eq. 'T_force') THEN
            READ (array,fmt=*,iostat=ierror) TGM_force(1)
            IF(ierror/=0) CALL stop_local('T_force value missing!')
            IF(TGM_force(1)<3600.OR.TGM_force(1)>7400)&
            & CALL stop_local('T_force must be 3600<T_force<7400 ') 
          ELSEIF (keyword .eq. 'G_force') THEN
            READ (array,fmt=*,iostat=ierror) TGM_force(2)
            IF(ierror/=0) CALL stop_local('G_force value missing!')
            IF(TGM_force(2)<0.2.OR.TGM_force(2)>5.0)&
            & CALL stop_local('G_force must be 0.2<G_force<5.4 ')
          ELSEIF (keyword .eq. 'error_est') THEN
            error_est=.TRUE.
          ELSEIF (keyword .eq. 'no_norm') THEN
            flag_norm=.FALSE.
          ELSEIF (keyword .eq. 'alpha') THEN
            flag_alpha=.TRUE.
          ELSEIF (keyword .eq. 'ABD_loop') THEN
            flag_ABD_loop=.TRUE.
          ELSEIF (keyword .eq. 'Salaris_MH') THEN
            flag_ABD_loop=.TRUE.
            flag_salaris_MH=.TRUE.
          ELSEIF (keyword .eq. 'norm_rad') THEN
            READ (array,fmt=*,iostat=ierror) norm_rad
            IF(ierror/=0) CALL stop_local('norm_rad missing!')
          ELSEIF (keyword .eq. 'null_value') THEN
            READ (array,fmt=*,iostat=ierror) null_val
            IF(ierror/=0) CALL stop_local('null_value is missing!')
            IF(.NOT.(null_val=='-9.99'.OR.null_val=='null'.OR.null_val=='NaN')) CALL stop_local('null_value is wrong!')
          ELSE
           CALL stop_local('There is an unrecognized keyword')
          END IF
        END DO
        CLOSE(unit=10)        

        !if not all the necessary keywords are present, then stop
        IF(.NOT.ALL(flag_keys)) CALL stop_local('one necessary keyword is missing,')


        !if alpha flag is given, then ele2write must contains
        !only the "element" 94 (alpha). All the other elements
        !will be considered as they were iron.
        IF(flag_alpha) THEN
         flag_ele2write=.TRUE.
         n_ele=1
         ele(1)=94
         !allocate ele2write with one place more to put Fe in the first position
         ALLOCATE(ele2write(2), STAT = AllocateStatus)
         IF (AllocateStatus /= 0) CALL stop_local('Not enough memory to allocate ele2write,')
         ele2write(1)=26
         ele2write(2)=94
        ELSE

         !if the keyword ele2write is found, then allocate
         IF(flag_ele2write) THEN
          DO i=1,n_ele
           WHERE(ele(i)==ele_allowed) mask_2write=.TRUE.
          END DO
          n_ele=INT(count(mask_2write),1)
          !allocate ele2write with one place more to put Fe in the first position
          ALLOCATE(ele2write(n_ele+1), STAT = AllocateStatus)
          IF (AllocateStatus /= 0) CALL stop_local('Not enough memory to allocate ele2write,')
          !put Fe in the first position
          ele2write(1)=26
          !now all the other elements requested by the user
           ele2write(2:)=PACK(ele_allowed,mask_2write)
         ELSE IF(.NOT.flag_ele2write) THEN!else give the default values
          n_ele=SIZE(ele_allowed)
          !allocate ele2write with one place more to put Fe in the first position
          ALLOCATE(ele2write(n_ele+1), STAT = AllocateStatus)
          IF (AllocateStatus /= 0) CALL stop_local('Not enough memory to allocate ele2write,')
          ele2write(1)=26
          ele2write(2:)=ele_allowed
          ELSE IF(flag_ele2write) THEN
          DO i=1,n_ele
           WHERE(ele(i)==ele_all) mask_any_ele=.TRUE.
          END DO
          n_ele=INT(count(mask_any_ele),1)
          !allocate ele2write with one place more to put Fe in the first position
          ALLOCATE(ele2write(n_ele+1), STAT = AllocateStatus)
          IF (AllocateStatus /= 0) CALL stop_local('Not enough memory to allocate ele2write,')
          !put Fe in the first position
          ele2write(1)=26
          !now all the other elements requested by the user
           ele2write(2:)=PACK(ele_all,mask_any_ele)
         END IF
        END IF

        END SUBROUTINE read_space_pars
!#########################
        SUBROUTINE stop_local(text)
        CHARACTER(LEN=*) :: text
!
        write(*,*)  text, ', SP_Ace stops!'
        OPEN(unit=10,file='space_msg.txt',action='WRITE')
        WRITE(unit=10,fmt=*) '### ', obs_sp_file
        WRITE(unit=10,fmt=*) '    ', text, ', SP_Ace stops!'
        CLOSE(unit=10)
        STOP 1

        END SUBROUTINE stop_local
!#########################

END MODULE space_pars

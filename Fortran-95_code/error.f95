!    This module contains routines to write null results when
!    SP_ace cannot converge to a meaningful result.
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


MODULE error
        USE num_type
        USE share, ONLY: obs_sp_file
        IMPLICIT NONE
        
        CONTAINS
        
!#########################
        SUBROUTINE error_msg(conv,text)
        
        USE space_pars, ONLY: null_val
        USE share, ONLY: ele2write,space_params_file
        USE data_lib, ONLY: ELE_symb
        IMPLICIT NONE
        INTEGER(I1B), INTENT(IN) :: conv
        INTEGER(I2B) :: i, ele_int
        CHARACTER(LEN=*) :: text
        CHARACTER(LEN=7) :: null
        CHARACTER(LEN=7) :: TGMc,lo_TGMc,up_TGMc
        CHARACTER(LEN=7) :: ABDc,lo_ABDc,up_ABDc
        CHARACTER(LEN=1500) :: values_TGM, values_ABD, values_lab
        CHARACTER(LEN=1500) :: line, header
        CHARACTER(LEN=30) :: string
        CHARACTER(LEN=120) ::  file_out_str

        write(*,*)  text, ' SP_Ace exit with no results!'

        !give the null value of the user
        IF(null_val.EQ.'NaN') THEN
        null='  NaN'
        ELSE IF(null_val.EQ.'null') THEN
        null=' null'
        ELSE IF (null_val.EQ.'-9.99') THEN
        null='-9.99'
        END IF

        TGMc=null
        lo_TGMc=null
        up_TGMc=null
        ABDc=null
        lo_ABDc=null
        up_ABDc=null

            write(line,fmt=*) ''
            write(header,fmt=*) ''

            !add convergence
            header=TRIM(header) // ' conv'
            write(values_ABD,fmt='(TR3,I2)') conv
            line=TRIM(line) // TRIM(values_ABD)
            !add radial velocity
          header=TRIM(header) // '      RV'
            write(values_ABD,fmt='(TR3,A5)') TGMc
            line=TRIM(line) // TRIM(values_ABD)
            !add fwhm
          header=TRIM(header) // '    FWHM'
            write(values_ABD,fmt='(TR3,A5)') TGMc
            line=TRIM(line) // TRIM(values_ABD)
            !add S/N
          header=TRIM(header) // '     S/N'
            write(values_ABD,fmt='(TR3,A5)') null
            line=TRIM(line) // TRIM(values_ABD)
            !add chisq
          header=TRIM(header) // '    chisq'
             write(values_ABD,fmt='(TR4,A5)') null
            line=TRIM(line) // TRIM(values_ABD)

            header=TRIM(header) // '  Teff' // '   T_l' // '   T_h'
            write(values_TGM,fmt='(TR1,A5,TR1,A5,TR1,A5)') TGMc,lo_TGMc,up_TGMc
            line=TRIM(line) // TRIM(values_TGM)
            header=TRIM(header) // '   logg' // '   L_l' // '   L_h'
            write(values_TGM,fmt='(TR2,A5,TR1,A5,TR1,A5)') TGMc,lo_TGMc,up_TGMc
            line=TRIM(line) // TRIM(values_TGM)
            header=TRIM(header) // '    MH' // '  MH_l' // '  MH_h'
            write(values_TGM,fmt='(TR1,A5,TR1,A5,TR1,A5)') TGMc,lo_TGMc,up_TGMc
            line=TRIM(line) // TRIM(values_TGM)



            !now write the elements
            DO i=1,INT(SIZE(ele2write,1),I1B)
              !define the position of the element ele2write(i) in the array ele2meas, if any
               ele_int = INT(ele2write(i),I2B)
               string = ELE_symb(ele_int) // '_l  ' // ELE_symb(ele_int) // '_h '&
                        & // ELE_symb(ele_int) // '_N' 
               write(values_lab,fmt='(TR3,A3,TR2,A19)') ELE_symb(INT(ele_int)), string
               header=TRIM(header) // TRIM(values_lab)
               write(values_ABD,fmt='(TR1,A6,TR1,A6,TR1,A6,TR1,I4)') ABDc,lo_ABDc,up_ABDc,0
               line=TRIM(line) // TRIM(values_ABD)
            END DO        

            !set the name of the output
            file_out_str=space_params_file(1:LEN_TRIM(space_params_file)-4) // '_TGM_ABD.dat'
            !write it
            OPEN(unit=10,file=file_out_str,action='WRITE')!,position='APPEND')        
             write(10,fmt=*) TRIM(header)
             write(10,fmt=*) TRIM(line)
            CLOSE(10)

            !here the program ends!
            write(*,*) 'SP_Ace exits with no solution.'        

        STOP 1
        
        END SUBROUTINE error_msg
!#########################
        SUBROUTINE stop_msg(text)
        IMPLICIT NONE
        CHARACTER(LEN=*) :: text
        CHARACTER(LEN=160) :: line

        write(*,*)  text, ', SP_Ace stops!'
        line=TRIM(obs_sp_file) // ' ' // TRIM(text) // ', SP_Ace stops!'
        OPEN(unit=10,file='space_msg.txt',action='WRITE')
        WRITE(unit=10,fmt=*) line
        CLOSE(unit=10)
        STOP 1

        END SUBROUTINE stop_msg
!#########################

END MODULE error
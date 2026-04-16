!    This routine write the final results.
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


        SUBROUTINE write_res(flag,chisq,conv)
        USE num_type
        USE share, ONLY: TGM, ABD, ele_ll,dim_ll,select_ll_mask,&
        &sn,w_sp,f_sp,f_sp_norm,f_model,cont,&
        &dimsp,ele2meas,weights,ele2write,n_weig,&
        &up_ABD,lo_ABD,sn_var,ew,wave_ll,wave_rej,&
        &rad_rej,dim_rej,sigma,TGM_mask,ABD_mask,space_params_file
        USE uncertains2, ONLY: up_TGM,lo_TGM
        USE space_pars, ONLY: null_val,flag_alpha
        USE data_lib, ONLY: temp_gridL,logg_gridL,met_gridL,ELE_symb
        IMPLICIT NONE        
        INTEGER(I2B) :: i,N_lin,constraints,ele_int
        INTEGER(I1B) :: k,j
        INTEGER(I4B) :: ii
        LOGICAL, DIMENSION(dim_ll) :: mask 
        CHARACTER(LEN=1500) :: values_TGM, values_ABD, values_lab
        CHARACTER(LEN=1500) :: line, header
        LOGICAL, INTENT(IN) :: flag
        REAL(DP), INTENT(IN) :: chisq
        INTEGER(I1B), INTENT(IN) :: conv
        REAL(DP) :: rad
        CHARACTER(LEN=5) :: null
        CHARACTER(LEN=7),DIMENSION(8) :: TGMc
        CHARACTER(LEN=7),DIMENSION(3) :: lo_TGMc,up_TGMc
        CHARACTER(LEN=7),DIMENSION(SIZE(ABD)) :: ABDc,lo_ABDc,up_ABDc
        CHARACTER(LEN=30) :: string
        CHARACTER(LEN=120) ::  file_out_str

        !before to compute the number of lines used, remove the rejected lines from select_ll_lines
        DO i=1,dim_rej
         rad=max(3*sigma,rad_rej(i))
         !plot the rejected lines that are inside a rejected interval
         WHERE(wave_ll>=(wave_rej(i)-rad_rej(i)).AND.wave_ll<=(wave_rej(i)+rad_rej(i))) select_ll_mask=.FALSE.
        END DO

        !give the null value of the user
        IF(null_val.EQ.'NaN') THEN
        null='  NaN'
        ELSE IF(null_val.EQ.'null') THEN
        null=' null'
        ELSE IF (null_val.EQ.'-9.99') THEN
        null='-9.99'
        END IF

        IF(flag) THEN
         !compute the semi-errors for TGM
         write(lo_TGMc(1),fmt='(I5)') INT(lo_TGM(1),I2B)
         write(up_TGMc(1),fmt='(I5)') INT(up_TGM(1),I2B)
         write(lo_TGMc(2),fmt='(F5.2)') lo_TGM(2)
         write(up_TGMc(2),fmt='(F5.2)') up_TGM(2)
         write(lo_TGMc(3),fmt='(F5.2)') lo_TGM(3)
         write(up_TGMc(3),fmt='(F5.2)') up_TGM(3)
        ELSE !if flag is false
         lo_TGMc=null
         up_TGMc=null
         up_ABDc=null
         lo_ABDc=null
        END IF

        !if the errors have no boudaries inside the grid, then give null value
        IF(lo_TGM(1)<temp_gridL(1)) lo_TGMc(1)=null
        IF(up_TGM(1)>temp_gridL(SIZE(temp_gridL))) up_TGMc(1)=null
        IF(lo_TGM(2)<logg_gridL(1)) lo_TGMc(2)=null
        IF(up_TGM(2)>logg_gridL(SIZE(logg_gridL))) up_TGMc(2)=null
        IF(lo_TGM(3)<met_gridL(1)) lo_TGMc(3)=null
        IF(up_TGM(3)>met_gridL(SIZE(met_gridL))) up_TGMc(3)=null

        !write the results of TGM in characters
        !Teff
         write(TGMc(1),fmt='(I5)') INT(TGM(1),I2B)
        !log g
         write(TGMc(2),fmt='(F5.2)') TGM(2)
        !met
         write(TGMc(3),fmt='(F5.2)') TGM(3)
        !fwhm
         write(TGMc(4),fmt='(F6.2)') ABS(TGM(4))*2.35
        !RV

         write(TGMc(5),fmt='(F7.1)') TGM(5)

        !now do it for ABD
        DO j=1,INT(SIZE(ABD,1),I1B)
         IF(ABD(j)>=-0.5.AND.ABD(j)<=0.7.AND.ABD_mask(j)) THEN
           write(ABDc(j),fmt='(F5.2)') ABD(j)+TGM(3)
           write(up_ABDc(j),fmt='(F5.2)') up_ABD(j)+TGM(3)
           write(lo_ABDc(j),fmt='(F5.2)') lo_ABD(j)+TGM(3)
!         write(*,*) j,TGM(3),ABD(j),up_ABD(j),lo_ABD(j)
         ELSE
        !if the enhancement is >1. or <-1, give the null value
          ABDc(j)=null
          up_ABDc(j)=null
          lo_ABDc(j)=null
         END IF
        END DO

        WHERE(lo_ABD<-0.5) lo_ABDc=null
        WHERE(up_ABD>0.7) up_ABDc=null

          !write the first line in the result file
            write(line,fmt=*) ''
            write(header,fmt=*) ''

            !add convergence
            header=TRIM(header) // '  conv'
            write(values_ABD,fmt='(TR3,I3)') conv
            line=TRIM(line) // TRIM(values_ABD)
            !add radial velocity
          header=TRIM(header) // '        RV'
            write(values_ABD,fmt='(TR3,A7)') TGMc(5)
            line=TRIM(line) // TRIM(values_ABD)
            !add fwhm
          header=TRIM(header) // '    FWHM'
            write(values_ABD,fmt='(TR2,A7)') TGMc(4)
            line=TRIM(line) // TRIM(values_ABD)
            !add S/N
          header=TRIM(header) // '     S/N'
            write(values_ABD,fmt='(TR3,F5.1)') sn
            line=TRIM(line) // TRIM(values_ABD)
            !add chisq
          constraints=INT(COUNT(TGM_mask),I2B)+INT(COUNT(ABD_mask),I2B)
          header=TRIM(header) // '    chisq'
            write(values_ABD,fmt='(TR3,F6.2)') chisq/(n_weig-constraints)
            line=TRIM(line) // TRIM(values_ABD)

            header=TRIM(header) // '  Teff' // '   T_l' // '   T_h'
            write(values_TGM,fmt='(TR1,A5,TR1,A5,TR1,A5)') TGMc(1),lo_TGMc(1),up_TGMc(1)
            line=TRIM(line) // TRIM(values_TGM)
            header=TRIM(header) // '   logg' // '   L_l' // '   L_h'
            write(values_TGM,fmt='(TR2,A5,TR1,A5,TR1,A5)') TGMc(2),lo_TGMc(2),up_TGMc(2)
            line=TRIM(line) // TRIM(values_TGM)
            header=TRIM(header) // '    MH' // '  MH_l' // '  MH_h'
            write(values_TGM,fmt='(TR1,A5,TR1,A5,TR1,A5)') TGMc(3),lo_TGMc(3),up_TGMc(3)
            line=TRIM(line) // TRIM(values_TGM)

            !now write the elements
            IF(flag_alpha) THEN
             !write "metals"
             mask=.FALSE.
             WHERE(INT(ele_ll)==26.AND.select_ll_mask) mask=.TRUE.
             string = ELE_symb(93) // '_l  ' // ELE_symb(93) // '_h ' // ELE_symb(93) // '_N' 
             write(values_lab,fmt='(TR3,A3,TR2,A19)') ELE_symb(93), string
             header=TRIM(header) // TRIM(values_lab)
             N_lin=INT(COUNT(mask),I2B)
             IF(N_lin>0) THEN
              j=INT(MINLOC(ABS(ele2meas-26),1),I1B)
              write(values_ABD,fmt='(TR1,A6,TR1,A6,TR1,A6,TR1,I4)') ABDc(j),lo_ABDc(j),up_ABDc(j),N_lin
             ELSE
              write(values_ABD,fmt='(TR1,A5,TR1,A6,TR1,A6,TR1,I4)') null,null,null,0              
             END IF
             line=TRIM(line) // TRIM(values_ABD)
             !write alpha
             mask=.FALSE.
             WHERE(INT(ele_ll)==94.AND.select_ll_mask) mask=.TRUE.
             string = ELE_symb(94) // '_l  ' // ELE_symb(94) // '_h ' // ELE_symb(94) // '_N' 
             write(values_lab,fmt='(TR3,A3,TR2,A19)') ELE_symb(94), TRIM(string)
             header=TRIM(header) // TRIM(values_lab)
             N_lin=INT(COUNT(mask),I2B)
             IF(N_lin>0) THEN
              j=INT(MINLOC(ABS(ele2meas-94),1),I1B)
              write(values_ABD,fmt='(TR1,A6,TR1,A6,TR1,A6,TR1,I4)') ABDc(j),lo_ABDc(j),up_ABDc(j),N_lin
             ELSE
              write(values_ABD,fmt='(TR1,A6,TR1,A6,TR1,A6,TR1,I4)') null,null,null,0              
             END IF
             line=TRIM(line) // TRIM(values_ABD)

            ELSE !if flag_alpha=.FALSE. then

            DO k=1,INT(SIZE(ele2write,1),I1B)
              !define the position of the element ele2meas(k) in the array ele2write, if any
              IF(ANY(ele2meas==ele2write(k))) THEN
               mask=.FALSE.
               j=INT(MINLOC(ABS(ele2meas-ele2write(k)),1),I1B)
               WHERE(INT(ele_ll)==ele2meas(j).AND.select_ll_mask) mask=.TRUE.
               ele_int = INT(ele2meas(j),I2B)
               string = ELE_symb(ele_int) // '_l  ' // ELE_symb(ele_int) // '_h '&
                        & // ELE_symb(ele_int) // '_N' 
                write(values_lab,fmt='(TR3,A3,TR2,A19)') ELE_symb(INT(ele_int)), string
               header=TRIM(header) // TRIM(values_lab)
               N_lin=INT(COUNT(mask),I2B)
               write(values_ABD,fmt='(TR1,A6,TR1,A6,TR1,A6,TR1,I4)') ABDc(j),lo_ABDc(j),up_ABDc(j),N_lin
               line=TRIM(line) // TRIM(values_ABD)
              ELSE
                ele_int = INT(ele2write(k),I2B)
                string = ELE_symb(ele_int) // '_l  ' // ELE_symb(ele_int) // '_h '&
                        & // ELE_symb(ele_int) // '_N' 
                write(values_lab,fmt='(TR3,A3,TR2,A19)') ELE_symb(INT(ele_int)), string
               header=TRIM(header) // TRIM(values_lab)
               write(values_ABD,fmt='(TR1,A6,TR1,A6,TR1,A6,TR1,I4)') null,null,null,0
               line=TRIM(line) // TRIM(values_ABD)
              END IF
            END DO        
           END IF

            !set the name of the output
            file_out_str=space_params_file(1:LEN_TRIM(space_params_file)-4) // '_TGM_ABD.dat'
            !write it
            OPEN(unit=10,file=TRIM(file_out_str),action='WRITE')
             write(10,fmt=*) TRIM(header)
             write(10,fmt=*) TRIM(line)
            CLOSE(10)

            !set the name of the output
            file_out_str=space_params_file(1:LEN_TRIM(space_params_file)-4) // '_ew_meas.dat'
            !now write the files with the EWs of the measured lines
            OPEN(unit=10,file=TRIM(file_out_str),action='WRITE')        
            DO i=1,dim_ll
             IF(select_ll_mask(i)) THEN
              write(10,fmt='(F8.3,2X,F5.1,2X,F6.1)') wave_ll(i), ele_ll(i), ew(i)*1000.
             END IF
            END DO
            CLOSE(10)        

          !set the name of the output
          file_out_str=space_params_file(1:LEN_TRIM(space_params_file)-4) // '_model.dat'
          !write it
          OPEN(unit=9,file=TRIM(file_out_str),action='WRITE')        
          write(9,fmt="(F9.3,TR1,F8.5,TR1,F8.5,TR1,F8.5,TR1,F8.5,TR1,F4.2,TR1,I4)")&
          & (w_sp(ii),f_sp(ii),f_sp_norm(ii),f_model(ii),cont(ii),weights(ii),NINT(sn_var(ii)), ii=1,dimsp)
          CLOSE(9)

        
        END SUBROUTINE write_res 

!    This routine fits the spectral continuum.
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


        SUBROUTINE fit_cont(f_sp,f_model,cont)
        USE num_type
        USE stats
        USE share, ONLY: dimsp,weights,rad_pix
        IMPLICIT NONE

        REAL(DP) :: avg_f,avg_r,avg_m,var,sig,sig_r,avg_mm,sig_mm,set,weig
        REAL(DP),DIMENSION(dimsp), INTENT(IN) :: f_sp,f_model
        REAL(DP),DIMENSION(dimsp) :: resid
        REAL(DP),DIMENSION(dimsp), INTENT(OUT) :: cont
        INTEGER(I4B) :: i,iinf,isup,dimsp_loop
        LOGICAL,DIMENSION(dimsp) :: mask

        resid=f_sp-f_model
 

        IF (rad_pix(1)<dimsp) THEN
         dimsp_loop=dimsp
        ELSE
         dimsp_loop=1
        END IF

       DO i=1,dimsp_loop

        IF(dimsp_loop>1) THEN
        weig=weights(i)
        ELSE
        weig=1.
        END IF
        
        IF(weig>1e-3) THEN
           mask=.FALSE.

           isup=MIN(i+rad_pix(i),dimsp);
           iinf=MAX(i-rad_pix(i),1_I4B);
           WHERE(weights(iinf:isup)>1e-3)
            mask(iinf:isup)=.TRUE.
           END WHERE

           IF(ANY(mask(iinf:isup))) THEN
            CALL avg_std_(f_sp(iinf:isup),mask(iinf:isup),avg_f,var,sig)
            CALL avg_std_(resid(iinf:isup),mask(iinf:isup),avg_r,var,sig_r)
            !pseudo-continuum level of the model
            CALL avg_std_(f_model(iinf:isup),mask(iinf:isup),avg_mm,var,sig_mm)
            !clip 2 sigma low
            WHERE(f_sp(iinf:isup)<=(avg_f-2.0_dp*sig_r)) mask(iinf:isup)=.FALSE. 
            !clip 4 sigma high
            WHERE(f_sp(iinf:isup)>=(1.0+4.0_dp*sig_r)) mask(iinf:isup)=.FALSE.
           END IF !IF(ANY(mask(iinf:isup))) THEN

           IF(COUNT(mask)>1) THEN
            CALL avg_std_(f_sp(iinf:isup),mask(iinf:isup),avg_f,var,sig)
            CALL avg_std_(f_model(iinf:isup),mask(iinf:isup),avg_m,var,sig)
            CALL avg_std_(resid(iinf:isup),mask(iinf:isup),avg_r,var,sig_r)
           ELSE
            avg_f=1.0_dp
            avg_m=1.0_dp
            sig_r=0.0_dp
           END IF !IF(COUNT(mask)>1) THEN

          set=4.*(1.-exp((1.-avg_mm)**3))
           cont(i)=1.0_dp+(avg_f-avg_m-set);

        ELSE !if weig<1e-3 then

           cont(i)=1.0_dp
          
        END IF

       END DO

        !if the radius is bigger than the extention of
        ! the spectrum, then cont(1) is valid for
        !the whole spectrum.
        IF (rad_pix(1)>=dimsp) THEN
         cont=cont(1)
        END IF



        END SUBROUTINE fit_cont
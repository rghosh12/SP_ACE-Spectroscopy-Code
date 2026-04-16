!    This module contains the necessary interfaces of some routines.
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


MODULE interfaces

        INTERFACE
                SUBROUTINE fit_cont(f_sp,f_model,cont)
                USE num_type
                USE stats
                USE share, ONLY: dimsp,weights,rad_pix
                IMPLICIT NONE
                REAL(DP) :: avg_f,avg_r,avg_m,var,sig,sig_r,avg_mm,sig_mm,set,weig
                REAL(DP),DIMENSION(dimsp), INTENT(IN) :: f_sp,f_model
                REAL(DP),DIMENSION(dimsp) :: resid
                REAL(DP),DIMENSION(dimsp), INTENT(OUT) :: cont
                INTEGER(I2B) :: i,iinf,isup,dimsp_loop
                LOGICAL,DIMENSION(dimsp) :: mask
                END SUBROUTINE fit_cont
        END INTERFACE
        INTERFACE
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
        INTEGER(I2B) :: i,N_lin
        INTEGER(I1B) :: k,j
        LOGICAL, DIMENSION(dim_ll) :: mask 
        CHARACTER(LEN=1500) :: values_TGM, values_ABD, values_lab
        CHARACTER(LEN=1500) :: line, header
        LOGICAL, INTENT(IN) :: flag
        REAL(DP), INTENT(IN) :: chisq
        INTEGER(I1B), INTENT(IN) :: conv
        REAL(DP) :: rad
        CHARACTER(LEN=5) :: null
        CHARACTER(LEN=5),DIMENSION(5) :: TGMc,lo_TGMc,up_TGMc
        CHARACTER(LEN=5),DIMENSION(SIZE(ABD)) :: ABDc,lo_ABDc,up_ABDc
        CHARACTER(LEN=120) ::  file_out_str

        END SUBROUTINE write_res
        END INTERFACE

END MODULE interfaces
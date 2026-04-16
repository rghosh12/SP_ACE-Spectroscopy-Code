!    
!    This module contains variable and arrays that can be shared among the
!    different modules/routines employed.
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

MODULE share
        USE num_type
        IMPLICIT NONE
        CHARACTER(120) :: obs_sp_file
        CHARACTER(LEN=120) :: space_params_file
        REAL(DP), DIMENSION(:), ALLOCATABLE :: w_sp, f_sp, f_sp_norm, f_model
        REAL(DP), DIMENSION(:), ALLOCATABLE :: cont0, cont, weights, f_discrep
        REAL(DP), DIMENSION(:), ALLOCATABLE :: wave_ll,ele_ll,ew,disp,Ex_inf
        REAL(DP), DIMENSION(:,:), ALLOCATABLE :: X_abd
        INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: wave_center_ll
        INTEGER(I2B), DIMENSION(:), ALLOCATABLE :: flag_lines
        REAL(DP), DIMENSION(:), ALLOCATABLE :: wave_rej,rad_rej,sn_var,sig_noise
        REAL(DP), DIMENSION(:,:), ALLOCATABLE :: coeff_4deg,coeff_4deg_quick
        INTEGER(I4B) :: dimsp,n_weig
        INTEGER(I2B) :: dim_ll,dim_rej,n_ele_symb
        INTEGER(I2B) :: dim_ele,dim_ele_dy
        INTEGER(I2B), DIMENSION(:), ALLOCATABLE :: ele2meas,ele2write
        INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: rad_pix
        REAL(DP), DIMENSION(6) :: TGM
        REAL(DP), DIMENSION(:), ALLOCATABLE :: ABD,up_ABD,lo_ABD,ABD_old,residABD
        REAL(DP) :: sn,sigma,fwhm,norm_rad
        REAL(DP) :: temp_infS, temp_supS, logg_infS, logg_supS, met_infS, met_supS
        REAL(DP), DIMENSION(3) :: TGM_prox
        LOGICAL, DIMENSION(:), ALLOCATABLE :: select_ll_mask,cosmic_mask
        LOGICAL, DIMENSION(:), ALLOCATABLE :: ABD_mask,write_ABD_mask,alpha_mask
        LOGICAL, DIMENSION(6) :: TGM_mask
END MODULE share
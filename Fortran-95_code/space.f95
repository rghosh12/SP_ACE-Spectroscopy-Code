!    SP_Ace derives stellar parameters, such as gravity, temperature, and element 
!    abundances from optical stellar spectra, assuming Local Thermodynamic 
!    Equilibrium (LTE) and 1D stellar atmosphere models.
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

PROGRAM space
        USE num_type
        USE error
        USE space_pars, ONLY: read_space_pars,sn_flag,&
         &TGM_force,error_est,flag_norm,flag_ABD_loop,rv_ini,&
         &flag_salaris_MH
        USE share
        USE utils
        USE make_model
        USE read_GCOG
        USE func_poly
        USE interfaces, ONLY: fit_cont,write_res
        USE minimize
        USE uncertains2, ONLY:TGM_errors,TGM_errors_null
        IMPLICIT NONE
        INTEGER(I2B) :: i,j,infoTGM,infoABD
        INTEGER(I1B) :: AllocateStatus,conv
        REAL(DP),DIMENSION(5) :: temp_gridS,logg_gridS,met_gridS
        REAL(DP), DIMENSION(3) :: TGM1,TGM_old, residTGM
        REAL(DP), DIMENSION(4) :: pars
        REAL(DP), DIMENSION(:), ALLOCATABLE :: TGMx,ABDx,TGMxx
        REAL(DP) :: chisq,alpha_mean
        REAL(DP) :: FeH,aFe,M_salaris,diff_M
        LOGICAL :: flag_lim,flag_limS,flag_move
        CHARACTER(LEN=5),DIMENSION(3) :: TGMc

        !write the SP_Ace version
        write(*,*) 'SP_Ace version v1.4'
        write(*,*) 'Copyright (C) 2020 Corrado Boeche'
        write(*,*) 'This program comes with ABSOLUTELY NO WARRANTY.'
        write(*,*) 'This is free software, and you are welcome to redistribute it'
        write(*,*) 'under certain conditions; see <http://www.gnu.org/licenses/> for details.'

        !initialize variables
        TGM_prox=REAL((/-99.,-99.,-99./),SP)
        flag_lim=.FALSE.
        flag_limS=.FALSE.
        flag_move=.FALSE.

        !retrieve the name of the parameters file from the command line
        CALL getarg(1, space_params_file)
        space_params_file = TRIM(space_params_file)
        IF (LEN_TRIM(space_params_file)<1) THEN
          space_params_file='space.par'
        END IF

        !read the parameters file
        CALL read_space_pars(space_params_file)
        !load the spectrum, the line lists
        CALL read_files(obs_sp_file,llist,llist_rej)
        !----------------------
        !allocate the mask for the llist that will be use for the measures
        ALLOCATE(select_ll_mask(dim_ll), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate select_ll_mask')

        !allocate the array sig_noise
        ALLOCATE(sig_noise(dimsp), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate sig_noise')
        !allocate the array f_discrep
        ALLOCATE(f_discrep(dimsp), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate f_discrep')

        ! allocate the array ele2meas which contains the list of elements
        ! present in the line list. It also allocate ABD, ADB_old,residABD
        ! ABD_mask. It outputs dim_ele
        CALL alloc_ABD(ele_ll,dim_ll,dim_ele,n_ele_symb)

        ALLOCATE(ABDx(dim_ele), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate ABDx')
        ABDx=PACK(ABD,ABD_mask)

        !initialize TGM_mask
        TGM_mask=.TRUE.
        !initialize conv=0 by default
        conv=0
        !initialize the sn
        sn = sn_ratio

        !give a starting point
        TGM(1:3)=(/5000.,3.0,-0.4/)
        !initialize Teff and gravity to the user's input, if any
        IF (TGM_force(1)>-9.) THEN
          TGM(1)=TGM_force(1)
          TGM_mask(1)=.FALSE.
        END IF  
        IF (TGM_force(2)>-9.) THEN
          TGM(2)=TGM_force(2)
          TGM_mask(2)=.FALSE.
        END IF  
        !initialize other TGM variables
        TGM(4)=sigma
        TGM(5)=rv_ini
        TGM(6)=1.0
        !if the user do not want the continuum normalization,
        !switch off the variable that control the continuum
        IF(.NOT.flag_norm) THEN
         TGM_mask(6)=.FALSE.
        END IF

         !allocate the array TGMx to be used as shorted array when TGM_force is employed
         !and fed to the lmdif1 routine. TGMx must be long as much as the number of .TRUE. in TGM_mask
         ALLOCATE(TGMx(COUNT(TGM_mask)), STAT = AllocateStatus)
         IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate TGMx')

        !initialize the internal normalized spectrum
        f_sp_norm=f_sp
        cont0=1.0_dp
        !prepare for the first TGM estimation
        CALL load_GCOG_4deg_quick
        weights=1.0_dp
        sig_noise=(1.0_dp/sn_var)/weights
        select_ll_mask=.TRUE.
        TGM_old=TGM(1:3)
        TGM1=TGM_old

        !start normalization loop that serves as a first guess too
        CALL update_gridS_ll
        normal_loop: DO j=1,30
         CALL make_model_TGM_quick(f_model,TGM)
         CALL normalize_pars(TGM, TGMx, TGM_mask)
         CALL lmdif1(chi_TGM_Q,dimsp,INT(SIZE(TGMx),I4B),TGMx,f_discrep,1E-3,infoTGM)
         CALL denormalize_pars(TGM, TGMx, TGM_mask)
         !the following check if the values obscillate.
         !if the values obscillate, then keep result and continue the analysis
         IF(MOD(j,2_I2B)==0) THEN
          IF(ABS(TGM1(1)-TGM(1))<1..AND.ALL(ABS(TGM1(2:3)-TGM(2:3))<0.01)) EXIT
          TGM1=TGM(1:3)
         END IF         
         write(TGMc(1),fmt='(I5)') INT(TGM(1))
         write(TGMc(2),fmt='(F5.2)') TGM(2)
         write(TGMc(3),fmt='(F5.2)') TGM(3)
         write(*,*) 'TGM norm. loop',j,' Teff=',TGMc(1),' log g=',TGMc(2),' [M/H]=',TGMc(3)
         !, 'sigma=', TGM(4)!, ' info=',infoTGM!,'cont=',TGM(6)

         !make the model
         sigma=TGM(4)
         CALL make_model_TGM_quick(f_model,TGM)
         CALL update_gridS_ll

         !now re-normalize the observed spectrum
         IF(flag_norm) THEN
          CALL fit_cont(f_sp,f_model,cont);
          cont=(cont+cont0)/2.0_dp
          f_sp_norm=f_sp/TGM(6)/cont
          cont0=cont
         END IF

          !compute the chisq
          chisq=SUM(((f_sp_norm-f_model)/sig_noise)**2)
          !update some variables
          residTGM=ABS(TGM_old-TGM(1:3))
          TGM_old=TGM(1:3)

         !if the user does not want the re-normalization, then exit this loop 
         IF(.NOT.flag_norm) EXIT 
         !check the convergence
         IF((residTGM(1)<10..AND.residTGM(2)<0.02.AND.residTGM(3)<0.02).AND.j>1) EXIT

        END DO normal_loop
        
        IF(j>=30) THEN
!         CALL error_msg(1_I1B,'norm loop did not converge')
          ! if it does not converge in 30 iterations, continue with flag conv=1
          conv=1
        END IF

        !final normalization of the observed spectrum
        f_sp_norm=f_sp/TGM(6)/cont

        CALL update_gridS_ll

        !allocate a new array that does not contain the variables FWHM and RV
        ALLOCATE(TGMxx(COUNT(TGM_mask)), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate TGMxx')


!#########################
!###### start iterations outer loop
!#########################
        outer: DO i=1,30
        !##############
        !search for TGM
        !##############

        CALL normalize_pars(TGM, TGMxx, TGM_mask)
        CALL lmdif1(chi_TGM,dimsp,INT(SIZE(TGMxx),I4B),TGMxx,f_discrep,1E-3,infoTGM)
        CALL denormalize_pars(TGM, TGMxx, TGM_mask)
        CALL make_model_TGM(f_model,TGM)
        chisq=SUM(((f_sp_norm-f_model)/sig_noise)**2)
        !write to the standard output the temporary solution found
        write(TGMc(1),fmt='(I5)') INT(TGM(1))
        write(TGMc(2),fmt='(F5.2)') TGM(2)
        write(TGMc(3),fmt='(F5.2)') TGM(3)
        write(*,*) 'TGM outer loop',i,' Teff=',TGMc(1),' log g=',TGMc(2),' [M/H]=',TGMc(3)!, ' info=',infoTGM

        !update the value sigma
        sigma=ABS(TGM(4))
        !find the closest point in the sub-grid and verify if TGM is inside the grid
        CALL find_proxS(TGM,temp_gridS,logg_gridS,met_gridS,TGM_prox,&
        &flag_lim,flag_limS,flag_move)
        !if flag_lim=.TRUE. means that TGM is out of the grid, then exit
        IF(flag_lim) THEN
           CALL error_msg(2_I1B,'stellar parameters out of the limits')
        END IF
        !if flag_limS=.TRUE. means that TGM is out of the 
        !gridS, then update the GCOG, the linelist and 
        !cycle the outer loop
        IF(flag_limS) THEN
           CALL update_gridS_ll
           TGM_old=TGM(1:3)
           CYCLE
        END IF
        IF(flag_move) THEN
           CALL update_gridS_ll
           TGM_old=TGM(1:3)
           CYCLE
        END IF


         !update some variables
         sigma=ABS(TGM(4))
         residTGM=ABS(TGM_old-TGM(1:3))
         TGM_old=TGM(1:3)
         !check the convergence of the outer loop
         IF((residTGM(1)>20..OR.residTGM(2)>0.05.OR.residTGM(3)>0.05).AND.i<30) THEN
          CYCLE
         ELSE
          
         !check if some of the parameters are meaningful.
         !check the FWHM: if it is too large, quit
         IF((sigma*2.35)>10.) THEN
          CALL error_msg(3_I1B,'FWHM does not converge')
         END IF
         !check the abs(RV): if it is too large, quit
         IF(ABS(TGM(5)-rv_ini)*w_sp(1)/300000.>sigma*2.35) THEN
          CALL error_msg(4_I1B,'RV too far from the initial one (beyond 1FWHM)')
         END IF

          !##################
          !now search for ABD
          !##################
        write(*,*) 'ABD loop'
        !assign to the dummy ABDx array the same values of ABD
        ABDx=ABD
        CALL lmdif1(chi_ABD,dimsp,INT(SIZE(ABDx),I4B),ABDx,f_discrep,1E-3,infoABD)
        !now assign to ABD the new values of ABDx
        ABD=UNPACK(ABDx,ABD_mask,ABD)
        CALL make_model_ABD(f_model,ABD)

          IF(flag_ABD_loop) THEN
           !compute the differences between the old and new ABD estimations
           residABD=ABS(ABD-ABD_old)
           !if the differences are large, then re-do the loop 
           IF(ANY((PACK(residABD,write_ABD_mask))>0.05).AND.i<30) THEN
           !Now substitute the Fe abundance (if the internal loop cycled at least one time) in TMG(3)  
           !if this is the first loop for the internal loop, then ABD=0 and TGM does not change
           !if flag_salaris_MH=.TRUE., put in TMG(3) the Salaris metallicity
           IF(flag_salaris_MH) THEN
             IF(flag_alpha) THEN
               alpha_mean = ABD(INT(MINLOC(ABS(ele2meas-94),1),I2B))
             ELSE
               alpha_mean = compute_aFe(ABD, ele2meas, ABD_mask) 
             END IF
             FeH=TGM(3)+ABD(1)
             aFe=(TGM(3)+alpha_mean)-FeH
             M_salaris=FeH+LOG10(0.638*10**aFe+0.362)
             diff_M=TGM(3)-M_salaris
!            write(*,*) 'TGM(3),Fe,aFe,M salaris,diff',TGM(3),FeH,aFe,M_salaris,diff_M
             WHERE(ABD_mask) ABD=ABD+diff_M
             TGM(3)=M_salaris
           ELSE
             TGM(3)=TGM(3)+ABD(1)
             WHERE(ABD_mask) ABD=ABD-ABD(1)
           END IF

           !update the ABD_old
           ABD_old=ABD

           !do not allow ABD to go beyond the grid
           WHERE((ABD>0.7.OR.ABD<-0.5).AND.ABD_mask) 
            ABD=0._dp
            ABD_mask=.FALSE.
           END WHERE

            CYCLE
           END IF
          END IF

          EXIT
          END IF

        !check the convergence and flags of the lmdif1 routines
        IF(infoTGM==0.OR.infoABD==0) THEN
         CALL error_msg(5_I1B,'improper input parameters in lmdif1)')
        ELSE IF (infoTGM==5) THEN
         CALL error_msg(6_I1B,'TGM minim routine exceeded the max N of iterations')
         EXIT
        ELSE IF (infoABD==5) THEN
         CALL error_msg(7_I1B,'ABD minim routine exceeded the max N of iterations')
         EXIT
        ELSE
         EXIT
        END IF

        END DO outer
!######################
        !if the outer loop has more than 30 loops, it means
        !it could not converge, then quit
        IF(i>=30) THEN
         !CALL error_msg(8_I1B,'TGM outer loop did not converge')
         conv=1
        END IF

        !make the final model and compute the chi square
        CALL make_model_TGM(f_model,TGM)
        chisq=SUM(((f_sp_norm-f_model)/sig_noise)**2)

         !compute the overall S/N
         IF(sn_flag) THEN
          CALL new_sn(sn)
         END IF

        !compute the errors, if requested by the user
        IF(error_est) THEN
         write(*,*) 'errors estimation in progress'
         CALL TGM_errors(chisq)
         !if you decomment the following line, the
         !chisq surface is output in the file space_chisq.dat
         !CALL write_chi2_grid(chisq)
!         CALL TGM_errors_poly(chisq)
         ELSE
         CALL TGM_errors_null
        END IF

        !write the results and exit
        CALL make_model_TGM(f_model,TGM)
        CALL write_res(error_est,chisq,conv)
        !here the program ends!
        write(*,*) 'SP_Ace successfully exits.'

        CONTAINS

!#######################################
        SUBROUTINE update_gridS_ll
        IMPLICIT NONE
        LOGICAL,DIMENSION(dimsp) :: mask1

        !make a new subgrids
        CALL find_prox(TGM,TGM_prox)
        CALL make_gridS(TGM_prox,temp_gridS,logg_gridS,met_gridS)

        ! allocate and upload the table coeff_4deg which contains the coefficient
        ! to construct the polynomial GCOG
        CALL load_GCOG_4deg(TGM_prox(1),TGM_prox(2),TGM_prox(3))

        !remove the cosmic rays with a 4sigma clipping
        CALL cosmic_rej


        IF(sn_flag) THEN
         !compute the REAL(sn), which is the spectrum overall S/N
         CALL new_sn(sn)
         !compute the array sn_var, which is the S/N pixel by pixel
         CALL find_sn_var(f_sp_norm,f_model,sn_var)
        END IF

        !select the lines strong enough to be measurable in this
        !region of the parameters space
        pars=(/TGM(1),TGM(2),TGM(3),0._dp/)
        CALL select_lines(ew_poly4,coeff_4deg,pars,.TRUE.)
        sig_noise=(1.0_dp/sn_var)/weights

        !give to n_weig the number of pixels having weight>0.01
        mask1=.TRUE.
        WHERE(weights<0.01) mask1=.FALSE.
        n_weig=INT(COUNT(mask1),I4B)
        IF(n_weig<1) THEN
         n_weig=dimsp
        END IF

        END SUBROUTINE update_gridS_ll
!########################################
END PROGRAM space

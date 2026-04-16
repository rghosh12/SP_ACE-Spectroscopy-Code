!    This module contains several routines for different uses.
!    See below for a short description of each routine. 
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


MODULE utils
        USE num_type
        USE error
        USE data_lib
        USE share, ONLY: ABD_mask,write_ABD_mask
        IMPLICIT NONE


        CONTAINS

!########################
        !this routine set the cosmic_mask
        SUBROUTINE cosmic_rej
        USE share, ONLY: f_sp,f_model,weights,cosmic_mask,dimsp
        USE stats
        IMPLICIT NONE
        REAL(DP) :: avg_r,var,sig_r
        REAL(DP),DIMENSION(dimsp) :: resid
        LOGICAL,DIMENSION(dimsp) :: mask
       
        mask=.FALSE. 
        resid=f_sp-f_model
        WHERE(weights>1e-3_dp) mask=.TRUE.
        CALL avg_std_(resid,mask,avg_r,var,sig_r)
        WHERE(f_sp>(1._dp+4.0_dp*sig_r)) cosmic_mask=.TRUE.
        WHERE(f_sp<=0.01) cosmic_mask=.TRUE.
        
        END SUBROUTINE cosmic_rej
!########################
        !this routine computes the overall S/N
        SUBROUTINE new_sn(sn)
        USE share, ONLY: f_sp_norm, f_model,dimsp,weights
        USE read_sp_ll, ONLY: add_pix_sp
        USE stats
        IMPLICIT NONE
        REAL(DP), INTENT(INOUT) :: sn
        REAL(DP), DIMENSION(dimsp) :: resid
        REAL(DP) :: mean,var,sdev1,sdev2
        LOGICAL,DIMENSION(dimsp) :: mask1,mask2,mask3

        !initialize
        mask1=.TRUE.
        mask2=.FALSE.
        mask3=.FALSE.
        
        resid=f_sp_norm-f_model
        WHERE(weights<0.001) mask1=.FALSE.
        !do not consider the extreme of the spectrum where the flux=1. artificially
        mask1(1:add_pix_sp)=.FALSE.
        mask1(dimsp-add_pix_sp:dimsp)=.FALSE.

        CALL avg_std_(resid,mask1,mean,var,sdev1)
        WHERE(ABS(resid-mean)<(3.*sdev1)) mask2=.TRUE.
        WHERE(mask1.AND.mask2) mask3=.TRUE.

        IF(ALL(.NOT.mask3)) THEN
         sn=-1
        ELSE
         CALL avg_std_(resid,mask3,mean,var,sdev2)
         !for some bad spectra resid=0.00 for many pixels, then
         IF(sdev2<1e-4) THEN 
          sdev2=10.0_dp
         END IF
         sn=1/sdev2
        END IF
        END SUBROUTINE new_sn

!########################
        !this routine computes the S/N in the neighbouring of each pixel
        SUBROUTINE find_sn_var(f_sp_norm,f_model,sn_var)
        USE stats
        USE share, ONLY: dimsp,weights
        IMPLICIT NONE

        REAL(DP) :: avg_r,var,sig_r
        REAL(DP),DIMENSION(dimsp), INTENT(IN) :: f_sp_norm,f_model
        REAL(DP),DIMENSION(dimsp) :: resid
        REAL(DP),DIMENSION(dimsp), INTENT(OUT) :: sn_var
        INTEGER(I4B) :: i,iinf,isup,int_pix
        LOGICAL,DIMENSION(dimsp) :: mask

        resid=f_sp_norm-f_model
        int_pix=25_I4B

       DO i=1,dimsp


        IF(weights(i)>1e-3) THEN
           mask=.FALSE.

           isup=MIN(i+int_pix,dimsp);
           iinf=MAX(i-int_pix,1_I4B);
           WHERE(weights(iinf:isup)>1e-3)
            mask(iinf:isup)=.TRUE.
           END WHERE
           
           CALL avg_std_(resid(iinf:isup),mask(iinf:isup),avg_r,var,sig_r)
           WHERE(ABS(resid(iinf:isup)-avg_r)>(3*sig_r)) mask(iinf:isup)=.FALSE.

           IF(COUNT(mask(iinf:isup))>2) THEN
            CALL avg_std_(resid(iinf:isup),mask(iinf:isup),avg_r,var,sig_r)
            !for some bad spectra resid=0.0, then
            IF(sig_r<1e-4) THEN 
             sig_r=10.0_dp
            END IF
            sn_var(i)=1.0_dp/sig_r
           ELSE
            sn_var(i)=0.1_dp
           END IF

        ELSE
         sn_var(i)=0.1_dp
        END IF

       END DO

       END SUBROUTINE find_sn_var

!########################
!### this routine finds the closest point of the grid to the stellar parameters given in TGM_
        SUBROUTINE find_prox(TGM_,TGM_prox)
        IMPLICIT NONE
        REAL(DP), DIMENSION(:), INTENT(IN) :: TGM_
        REAL(DP), DIMENSION(3), INTENT(OUT) :: TGM_prox
        REAL(DP), DIMENSION(1) :: prox_dummy

        prox_dummy=temp_grid(MINLOC(ABS(temp_grid-TGM_(1))))
        TGM_prox(1)=prox_dummy(1)
        prox_dummy=logg_grid(MINLOC(ABS(logg_grid-TGM_(2))))
        TGM_prox(2)=prox_dummy(1)
        prox_dummy=met_grid(MINLOC(ABS(met_grid-TGM_(3))))
        TGM_prox(3)=prox_dummy(1)

        IF(TGM_prox(1)>5600.AND.TGM_prox(2)<1.4) THEN
         IF((TGM_prox(1)-5600)/200.>(1.4-TGM_prox(2))/0.4) THEN
          TGM_prox(2)=1.4_dp
         ELSE
           TGM_prox(1)=5600._dp
         END IF
        END IF
        END SUBROUTINE find_prox
!##############################
        !this routine, given one grid point, creates the small grid "gridS"
        SUBROUTINE make_gridS(TGM_prox,temp_gridS,logg_gridS,met_gridS)
        IMPLICIT NONE
        REAL(DP), DIMENSION(3), INTENT(IN) :: TGM_prox
        REAL(DP), DIMENSION(5),INTENT(OUT) :: temp_gridS,logg_gridS,met_gridS

        INTEGER(I1B) :: inf,sup,pos
        
        pos=INT(MINLOC(ABS(temp_gridL-TGM_prox(1)),1),I1B)
        inf=pos-2_I1B
        sup=pos+2_I1B
        temp_gridS=temp_gridL(inf:sup)        
        pos=INT(MINLOC(ABS(logg_gridL-TGM_prox(2)),1),I1B)
        inf=pos-2_I1B
        sup=pos+2_I1B
        logg_gridS=logg_gridL(inf:sup)
        pos=INT(MINLOC(ABS(met_gridL-TGM_prox(3)),1),I1B)
        inf=pos-2_I1B
        sup=pos+2_I1B
        met_gridS=met_gridL(inf:sup)

        END SUBROUTINE make_gridS
!########################
        !this routine finds the closest grid point to the present
        !TGM estimation and set the three flags flag_lim,flag_limS,
        !and flag_move. The first and the second ones report when the TGM estimation is beyond 
        !the border of the large grid and the small grid, respectively.
        !The third one indicate if the estimation is beyond 1.5 steps from the central point.
        SUBROUTINE find_proxS(TGM_,temp_gridS,logg_gridS,met_gridS,&
          TGM_prox,flag_lim,flag_limS,flag_move)
        IMPLICIT NONE
        REAL(DP), DIMENSION(:), INTENT(IN) :: TGM_
        REAL(DP), DIMENSION(5),INTENT(IN) :: temp_gridS,logg_gridS,met_gridS
        REAL(DP), DIMENSION(3),INTENT(INOUT) :: TGM_prox
        LOGICAL,INTENT(OUT) :: flag_lim,flag_limS,flag_move
        REAL(DP), DIMENSION(3) :: TGM_dummy
        !initialize
        flag_lim=.FALSE.
        flag_limS=.FALSE.
        flag_move=.FALSE.
        
        !remember that: gridL means Large grid
        !               grid means the intermediate grid
        !               gridS means Small grid
        !see in data_lab their values

        TGM_dummy(1)=temp_gridS(MINLOC(ABS(temp_gridS-TGM_(1)),1))
        TGM_dummy(2)=logg_gridS(MINLOC(ABS(logg_gridS-TGM_(2)),1))
        TGM_dummy(3)=met_gridS(MINLOC(ABS(met_gridS-TGM_(3)),1))

          CALL find_prox(TGM_dummy,TGM_prox)

        !if TGM_dummy is on the border of the gridL, then exit the program
        IF((TGM_dummy(1)==temp_gridL(1).AND.TGM_(1)<=temp_gridL(1))&
          .OR.(TGM_dummy(1)==temp_gridL(SIZE(temp_gridL)).AND.TGM_(1)>=temp_gridL(SIZE(temp_gridL)))&
          .OR.(TGM_dummy(2)==logg_gridL(1).AND.TGM_(2)<=logg_gridL(1))&
          .OR.(TGM_dummy(2)==logg_gridL(SIZE(logg_gridL)).AND.TGM_(2)>=logg_gridL(SIZE(logg_gridL)))&
          .OR.(TGM_dummy(3)==met_gridL(1).AND.TGM_(3)<=met_gridL(1))&
          .OR.(TGM_dummy(3)==met_gridL(SIZE(met_gridL)).AND.TGM_(3)>=met_gridL(SIZE(met_gridL)))) THEN
          flag_lim=.TRUE.
          RETURN
        END IF


        !if the TGM point is out of the gridS then flag_limS=.TRUE.
        !(which exit the inner loop in the main program)
        IF(TGM_(1)<=temp_gridS(1).OR.TGM_(1)>=temp_gridS(5)&
          .OR.TGM_(2)<=logg_gridS(1).OR.TGM_(2)>=logg_gridS(5)&
          .OR.TGM_(3)==met_gridS(1).OR.TGM_(3)>=met_gridS(5)) THEN
          flag_limS=.TRUE.
          !TGM_dummy can be out of the grid (because gridS can be as well)
          !then take the closest point of grid (we don't want TGM_prox to be out
          ! of the grid)

          RETURN
        END IF



        !if the TGM point moved more than 1 point of the gridS, then flag_move=.TRUE.
        !(which call the routine update_gridS_ll in the main program)
        IF(TGM_(1)<=temp_gridS(2).OR.TGM_(1)>=temp_gridS(4)&
        &.OR.TGM_(2)<=logg_gridS(2).OR.TGM_(2)>=logg_gridS(4)&
        .OR.TGM_(3)<=met_gridS(2).OR.TGM_(3)>=met_gridS(4)) THEN
           flag_move=.TRUE.
        ELSE 
         flag_move=.FALSE.
        END IF

        !if TGM_dummy is outside of the grid, it cannot move forward, then
        IF(TGM_dummy(1)/=TGM_prox(1).OR.TGM_dummy(2)/=TGM_prox(2)&
          &.OR.TGM_dummy(3)/=TGM_prox(3)) THEN
          flag_move=.FALSE.
        END IF
          


        END SUBROUTINE find_proxS
!########################
        !this routine allocates the array ABD and other arrays
        !related to it.
        SUBROUTINE alloc_ABD(ele_in,dim_in,dim_ele,n_ele_symb)
        USE error
        USE share, ONLY: ele2meas,ABD,ele2write,up_ABD,lo_ABD,&
         & ABD_old,residABD,X_abd,alpha_mask
        IMPLICIT NONE
        INTEGER(I2B), INTENT(IN) :: dim_in,n_ele_symb
        INTEGER(I2B), INTENT(OUT) :: dim_ele
        INTEGER(I2B) :: i, j, AllocateStatus
        REAL(DP), DIMENSION(dim_in), INTENT(IN) :: ele_in
        LOGICAL, DIMENSION(dim_in) :: mask
        LOGICAL, DIMENSION(n_ele_symb) :: mask_ele
        !initialize variables
        mask_ele=.FALSE.
        mask=.FALSE.

        DO i=2,n_ele_symb
          WHERE(ABS(ele_in-i)<0.3) mask=.TRUE.
          IF(ANY(mask)) mask_ele(i)=.TRUE.
          mask=.FALSE.
        END DO

        !if there is no Fe lines the program stops
        IF(mask_ele(26).NEQV..TRUE.) CALL stop_msg('no iron lines!')

        dim_ele=INT(COUNT(mask_ele),I1B)

        ALLOCATE(ABD(dim_ele), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate ABD')
        ALLOCATE(lo_ABD(dim_ele), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate lo_ABD')
        ALLOCATE(up_ABD(dim_ele), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate up_ABD')
        ALLOCATE(ABD_old(dim_ele), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate ABD_old')
        ALLOCATE(residABD(dim_ele), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate residABD')
        ALLOCATE(ABD_mask(dim_ele), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate ABD_mask')
        ALLOCATE(alpha_mask(dim_ele), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate alpha_mask')
        ALLOCATE(write_ABD_mask(dim_ele), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate write_ABD_mask')
        ALLOCATE(X_abd(70,dim_ele), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate X_abd')
        !initialize ABD_mask
        ABD_mask=.TRUE.
        alpha_mask=.FALSE.
        
        !initialize ABD (it must come before quick estimation)
        ABD=0.0
        !initializing at -0.1 instead of 0.0 allows 2 abundances loops instead of one for
        !those spectra with met=0.0
        ABD_old=-0.1


        ALLOCATE(ele2meas(dim_ele), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) CALL stop_msg('Not enough memory to allocate ele2meas')
        ele2meas=0


        !allocate Fe in the first position of ele2meas
        ele2meas(1)=26
        !now allocate the other elements
        j=2
        DO i=2,n_ele_symb
          IF(mask_ele(i).AND.i/=26) THEN
          ele2meas(j)=i        
          j=j+1_I1B
          END IF
        END DO

        !check if ele2write has been allocate. if not, allocate it
        !it has to be equal to ele2meas (this happen only when the user
        !did not put the keywork 'ELE2WRITE' in the parameter file
        IF(.NOT.ALLOCATED(ele2write)) THEN
         ALLOCATE(ele2write(dim_ele), STAT = AllocateStatus)
         IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate ele2write')
         ele2write=ele2meas
        END IF

        END SUBROUTINE alloc_ABD
!##############################
        REAL(DP) FUNCTION compute_aFe(ABD, ele2meas, ABD_mask)
        USE stats, ONLY: avg_std_mask
        USE share, ONLY: alpha_mask
        IMPLICIT NONE
        REAL(DP), DIMENSION(:), INTENT(IN) :: ABD
        INTEGER(I2B), DIMENSION(:), INTENT(IN) :: ele2meas
        LOGICAL, DIMENSION(:), INTENT(IN) :: ABD_mask
        REAL(DP) :: alpha_mean,var,sdev

        alpha_mask=.FALSE.

        WHERE((ele2meas==6.OR.ele2meas==8.OR.ele2meas==12.OR.ele2meas==14.OR.ele2meas==20.OR.ele2meas==22)&
        &.AND.ABD_mask) alpha_mask=.TRUE.
        IF (COUNT(alpha_mask)>0) THEN 
           CALL avg_std_mask(ABD,alpha_mask,alpha_mean,var,sdev)
        ELSE
           alpha_mean = 0.
        END IF

        compute_aFe = alpha_mean

        END FUNCTION compute_aFe
!##############################
        !this subroutine creates the logical arrays select_ll_mask
        !with dimension dim_ll to select only the lines
        !that will be measured. Besides, it gives weights=1e-6 to the
        !rejected pixels, included those pixels that are recognised as 
        !cosmic rays by the cosmic_rej subroutine above
        SUBROUTINE select_lines(fun,coeff,pars,flag)
        USE share, ONLY: wave_ll,ele_ll,dim_ll,dim_rej, &
        & wave_rej,rad_rej,select_ll_mask,w_sp,weights,&
        &ew,TGM,ABD,ele2meas,ele2write,cosmic_mask,dimsp,&
        &sn_var!,Ex_inf
        USE space_pars, ONLY: N_w_int,w_inf,w_sup
        USE read_sp_ll, ONLY: add_pix_sp
        USE func_poly, ONLY: poly4_transform

        IMPLICIT NONE
        REAL(DP), EXTERNAL :: fun
        INTEGER(I2B) :: i
        INTEGER(I4B) :: line_pos
        REAL(DP), DIMENSION(4), INTENT(IN) :: pars
        REAL(DP), DIMENSION(:,:), INTENT(IN) :: coeff
        REAL(DP) :: ew_neigh,rad,c,sigma
        REAL(DP), DIMENSION(70) :: X
        LOGICAL, INTENT(IN) :: flag
        LOGICAL, DIMENSION(dim_ll) :: mask_ew,select_ll_mask2
        

        !initialize variables
        !this mask is used to select all the lines used plus the unused inside the rejected intervals
        !the latter will have very low weights
        select_ll_mask=.FALSE.
        !this mask select the line used for measurements only. It will not be used for plots
        !but only to counts which elements are involved in the measurements
        select_ll_mask2=.FALSE.
        ew_neigh=0.
        weights=1.0_dp
        mask_ew=.FALSE.
        c=1._dp+TGM(5)/299792._dp
        sigma=TGM(4)

        !compute the EW of the lines from the polynomial GCOG
        CALL poly4_transform(X,ABD,pars)
        DO i=1,dim_ll
         ew(i)=fun(X,ele_ll(i),coeff(:,i))
        END DO

        !create the ew_mask to sum up the EW of the lines on a radius=sigma from the i-th line
        DO i=1,dim_ll
          mask_ew=.FALSE.
          WHERE(wave_ll>(wave_ll(i)-sigma).AND.wave_ll<(wave_ll(i)+sigma)) mask_ew=.TRUE.
          ew_neigh=SUM(ew,mask_ew)
          !find the pixel at the center of the line
          line_pos=INT(MINLOC(ABS(wave_ll(i)-w_sp),1),I4B) 
          !now create the mask to select only lines for which the sum of the EW of the neighborhood
          ! is greater than a given value
          IF ((ew_neigh/2.5/sigma)>(1./sn_var(line_pos))) THEN
           WHERE((ABS(wave_ll(i)-wave_ll)<3*sigma)&
                 &.AND.(ew>0.1/sn_var(line_pos))) 
                 select_ll_mask=.TRUE.
                 select_ll_mask2=.TRUE.
           END WHERE
          END IF
        END DO


        !set low weights around the rejected lines
        DO i=1,dim_rej
         !determine the radius over which the lines must be rejected
         rad=max(3*sigma,rad_rej(i))
         WHERE(wave_ll>=(wave_rej(i)-rad_rej(i)).AND.wave_ll<=(wave_rej(i)+rad_rej(i))) select_ll_mask=.TRUE.
         !set the weights
         WHERE(ABS(w_sp-wave_rej(i)*c)<=rad) weights=1e-6_dp
        END DO

        !set low weights on the 4 Angstroms tails added to the
        !blue and red extremes of the spectrum 
        WHERE(w_sp<w_sp(add_pix_sp)) weights=1e-6
        WHERE(w_sp>w_sp(dimsp-add_pix_sp)) weights=1e-6
        !set low weights on the rejected intervals, if any
        IF(N_w_int>1) THEN
         DO i=1,INT(N_w_int,I2B)-1_I2B
          WHERE(w_sp>w_sup(i).AND.w_sp<w_inf(i+1)) weights=1e-6
         END DO
        END IF
         
        !set weight=0 where there are cosmic rays
        WHERE(cosmic_mask) weights=1e-6_dp


        !set higher weights of Fe lines
!        DO i=1,dim_ll
!          IF (Ex_inf(i)<20000.) THEN
!              WHERE(w_sp>(wave_ll(i)-sigma).AND.w_sp<(wave_ll(i)+sigma)) weights=5_dp
!          END IF
!        END DO
        
        
        IF(flag) THEN
        ! set the ABD_mask
        ABD_mask=.FALSE.
        !now give the right mask to the other elements
        DO i=1,INT(SIZE(ele2meas,1),I2B)
        mask_ew=.FALSE. !now I use mask_ew as dummy mask
        !I used select_ll_mask2 because it select only the lines employed for measurements
        WHERE(ABS(ele_ll-ele2meas(i))<0.3.AND.select_ll_mask2) mask_ew=.TRUE.
        IF((COUNT(mask_ew)>0).AND.ABD(i)<0.7.AND.ABD(i)>-0.5) THEN
         ABD_mask(i)=.TRUE.
        ELSE
         ABD_mask(i)=.FALSE.
         ABD(i)=0.
        END IF
        END DO
        
        
        write_ABD_mask=.FALSE.
        DO i=1,INT(SIZE(ele2write,1),I2B)
         WHERE (ele2meas==ele2write(i)) write_ABD_mask=.TRUE.
        END DO
        END IF

      END SUBROUTINE select_lines
!#########################
        SUBROUTINE normalize_pars(tgm, tgmx, tgm_mask)
        IMPLICIT NONE
        REAL(DP), DIMENSION(:),INTENT(IN) :: tgm
        REAL(DP), DIMENSION(:),INTENT(INOUT) :: tgmx
        LOGICAL, DIMENSION(:),INTENT(IN) :: tgm_mask
        REAL(DP), DIMENSION(6) :: tgm_local

        tgm_local=tgm
        !normalize Teff and log g
        tgm_local(1)=tgm_local(1)/5000.
        tgm_local(2)=tgm_local(2)/3.
        !pack it again for the minimization routine
        tgmx=PACK(tgm_local,tgm_mask)

        END SUBROUTINE normalize_pars
!#########################
        SUBROUTINE denormalize_pars(tgm_local, tgmx, tgm_mask)
        USE share, ONLY: TGM
        IMPLICIT NONE
        REAL(DP), DIMENSION(:),INTENT(INOUT) :: tgm_local
        REAL(DP), DIMENSION(:),INTENT(IN) :: tgmx
        LOGICAL, DIMENSION(:),INTENT(IN) :: tgm_mask

        !unpack
        tgm_local=UNPACK(tgmx,tgm_mask,TGM)

        IF (tgm_mask(1)) THEN
            tgm_local(1)=tgm_local(1)*5000.
        END IF
        IF (tgm_mask(2)) THEN
            tgm_local(2)=tgm_local(2)*3.
        END IF


        END SUBROUTINE denormalize_pars
!#####################################
END MODULE utils

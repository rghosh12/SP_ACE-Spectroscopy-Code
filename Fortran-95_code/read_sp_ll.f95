!    This module contains the routine that upload the spectrum and
!    the absorption lines expected to be in the wavelength interval
!    given by the user and allocate some fundamental arrays.
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


MODULE read_sp_ll
        USE num_type
        USE space_pars, ONLY: llist,llist_rej,&
           &w_inf,w_sup,N_w_int,flag_rej,flag_alpha,norm_rad,sn_ratio
        USE data_lib, ONLY: w_rej_op,r_rej_op,&
            &w_rej_nlte,r_rej_nlte,w_rej_unknown,r_rej_unknown,w_rej_bad,r_rej_bad
        USE share
        USE error
        IMPLICIT NONE
        LOGICAL, DIMENSION(:), ALLOCATABLE :: wave_mask,llist_mask,llist_mask_w,llist_mask_ll
        INTEGER(I4B) :: add_pix_sp,add_pix_sp_l,add_pix_sp_u

        CONTAINS
        
        SUBROUTINE read_files(name1,name2,name4)        
        CHARACTER(len=120) :: name1,name2,name4
        INTEGER(I1B) :: AllocateStatus
        INTEGER(I4B) :: ierror
        INTEGER(I4B) :: i,dims,sp_dim_trim,ll_dim_trim,w_center_line
        INTEGER(I4B), PARAMETER :: MAX_SIZE_SP=128000
        INTEGER(I4B) :: dim_trim1,dim_trim2,dim_trim3,dim_trim
        REAL(DP), DIMENSION(MAX_SIZE_SP) :: wave,flux
        REAL(DP), DIMENSION(MAX_SIZE_SP) :: w_ll,e_ll,ex_ll
        REAL(DP), DIMENSION(MAX_SIZE_SP) :: col1,col2
        REAL(DP) :: wave_lbound,wave_ubound
        LOGICAL, DIMENSION(:), ALLOCATABLE :: mask

        !initialize variables
        wave=0;flux=1;w_ll=0;e_ll=0;ex_ll=0;i=0
!#############################

!######### read the spectrum        
        write(*,*) 'spectrum: ',name1
        OPEN(unit=10,file=name1,status='OLD',action='READ',iostat=ierror)
        IF(ierror/=0) CALL error_msg(9_I1B,'I cannot open the spectrum! (maybe wrong file?)')

        sp_dim_trim=1_I4B
        DO
          IF(sp_dim_trim>MAX_SIZE_SP-1) THEN
           write(*,*) 'Warning: the spectrum exceeds ',MAX_SIZE_SP,' pixels!'
           EXIT
          END IF
          READ(unit=10,fmt=*,iostat=ierror) wave(sp_dim_trim),flux(sp_dim_trim)
          IF(ierror/=0) EXIT
             IF((wave(sp_dim_trim)>=6860.).AND.(wave(sp_dim_trim)<=8400.)) THEN
              CYCLE
             END IF
            DO i=1,N_w_int
             IF((wave(sp_dim_trim)>=w_inf(i)).AND.(wave(sp_dim_trim)<=w_sup(i))) THEN
              sp_dim_trim=sp_dim_trim+1_I4B
             ELSEIF (wave(sp_dim_trim)>w_sup(i).AND.wave(sp_dim_trim)<(MIN(w_sup(N_w_int),w_sup(i)+4.0_dp))) THEN
              flux(sp_dim_trim)=1.0_dp
              sp_dim_trim=sp_dim_trim+1_I4B
             ELSEIF (wave(sp_dim_trim)>(MAX(w_inf(1),w_inf(i)-4.0_dp)).AND.wave(sp_dim_trim)<w_inf(i)) THEN
              flux(sp_dim_trim)=1.0_dp
              sp_dim_trim=sp_dim_trim+1_I4B
             END IF
            END DO
        END DO
        CLOSE(unit=10)

        sp_dim_trim=sp_dim_trim-1_I4B
        IF(sp_dim_trim<1) CALL stop_msg('No readable spectrum or wrong wave_lims.')
        !compute how many pixels to add to the left and right of the spectrum
        !in order to have 4 angstrom more each side
        !use_add_pix_sp_l and add_pix_sp_u because they can be different
        !when the spectrum has variable dispersion 
        add_pix_sp_l=NINT(4./(wave(2)-wave(1)),I4B)
        add_pix_sp_u=NINT(4./(wave(sp_dim_trim)-wave(sp_dim_trim-1)),I4B)
        add_pix_sp=MAX(add_pix_sp_l,add_pix_sp_u,5_I4B)
!###############################

!######### read the linelist
        OPEN(unit=10,file=name2,status='OLD',action='READ',iostat=ierror)
        IF(ierror/=0) CALL stop_msg('I cannot open the line list')

        ll_dim_trim=1_I4B
        DO
          READ(unit=10,fmt=*,iostat=ierror) w_ll(ll_dim_trim),e_ll(ll_dim_trim), ex_ll(i)
          IF(ierror/=0) THEN
           ll_dim_trim=ll_dim_trim-1_I4B
           EXIT
          END IF
            ll_dim_trim=ll_dim_trim+1_I4B
        END DO
        CLOSE(unit=10)




!######### select the upper and lower limits of the spectrum wavelength range
!######### according to the limits given by the user and the spectrum

        ALLOCATE(wave_mask(sp_dim_trim), STAT = AllocateStatus)
           IF (AllocateStatus /= 0) CALL stop_msg('Not enough memory to allocate wave_mask')
        wave_mask=.FALSE.

        wave_lbound=MAXVAL((/w_inf(1),w_ll(1),wave(1)/))        
        wave_ubound=MINVAL((/w_sup(N_w_int),w_ll(ll_dim_trim)/))
        WHERE(wave(1:sp_dim_trim)>wave_lbound.AND.wave(1:sp_dim_trim)<wave_ubound) wave_mask=.TRUE.
        dimsp=INT(COUNT(wave_mask),I4B)


        ALLOCATE(w_sp(dimsp+add_pix_sp*2), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate w_sp')
        ALLOCATE(f_sp(dimsp+add_pix_sp*2), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate f_sp')
        ALLOCATE(f_sp_norm(dimsp+add_pix_sp*2), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate f_cont')
        ALLOCATE(f_model(dimsp+add_pix_sp*2), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate f_model')
        ALLOCATE(cont(dimsp+add_pix_sp*2), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate cont')
        ALLOCATE(cont0(dimsp+add_pix_sp*2), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate cont0')
        ALLOCATE(weights(dimsp+add_pix_sp*2), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate weights')
        ALLOCATE(cosmic_mask(dimsp+add_pix_sp*2), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate cosmic_mask')
        ALLOCATE(sn_var(dimsp+add_pix_sp*2), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate sn_var')


        !initialize w_sp
        w_sp(add_pix_sp+1:dimsp+add_pix_sp)=PACK(wave(1:sp_dim_trim),wave_mask)
        !add the right values to the left and right of the spectrum
        DO i=1,add_pix_sp
        w_sp(add_pix_sp+1-i)=w_sp(add_pix_sp+1)-(w_sp(add_pix_sp+2)-w_sp(add_pix_sp+1))*i
        w_sp(dimsp+add_pix_sp+i)=w_sp(dimsp+add_pix_sp)+(w_sp(dimsp+add_pix_sp)-w_sp(dimsp+add_pix_sp-1))*i
        END DO
        !initialize f_sp
        f_sp=1.0_dp
        f_sp(add_pix_sp+1:dimsp+add_pix_sp)=PACK(flux(1:sp_dim_trim),wave_mask)
        !initialize the other arrays
        f_sp_norm=1.0_dp
        f_model=1.0_dp
        cont=1.0_dp
        weights=1.0_dp
        cosmic_mask=.FALSE.
        sn_var=sn_ratio

        !re-define the spectrum dimension
        dimsp=INT(SIZE(w_sp,1),I4B)
        !define rad_pix to be used in fit_cont.f95
        ALLOCATE(rad_pix(dimsp+add_pix_sp*2), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate rad_pix')

        IF(norm_rad<(w_sp(dimsp)-w_sp(1))) THEN
         IF(norm_rad<5.) norm_rad=5.
         rad_pix(2:dimsp)=NINT(norm_rad/(w_sp(2:dimsp)-w_sp(1:dimsp-1)),I4B)
         rad_pix(1)=rad_pix(2)
        ELSE
         rad_pix=dimsp
        END IF
        !the following 7 lines fix the value at the pixels corresponding to w_inf(i)-4.0,
        !where, if N_w_int>1, the value rad_pix would be wrong because w_sp(2:dimsp)-w_sp(1:dimsp-1) would
        !be too large (because it measures the excluded interval)
        IF(N_w_int>1) THEN
        DO i=2,N_w_int
         w_center_line=INT(MINLOC(ABS(w_sp-(w_inf(i)-4.0)),1),I4B)
         rad_pix(w_center_line)=rad_pix(w_center_line+1)
        END DO
        END IF
!######### select the upper and lower limits of the line list
!######### according to the limits given by the user and the spectrum
        ALLOCATE(llist_mask(ll_dim_trim), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate llist_mask')
        llist_mask=.FALSE.
        ALLOCATE(llist_mask_w(ll_dim_trim), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate llist_mask_w')
        llist_mask=.FALSE.
        ALLOCATE(llist_mask_ll(ll_dim_trim), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate llist_mask_ll')
        llist_mask=.FALSE.


        DO i=1,N_w_int
          WHERE(w_ll(1:ll_dim_trim)>w_inf(i).AND.w_ll(1:ll_dim_trim)<w_sup(i)) llist_mask_w=.TRUE.
          WHERE(w_ll(1:ll_dim_trim)>wave(1).AND.w_ll(1:ll_dim_trim)<wave(sp_dim_trim)) llist_mask_ll=.TRUE.
        END DO
        WHERE(llist_mask_w.AND.llist_mask_ll) llist_mask=.TRUE.

        dim_ll=INT(COUNT(llist_mask),I2B)

        ALLOCATE(wave_ll(dim_ll), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate wave_ll')
        ALLOCATE(wave_center_ll(dim_ll), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate wave_center_ll')
        ALLOCATE(ele_ll(dim_ll), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate ele_ll')
        ALLOCATE(Ex_inf(dim_ll), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate Ex_inf')
        ALLOCATE(ew(dim_ll), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate ew')
        ALLOCATE(disp(dim_ll), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) CALL stop_msg('Not enough memory to allocate disp')
        ALLOCATE(coeff_4deg(70,dim_ll), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) CALL stop_msg('Not enough memory to allocate coeff_4deg')
        ALLOCATE(coeff_4deg_quick(84,dim_ll), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) CALL stop_msg('Not enough memory to allocate coeff_4deg_quick')

        wave_ll=PACK(w_ll(1:ll_dim_trim),llist_mask)
        ele_ll=PACK(e_ll(1:ll_dim_trim),llist_mask)
        Ex_inf=PACK(ex_ll(1:ll_dim_trim),llist_mask)

        !define dispersion
        DO i=1,dim_ll
         w_center_line=INT(MINLOC(ABS(w_sp-wave_ll(i)),1),I4B)
         disp(i)=w_sp(w_center_line+1)-w_sp(w_center_line)
         wave_center_ll(i)=w_center_line
        END DO

        !for the molecule use the dummy atomic values
        WHERE(ele_ll==106.0) ele_ll=95.0
        WHERE(ele_ll==107.0) ele_ll=96.0
        WHERE(ele_ll==112.0) ele_ll=97.0
        WHERE(ele_ll==114.0) ele_ll=98.0
        WHERE(ele_ll==606.0) ele_ll=99.0
        WHERE(ele_ll==607.0) ele_ll=100.0
        WHERE(ele_ll==814.0) ele_ll=101.0
        !if flag_alpha=true, then the alpha elements Mg, Si, Ca, Ti 
        !become one element (fictitious atomic number=94), 
        !while all the other elements become Fe
        IF(flag_alpha) THEN
         WHERE(ABS(ele_ll-12)<0.3.OR.&
            &ABS(ele_ll-14)<0.3.OR.ABS(ele_ll-20)<0.3.OR.&
            &ABS(ele_ll-22)<0.3) ele_ll=94.0
         WHERE(ABS(ele_ll-6)>0.3.AND.ABS(ele_ll-7)>0.3.AND.ABS(ele_ll-8)>0.3.AND.&
         &ABS(ele_ll-94)>0.3.AND.ele_ll>94.0) ele_ll=26.0
        END IF
!###########################################

!######### here allocate the array flag_lines to label the lines 

        ALLOCATE(flag_lines(dim_ll), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate flag_lines')

        !initialize. All the "normal" lines are labelled with flag_lines=0
        flag_lines=0
        !the Ha and Hb lines are labelled with 1
        WHERE(ele_ll==1.0.and.wave_ll<6563.) 
           flag_lines=1
        END WHERE
        !the Na lines are labelled with 2
        WHERE(ele_ll==11.0.AND.wave_ll>5889.AND.wave_ll<5896.) 
           flag_lines=2
        END WHERE

! this was experimental and not used now
!        WHERE(ele_ll==12.0.and.wave_ll<5184..AND.wave_ll>5160) 
!           flag_lines=2
!        END WHERE
!
!        WHERE(ele_ll==24.0.AND.(wave_ll==5204.498.OR.&
!                               &wave_ll==5206.023.OR.&
!                               &wave_ll==5208.409))
!           flag_lines=3
!        END WHERE

!######### prepare the rejected list taken from the llist in data_lib
        !put the lines rejected because bad opacity correction
        dim_trim=INT(SIZE(w_rej_op,1),I2B)
        col1(1:dim_trim)=w_rej_op
        col2(1:dim_trim)=r_rej_op
        dims=dim_trim
        !put the lines rejected because NLTE effect neglected
        dim_trim1=INT(SIZE(w_rej_nlte,1),I2B)
        col1(dims+1:dims+dim_trim1)=w_rej_nlte
        col2(dims+1:dims+dim_trim1)=r_rej_nlte
        dims=dims+dim_trim1
        !put the lines rejected because unknown
        dim_trim2=INT(SIZE(w_rej_unknown,1),I2B)
        col1(dims+1:dims+dim_trim2)=w_rej_unknown
        col2(dims+1:dims+dim_trim2)=r_rej_unknown    
        dims=dims+dim_trim2
        !put the lines rejected because they fit bad for unknown reasons
        dim_trim3=INT(SIZE(w_rej_bad,1),I2B)
        col1(dims+1:dims+dim_trim3)=w_rej_bad
        col2(dims+1:dims+dim_trim3)=r_rej_bad    
        dims=dims+dim_trim3
!######### read the user's list of lines to reject, if any
        IF(flag_rej) THEN
         OPEN(unit=10,file=name4,status='OLD',action='READ',iostat=ierror)
         IF(ierror/=0) CALL stop_msg('I cannot open the reject line list')

         dims=dims+1_I2B

         DO
           READ(unit=10,fmt=*,iostat=ierror) col1(dims),col2(dims)
           IF(ierror/=0) THEN
            dims=dims-1_I2B
            EXIT
           END IF
           dims=dims+1_I2B
         END DO
         CLOSE(unit=10)
        END IF

         !sort the array
         CALL sorting(col1(1:dims),col2(1:dims),dims)
!######### select the upper and lower limits of the rejected line list
!######### according to the limits given by the user and the spectrum

        ALLOCATE(mask(dims), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate mask(dim_rej)')
        mask=.FALSE.

        IF(N_w_int==1) THEN
         wave_lbound=MAXVAL((/w_inf(1),col1(1),wave(1)/))        
         wave_ubound=MINVAL((/w_sup(1),col1(dims),wave(sp_dim_trim)/))
         WHERE((col1(1:dims)+col2(1:dims))>=wave_lbound.AND.(col1(1:dims)-col2(1:dims))<=wave_ubound) mask=.TRUE.
        ELSE 
         DO i=1,N_w_int
          IF(i==1) THEN
           wave_lbound=MAXVAL((/w_inf(i),col1(1),wave(1)/))        
           wave_ubound=MINVAL((/w_sup(i),col1(dims),wave(sp_dim_trim)/))
           WHERE((col1(1:dims)+col2(1:dims))>=wave_lbound.AND.(col1(1:dims)-col2(1:dims))<=wave_ubound) mask=.TRUE.
          ELSEIF(i==N_w_int) THEN
           wave_lbound=MAXVAL((/w_inf(i),col1(1)/))        
           wave_ubound=MINVAL((/w_sup(i),col1(dims),wave(sp_dim_trim)/))
           WHERE((col1(1:dims)+col2(1:dims))>=wave_lbound.AND.(col1(1:dims)-col2(1:dims))<=wave_ubound) mask=.TRUE.
          ELSE
           wave_lbound=w_inf(i)        
           wave_ubound=w_sup(i)
           WHERE((col1(1:dims)+col2(1:dims))>=wave_lbound.AND.(col1(1:dims)-col2(1:dims))<=wave_ubound) mask=.TRUE.
          END IF
         END DO
        END IF
        dim_rej=INT(COUNT(mask),I2B)


        ALLOCATE(wave_rej(dim_rej), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate wave_rej')
        ALLOCATE(rad_rej(dim_rej), STAT = AllocateStatus)
        IF (AllocateStatus/=0) CALL stop_msg('Not enough memory to allocate rad_rej')

        wave_rej=PACK(col1(1:dims),mask)
        rad_rej=PACK(col2(1:dims),mask)

        DEALLOCATE(mask)
        END SUBROUTINE read_files        
!##################
        SUBROUTINE sorting(x,x1,dim)
        IMPLICIT NONE
        INTEGER (I4B), INTENT(IN) :: dim
        REAL(DP), DIMENSION(dim), INTENT(INOUT) :: x,x1
        REAL(DP), DIMENSION(dim) :: x_sort,x1_sort
        INTEGER (I4B) :: a,i
        LOGICAL, DIMENSION(dim) :: mask

        mask=.TRUE.

        DO i=1,INT(SIZE(x),I2B)
         a=INT(MINLOC(x,1,mask),I2B)
         mask(a)=.FALSE.
         x_sort(i)=x(a)
         x1_sort(i)=x1(a)
        END DO

        x=x_sort
        x1=x1_sort

        END SUBROUTINE sorting
!##################

END MODULE read_sp_ll

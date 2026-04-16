!    
!    This module contains subroutines that estimate uncertains
!    in the stellar parameters and chemical abundances.
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

MODULE uncertains2
        USE num_type
        USE share, ONLY: up_ABD,lo_ABD,ele2meas
        USE utils, ONLY: ABD_mask,write_ABD_mask
        IMPLICIT NONE
        REAL (DP), DIMENSION(3) :: up_TGM,lo_TGM
        INTEGER(I1B),PARAMETER :: grado=2
        INTEGER(I2B),PARAMETER :: dim_coeff=10! 10 for grado=2, 20 for grado=3, 35 for grado=4
        REAL (DP), DIMENSION(3), PARAMETER :: prob=(/3.53,6.25,8.02/)!for 3 deg of freedom, correspond to %68.3,%90,%95.4
        REAL(DP), DIMENSION(27,10) :: matrixA,matrixA_orig
        REAL(DP), DIMENSION(27,1) :: vecB,vecB_orig
        REAL(DP), DIMENSION(27,3) :: TGM_test,TGM_test_orig
        REAL (DP) :: chisq_target
        INTEGER (I2B) :: pos=1,n_eq
        INTEGER(I1B) :: direc
        INTEGER (I2B),SAVE :: count_point
        LOGICAL :: flag_TGM=.TRUE.,flag_write

      INTERFACE LA_ILAENV

      FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
         INTEGER :: ILAENV
         CHARACTER(LEN=*), INTENT(IN) :: NAME, OPTS
         INTEGER, INTENT(IN) :: ISPEC, N1, N2, N3, N4
      END FUNCTION ILAENV

      END INTERFACE

      INTERFACE LA_GELSS
       SUBROUTINE DGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,   &
     &                    WORK, LWORK, INFO )
!         USE LA_PRECISION, ONLY: WP => DP
         USE num_type
         INTEGER, INTENT(IN) :: NRHS, M, N, LDA, LDB, LWORK
         INTEGER, INTENT(OUT) :: INFO, RANK
         REAL(DP), INTENT(IN) :: RCOND
         REAL(DP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         REAL(DP), INTENT(OUT) :: S(*)
         REAL(DP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DGELSS
       END INTERFACE

        CONTAINS
!###############################
SUBROUTINE normalize3D(tgm)
        IMPLICIT NONE
        REAL(DP), DIMENSION(:),INTENT(INOUT) :: tgm

        tgm(1)=tgm(1)/5000.
        tgm(2)=tgm(2)/3.
END SUBROUTINE normalize3D
!###############################
SUBROUTINE denormalize3D(tgm)
        IMPLICIT NONE
        REAL(DP), DIMENSION(:),INTENT(INOUT) :: tgm

        tgm(1)=tgm(1)*5000.
        tgm(2)=tgm(2)*3.
END SUBROUTINE denormalize3D

!###############################
SUBROUTINE eval_poly(tgm,coeff,chisq_out)
       IMPLICIT NONE
       REAL(DP), DIMENSION(3), INTENT(IN) :: tgm
       REAL(DP), DIMENSION(dim_coeff),INTENT(IN) :: coeff
       REAL(DP), INTENT(OUT) :: chisq_out
       INTEGER(I1B) :: i1,i2,i3,count
       REAL(DP) :: a,b,c,chisq_poly

       chisq_poly=0
       count=0

        DO i1=0,grado
         a=(tgm(1)**i1)
         DO i2=0,grado-i1
          b=(tgm(2)**i2)
          DO i3=0,grado-i2-i1
           c=(tgm(3)**i3)
             count=count+1_I1B
             chisq_poly=chisq_poly+coeff(count)*a*b*c
          END DO
         END DO
        END DO                                 


        chisq_out=chisq_poly

END SUBROUTINE eval_poly
!###############################
SUBROUTINE coeff_poly_find(TGM_ini,step,coeff)
        IMPLICIT NONE
        REAL(DP), DIMENSION(3),INTENT(IN) :: TGM_ini,step
        REAL(DP), DIMENSION(dim_coeff),INTENT(INOUT) :: coeff
        REAL(DP), DIMENSION(n_eq,10) :: matrixA_int
        REAL(DP), DIMENSION(n_eq,1) :: vecB_int
        REAL(DP), DIMENSION(n_eq,3) :: tgm_grid
        REAL(DP), DIMENSION(dim_coeff) :: s
        REAL(DP) :: rcond=1e-18
        INTEGER :: info,rank

        vecB_int=0.
        !this routine sets tgm_grid and also vecB
        CALL define_grid(TGM_ini,step,tgm_grid)
        TGM_test(1:n_eq,:)=tgm_grid        
        CALL prepare_matrix(1_I2B,n_eq)
        !fit the chi2 with a quadratic function
        matrixA_int=matrixA(1:n_eq,:)
!        write(*,*) 'vecB',vecB(1:n_eq,:)
        vecB_int=vecB(1:n_eq,:)
        CALL DGELSS_F95( matrixA_int, vecB_int, rank,s,rcond,info )
        coeff=vecB_int(1:dim_coeff,1)
!        IF(info/=0) write(*,*) 'error',coeff
!        IF(info/=0) READ'(I1)', info
!        write(*,*) 'rank',rank        

END SUBROUTINE coeff_poly_find
!###############################
SUBROUTINE check_up_limits(step_loc,TGM_loc,flag_out_limit)
        USE data_lib, ONLY: temp_gridL,logg_gridL,met_gridL
        IMPLICIT NONE
        REAL(DP), DIMENSION(3), INTENT(INOUT) :: TGM_loc
        REAL(DP), DIMENSION(3), INTENT(IN) :: step_loc
        LOGICAL, INTENT(OUT) :: flag_out_limit

        CALL denormalize3D(TGM_loc)

        flag_out_limit=.FALSE.

        IF(TGM_loc(1)>temp_gridL(SIZE(temp_gridL))) THEN
         TGM_loc(1)=MIN(TGM_loc(1),temp_gridL(SIZE(temp_gridL)))
         TGM_loc(1)=TGM_loc(1)-step_loc(1)
         flag_out_limit=.TRUE.
        END IF

        IF(TGM_loc(2)>logg_gridL(SIZE(logg_gridL))) THEN
         TGM_loc(2)=MIN(TGM_loc(2),logg_gridL(SIZE(logg_gridL)))
         TGM_loc(2)=TGM_loc(2)-step_loc(2)
         flag_out_limit=.TRUE.
        END IF

        IF(flag_TGM) THEN
         IF(TGM_loc(3)>met_gridL(SIZE(met_gridL))) THEN
          TGM_loc(3)=MIN(TGM_loc(3),met_gridL(SIZE(met_gridL)))
          TGM_loc(3)=TGM_loc(3)-step_loc(3)
         flag_out_limit=.TRUE.
         END IF
        ELSE
         IF(TGM_loc(3)>0.8) THEN
          TGM_loc(3)=MIN(TGM_loc(3),0.8)
          TGM_loc(3)=TGM_loc(3)-step_loc(3)
         flag_out_limit=.TRUE.
         END IF
        END IF

        CALL normalize3D(TGM_loc)
        
END SUBROUTINE check_up_limits
!###############################
SUBROUTINE check_lo_limits(step_loc,TGM_loc,flag_out_limit)
        USE data_lib, ONLY: temp_gridL,logg_gridL,met_gridL
        IMPLICIT NONE
        REAL(DP), DIMENSION(3), INTENT(INOUT) :: TGM_loc
        REAL(DP), DIMENSION(3), INTENT(IN) :: step_loc
        LOGICAL, INTENT(OUT) :: flag_out_limit

        CALL denormalize3D(TGM_loc)

        flag_out_limit=.FALSE.

        IF(TGM_loc(1)<temp_gridL(1)) THEN
         TGM_loc(1)=MAX(TGM_loc(1),temp_gridL(1))
         TGM_loc(1)=TGM_loc(1)+step_loc(1)
         flag_out_limit=.TRUE.
        END IF

        IF(TGM_loc(2)<logg_gridL(1)) THEN
         TGM_loc(2)=MAX(TGM_loc(2),logg_gridL(1))
         TGM_loc(2)=TGM_loc(2)+step_loc(2)
         flag_out_limit=.TRUE.
        END IF

        IF(flag_TGM) THEN
         IF(TGM_loc(3)<met_gridL(1)) THEN
          TGM_loc(3)=MAX(TGM_loc(3),met_gridL(1))
          TGM_loc(3)=TGM_loc(3)+step_loc(3)
         flag_out_limit=.TRUE.
         END IF
        ELSE
         IF(TGM_loc(3)<-0.4) THEN  
          TGM_loc(3)=MAX(TGM_loc(3),-0.4)
          TGM_loc(3)=TGM_loc(3)+step_loc(3)
         flag_out_limit=.TRUE.
         END IF
        END IF

        CALL normalize3D(TGM_loc)

END SUBROUTINE check_lo_limits
!###############################
!SUBROUTINE add_test_points(TGM_ini,step,coeff)
!        IMPLICIT NONE
!        REAL(DP), DIMENSION(3),INTENT(IN) :: TGM_ini,step
!        REAL(DP), DIMENSION(10),INTENT(INOUT) :: coeff
!        REAL(DP), DIMENSION(500,10) :: matrixA_int
!        REAL(DP), DIMENSION(500,1) :: vecB_int
!        INTEGER (I1B) :: i,j,k
!        INTEGER (I2B) :: count
!        INTEGER(I2B), SAVE :: inf_eq,sup_eq,max_eq
!        REAL(DP), DIMENSION(dim_coeff) :: s
!        REAL(DP) :: rcond=1e-18
!        INTEGER :: info,rank
!
!        max_eq=81
!
!        IF(n_eq>=27.AND.n_eq<max_eq) THEN
!         count=n_eq        
!         inf_eq=n_eq+1_I2B
!         sup_eq=n_eq+27_I2B
!         n_eq=sup_eq
!        count=inf_eq-1_I2B
!        ELSE IF(n_eq==max_eq.AND.sup_eq==max_eq) THEN
!         inf_eq=27_I2B
!         sup_eq=inf_eq+27_I2B
!        count=inf_eq-1_I2B
!        ELSE IF(n_eq==max_eq.AND.sup_eq<max_eq) THEN
!         inf_eq=inf_eq+1_I2B       
!         sup_eq=sup_eq+27_I2B
!        count=inf_eq-1_I2B
!        END IF
! write(*,*) 'n_eq ini',n_eq,inf_eq,sup_eq
!
!        DO k=-1,1
!        DO j=-1,1
!        DO i=-1,1
!         count=count+1_I1B
!         TGM_test(count,1)=TGM_ini(1)+step(1)*i
!         TGM_test(count,2)=TGM_ini(2)+step(2)*j
!         TGM_test(count,3)=TGM_ini(3)+step(3)*k
!         IF (flag_TGM) CALL chisq_tgm_e(TGM_test(count,1:3),vecB(count,1))
!!         write(*,*) 'add', count,TGM_test(count,1:3),vecB(count,1)
!         IF (.NOT.flag_TGM) CALL chisq_abd_e(TGM_test(count,1:3),vecB(count,1))
!        ENDDO
!        ENDDO
!        ENDDO
!
!
!        inf_eq=count
!        CALL prepare_matrix(inf_eq,sup_eq)
!        !fit the chi2 with a quadratic function
!!        write(*,*) 'n_eq',n_eq,inf_eq,sup_eq
!        matrixA_int(1:n_eq,:)=matrixA(1:n_eq,:)
!        vecB_int(1:n_eq,1)=vecB(1:n_eq,1)
!        CALL DGELSS_F95( matrixA_int(1:n_eq,:), vecB_int(1:n_eq,:), rank,s,rcond,info )
!        coeff=vecB_int(1:dim_coeff,1)
!!        IF(info/=0) write(*,*) 'error',coeff
!!        IF(info/=0) READ'(I1)', info
!!        write(*,*) 'rank',rank
!
!END SUBROUTINE add_test_points
!###############################
!SUBROUTINE add_one_point(TGM_ini,coeff)
!        IMPLICIT NONE
!        REAL(DP), DIMENSION(3),INTENT(IN) :: TGM_ini
!        REAL(DP), DIMENSION(10),INTENT(INOUT) :: coeff
!        REAL(DP), DIMENSION(500,10) :: matrixA_int
!        REAL(DP), DIMENSION(500,1) :: vecB_int
!        INTEGER (I2B),SAVE :: count
!        REAL(DP), DIMENSION(dim_coeff) :: s
!        REAL(DP) :: rcond=1e-60
!        INTEGER :: info,rank
!
!         !put a new point at the beginning of the array TGM_test
!         count=n_eq
!         count=count+1_I2B
!         n_eq=count
!         write(*,*) 'n_eq',n_eq
!
!         TGM_test(count,1)=TGM_ini(1)
!         TGM_test(count,2)=TGM_ini(2)
!         TGM_test(count,3)=TGM_ini(3)
!
!         IF (flag_TGM) CALL chisq_tgm_e(TGM_test(count,1:3),vecB(count,1))
!         IF (.NOT.flag_TGM) CALL chisq_abd_e(TGM_test(count,1:3),vecB(count,1))
!
!
!        CALL prepare_matrix(count,count)
!        !fit the chi2 with a quadratic function
!!        write(*,*) 'n_eq',n_eq,inf_eq,sup_eq
!        matrixA_int(1:n_eq,:)=matrixA(1:n_eq,:)
!        vecB_int(1:n_eq,1)=vecB(1:n_eq,1)
!        CALL DGELSS_F95( matrixA_int(1:n_eq,:), vecB_int(1:n_eq,:), rank,s,rcond,info )
!        coeff=vecB_int(1:dim_coeff,1)
!!        IF(info/=0) write(*,*) 'error'
!!        IF(info/=0) READ'(I1)', info
!!        write(*,*) 'rank',rank
!
!END SUBROUTINE add_one_point
!####################################
SUBROUTINE stick_one_point(TGM_ini,coeff)
        IMPLICIT NONE
        REAL(DP), DIMENSION(3),INTENT(IN) :: TGM_ini
        REAL(DP), DIMENSION(10),INTENT(INOUT) :: coeff
        REAL(DP), DIMENSION(n_eq,10) :: matrixA_int
        REAL(DP), DIMENSION(n_eq,1) :: vecB_int

        REAL(DP), DIMENSION(dim_coeff) :: s
        REAL(DP) :: rcond=1e-60
        INTEGER :: info,rank

         !put a new point at the beginning of the array TGM_test
         count_point=count_point+1_I2B

         TGM_test(count_point,1)=TGM_ini(1)
         TGM_test(count_point,2)=TGM_ini(2)
         TGM_test(count_point,3)=TGM_ini(3)

         IF (flag_TGM) CALL chisq_tgm_e(TGM_test(count_point,1:3),vecB(count_point,1))
         IF (.NOT.flag_TGM) CALL chisq_abd_e(TGM_test(count_point,1:3),vecB(count_point,1))


        CALL prepare_matrix(count_point,count_point)
        !fit the chi2 with a quadratic function
!        write(*,*) 'n_eq',n_eq,inf_eq,sup_eq
        matrixA_int(1:n_eq,:)=matrixA(1:n_eq,:)
        vecB_int(1:n_eq,1)=vecB(1:n_eq,1)
        CALL DGELSS_F95( matrixA_int(1:n_eq,:), vecB_int(1:n_eq,:), rank,s,rcond,info )
        coeff=vecB_int(1:dim_coeff,1)
!        IF(info/=0) write(*,*) 'error'
!        IF(info/=0) READ'(I1)', info
!        write(*,*) 'rank',rank

END SUBROUTINE stick_one_point
!####################################
SUBROUTINE TGM_errors (chisq_best)
        USE share, ONLY: TGM,ABD
        IMPLICIT NONE
        REAL(DP), INTENT(IN) :: chisq_best
        REAL(DP), DIMENSION(dim_coeff) :: coeff
        REAL(DP), DIMENSION(3) :: step,TGM_temp
        REAL(DP) :: chisq_tol
        INTEGER(I1B) :: i

        chisq_target=chisq_best+prob(1)
        flag_write=.FALSE.


        step=(/5.,0.02,0.01/)
        CALL normalize3D(step)
        chisq_tol=0.2

!########### start TGM errors estmates #########
!###############################################
        !if flag_TGM=true we estimate TGM, otherwise ABD
        flag_TGM=.TRUE.

        !initialize TGM_temp
        TGM_temp=TGM(1:3)
        CALL normalize3D(TGM_temp)
        n_eq=27
        CALL coeff_poly_find(TGM_temp,step,coeff)

!       CALL write_poly_surf(TGM(1:3),chisq_best,coeff)
        !estimate the poly up_TGM errors
        pos=1
        count_point=0
        CALL est_errors_p(coeff,chisq_best,chisq_tol,step,TGM_temp)
        pos=2
        count_point=0
        CALL est_errors_p(coeff,chisq_best,chisq_tol,step,TGM_temp)
        pos=3
        count_point=0
        CALL est_errors_p(coeff,chisq_best,chisq_tol,step,TGM_temp)

        !initialize the coeff, TGM_temp and n_eq at the best_chisq
        TGM_temp=TGM(1:3)
        CALL normalize3D(TGM_temp)
        n_eq=27
        CALL coeff_poly_find(TGM_temp,-step,coeff)
        !estimate the poly lo_TGM errors
        pos=1
        count_point=0
        CALL est_errors_m(coeff,chisq_best,chisq_tol,step,TGM_temp)
        pos=2
        count_point=0
        CALL est_errors_m(coeff,chisq_best,chisq_tol,step,TGM_temp)
        pos=3
        count_point=0
        CALL est_errors_m(coeff,chisq_best,chisq_tol,step,TGM_temp)

!########### start ABD errors estimates #########
!###############################################

!        write(*,*) 'compute the ABD poly errors'
        !if flag_TGM=true we estimate TGM, otherwise ABD
        flag_TGM=.FALSE.

        DO i=1,INT(SIZE(ABD,1),I1B)
        step=(/0.001,0.03,0.01/)
        pos=i
!        write(*,*) i,ele2meas(i)
         IF(write_ABD_mask(i).AND.ABD_mask(i).AND.ABD(i)>-0.599.AND.ABD(i)<0.799) THEN
        !initialize TGM_temp
          TGM_temp(1:2)=TGM(1:2)
          TGM_temp(3)=ABD(pos)
          CALL normalize3D(TGM_temp)
          !fit the chi2 and estimates the coeff
          n_eq=27
          CALL coeff_poly_find(TGM_temp,step,coeff)
          !estimate the poly up_ABD(pos) errors
          count_point=0
          CALL est_errors_p(coeff,chisq_best,chisq_tol,step,TGM_temp)
          !estimate the poly lo_ABD(pos) errors
          CALL coeff_poly_find(TGM_temp,-step,coeff)
          count_point=0
          CALL est_errors_m(coeff,chisq_best,chisq_tol,step,TGM_temp)
!          write(*,*) pos,lo_ABD(pos),ABD(pos),up_ABD(pos)
         ELSE
          up_ABD(pos)=10.
          lo_ABD(pos)=-10.
         END IF
        END DO

END SUBROUTINE TGM_errors
!###############################
SUBROUTINE prepare_matrix(m,n)
       IMPLICIT NONE
       INTEGER(I2B), INTENT(IN) :: n,m
       INTEGER(I2B) :: i
       INTEGER(I1B) :: i1,i2,i3,count
       REAL(DP) :: a,b,c,chisq_out

        DO i=m,n
         count=0

         DO i1=0,grado
          a=(TGM_test(i,1)**i1)
          DO i2=0,grado-i1
           b=(TGM_test(i,2)**i2)
           DO i3=0,grado-i2-i1
           c=(TGM_test(i,3)**i3)
            count=count+1_I1B
            matrixA(i,count)=a*b*c
           END DO
          END DO
         END DO                                 
         IF(flag_TGM) THEN
          CALL chisq_tgm_e(TGM_test(i,:),chisq_out)
         ELSE
          CALL chisq_abd_e(TGM_test(i,:),chisq_out)
         END IF
         vecB(i,1)=chisq_out
        END DO

END SUBROUTINE prepare_matrix
!############################
SUBROUTINE est_errors_p(coeff,chisq_best,chisq_tol,step,TGM_p)
        IMPLICIT NONE
        REAL(DP), DIMENSION(dim_coeff),INTENT(IN) :: coeff
        REAL(DP), INTENT(IN) :: chisq_best,chisq_tol
        REAL(DP), DIMENSION(3),INTENT(IN) :: step
        REAL(DP), DIMENSION(3),INTENT(IN) :: TGM_p
        REAL(DP) :: chisq,chisq_poly,chisq_best_loc
        REAL(DP), DIMENSION(dim_coeff) :: coeff_loc
        REAL(DP), DIMENSION(3) :: TGM_loc,step_loc,TGM_best_loc
        LOGICAL :: flag_out_limit,flag_solved,flag_D1
        INTEGER(I1B) :: i,j
        INTEGER(I2B) :: pos_loc

         flag_out_limit=.FALSE.
         flag_solved=.FALSE.
         coeff_loc=coeff
         TGM_loc=TGM_p
         step_loc=step

         TGM_best_loc=TGM_p
         chisq_best_loc=chisq_best
         IF(flag_TGM) THEN
         pos_loc=pos
         ELSE
         pos_loc=3
         END IF

        DO i=1,20
!###     if flag_TGM=true, we evaluate TGM
          CALL find_extremes(coeff_loc,chisq_best,.TRUE.,TGM_loc,flag_D1)
          IF(flag_D1) THEN
           !if the ellopsoid fitting does not converge, change one point
           CALL stick_one_point(TGM_loc,coeff_loc)
           CYCLE
          END IF
          !if the result is lower than the best result, then cycle
          IF(TGM_loc(pos_loc)<=TGM_p(pos_loc)) EXIT
          !if TGM_loc is out of the parameter limits, then move it inside them  
!          write(*,*) 'out extremes',pos,TGM_loc
          CALL check_up_limits(step_loc,TGM_loc,flag_out_limit)

          !compute the chisq and chisq_poly
          IF(flag_TGM) THEN
           CALL chisq_tgm_e(TGM_loc,chisq)
          ELSE
           CALL chisq_abd_e(TGM_loc,chisq)
          END IF
          CALL eval_poly(TGM_loc,coeff_loc,chisq_poly)
!         write(*,*) 'TGM_loc pos=',flag_out_limit,pos,TGM_loc,chisq,chisq_poly,chisq_target
!##############################################
!############### if chisq is too small
         IF(chisq<chisq_target-chisq_tol) THEN
         IF(chisq_poly<chisq_target-chisq_tol.AND.flag_out_limit) THEN
!          write(*,*) 'exit',flag_out_limit
          EXIT
         END IF

          !re-compute the coeff around TGM_p
!          CALL denormalize3D(TGM_loc)
!          write(*,*) 'too low up_TGM pos=     ',pos,TGM_loc,chisq,chisq_poly,chisq_best+prob(1)
!          CALL normalize3D(TGM_loc)
!          CALL coeff_poly_find(TGM_loc,step_loc,coeff_loc)
          IF(flag_write) THEN
          CALL write_poly_surf(TGM_p,chisq,coeff_loc)
          WRITE(*,*) 'written'
          READ'(I1)', j
          END IF
!          CALL add_test_points(TGM_loc,step_loc,coeff_loc)
           CALL stick_one_point(TGM_loc,coeff_loc)
           IF(flag_TGM) THEN
            IF(TGM_loc(pos)>TGM_best_loc(pos)) TGM_best_loc=TGM_loc 
           ELSE
            IF(TGM_loc(3)>TGM_best_loc(3)) TGM_best_loc=TGM_loc            
           END IF
!##############################################
!############### if chisq is too large
         ELSE IF(chisq>chisq_target+chisq_tol) THEN
!          CALL denormalize3D(TGM_loc)
!          write(*,*) 'too high up_TGM pos=',pos,TGM_loc,chisq,chisq_poly,chisq_best+prob(1)
!          CALL normalize3D(TGM_loc)

          IF(flag_write) THEN
          CALL write_poly_surf(TGM_p,chisq,coeff_loc)
          WRITE(*,*) 'written'
          READ'(I1)', j
          END IF
           CALL stick_one_point(TGM_loc,coeff_loc)
         ELSE
!          CALL denormalize3D(TGM_loc)
!          write(*,*) '**** right up_TGM pos=',pos,TGM_loc,chisq,chisq_poly,chisq_best+prob(1)
!          CALL normalize3D(TGM_loc)
          flag_solved=.TRUE.
          EXIT
         END IF
        END DO

!#########################################
!### write the results on lo_TGM and up_TGM
!### if we measure TGM
     IF(flag_solved) THEN
         CALL denormalize3D(TGM_loc)
         IF(flag_TGM) THEN
           up_TGM(pos)=TGM_loc(pos)
         END IF
!### if we measure ABD
         IF(.NOT.flag_TGM) THEN
            up_ABD(pos)=TGM_loc(3)
         END IF
     ELSE

      IF(flag_TGM) THEN
       up_TGM(pos)=10000
      ELSE
      up_ABD(pos)=10000
      END IF
     END IF

END SUBROUTINE est_errors_p
!############################
SUBROUTINE est_errors_m(coeff,chisq_best,chisq_tol,step,TGM_m)
        IMPLICIT NONE
        REAL(DP), DIMENSION(dim_coeff),INTENT(IN) :: coeff
        REAL(DP), INTENT(IN) :: chisq_best,chisq_tol
        REAL(DP), DIMENSION(3),INTENT(IN) :: step
        REAL(DP), DIMENSION(3),INTENT(IN) :: TGM_m
        REAL(DP) :: chisq,chisq_poly,chisq_best_loc
        REAL(DP), DIMENSION(dim_coeff) :: coeff_loc
        REAL(DP), DIMENSION(3) :: TGM_loc,step_loc,TGM_best_loc
        LOGICAL :: flag_out_limit,flag_solved,flag_D1
        INTEGER(I1B) :: i,j
        INTEGER(I2B) :: pos_loc

         flag_out_limit=.FALSE.
         flag_solved=.FALSE.
         coeff_loc=coeff
         TGM_loc=TGM_m
         step_loc=step
         TGM_best_loc=TGM_m
         chisq_best_loc=chisq_best
         IF(flag_TGM) THEN
         pos_loc=pos
         ELSE
         pos_loc=3
         END IF

        DO i=1,20
!###     if flag_TGM=true, we evaluate TGM
          CALL find_extremes(coeff_loc,chisq_best,.FALSE.,TGM_loc,flag_D1)
          IF(flag_D1) THEN
           CALL stick_one_point(TGM_loc,coeff_loc)
           CYCLE
          END IF
          !if TGM_loc is out of the parameter limits, then move it inside them  
          CALL check_lo_limits(step_loc,TGM_loc,flag_out_limit)
          !if the result is higher than the best result, then cycle
          IF(TGM_loc(pos_loc)>=TGM_m(pos_loc)) EXIT
          !compute the chisq and chisq_poly
          IF(flag_TGM) THEN
           CALL chisq_tgm_e(TGM_loc,chisq)
          ELSE
           CALL chisq_abd_e(TGM_loc,chisq)
          END IF
          CALL eval_poly(TGM_loc,coeff_loc,chisq_poly)

!         write(*,*) 'TGM_loc pos=',flag_out_limit,pos,TGM_loc,chisq,chisq_poly,chisq_best+prob(1)
!##############################################
!############### if chisq is too small
         IF(chisq<chisq_best+prob(1)-chisq_tol) THEN
         IF(chisq_poly<chisq_best+prob(1)-chisq_tol.AND.flag_out_limit) THEN
          EXIT
         END IF
          !re-compute the coeff around TGM_p
!          CALL denormalize3D(TGM_loc)
!          write(*,*) 'too low lo_TGM pos=     ',pos,TGM_loc,chisq,chisq_poly,chisq_best+prob(1)
!          CALL normalize3D(TGM_loc)
         CALL stick_one_point(TGM_loc,coeff_loc)
          IF(flag_write) THEN
          CALL write_poly_surf(TGM_m,chisq,coeff_loc)
          WRITE(*,*) 'written'
          READ'(I1)', j
          END IF
           IF(flag_TGM) THEN
            IF(TGM_loc(pos)<TGM_best_loc(pos)) TGM_best_loc=TGM_loc 
           ELSE
            IF(TGM_loc(3)<TGM_best_loc(3)) TGM_best_loc=TGM_loc            
           END IF
!##############################################
!############### if chisq is too large
         ELSE IF(chisq>chisq_best+prob(1)+chisq_tol) THEN
!          CALL denormalize3D(TGM_loc)
!          write(*,*) 'too high lo_TGM pos=',pos,TGM_loc,chisq,chisq_poly,chisq_best+prob(1)
!          CALL normalize3D(TGM_loc)

          CALL stick_one_point(TGM_loc,coeff_loc)
          IF(flag_write) THEN
           CALL write_poly_surf(TGM_m,chisq,coeff_loc)
           WRITE(*,*) 'written'
           READ'(I1)', j
          END IF
         ELSE
!          write(*,*) '***** right lo_TGM pos=',pos,TGM_loc,chisq,chisq_poly,chisq_best+prob(1)
          flag_solved=.TRUE.
          EXIT
         END IF
        END DO

!#########################################
!### write the results on lo_TGM and up_TGM
!### if we measure TGM
     IF(flag_solved) THEN
         CALL denormalize3D(TGM_loc)
         IF(flag_TGM) THEN
           lo_TGM(pos)=TGM_loc(pos)
         END IF
!### if we measure ABD
         IF(.NOT.flag_TGM) THEN
            lo_ABD(pos)=TGM_loc(3)
         END IF
     ELSE
!       write(*,*) pos,'not solved'
      IF(flag_TGM) THEN
       lo_TGM(pos)=-10000
      ELSE
      lo_ABD(pos)=-10000
      END IF
     END IF


END SUBROUTINE est_errors_m
!############################
SUBROUTINE find_extremes(coeff,chisq_best,flag_up_low,TGM_loc,flag_D1)
!       USE data_lib, ONLY: temp_gridL,logg_gridL,met_gridL
       IMPLICIT NONE
       REAL(DP), DIMENSION(:),INTENT(IN) :: coeff
       REAL(DP), INTENT(IN) :: chisq_best
       REAL(DP), DIMENSION(3), INTENT(OUT) :: TGM_loc
       LOGICAL, INTENT(IN) :: flag_up_low
       LOGICAL, INTENT(OUT) :: flag_D1
       REAL(DP) :: a,b,c,d,e,f,g,h,i,j,A1,B1,C1,D1
       REAL(DP) :: x_p,x_m,y_p,y_m,z_p,z_m
       REAL(DP) :: cond1,cond2,cond3
       REAL(DP), DIMENSION(4,3) :: mat, mat_

!Here we solve the quadratic function 
!like ax^2+by^2+cz^2+dxy+exz+fyz+gx+hy+iz+j=chisq_best+prob(1) 
!which fits the chi^2 surface. First, I re-order it as a 2nd degree equation in z.
!Then I impose b^2-4ac=0, obtaining a 2nd degree equation in y.
!I impose again b^2-4ac=0, obtainig a 2nd degree equation in x.
!Solving this last equation I obtain the extreme x1 and x2 for which two
!planes parallel to yz are tangent to the quadratic function.
!By back-substituting x I obtain the y and z of these tangent points.
!Here I compute the analytic solution.

!### here I give the order of the coefficients of the quadratic function
!### in a matrix, where x=Teff, y=logg, z=[M/H]
!### the matrix is at first given to solve the max amn min of Teff.
!### then, by swapping some terms I can solve max and min for the 
!### other parameters (because the quadratic does not change shape 
!### by swapping x with y or z)
 mat(1,1)=coeff(10)!a, Teff^2
 mat(2,1)=coeff(9)!d, Teff*logg
 mat(3,1)=coeff(8)!e, Teff*met
 mat(4,1)=coeff(7)!g, Teff
 mat(1,2)=coeff(6)!b, logg^2
 mat(2,2)=coeff(9)!d, Teff*logg
 mat(3,2)=coeff(5)!f, logg*met
 mat(4,2)=coeff(4)!h, logg
 mat(1,3)=coeff(3)!c, met^2
 mat(2,3)=coeff(5)!f, logg*met
 mat(3,3)=coeff(8)!e, Teff*met
 mat(4,3)=coeff(2)!i, met

!initialize mat_
mat_=0.
!swap the terms as a function of pos
IF(pos==1) mat_=mat
IF(pos==2) THEN
 mat_(:,1)=mat(:,2)
 mat_(:,2)=mat(:,1)
 mat_(:,3)=mat(:,3)
 mat_(2,3)=mat(3,3)
 mat_(3,3)=mat(2,3)
END IF
IF(pos==3.OR..NOT.flag_TGM) THEN
 mat_(:,1)=mat(:,3)
 mat_(:,2)=mat(:,2)
 mat_(:,3)=mat(:,1)
 mat_(2,2)=mat(3,2)
 mat_(3,2)=mat(2,2)
END IF

!### assign the matrix terms to the coefficients named like
!### the quadratic function above
 j=coeff(1)-(chisq_best+prob(1))! known term
 i=mat_(4,3)!met
 c=mat_(1,3)!met^2
 h=mat_(4,2)!logg
 f=mat_(3,2)!logg*met
 b=mat_(1,2)!logg^2 
 g=mat_(4,1)!Teff
 e=mat_(3,1)!Teff*met
 d=mat_(2,2)!Teff*logg
 a=mat_(1,1)!Teff^2

 
!the 2nd degree equation in x can be written as A1x^2*B1x*C1=0 where A1, B1, C1 are
!the coefficients which solution is
 A1=-c*d*e*f+c**2*d**2+a*c*f**2+b*c*e**2-4*a*b*c**2
 B1=-c*d*f*i-c*e*f*h+2*c**2*d*h+2*b*c*e*i+c*f**2*g-4*b*c**2*g
 C1=-c*f*h*i+c**2*h**2+c*f**2*j+b*c*i**2-4*b*c**2*j
! define D1 as the argument of the sqrt
 D1=B1**2-4._dp*A1*C1
!if cond1<0, the plane z=0 it is a ellipse
 cond1=d**2-4*a*b
!if cond2<0, the plane y=0 it is a ellipse
 cond2=e**2-4*a*c
!if cond3<0, the plane x=0 it is a ellipse
 cond3=f**2-4*b*c
 IF(D1>0..AND.cond1<0.AND.cond2<0.AND.cond3<0) THEN 
  flag_D1=.FALSE.

!remember that A1 results negative. This is why
!(-B1+sqrt(D1))/(2*A1) is equal to x_m and not to x_p.
!The rest is similar.
 x_m=(-B1+sqrt(D1))/(2*A1)
 x_p=(-B1-sqrt(D1))/(2*A1)
 y_m=-(2*e*f*x_m+2*f*i-4*c*d*x_m-4*c*h)/(2*(f**2-4*b*c))
 y_p=-(2*e*f*x_p+2*f*i-4*c*d*x_p-4*c*h)/(2*(f**2-4*b*c))
 z_m=-(e*x_m+f*y_m+i)/(2*c)
 z_p=-(e*x_p+f*y_p+i)/(2*c)

 
!### if we compute the sup errors of TGM
 IF(flag_TGM.AND.flag_up_low) THEN
  IF(pos==1) THEN
   TGM_loc=(/x_p,y_p,z_p/)
  END IF
  IF(pos==2) THEN
   TGM_loc=(/y_p,x_p,z_p/)
  END IF
  IF(pos==3) THEN
   TGM_loc=(/z_p,y_p,x_p/)
  END IF
 END IF
!### if we compute the low errors of TGM
 IF(flag_TGM.AND..NOT.flag_up_low) THEN
  IF(pos==1) THEN
   TGM_loc=(/x_m,y_m,z_m/)
  END IF
  IF(pos==2) THEN
   TGM_loc=(/y_m,x_m,z_m/)
  END IF
  IF(pos==3) THEN
   TGM_loc=(/z_m,y_m,x_m/)
  END IF
 END IF
!### if we compute the sup errors of ABD
 IF(flag_up_low.AND..NOT.flag_TGM) THEN
   TGM_loc=(/z_p,y_p,x_p/)
 END IF

!### if we compute the inf errors of ABD
 IF(.NOT.flag_up_low.AND..NOT.flag_TGM) THEN
   TGM_loc=(/z_m,y_m,x_m/)
 END IF
 ELSE
!   IF(cond1<0.AND.cond2<0.AND.cond3<0) THEN 
!     write(*,*) 'this is not an ellipse'
!   ELSE
!     write(*,*) 'D1 negative'
!   END IF
 flag_D1=.TRUE.
END IF
END SUBROUTINE find_extremes

!############################
SUBROUTINE chisq_tgm_e(x,chisq)
       USE make_model, ONLY: make_model_TGM
       USE share, ONLY: dimsp,TGM,f_sp_norm,sig_noise
       IMPLICIT NONE
       REAL (DP), DIMENSION(3), INTENT(IN) :: x 
       REAL (DP), INTENT(OUT) :: chisq 
       REAL (DP), DIMENSION(dimsp) :: model
       REAL (DP), DIMENSION(5) :: xTGM 

       xTGM(1:3)=x(1:3)
       xTGM(4:5)=TGM(4:5)
       CALL denormalize3D(xTGM)

       CALL make_model_TGM(model,xTGM)
       chisq=SUM(((f_sp_norm-model)/sig_noise)**2)

END SUBROUTINE chisq_tgm_e
!############################
SUBROUTINE chisq_abd_e(x,chisq)
       USE make_model, ONLY: make_model_ABDerr
       USE share, ONLY: dimsp,dim_ele,TGM,ABD,f_sp_norm,sig_noise
       IMPLICIT NONE
       REAL (DP), DIMENSION(3), INTENT(IN) :: x 
       REAL (DP), INTENT(OUT) :: chisq 
       REAL (DP), DIMENSION(dimsp) :: model
       REAL (DP), DIMENSION(5) :: parTGM 
       REAL (DP), DIMENSION(dim_ele) :: parABD

       parTGM(1:2)=x(1:2)
       parTGM(3:5)=TGM(3:5)
       CALL denormalize3D(parTGM)

       parABD=ABD
       parABD(pos)=x(3)

       CALL make_model_ABDerr(model,parTGM,parABD)
       chisq=SUM(((f_sp_norm-model)/sig_noise)**2)

END SUBROUTINE chisq_abd_e
!############################
      INTEGER FUNCTION LA_WS_GELSS( VER, M, N, NRHS )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!   .. USE STATEMENTS ..
!corrado      USE F77_LAPACK, ONLY: ILAENV_F77 => LA_ILAENV
!     .. IMPLICIT STATEMENT ..
      IMPLICIT NONE
!     .. SCALAR ARGUMENTS ..
      CHARACTER(LEN=1), INTENT(IN) :: VER
      INTEGER, INTENT(IN) :: M, N, NRHS
!     .. PARAMETERS ..
      CHARACTER(LEN=5), PARAMETER :: NAME1='GELSS', NAME2='GEQRF', NAME3='ORMQR', NAME4='GEBRD', &
                                     NAME5='ORMBR', NAME6='ORGBR', NAME7='GELQF', NAME8='ORMLQ'
!     .. LOCAL SCALARS ..
      INTEGER :: MNTHR, MINWRK, MAXWRK, MM, BDSPAC
!     .. INTRINSIC FUNCTIONS ..
      INTRINSIC MAX
!     .. EXECUTABLE STATEMENTS ..
!Corrado      MNTHR = ILAENV_F77( 6, VER//NAME1, ' ', M, N, NRHS, -1 )
      MNTHR = LA_ILAENV( 6, VER//NAME1, ' ', M, N, NRHS, -1 )
!
!     COMPUTE WORKSPACE
!      (NOTE: COMMENTS IN THE CODE BEGINNING "Workspace:" DESCRIBE THE
!       MINIMAL AMOUNT OF WORKSPACE NEEDED AT THAT POINT IN THE CODE,
!       AS WELL AS THE PREFERRED AMOUNT FOR GOOD PERFORMANCE.
!       NB REFERS TO THE OPTIMAL BLOCK SIZE FOR THE IMMEDIATELY
!       FOLLOWING SUBROUTINE, AS RETURNED BY ILAENV.)
!
      MINWRK = 1
      MAXWRK = 0
      MM = M
      IF( M.GE.N .AND. M.GE.MNTHR ) THEN
!
!        PATH 1A - OVERDETERMINED, WITH MANY MORE ROWS THAN COLUMNS
!
         MM = N
         MAXWRK = MAX(MAXWRK,N+N*LA_ILAENV(1,VER//NAME2,' ',M,N,-1,-1))
         MAXWRK = MAX( MAXWRK, N+NRHS*                                  &
     &            LA_ILAENV( 1, VER//NAME3, 'LT', M, NRHS, N, -1 ) )
      END IF
      IF( M.GE.N ) THEN
!
!        PATH 1 - OVERDETERMINED OR EXACTLY DETERMINED
!
!        COMPUTE WORKSPACE NEEDE FOR DBDSQR
!
         BDSPAC = MAX( 1, 5*N-4 )
         MAXWRK = MAX( MAXWRK, 3*N+( MM+N )*                            &
     &            LA_ILAENV( 1, VER//NAME4, ' ', MM, N, -1, -1 ) )
         MAXWRK = MAX( MAXWRK, 3*N+NRHS*                                &
     &            LA_ILAENV( 1, VER//NAME5, 'QLT', MM, NRHS, N, -1 ) )
         MAXWRK = MAX( MAXWRK, 3*N+( N-1 )*                             &
     &            LA_ILAENV( 1, VER//NAME6, 'P', N, N, N, -1 ) )
         MAXWRK = MAX( MAXWRK, BDSPAC )
         MAXWRK = MAX( MAXWRK, N*NRHS )
         MINWRK = MAX( 3*N+MM, 3*N+NRHS, BDSPAC )
         MAXWRK = MAX( MINWRK, MAXWRK )
      END IF
      IF( N.GT.M ) THEN
!
!        COMPUTE WORKSPACE NEEDE FOR DBDSQR
!
         BDSPAC = MAX( 1, 5*M-4 )
         MINWRK = MAX( 3*M+NRHS, 3*M+N, BDSPAC )
         IF( N.GE.MNTHR ) THEN
!
!           PATH 2A - UNDERDETERMINED, WITH MANY MORE COLUMNS
!           THAN ROWS
!
            MAXWRK = M + M*LA_ILAENV( 1,VER//NAME7,' ',M,N,-1,-1 )
            MAXWRK = MAX( MAXWRK, M*M+4*M+2*M*                          &
     &               LA_ILAENV( 1, VER//NAME4, ' ', M, M, -1, -1 ) )
            MAXWRK = MAX( MAXWRK, M*M+4*M+NRHS*                         &
     &               LA_ILAENV( 1,VER//NAME5,'QLT',M,NRHS,M,-1 ) )
            MAXWRK = MAX( MAXWRK, M*M+4*M+( M-1 )*                      &
     &               LA_ILAENV( 1, VER//NAME6, 'P', M, M, M, -1 ) )
            MAXWRK = MAX( MAXWRK, M*M+M+BDSPAC )
            IF( NRHS.GT.1 ) THEN
               MAXWRK = MAX( MAXWRK, M*M+M+M*NRHS )
            ELSE
               MAXWRK = MAX( MAXWRK, M*M+2*M )
            END IF
            MAXWRK = MAX( MAXWRK, M+NRHS*                               &
     &               LA_ILAENV( 1, VER//NAME8, 'LT', N, NRHS, M, -1 ) )
         ELSE
!
!           PATH 2 - UNDERDETERMINED
!
            MAXWRK = 3*M+(N+M)*LA_ILAENV(1,VER//NAME4,' ',M,N,-1,-1)
            MAXWRK = MAX( MAXWRK, 3*M+NRHS*                             &
     &               LA_ILAENV( 1,VER//NAME5,'QLT',M,NRHS,M,-1 ) )
            MAXWRK = MAX( MAXWRK, 3*M+M*                                &
     &               LA_ILAENV( 1, VER//NAME6, 'P', M, N, M, -1 ) )
            MAXWRK = MAX( MAXWRK, BDSPAC )
            MAXWRK = MAX( MAXWRK, N*NRHS )
         END IF
      END IF
      LA_WS_GELSS = MAX( MINWRK, MAXWRK )
      END FUNCTION LA_WS_GELSS

!####################################################
 SUBROUTINE DGELSS_F95( A, B, RANK, S, RCOND, INFO )
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!   .. USE STATEMENTS ..
!    USE LA_PRECISION, ONLY: WP => DP
!    USE LA_AUXMOD, ONLY: ERINFO, LA_WS_GELSS
!corrado    USE F77_LAPACK, ONLY: GELSS_F77 => LA_GELSS
!   .. IMPLICIT STATEMENT ..
    IMPLICIT NONE
!   .. SCALAR ARGUMENTS ..
    INTEGER, INTENT(OUT), OPTIONAL :: RANK
    INTEGER, INTENT(OUT), OPTIONAL :: INFO
    REAL(DP), INTENT(IN), OPTIONAL :: RCOND
!   .. ARRAY ARGUMENTS ..
    REAL(DP), INTENT(INOUT) :: A(:,:), B(:,:)
    REAL(DP), INTENT(OUT), OPTIONAL, TARGET :: S(:)
!----------------------------------------------------------------------
! 
! Purpose
! ======= 
!
!       LA_GELSS and LA_GELSD compute the minimum-norm least squares 
! solution to one or more real or complex linear systems A*x = b using
! the singular value decomposition of A. Matrix A is rectangular and may
! be rank-deficient. The vectors b and corresponding solution vectors x
! are the columns of matrices denoted B and X , respectively.
!       The effective rank of A is determined by treating as zero those 
! singular values which are less than RCOND times the largest singular 
! value. In addition to X , the routines also return the right singular
! vectors and, optionally, the rank and singular values of A.
!       LA_GELSD combines the singular value decomposition with a divide
! and conquer technique. For large matrices it is often much faster than 
! LA_GELSS but uses more workspace.
! 
! ==========
! 
!        SUBROUTINE LA_GELSS / LA_GELSD( A, B, RANK=rank, S=s, &
!                                          RCOND=rcond, INFO=info )
!             <type>(<wp>), INTENT( INOUT ) :: A( :, : ), <rhs>
!             INTEGER, INTENT(OUT), OPTIONAL :: RANK
!             REAL(<wp>), INTENT(OUT), OPTIONAL :: S(:)
!             REAL(<wp>), INTENT(IN), OPTIONAL :: RCOND
!             INTEGER, INTENT(OUT), OPTIONAL :: INFO
!        where
!             <type> ::= REAL | COMPLEX
!             <wp>   ::= KIND(1.0) | KIND(1.0D0)
!             <rhs>  ::= B(:,:) | B(:)
! 
! Arguments
! =========
! 
! A       (input/output) REAL or COMPLEX array, shape (:,:).
!         On entry, the matrix A.
!         On exit, the first min(size(A,1), size(A,2)) rows of A are 
!         overwritten with its right singular vectors, stored rowwise.
! B       (input/output) REAL or COMPLEX array, shape (:,:) with 
!         size(B,1) = max(size(A,1), size(A,2)) or shape (:) with 
!         size(B) = max(size(A,1), size(A,2)).
!         On entry, the matrix B.
!         On exit, the solution matrix X .
!         If size(A,1) >= size(A,2) and RANK = size(A,2), the residual 
!         sum-of-squares for the solution in a column of B is given by 
!         the sum of squares of elements in rows size(A,2)+1:size(A,1)
!         of that column.
! RANK    Optional (output) INTEGER.
!         The effective rank of A, i.e., the number of singular values 
!         of A which are greater than the product RCOND*sigma1 , where
!  	  sigma1 is the greatest singular value.
! S       Optional (output) REAL array, shape (:) with size(S) = 
!         min(size(A,1), size(A,2)).
!         The singular values of A in decreasing order.
!         The condition number of A in the 2-norm is
! 	  K2(A)= sigma1/sigma(min(size(A,1),size(A,2)) .
! RCOND   Optional (input) REAL.
!         RCOND is used to determine the effective rank of A.
!         Singular values sigma(i)<=RCOND*sigma1  are treated as zero.
!         Default value: 10*max(size(A,1), size(A,2))*EPSILON(1.0_<wp>),
!         where <wp> is the working precision.
! INFO    Optional (output) INTEGER.
!         = 0: successful exit.
!         < 0: if INFO = -i, the i-th argument had an illegal value.
!         > 0: the algorithm for computing the SVD failed to converge; 
! 	  if INFO = i,i off-diagonal elements of an intermediate 
! 	  bidiagonal form did not converge to zero.
!         If INFO is not present and an error occurs, then the program 
! 	  is terminated with an error message.
!----------------------------------------------------------------------
!   .. PARAMETERS ..
    CHARACTER(LEN=8), PARAMETER :: SRNAME = 'LA_GELSS'
    CHARACTER(LEN=1), PARAMETER :: VER = 'D'
!   .. LOCAL SCALARS ..
    INTEGER :: LINFO, ISTAT, ISTAT1, LWORK, N, M, MN, NRHS, LRANK, SS
    REAL(DP) :: LRCOND
!   .. LOCAL POINTERS ..
    REAL(DP), POINTER :: WORK(:), LS(:)
!   .. INTRINSIC FUNCTIONS ..
    INTRINSIC SIZE, PRESENT, MAX, MIN, EPSILON
!   .. EXECUTABLE STATEMENTS ..
    LINFO = 0; ISTAT = 0; M = SIZE(A,1); N = SIZE(A,2); NRHS = SIZE(B,2)
    MN = MIN(M,N)
    IF( PRESENT(RCOND) )THEN; LRCOND = RCOND; ELSE
        LRCOND = 100*EPSILON(1.0_DP) ; ENDIF
    IF( PRESENT(S) )THEN; SS = SIZE(S); ELSE; SS =MN; ENDIF
!   .. TEST THE ARGUMENTS
    IF( M < 0 .OR. N < 0 ) THEN; LINFO = -1
    ELSE IF( SIZE( B, 1 ) /= MAX(1,M,N) .OR. NRHS < 0 ) THEN; LINFO = -2
    ELSE IF( SS /= MN ) THEN; LINFO = -4
    ELSE IF( LRCOND <= 0.0_DP ) THEN; LINFO = -5
    ELSE
       IF( PRESENT(S) )THEN; LS => S
       ELSE; ALLOCATE( LS(MN), STAT = ISTAT ); END IF
       IF( ISTAT == 0 ) THEN
          LWORK = LA_WS_GELSS( VER, M, N, NRHS )
          ALLOCATE( WORK(LWORK), STAT = ISTAT )
          IF( ISTAT /= 0 ) THEN
             DEALLOCATE( WORK, STAT=ISTAT1 )
             LWORK = MAX( 1, 3*MIN(M,N) + MAX( 2*MIN(M,N), MAX(M,N), NRHS ) )
             ALLOCATE( WORK(LWORK), STAT = ISTAT )
             IF( ISTAT /= 0 ) CALL ERINFO( -200, SRNAME, LINFO )
          END IF
       END IF
       IF ( ISTAT == 0 ) THEN
          CALL LA_GELSS( M, N, NRHS, A, MAX(1,M), B, MAX(1,M,N), &
                         LS, LRCOND, LRANK, WORK, LWORK, LINFO )
       ELSE; LINFO = -100; END IF
       IF( PRESENT(RANK) ) RANK = LRANK
       DEALLOCATE(WORK, STAT = ISTAT1 )
    END IF
    CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )
 END SUBROUTINE DGELSS_F95


!###########################################
      SUBROUTINE ERINFO(LINFO, SRNAME, INFO, ISTAT)
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
!  .. IMPLICIT STATEMENT ..
         IMPLICIT NONE
!  .. SCALAR ARGUMENTS ..
         CHARACTER( LEN = * ), INTENT(IN)              :: SRNAME
         INTEGER             , INTENT(IN)              :: LINFO
         INTEGER             , INTENT(OUT), OPTIONAL   :: INFO
         INTEGER             , INTENT(IN), OPTIONAL    :: ISTAT
!  .. EXECUTABLE STATEMENTS ..
!         IF( ( LINFO < 0 .AND. LINFO > -200 ) .OR.                     &
!    &       ( LINFO > 0 .AND. .NOT.PRESENT(INFO) ) )THEN
      IF( ( ( LINFO < 0 .AND. LINFO > -200 ) .OR. LINFO > 0 )           &
     &           .AND. .NOT.PRESENT(INFO) )THEN
        WRITE (*,*) 'Program terminated in LAPACK95 subroutine ',SRNAME
        WRITE (*,*) 'Error indicator, INFO = ',LINFO
        IF( PRESENT(ISTAT) )THEN
          IF( ISTAT /= 0 ) THEN
            IF( LINFO == -100 )THEN
              WRITE (*,*) 'The statement ALLOCATE causes STATUS = ',    &
     &                    ISTAT
            ELSE
              WRITE (*,*) 'LINFO = ', LINFO, ' not expected'
            END IF
          END IF   
        END IF
        STOP
         ELSE IF( LINFO <= -200 ) THEN
           WRITE(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++'
           WRITE(*,*) '*** WARNING, INFO = ', LINFO, ' WARNING ***'
           IF( LINFO == -200 )THEN
             WRITE(*,*)                                                 &
     &        'Could not allocate sufficient workspace for the optimum'
             WRITE(*,*)                                                 &
     &        'blocksize, hence the routine may not have performed as'
             WRITE(*,*) 'efficiently as possible'
         ELSE
           WRITE(*,*) 'Unexpected warning'
         END IF
           WRITE(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++'
        END IF
        IF( PRESENT(INFO) ) THEN
          INFO = LINFO
        END IF
      END SUBROUTINE ERINFO
!############################
SUBROUTINE write_poly_surf(TGM_in,chisq_best,coeff)
        IMPLICIT NONE
        REAL(DP), INTENT(IN) :: chisq_best
        REAL(DP), DIMENSION(dim_coeff),INTENT(IN) :: coeff
        REAL(DP), DIMENSION(3),INTENT(IN) :: TGM_in
        REAL(DP), DIMENSION(3) :: TGM_,TGM_test_loc,step
        REAL(DP) :: chisq,chisq_poly
        INTEGER (I1B) :: i,j,k

        TGM_=TGM_in
        CALL denormalize3D(TGM_)

        OPEN(unit=10,file='space_chisq.dat',action='WRITE')        
        OPEN(unit=11,file='space_chisq_poly.dat',action='WRITE')        
        write(11,fmt='(F8.3,2X,F8.3,2X,F8.3,2X,F12.4)') TGM_(1),TGM_(2),TGM_(3),chisq_best
        write(10,fmt='(F8.3,2X,F8.3,2X,F8.3,2X,F12.4)') TGM_(1),TGM_(2),TGM_(3),chisq_best

        step=(/20.0,0.05,0.02/)
        CALL normalize3D(step)
        DO k=-10,10
        DO j=-10,10
        DO i=-10,10
        
         TGM_test_loc(1)=TGM_in(1)+step(1)*i
         TGM_test_loc(2)=TGM_in(2)+step(2)*j
         TGM_test_loc(3)=TGM_in(3)+step(3)*k

         CALL eval_poly(TGM_test_loc,coeff,chisq_poly)
         IF(flag_TGM) CALL chisq_tgm_e(TGM_test_loc,chisq)
         IF(.NOT.flag_TGM) CALL chisq_abd_e(TGM_test_loc,chisq)
!         CALL denormalize3D(TGM_test_loc)
!         write(*,*) TGM_test_loc,chisq_poly
         CALL denormalize3D(TGM_test_loc)
         write(11,fmt='(F8.3,2X,F8.3,2X,F8.3,2X,F12.4)') TGM_test_loc(1),TGM_test_loc(2),TGM_test_loc(3),chisq_poly
         write(10,fmt='(F8.3,2X,F8.3,2X,F8.3,2X,F12.4)') TGM_test_loc(1),TGM_test_loc(2),TGM_test_loc(3),chisq
        ENDDO
        ENDDO
        ENDDO

        CLOSE(11)
        CLOSE(10)
END SUBROUTINE write_poly_surf
!###############################
SUBROUTINE define_grid(TGM_ini,step,tgm_grid)
        USE data_lib, ONLY: temp_gridL,logg_gridL,met_gridL
        IMPLICIT NONE
        REAL(DP), DIMENSION(3),INTENT(IN) :: TGM_ini,step
        REAL(DP), DIMENSION(n_eq,3),INTENT(OUT) :: tgm_grid
        REAL(DP), DIMENSION(3) :: TGM_,step_,step_loc,tgm_out,rnd
        LOGICAl :: flag_found
        REAL(DP) :: chi
        INTEGER (I1B) :: i,j,k,l
        INTEGER (I2B) :: count

       TGM_=TGM_ini
       step_=step
       CALL denormalize3D(TGM_)
       CALL denormalize3D(step_)       

!       write(*,*) 'define grid TGM_, step_',TGM_,step_

        IF(TGM_(1)>temp_gridL(SIZE(temp_gridL))) THEN
         TGM_(1)=MIN(TGM_(1),temp_gridL(SIZE(temp_gridL)))
         TGM_(1)=TGM_(1)-1.1*step_(1)
        END IF

        IF(TGM_(2)>logg_gridL(SIZE(logg_gridL))) THEN
         TGM_(2)=MIN(TGM_(2),logg_gridL(SIZE(logg_gridL)))
         TGM_(2)=TGM_(2)-1.1*step_(2)
        END IF

        IF(flag_TGM) THEN
         IF(TGM_(3)>met_gridL(SIZE(met_gridL))) THEN
          TGM_(3)=MIN(TGM_(3),met_gridL(SIZE(met_gridL)))
          TGM_(3)=TGM_(3)-1.1*step_(3)
         END IF
        ELSE
         IF(TGM_(3)>0.8) THEN
          TGM_(3)=MIN(TGM_(3),0.8)
          TGM_(3)=TGM_(3)-1.1*step_(3)
         END IF
        END IF

        IF(TGM_(1)<temp_gridL(1)) THEN
         TGM_(1)=MAX(TGM_(1),temp_gridL(1))
         TGM_(1)=TGM_(1)+1.1*step_(1)
        END IF

        IF(TGM_(2)<logg_gridL(1)) THEN
         TGM_(2)=MAX(TGM_(2),logg_gridL(1))
         TGM_(2)=TGM_(2)+1.1*step_(2)
        END IF

        IF(flag_TGM) THEN
         IF(TGM_(3)<met_gridL(1)) THEN
          TGM_(3)=MAX(TGM_(3),met_gridL(1))
          TGM_(3)=TGM_(3)+1.1*step_(3)
         END IF
        ELSE
         IF(TGM_(3)<-0.4) THEN  
          TGM_(3)=MAX(TGM_(3),-0.4)
          TGM_(3)=TGM_(3)+1.1*step_(3)
         END IF
        END IF

        CALL normalize3D(TGM_)
        CALL normalize3D(step_)
        !initialize count
        count=0

        DO k=-1,1
        DO j=-1,1
        DO i=-1,1
 
         IF(.NOT.(k==0.AND.j==0.AND.i==0)) THEN

         count=count+1_I1B
         tgm_grid(count,1)=TGM_(1)+step_(1)*i
         tgm_grid(count,2)=TGM_(2)+step_(2)*j
         tgm_grid(count,3)=TGM_(3)+step_(3)*k

         CALL find_chisq_border(TGM_,tgm_grid(count,:),tgm_out,chi,flag_found)
         ! if the border is found, then
         IF(flag_found) THEN
          tgm_grid(count,:)=tgm_out
!          write(*,*) 'output',tgm_out,chi,chisq_target,count
         ELSE
         !if the border is not found, chose a random direction
         ! and look for it until 100 times
          DO l=1,100
           CALL RANDOM_NUMBER(rnd)
           step_loc=step_*(rnd*2-1)
           tgm_grid(count,:)=TGM_+step_loc
           CALL find_chisq_border(TGM_,tgm_grid(count,:),tgm_out,chi,flag_found)
           IF(flag_found) THEN
            tgm_grid(count,:)=tgm_out
            write(*,*) 'output repeated',tgm_out,chi,chisq_target,count
            EXIT
           END IF
          END DO
          !if after 100 times the border has not been found, assign the point at minimum chisq 
          IF(l>=100) THEN
           tgm_grid(count,:)=TGM_
          END IF
         END IF
!         weight_test(count)=1_dp
!  write(*,*) count,tgm_grid(count,1:3)
        END IF
        ENDDO
        ENDDO
        ENDDO

         count=count+1_I2B
         tgm_grid(count,:)=TGM_

END SUBROUTINE define_grid
!############################
        SUBROUTINE TGM_errors_null
        USE utils, ONLY: write_ABD_mask
        USE share, ONLY: ABD
        IMPLICIT NONE
        INTEGER(I1B) :: i

        DO i=1,3
         up_TGM(i)=9999._dp
         lo_TGM(i)=-9999._dp
        END DO
        DO i=1,INT(SIZE(ABD,1),I1B)
         IF(write_ABD_mask(i)) THEN
         up_ABD(i)=9.99_dp
         lo_ABD(i)=-9.99_dp
         END IF
        END DO

        END SUBROUTINE TGM_errors_null
!###########################################
SUBROUTINE find_chisq_border(point1,point2,tgm_out,chi,flag_found)
  IMPLICIT NONE
  REAL(DP), DIMENSION(3),INTENT(IN) :: point1,point2
  REAL(DP), DIMENSION(3),INTENT(OUT) :: tgm_out
  REAL(DP), INTENT(OUT) :: chi
  LOGICAL, INTENT(OUT) :: flag_found
  REAL(DP), DIMENSION(3) :: v_comp,tgm_loc,vec,t
  REAL(DP) :: t_p,t_m
  INTEGER (I1B) :: i
  INTEGER (I1B),DIMENSION(1) :: a
  LOGICAL :: flag_solved

  v_comp=point2-point1
  t=(/0,1,2/)
  CALL set_t_vec(point1,v_comp,t,vec)


  chi=-1.
  DO i=1,10 
!    write(*,*) 'input find_zeros',t,vec
    CALL find_zeros(t,vec,t_p,t_m,flag_solved)
    IF(.NOT.flag_solved) THEN
!      write(*,*) 'find_zeros not solved',i
      t=t+(/0,1,2/)
!      write(*,*) 't*2',t
      IF(ANY(t>100)) THEN
        flag_found=.TRUE.
        t_p=t(2)
        RETURN
      END IF
      CALL set_t_vec(point1,v_comp,t,vec)
      CYCLE
    END IF
!    write(*,*) 't_p,t_m,flag_solved',t_p,t_m,flag_solved
    tgm_loc=point1+v_comp*t_p
!    write(*,*) 'point1,tgm_loc',i,point1,tgm_loc,vec
    IF(flag_TGM) THEN
     CALL chisq_tgm_e(tgm_loc,chi)
    ELSE
     CALL chisq_abd_e(tgm_loc,chi)         
    END IF     

    IF (ABS(chi-chisq_target)<0.05) THEN
!     write(*,*) 'converged,chi,chisq_target',tgm_loc,chi,chisq_target
      EXIT
    END IF
    a=INT(MAXLOC(ABS(vec-chi)),I1B)
    vec(a(1))=chi
    t(a(1))=t_p
!    write(*,*) 'chi,chisq_target,vec,t',chi,chisq_target,vec,t
  END DO

  tgm_out=tgm_loc
  flag_found=.TRUE.

END SUBROUTINE find_chisq_border
!############################
!set the arrays t and vec
SUBROUTINE set_t_vec(point1,v_comp,t,vec)
  REAL(DP), DIMENSION(3),INTENT(IN) :: point1,v_comp,t
  REAL(DP), DIMENSION(3),INTENT(OUT) :: vec
  REAL(DP), DIMENSION(3) :: tgm_loc
  INTEGER (I1B) :: i

  DO i=1,3
    tgm_loc=point1+v_comp*t(i)
    IF(flag_TGM) THEN
     CALL chisq_tgm_e(tgm_loc,vec(i))
    ELSE
     CALL chisq_abd_e(tgm_loc,vec(i))         
    END IF
  END DO

END SUBROUTINE set_t_vec
!############################
!find a parabola through the three points and find where is equal to zero
SUBROUTINE find_zeros(t,vec,t_p,t_m,flag_solved)
  LOGICAL, INTENT(OUT) :: flag_solved
  REAL(DP), INTENT(OUT) :: t_p,t_m
  REAL(DP), DIMENSION(3),INTENT(IN) :: t,vec
  REAL(DP), DIMENSION(3,3) :: mat,mat1
  REAL(DP) :: D,Da,Db,Dc,a,b,c,arg_sq,t1_p,t1_m
  INTEGER (I2B) :: i
    DO i=1,3
   mat(i,1)=t(i)**2
   mat(i,2)=t(i)
   mat(i,3)=1._dp
  END DO


  D=det(mat)
  mat1=mat
  mat1(:,1)=vec
  Da=det(mat1)
  mat1=mat
  mat1(:,2)=vec
  Db=det(mat1)
  mat1=mat
  mat1(:,3)=vec
  Dc=det(mat1)

  a=Da/D
  b=Db/D
  c=Dc/D-chisq_target
  arg_sq=b**2-4*a*c
  IF(arg_sq<0.) THEN
   flag_solved=.FALSE.
   t_p=t(2)
   t_m=t(2)
!   write(*,*) 'arg_sq<0.',arg_sq,t,t_p,vec
  ELSE
   t1_p=(-b+sqrt(arg_sq))/(2*a)
   t1_m=(-b-sqrt(arg_sq))/(2*a)

   t_p=MAX(t1_p,t1_m)
   t_m=MIN(t1_p,t1_m)
   flag_solved=.TRUE.
  END IF

END SUBROUTINE find_zeros
!############################
REAL FUNCTION det(mat)
     IMPLICIT NONE
     REAL(DP), DIMENSION(3,3), INTENT(IN) :: mat
     
     det = mat(1,1)*(mat(2,2)*mat(3,3) - mat(3,2)*mat(2,3)) &
      & + mat(1,2)*(mat(3,1)*mat(2,3) - mat(2,1)*mat(3,3))  &
      & + mat(1,3)*(mat(2,1)*mat(3,2) - mat(3,1)*mat(2,2))

END FUNCTION det
!###########################################



END MODULE uncertains2
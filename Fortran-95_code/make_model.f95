!    This module contains the routines that construct the spectrum models.
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


MODULE make_model
        USE num_type
        CONTAINS
!#########################
        SUBROUTINE make_model_TGM_quick(model,params)
        USE share, ONLY: w_sp,wave_ll,select_ll_mask&
        &,dim_ll,wave_center_ll,coeff_4deg_quick,dimsp,ew
        USE func_poly, ONLY: poly6_transform,ew_poly6, voigt_logg
        IMPLICIT NONE
        
        REAL(DP) :: sig,c,gamg,gaml,vv
        INTEGER(I4B) :: i,j,w_center_line
        REAL(DP),DIMENSION(dimsp), INTENT(INOUT) :: model
        REAL(DP),DIMENSION(:) :: params
        REAL(DP),DIMENSION(4) :: st_pars
        REAL(DP), DIMENSION(84) :: X
        !initialize model
        model=1._dp

        IF(params(1)<3600) params(1)=3600._dp
        IF(params(1)>7400) params(1)=7400._dp
        IF(params(2)<0.1) params(2)=0.1_dp
        IF(params(2)>5.5) params(2)=5.5_dp
        IF(params(3)<-2.6) params(3)=-2.6_dp
        IF(params(3)>0.5) params(3)=0.5_dp
        IF(params(4)<0.01) params(4)=0.01_dp
        IF(ABS(params(4))>5.) params(4)=5.0_dp

        st_pars(1:3)=params(1:3)
        st_pars(4)=0.0001_dp        
        sig=ABS(params(4))
        gamg=2.0_dp*sig*1.17741002252_dp!### sqrt(2.*log(2.))=1.17741002252
        c=1.0_dp+params(5)/299792.0_dp

        CALL poly6_transform(X,st_pars)

        DO i=1,dim_ll
          IF(select_ll_mask(i)) THEN
            ew(i)=ew_poly6(X,coeff_4deg_quick(:,i))
            w_center_line=wave_center_ll(i)
            !compute the parameters used by the voigt function
            gaml=gammaL(i,ew(i),gamg,st_pars(2),st_pars(1))

            DO j=w_center_line-1_I4B,1,-1
             vv=voigt_logg(w_sp(j),wave_ll(i)*c,gamg,ew(i),gaml)
             IF(vv>1.E-4) THEN
              model(j)=model(j)-vv
             ELSE
              EXIT
             END IF
            END DO
            DO j=w_center_line,dimsp
             vv=voigt_logg(w_sp(j),wave_ll(i)*c,gamg,ew(i),gaml)
             IF(vv>1.E-4) THEN
              model(j)=model(j)-vv
             ELSE
              EXIT
             END IF
            END DO
          END IF
        END DO

        END SUBROUTINE make_model_TGM_quick
!#########################
        SUBROUTINE make_model_TGM(model,params)
        USE share, ONLY: ele_ll,w_sp,wave_ll,select_ll_mask&
        &,dim_ll,ABD,wave_center_ll,coeff_4deg,dimsp,ew
        USE func_poly, ONLY: poly4_transform,ew_poly4, voigt_logg
        IMPLICIT NONE
        
        REAL(DP) :: sig,c,gamg,gaml,vv
        INTEGER(I4B) :: i,j,w_center_line
        REAL(DP),DIMENSION(dimsp), INTENT(INOUT) :: model
        REAL(DP),DIMENSION(:) :: params
        REAL(DP),DIMENSION(4) :: st_pars
        REAL(DP), DIMENSION(70) :: X
        !initialize model
        model=1._dp

        IF(params(1)<3600) params(1)=3600._dp
        IF(params(1)>7400) params(1)=7400._dp
        IF(params(2)<0.1) params(2)=0.1_dp
        IF(params(2)>5.5) params(2)=5.5_dp
        IF(params(3)<-2.6) params(3)=-2.6_dp
        IF(params(3)>0.5) params(3)=0.5_dp
        IF(params(4)<0.01) params(4)=0.01_dp
        IF(ABS(params(4))>5.) params(4)=5.0_dp

        st_pars(1:3)=params(1:3)
        st_pars(4)=0.0001_dp        
        sig=ABS(params(4))
        gamg=2.0_dp*sig*1.17741002252_dp!### sqrt(2.*log(2.))=1.17741002252
        c=1.0_dp+params(5)/299792.0_dp

        CALL poly4_transform(X,ABD,st_pars)

        DO i=1,dim_ll
          IF(select_ll_mask(i)) THEN
            ew(i)=ew_poly4(X,ele_ll(i),coeff_4deg(:,i))
            w_center_line=wave_center_ll(i)
            !compute the parameters used by the voigt function
            gaml=gammaL(i,ew(i),gamG,st_pars(2),st_pars(1))

            DO j=w_center_line-1_I4B,1,-1
             vv=voigt_logg(w_sp(j),wave_ll(i)*c,gamg,ew(i),gaml)
             IF(vv>1.E-4) THEN
              model(j)=model(j)-vv
             ELSE
              EXIT
             END IF
            END DO
            DO j=w_center_line,dimsp
             vv=voigt_logg(w_sp(j),wave_ll(i)*c,gamg,ew(i),gaml)
             IF(vv>1.E-4) THEN
              model(j)=model(j)-vv
             ELSE
              EXIT
             END IF
            END DO
          END IF
        END DO
        END SUBROUTINE make_model_TGM
!#########################
        SUBROUTINE make_model_ABD(model,params)
        USE share, ONLY: ele_ll,TGM,w_sp,wave_ll,select_ll_mask&
        &,dim_ll,wave_center_ll,coeff_4deg,dimsp,ew
        USE func_poly, ONLY: poly4_transform,ew_poly4, voigt_logg
        IMPLICIT NONE
        
        REAL(DP) :: sig,c,gamg,gaml,vv
        INTEGER(I4B) :: i,j,w_center_line
        REAL(DP),DIMENSION(dimsp), INTENT(INOUT) :: model
        REAL(DP),DIMENSION(:) :: params
        REAL(DP),DIMENSION(4) :: st_pars
        REAL(DP), DIMENSION(70) :: X
        !initialize model
        model=1._dp

        st_pars(1:3)=TGM(1:3)
        st_pars(4)=0.0001_dp        
        sig=ABS(TGM(4))
        gamg=2.0_dp*sig*1.17741002252_dp!### sqrt(2.*log(2.))=1.17741002252
        c=1.0_dp+TGM(5)/299792.0_dp



        WHERE(params>0.8) params=0.8_dp
        WHERE(params<-0.6) params=-0.6_dp

        CALL poly4_transform(X,params,st_pars)

        DO i=1,dim_ll
          IF(select_ll_mask(i)) THEN
            ew(i)=ew_poly4(X,ele_ll(i),coeff_4deg(:,i))
            w_center_line=wave_center_ll(i)
            !compute the parameters used by the voigt function
            gaml=gammaL(i,ew(i),gamG,st_pars(2),st_pars(1))
            DO j=w_center_line-1_I4B,1,-1
             vv=voigt_logg(w_sp(j),wave_ll(i)*c,gamg,ew(i),gaml)
             IF(vv>1.E-4) THEN
              model(j)=model(j)-vv
             ELSE
              EXIT
             END IF
            END DO
            DO j=w_center_line,dimsp
             vv=voigt_logg(w_sp(j),wave_ll(i)*c,gamg,ew(i),gaml)
             IF(vv>1.E-4) THEN
              model(j)=model(j)-vv
             ELSE
              EXIT
             END IF
            END DO
          END IF
        END DO

        END SUBROUTINE make_model_ABD
!#########################
        SUBROUTINE make_model_ABDerr(model,parTGM,parABD)
        USE share, ONLY: ele_ll,w_sp,wave_ll,select_ll_mask&
        &,dim_ll,wave_center_ll,coeff_4deg,dimsp,ew
        USE func_poly, ONLY: poly4_transform,ew_poly4, voigt_logg
        IMPLICIT NONE
        
        REAL(DP) :: sig,c,gamg,gaml,vv
        INTEGER(I4B) :: i,j,w_center_line
        REAL(DP),DIMENSION(dimsp), INTENT(INOUT) :: model
        REAL(DP),DIMENSION(:) :: parTGM,parABD
        REAL(DP),DIMENSION(4) :: st_pars
        REAL(DP), DIMENSION(70) :: X
        !initialize model
        model=1._dp

        st_pars(1:3)=parTGM(1:3)
        st_pars(4)=0.0001_dp        
        sig=ABS(parTGM(4))
        gamg=2.0_dp*sig*1.17741002252_dp!### sqrt(2.*log(2.))=1.17741002252
        c=1._dp+parTGM(5)/299792._dp

        WHERE(parABD>0.8) parABD=0.8_dp
        WHERE(parABD<-0.6) parABD=-0.6_dp

        CALL poly4_transform(X,parABD,st_pars)

        DO i=1,dim_ll
          IF(select_ll_mask(i)) THEN
            ew(i)=ew_poly4(X,ele_ll(i),coeff_4deg(:,i))
            w_center_line=wave_center_ll(i)
            !compute the parameters used by the voigt function
            gaml=gammaL(i,ew(i),gamG,st_pars(2),st_pars(1))
            DO j=INT(w_center_line-1,I4B),1,-1
             vv=voigt_logg(w_sp(j),wave_ll(i)*c,gamg,ew(i),gaml)
             IF(vv>1.E-4) THEN
              model(j)=model(j)-vv
             ELSE
              EXIT
             END IF
            END DO
            DO j=w_center_line,dimsp
             vv=voigt_logg(w_sp(j),wave_ll(i)*c,gamg,ew(i),gaml)
             IF(vv>1.E-4) THEN
              model(j)=model(j)-vv
             ELSE
              EXIT
             END IF
            END DO
          END IF
        END DO

        END SUBROUTINE make_model_ABDerr
!#########################
!        REAL(DP) FUNCTION gammaL(i_line,x0,x1,x2)
!        USE data_lib, ONLY: gamL_coeff,gamL_coeff_H,gamL_coeff_Na
!        USE share, ONLY: flag_lines
!        IMPLICIT NONE
!        REAL(DP), INTENT(IN) :: x0,x1,x2 !correspond to ew, gamG, logg
!        INTEGER(I4B), INTENT(IN) :: i_line
!        REAL(DP), DIMENSION(10) :: transf
!        REAL(DP), DIMENSION(4) :: x
!        INTEGER(I1B) :: i, j, count
!
!        x(1) = 1._dp
!        x(2) = x0
!        x(3) = x1
!        x(4) = x2
!
!        count = 0
!        DO i=1,4
!          DO j=i,4
!            count = count + 1_I1B
!            transf(count) = x(i)*x(j)          
!          ENDDO
!        ENDDO
!
!        !if the line is Ha or Hb
!        IF (flag_lines(i_line)==1) THEN
!            gammaL = DOT_PRODUCT(transf,gamL_coeff_H)
!        !else if the line belong to the NaI doublet
!        ELSE IF (flag_lines(i_line)==2) THEN
!            gammaL = DOT_PRODUCT(transf,gamL_coeff_Na)
!        ELSE !if the line is a generic line
!            gammaL = DOT_PRODUCT(transf,gamL_coeff)
!        END IF
!
!        gammaL = MAX(gammaL,0.01_dp)
!
!
!
!        END FUNCTION gammaL
!#########################
        REAL(DP) FUNCTION gammaL(i_line,x0,x1,x2,x3)
        USE data_lib, ONLY: gamL_coeff,gamL_coeff_H,gamL_coeff_Na
        USE share, ONLY: flag_lines
        IMPLICIT NONE
        REAL(DP), INTENT(IN) :: x0,x1,x2,x3 !correspond to ew, gamG, logg, teff
        INTEGER(I4B), INTENT(IN) :: i_line
        REAL(DP), DIMENSION(15) :: transf
        REAL(DP), DIMENSION(5) :: x
        INTEGER(I1B) :: i, j, count

        x(1) = 1._dp
        x(2) = x0
        x(3) = x1
        x(4) = x2
        x(5) = x3

        count = 0
        DO i=1,5
          DO j=i,5
            count = count + 1_I1B
            transf(count) = x(i)*x(j)
          ENDDO
        ENDDO

        !if the line is Ha or Hb
        IF (flag_lines(i_line)==1) THEN
            gammaL = DOT_PRODUCT(transf,gamL_coeff)
        !else if the line belong to the NaI doublet
        ELSE IF (flag_lines(i_line)==2) THEN
            gammaL = DOT_PRODUCT(transf,gamL_coeff)
        ELSE !if the line is a generic line
            gammaL = DOT_PRODUCT(transf,gamL_coeff)
        END IF


        gammaL = MAX(gammaL,0.01_dp)



        END FUNCTION gammaL

!#########################
!        REAL(DP) FUNCTION gammaG(x0,sig)
!        USE data_lib, ONLY: gamG_coeff
!        IMPLICIT NONE
!        REAL(DP), INTENT(IN) :: x0,sig !x0 correspond to ew
!        REAL(DP), DIMENSION(6) :: transf
!        REAL(DP), DIMENSION(3) :: x
!        INTEGER(I1B) :: i, j, count
!
!        REAL(DP) :: gamma_loc
!
!        x(1) = 1._dp
!        x(2) = x0
!        x(3) = sig
!
!        count = 0
!        DO i=1,3
!          DO j=i,3
!            count = count + 1_I1B
!            transf(count) = x(i)*x(j)          
!          ENDDO
!        ENDDO
!
!        gamma_loc = DOT_PRODUCT(transf,gamG_coeff)
!
!        IF (gamma_loc>sig*(2.0_dp*1.17741002252_dp)) THEN
!            gammaG = gamma_loc
!        ELSE
!            gammaG = 2.0_dp*sig*1.17741002252_dp
!        END IF
!
!        END FUNCTION gammaG
!#########################
END MODULE make_model
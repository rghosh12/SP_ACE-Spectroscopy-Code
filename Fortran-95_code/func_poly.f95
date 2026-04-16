!    This module contains a routine to compute the EWs of the lines from their
!    polynomial GCOG, and a function for the voigt profile.
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


MODULE func_poly

        USE num_type
        CONTAINS
!################################################        
        REAL(DP) FUNCTION ew_poly4(X,ele,coeff)
        USE share, ONLY: ele2meas, X_abd
        IMPLICIT NONE
        REAL(DP), DIMENSION(70), INTENT(IN) :: X
        REAL(DP), INTENT(IN) :: ele
        REAL(DP), DIMENSION(70), INTENT(IN) :: coeff
        REAL(DP), DIMENSION(70) :: X_pars
        INTEGER(I2B) :: ele_pos

        ele_pos=INT(MINLOC(ABS(ele-ele2meas),1),I2B)
        X_pars = X*X_abd(:,ele_pos)
        !initialize variables
        ew_poly4=0.0_dp
        ew_poly4 = DOT_PRODUCT(X_pars,coeff)

        IF(ew_poly4<1e-6) ew_poly4=1e-6_dp

        END FUNCTION ew_poly4        
!####################################
        SUBROUTINE poly4_transform(X,ABD,pars)
        USE share, ONLY: dim_ele,X_abd,dim_ele
        IMPLICIT NONE
        REAL(DP), DIMENSION(4), INTENT(IN) :: pars
        REAL(DP), DIMENSION(dim_ele), INTENT(IN) :: ABD
        REAL(DP), DIMENSION(70), INTENT(OUT) :: X
        INTEGER(I2B), DIMENSION(70) :: X_exp
        INTEGER(I1B), PARAMETER :: grado=4
        INTEGER(I2B) :: i1,i2,i3,i4
        INTEGER(I1B) :: count
        REAL(DP) :: a,b,c

        !initialize variables
        X=0.0_dp
        X_exp=0
        count=0

        DO i1=0,grado
         a=(pars(1)**i1)
         DO i2=0,grado-i1
          b=a*(pars(2)**i2)
          DO i3=0,grado-i2-i1
           c=b*(pars(3)**i3)
           DO i4=0,grado-i3-i2-i1
             count=count+1_I1B
             X(count)=c
             X_exp(count)=i4
           END DO
          END DO
         END DO
        END DO                                 

        !set the X_abd matrix
        DO i1=1,dim_ele
          DO i2=1,70
            X_abd(i2,i1)=ABD(i1)**X_exp(i2)
          END DO
        END DO

        END SUBROUTINE poly4_transform
!#######################
        REAL(DP) FUNCTION ew_poly6(X,coeff)
        IMPLICIT NONE
        REAL(DP), DIMENSION(84), INTENT(IN) :: X
        REAL(DP), DIMENSION(84), INTENT(IN) :: coeff

        !initialize variables
        ew_poly6=0.0_dp
        ew_poly6 = DOT_PRODUCT(X,coeff)

        IF(ew_poly6<1e-6) ew_poly6=1e-6_dp

        END FUNCTION ew_poly6
!#####################
        SUBROUTINE poly6_transform(X,pars)
        IMPLICIT NONE
        REAL(DP), DIMENSION(4), INTENT(IN) :: pars
        REAL(DP), DIMENSION(84), INTENT(OUT) :: X
        INTEGER(I1B), PARAMETER :: grado=6
        INTEGER(I1B) :: i1,i2,i3
        INTEGER(I1B) :: count
        REAL(DP) :: a,b

        !initialize variables
        X=0.0_dp
        count=0

        DO i1=0,grado
         a=(pars(1)**i1)
         DO i2=0,grado-i1
          b=a*(pars(2)**i2)
          DO i3=0,grado-i2-i1
             count=count+1_I1B
             X(count)=b*(pars(3)**i3)
          END DO
         END DO
        END DO                                 


        END SUBROUTINE poly6_transform
!#####################
        REAL(DP) FUNCTION voigt_logg(w,mu,gamG,ew,gamL)
        USE num_type
        IMPLICIT NONE
        REAL(DP),INTENT(IN) :: w,mu,gamG,ew,gamL
        REAL(DP),PARAMETER :: sqrtln2=0.832554611_dp,sqrtpi=1.772453851_dp
        REAL(DP),PARAMETER :: PIxsqrtPI=5.568327996831707_dp ! =pi*sqrtpi
        REAL(DP),PARAMETER :: sqrtpixsqrtln2=1.4756646266356057_dp ! =sqrtpi*sqrtln2
        REAL(DP),PARAMETER :: sqrtpi2=0.28209479177387814_dp ! =0.5/sqrtpi
        REAL(DP),PARAMETER :: sqrtln2x2=1.6651092223153954_dp ! =2.0_dp*sqrtln2
        REAL(DP), DIMENSION(4), PARAMETER :: A=(/-1.2150,-1.3509,-1.2150,-1.3509/)
        REAL(DP), DIMENSION(4), PARAMETER :: B=(/1.2359,0.3786,-1.2359,-0.3786/)
        REAL(DP), DIMENSION(4), PARAMETER :: C=(/-0.3085,0.5906,-0.3085,0.5906/)
        REAL(DP), DIMENSION(4), PARAMETER :: D=(/0.0210,-1.1858,-0.0210,1.1858/)
        REAL(DP), DIMENSION(4) :: V
        REAL(DP) :: sigmaL,aL,X,Y
        REAL(DP), PARAMETER :: pi=3.14159265358979_dp

        sigmaL=gamL*sqrtpi2
        aL=ew/(sigmaL*PIxsqrtPI)
        X=(w-mu)*sqrtln2x2/gamG

        Y=gamL*sqrtln2/gamG

        V=(C*(Y-A)+D*(X-B))/((Y-A)**2+(X-B)**2)
        voigt_logg=SUM(V)*(gamL*aL*sqrtpixsqrtln2/gamG)
        voigt_logg=MAX(voigt_logg,1e-6_dp)
        
        END FUNCTION voigt_logg
!###################################################

END MODULE func_poly

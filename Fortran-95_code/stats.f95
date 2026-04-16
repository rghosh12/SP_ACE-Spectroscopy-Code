!    
!    This module contains subroutines to compute statistics like mean, 
!    variance, standard deviation.
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

MODULE stats
        USE num_type
        IMPLICIT NONE

        INTERFACE avg_std_
                MODULE PROCEDURE avg_std_mask,avg_std
        END INTERFACE

        CONTAINS

        SUBROUTINE avg_std_mask(dat,mask,mean,var,sdev)
        USE num_type
        IMPLICIT NONE

        REAL(DP), DIMENSION(:), INTENT(IN) :: dat
        REAL(DP), INTENT(out) :: mean, var, sdev
        INTEGER(I4B) :: n
        REAL(DP), DIMENSION(size(dat)) :: discr_sq,discr
        REAL(DP) :: sum_disc
        LOGICAL, DIMENSION(size(dat)) :: mask

        n=INT(COUNT(mask),I4B)
        mean=sum(dat,mask)/n
        discr=dat-mean
        sum_disc=sum(discr,mask)
        discr_sq=discr**2
        var=sum(discr_sq,mask)/n
        sdev=sqrt(var)

        END SUBROUTINE avg_std_mask

        SUBROUTINE avg_std(dat,mean,var,sdev)
        IMPLICIT NONE

        REAL(DP), DIMENSION(:), INTENT(IN) :: dat
        REAL(DP), INTENT(out) :: mean, var, sdev
        INTEGER(I4B) :: n
        REAL(DP), DIMENSION(size(dat)) :: discr_sq,discr
        REAL(DP) :: sum_disc

        n=INT(SIZE(dat,1),I4B)
        mean=sum(dat)/n
        discr=dat-mean
        sum_disc=sum(discr)
        discr_sq=discr**2
        var=sum(discr_sq)/n
        sdev=sqrt(var)

        END SUBROUTINE avg_std

END MODULE stats

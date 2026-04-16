!    This module contains the routines that minimize the chi square
!    between observed spectrum and models. All the routines depending 
!    from the lmdif1 routine has been rewritten in F95 style.

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


!######################################################
MODULE minimize
        USE num_type
        USE error
        CONTAINS

!###############################################
SUBROUTINE chi_TGM_Q(dimsp,dim_var,tgmx,discrep)
       USE make_model, ONLY: make_model_TGM_quick
       USE share, ONLY: f_sp,cont,sig_noise,TGM_mask
       USE utils, ONLY: denormalize_pars
       IMPLICIT NONE
       INTEGER(I4B), INTENT(IN) :: dimsp,dim_var
       REAL(DP), DIMENSION(dimsp), INTENT(OUT) :: discrep
       REAL(DP), DIMENSION(dimsp) :: model,sp_norm
       REAL(DP), DIMENSION(dim_var), INTENT(IN) :: tgmx
       REAL(DP),DIMENSION(6) :: tgm_local

        CALL denormalize_pars(tgm_local, tgmx, TGM_mask)
        sp_norm=f_sp/tgm_local(6)/cont
        CALL make_model_TGM_quick(model,tgm_local)
        discrep=(sp_norm-model)/sig_noise

END SUBROUTINE chi_TGM_Q
!#########################
SUBROUTINE chi_TGM(dimsp,dim_var,tgmx,discrep)
       USE make_model, ONLY: make_model_TGM
       USE share, ONLY: f_sp_norm,sig_noise,TGM_mask
       USE utils, ONLY: denormalize_pars
       IMPLICIT NONE
       INTEGER(I4B), INTENT(IN) :: dimsp,dim_var
       REAL(DP), DIMENSION(dimsp), INTENT(OUT) :: discrep
       REAL(DP), DIMENSION(dimsp) :: model
       REAL(DP), DIMENSION(dim_var), INTENT(IN) :: tgmx
       REAL (DP),DIMENSION(6) :: tgm_local


        CALL denormalize_pars(tgm_local, tgmx, TGM_mask)
        CALL make_model_TGM(model,tgm_local)
        discrep=(f_sp_norm-model)/sig_noise

END SUBROUTINE chi_TGM
!#########################
SUBROUTINE chi_ABD(dimsp,dim_var,abdx,discrep)
       USE make_model, ONLY: make_model_ABD
       USE share, ONLY: f_sp_norm,ABD,sig_noise,ABD_mask,dim_ele
       IMPLICIT NONE
       INTEGER(I4B), INTENT(IN) :: dimsp,dim_var
       REAL(DP), DIMENSION(dimsp), INTENT(OUT) :: discrep
       REAL(DP), DIMENSION(dimsp) :: model
       REAL(DP), DIMENSION(dim_var), INTENT(IN) :: abdx
       REAL (DP),DIMENSION(dim_ele) :: abd_in,penal
       LOGICAL, DIMENSION(dim_ele) :: mask_internal

       mask_internal=.FALSE.

       abd_in=UNPACK(abdx,ABD_mask,abd)
       WHERE(abd_in>0.8) mask_internal=.TRUE.
       WHERE(abd_in<-0.6) mask_internal=.TRUE.
       penal=ABS(abd_in)
       model=1_DP

       CALL make_model_ABD(model,abd_in)
       discrep=(f_sp_norm-model)/sig_noise

       !give a penalization to the abundances that are >0.8 or <-0.6
       IF(COUNT(mask_internal)>0) THEN
         !this give to the chi^2 a minimum at 2 and -2. Later the abundances
         !that converged at -2 will be rejected and assigned a null value
         where(mask_internal) penal=abs(penal-2)/1.E6
         discrep=discrep+(SIGN(1.,discrep)*SUM(penal,mask_internal)/dimsp)
       END IF

END SUBROUTINE chi_ABD
!#########################
subroutine qrsolv ( n, r, ldr, ipvt, diag, qtb, x, sdiag )

!*****************************************************************************80
!
!! QRSOLV solves a rectangular linear system A*x=b in the least squares sense.
!
!  Discussion:
!
!    Given an M by N matrix A, an N by N diagonal matrix D,
!    and an M-vector B, the problem is to determine an X which
!    solves the system
!
!      A*X = B
!      D*X = 0
!
!    in the least squares sense.
!
!    This subroutine completes the solution of the problem
!    if it is provided with the necessary information from the
!    QR factorization, with column pivoting, of A.  That is, if
!    Q*P = Q*R, where P is a permutation matrix, Q has orthogonal
!    columns, and R is an upper triangular matrix with diagonal
!    elements of nonincreasing magnitude, then QRSOLV expects
!    the full upper triangle of R, the permutation matrix p,
!    and the first N components of Q'*B.
!
!    The system is then equivalent to
!
!      R*Z = Q'*B
!      P'*D*P*Z = 0
!
!    where X = P*Z.  If this system does not have full rank,
!    then a least squares solution is obtained.  On output QRSOLV
!    also provides an upper triangular matrix S such that
!
!      P'*(A'*A + D*D)*P = S'*S.
!
!    S is computed within QRSOLV and may be of separate interest.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of R.
!
!    Input/output, real ( kind = 8 ) R(LDR,N), the N by N matrix.
!    On input the full upper triangle must contain the full upper triangle
!    of the matrix R.  On output the full upper triangle is unaltered, and
!    the strict lower triangle contains the strict upper triangle
!    (transposed) of the upper triangular matrix S.
!
!    Input, integer ( kind = 4 ) LDR, the leading dimension of R, which must be
!    at least N.
!
!    Input, integer ( kind = 4 ) IPVT(N), defines the permutation matrix P such 
!    that A*P = Q*R.  Column J of P is column IPVT(J) of the identity matrix.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal elements of the matrix D.
!
!    Input, real ( kind = 8 ) QTB(N), the first N elements of the vector Q'*B.
!
!    Output, real ( kind = 8 ) X(N), the least squares solution.
!
!    Output, real ( kind = 8 ) SDIAG(N), the diagonal elements of the upper
!    triangular matrix S.
!
  implicit none

  integer(I4B),INTENT(IN) :: ldr,n
  integer(I4B) :: i,j,k,l,nsing
  integer(I4B), DIMENSION(n), INTENT(IN) :: ipvt

  real (DP) :: c,cotan,qtbpj,s,sum2,t,temp
  real (DP), DIMENSION(n), INTENT(IN) :: diag,qtb
  real (DP), DIMENSION(n), INTENT(OUT) ::  x,sdiag
  real (DP), DIMENSION(n) :: wa
  real (DP), DIMENSION(ldr,n), INTENT(INOUT) :: r

!
!  Copy R and Q'*B to preserve input and initialize S.
!
!  In particular, save the diagonal elements of R in X.
!
  do j = 1, n
    r(j:n,j) = r(j,j:n)
    x(j) = r(j,j)
  end do

  wa(1:n) = qtb(1:n)
!
!  Eliminate the diagonal matrix D using a Givens rotation.
!
  do j = 1, n
!
!  Prepare the row of D to be eliminated, locating the
!  diagonal element using P from the QR factorization.
!
    l = ipvt(j)

    if ( diag(l) /= 0.0D+00 ) then

      sdiag(j:n) = 0.0D+00
      sdiag(j) = diag(l)
!
!  The transformations to eliminate the row of D
!  modify only a single element of Q'*B
!  beyond the first N, which is initially zero.
!
      qtbpj = 0.0D+00

      do k = j, n
!
!  Determine a Givens rotation which eliminates the
!  appropriate element in the current row of D.
!
        if ( sdiag(k) /= 0.0D+00 ) then

          if ( abs ( r(k,k) ) < abs ( sdiag(k) ) ) then
            cotan = r(k,k) / sdiag(k)
            s = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * cotan**2 )
            c = s * cotan
          else
            t = sdiag(k) / r(k,k)
            c = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * t**2 )
            s = c * t
          end if
!
!  Compute the modified diagonal element of R and
!  the modified element of (Q'*B,0).
!
          r(k,k) = c * r(k,k) + s * sdiag(k)
          temp = c * wa(k) + s * qtbpj
          qtbpj = - s * wa(k) + c * qtbpj
          wa(k) = temp
!
!  Accumulate the tranformation in the row of S.
!
          do i = k+1_I2B, n
            temp = c * r(i,k) + s * sdiag(i)
            sdiag(i) = - s * r(i,k) + c * sdiag(i)
            r(i,k) = temp
          end do

        end if

      end do

    end if
!
!  Store the diagonal element of S and restore
!  the corresponding diagonal element of R.
!
    sdiag(j) = r(j,j)
    r(j,j) = x(j)

  end do
!
!  Solve the triangular system for Z.  If the system is
!  singular, then obtain a least squares solution.
!
  nsing = n

  do j = 1, n

    if ( sdiag(j) == 0.0D+00 .and. nsing == n ) then
      nsing = j - 1_I2B
    end if

    if ( nsing < n ) then
      wa(j) = 0.0D+00
    end if

  end do

  do j = nsing, 1, -1
    sum2 = dot_product ( wa(j+1:nsing), r(j+1:nsing,j) )
    wa(j) = ( wa(j) - sum2 ) / sdiag(j)
  end do
!
!  Permute the components of Z back to components of X.
!
  do j = 1, n
    l = ipvt(j)
    x(l) = wa(j)
  end do

  return
end subroutine qrsolv
!#################################
subroutine qrfac ( m, n, a, lda, pivot, ipvt, lipvt, rdiag, acnorm )

!*****************************************************************************80
!
!! QRFAC computes a QR factorization using Householder transformations.
!
!  Discussion:
!
!    This subroutine uses Householder transformations with column
!    pivoting (optional) to compute a QR factorization of the
!    M by N matrix A.  That is, QRFAC determines an orthogonal
!    matrix Q, a permutation matrix P, and an upper trapezoidal
!    matrix R with diagonal elements of nonincreasing magnitude,
!    such that A*P = Q*R.  The Householder transformation for
!    column K, K = 1,2,...,min(M,N), is of the form
!
!      I - ( 1 / U(K) ) * U * U'
!
!    where U has zeros in the first K-1 positions.  The form of
!    this transformation and the method of pivoting first
!    appeared in the corresponding LINPACK subroutine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, real ( kind = 8 ) A(LDA,N), the M by N array.
!    On input, A contains the matrix for which the QR factorization is to
!    be computed.  On output, the strict upper trapezoidal part of A contains
!    the strict upper trapezoidal part of R, and the lower trapezoidal
!    part of A contains a factored form of Q (the non-trivial elements of
!    the U vectors described above).
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A, which must
!    be no less than M.
!
!    Input, logical PIVOT, is TRUE if column pivoting is to be carried out.
!
!    Output, integer ( kind = 4 ) IPVT(LIPVT), defines the permutation matrix P 
!    such that A*P = Q*R.  Column J of P is column IPVT(J) of the identity 
!    matrix.  If PIVOT is false, IPVT is not referenced.
!
!    Input, integer ( kind = 4 ) LIPVT, the dimension of IPVT, which should 
!    be N if pivoting is used.
!
!    Output, real ( kind = 8 ) RDIAG(N), contains the diagonal elements of R.
!
!    Output, real ( kind = 8 ) ACNORM(N), the norms of the corresponding
!    columns of the input matrix A.  If this information is not needed,
!    then ACNORM can coincide with RDIAG.
!
  implicit none

  integer (I4B), INTENT(IN) :: m,n,lda,lipvt


  real (DP), DIMENSION(lda,n), INTENT(INOUT) :: a
  real (DP), DIMENSION(n), INTENT(OUT) :: acnorm
  real (DP) :: epsmch, ajnorm
  integer (I4B) :: i4_temp,j,k,kmax,minmn
  integer (I4B), DIMENSION(lipvt), INTENT(OUT) :: ipvt
  logical, INTENT(IN) :: pivot
  real (DP), DIMENSION(m) :: r8_temp
  real (DP), DIMENSION(n) :: rdiag,wa
  real (DP) :: temp

  epsmch = epsilon ( epsmch )
!
!  Compute the initial column norms and initialize several arrays.
!
  do j = 1, n
    acnorm(j) = enorm ( m, a(1:m,j) )
  end do

  rdiag(1:n) = acnorm(1:n)
  wa(1:n) = acnorm(1:n)

  if ( pivot ) then
    do j = 1, n
      ipvt(j) = j
    end do
  end if
!
!  Reduce A to R with Householder transformations.
!
  minmn = min ( m, n )

  do j = 1, minmn
!
!  Bring the column of largest norm into the pivot position.
!
    if ( pivot ) then

      kmax = j

      do k = j, n
        if ( rdiag(k) > rdiag(kmax) ) then
          kmax = k
        end if
      end do

      if ( kmax /= j ) then

        r8_temp(1:m) = a(1:m,j)
        a(1:m,j)     = a(1:m,kmax)
        a(1:m,kmax)  = r8_temp(1:m)

        rdiag(kmax) = rdiag(j)
        wa(kmax) = wa(j)

        i4_temp    = ipvt(j)
        ipvt(j)    = ipvt(kmax)
        ipvt(kmax) = i4_temp

      end if

    end if
!
!  Compute the Householder transformation to reduce the
!  J-th column of A to a multiple of the J-th unit vector.
!
    ajnorm = enorm ( m-j+1_I2B, a(j,j) )

    if ( ajnorm /= 0.0D+00 ) then

      if ( a(j,j) < 0.0D+00 ) then
        ajnorm = -ajnorm
      end if

      a(j:m,j) = a(j:m,j) / ajnorm
      a(j,j) = a(j,j) + 1.0D+00
!
!  Apply the transformation to the remaining columns and update the norms.
!
      do k = j+1_I2B, n

        temp = dot_product ( a(j:m,j), a(j:m,k) ) / a(j,j)

        a(j:m,k) = a(j:m,k) - temp * a(j:m,j)

        if ( pivot .and. rdiag(k) /= 0.0D+00 ) then

          temp = a(j,k) / rdiag(k)
          rdiag(k) = rdiag(k) * sqrt ( max ( 0.0D+00, 1.0D+00-temp**2 ) )

          if ( 0.05D+00 * ( rdiag(k) / wa(k) )**2 <= epsmch ) then
            rdiag(k) = enorm ( m-j, a(j+1,k) )
            wa(k) = rdiag(k)
          end if

        end if

      end do

    end if

    rdiag(j) = -ajnorm

  end do

  return
end subroutine qrfac
!#################################
subroutine fdjac2 ( fcn, m, n, x, fvec, fjac, ldfjac, iflag, epsfcn)

!  Corrado comment: removed the iflag in the functions because useless for our purpose.

!*****************************************************************************80
!
!! FDJAC2 estimates an M by N jacobian matrix using forward differences.
!
!  Discussion:
!
!    This subroutine computes a forward-difference approximation
!    to the M by N jacobian matrix associated with a specified
!    problem of M functions in N variables.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions.  The routine should have the form:
!
!      subroutine fcn ( m, n, x, fvec, iflag )
!      integer ( kind = 4 ) n
!      real fvec(m)
!      integer ( kind = 4 ) iflag
!      real x(n)
!
!    The value of IFLAG should not be changed by FCN unless
!    the user wants to terminate execution of the routine.
!    In this case set IFLAG to a negative integer.
!
!
!
!    Input, integer ( kind = 4 ) M, is the number of functions.
!
!    Input, integer ( kind = 4 ) N, is the number of variables.  
!    N must not exceed M.
!
!    Input, real ( kind = 8 ) X(N), the point where the jacobian is evaluated.
!
!    Input, real ( kind = 8 ) FVEC(M), the functions evaluated at X.
!
!    Output, real ( kind = 8 ) FJAC(LDFJAC,N), the M by N approximate
!    jacobian matrix.
!
!    Input, integer ( kind = 4 ) LDFJAC, the leading dimension of FJAC, 
!    which must not be less than M.
!
!    Output, integer ( kind = 4 ) IFLAG, is an error flag returned by FCN.  
!    If FCN returns a nonzero value of IFLAG, then this routine returns 
!    immediately to the calling program, with the value of IFLAG.
!
!    Input, real ( kind = 8 ) EPSFCN, is used in determining a suitable
!    step length for the forward-difference approximation.  This approximation
!    assumes that the relative errors in the functions are of the order of
!    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
!    the relative errors in the functions are of the order of the machine
!    precision.
!
  implicit none

  integer (I4B), INTENT(IN) :: m,n,ldfjac
  integer (I2B), INTENT(OUT) :: iflag
  integer (I4B) :: j
  real (DP), INTENT(IN) :: epsfcn
  real (DP), DIMENSION(ldfjac,n), INTENT(OUT) :: fjac
  real (DP), DIMENSION(n) :: x
  real (DP), DIMENSION(m), INTENT(IN) :: fvec
  real (DP) :: eps,epsmch,h,temp
  real (DP), DIMENSION(m) :: wa
  external fcn

  epsmch = epsilon ( epsmch )

  eps = sqrt ( max ( epsfcn, epsmch ) )

  do j = 1, n

    temp = x(j)
    h = eps * abs ( temp )
    if ( h == 0.0D+00 ) then
      h = eps
    end if

    x(j) = temp + h
    call fcn ( m, n, x, wa)

    if ( iflag < 0 ) then
      exit
    end if

    x(j) = temp
    fjac(1:m,j) = ( wa(1:m) - fvec(1:m) ) / h

  end do

  return
end subroutine fdjac2
!#################################
real(DP) function enorm ( n, x )

!*****************************************************************************80
!
!! ENORM computes the Euclidean norm of a vector.
!
!  Discussion:
!
!    This is an extremely simplified version of the original ENORM
!    routine, which has been renamed to "ENORM2".
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the length of the vector.
!
!    Input, real ( kind = 8 ) X(N), the vector whose norm is desired.
!
!    Output, real ( kind = 8 ) ENORM, the Euclidean norm of the vector.
!
  implicit none

  integer (I4B), INTENT(IN) :: n
  real (DP), DIMENSION(n), INTENT(IN) :: x

  enorm = sqrt ( sum ( x**2 ))

  return
end function enorm
!#################################
subroutine lmdif ( fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, &
  diag, mode, factor, nprint, info, nfev, fjac, ldfjac, ipvt, qtf)

!Corrado comment: I did some changes to this routine. I removed the iflag in
!the functions because useless for our purpose.  

!*****************************************************************************80
!
!! LMDIF minimizes M functions in N variables by the Levenberg-Marquardt method.
!
!  Discussion:
!
!    LMDIF minimizes the sum of the squares of M nonlinear functions in
!    N variables by a modification of the Levenberg-Marquardt algorithm.
!    The user must provide a subroutine which calculates the functions.
!    The jacobian is then calculated by a forward-difference approximation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions.  The routine should have the form:
!
!      subroutine fcn ( m, n, x, fvec, iflag )
!      integer ( kind = 4 ) m
!      integer ( kind = 4 ) n
!
!      real fvec(m)
!      integer ( kind = 4 ) iflag
!      real x(n)
!
!    The value of IFLAG should not be changed by FCN unless
!    the user wants to terminate execution of the routine.
!    In this case set IFLAG to a negative integer.
!
!
!    Input, integer ( kind = 4 ) M, the number of functions.
!
!    Input, integer ( kind = 4 ) N, the number of variables.  
!    N must not exceed M.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real ( kind = 8 ) FVEC(M), the functions evaluated at the output X.
!
!    Input, real ( kind = 8 ) FTOL.  Termination occurs when both the actual
!    and predicted relative reductions in the sum of squares are at most FTOL.
!    Therefore, FTOL measures the relative error desired in the sum of
!    squares.  FTOL should be nonnegative.
!
!    Input, real ( kind = 8 ) XTOL.  Termination occurs when the relative error
!    between two consecutive iterates is at most XTOL.  Therefore, XTOL
!    measures the relative error desired in the approximate solution.  XTOL
!    should be nonnegative.
!
!    Input, real ( kind = 8 ) GTOL. termination occurs when the cosine of the
!    angle between FVEC and any column of the jacobian is at most GTOL in
!    absolute value.  Therefore, GTOL measures the orthogonality desired
!    between the function vector and the columns of the jacobian.  GTOL should
!    be nonnegative.
!
!    Input, integer ( kind = 4 ) MAXFEV.  Termination occurs when the number of
!    calls to FCN is at least MAXFEV by the end of an iteration.
!
!    Input, real ( kind = 8 ) EPSFCN, is used in determining a suitable step 
!    length for the forward-difference approximation.  This approximation 
!    assumes that the relative errors in the functions are of the order of 
!    EPSFCN.  If EPSFCN is less than the machine precision, it is assumed that
!    the relative errors in the functions are of the order of the machine
!    precision.
!
!    Input/output, real ( kind = 8 ) DIAG(N).  If MODE = 1, then DIAG is set
!    internally.  If MODE = 2, then DIAG must contain positive entries that
!    serve as multiplicative scale factors for the variables.
!
!    Input, integer ( kind = 4 ) MODE, scaling option.
!    1, variables will be scaled internally.
!    2, scaling is specified by the input DIAG vector.
!
!    Input, real ( kind = 8 ) FACTOR, determines the initial step bound.  
!    This bound is set to the product of FACTOR and the euclidean norm of
!    DIAG*X if nonzero, or else to FACTOR itself.  In most cases, FACTOR should 
!    lie in the interval (0.1, 100) with 100 the recommended value.
!
!    Input, integer ( kind = 4 ) NPRINT, enables controlled printing of iterates
!    if it is positive.  In this case, FCN is called with IFLAG = 0 at the
!    beginning of the first iteration and every NPRINT iterations thereafter
!    and immediately prior to return, with X and FVEC available
!    for printing.  If NPRINT is not positive, no special calls
!    of FCN with IFLAG = 0 are made.
!
!    Output, integer ( kind = 4 ) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See description
!    of FCN.  Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, both actual and predicted relative reductions in the sum of squares
!       are at most FTOL.
!    2, relative error between two consecutive iterates is at most XTOL.
!    3, conditions for INFO = 1 and INFO = 2 both hold.
!    4, the cosine of the angle between FVEC and any column of the jacobian
!       is at most GTOL in absolute value.
!    5, number of calls to FCN has reached or exceeded MAXFEV.
!    6, FTOL is too small.  No further reduction in the sum of squares
!       is possible.
!    7, XTOL is too small.  No further improvement in the approximate
!       solution X is possible.
!    8, GTOL is too small.  FVEC is orthogonal to the columns of the
!       jacobian to machine precision.
!
!    Output, integer ( kind = 4 ) NFEV, the number of calls to FCN.
!
!    Output, real ( kind = 8 ) FJAC(LDFJAC,N), an M by N array.  The upper
!    N by N submatrix of FJAC contains an upper triangular matrix R with
!    diagonal elements of nonincreasing magnitude such that
!
!      P' * ( JAC' * JAC ) * P = R' * R,
!
!    where P is a permutation matrix and JAC is the final calculated jacobian.
!    Column J of P is column IPVT(J) of the identity matrix.  The lower
!    trapezoidal part of FJAC contains information generated during
!    the computation of R.
!
!    Input, integer ( kind = 4 ) LDFJAC, the leading dimension of FJAC.
!    LDFJAC must be at least M.
!
!    Output, integer ( kind = 4 ) IPVT(N), defines a permutation matrix P such
!    that JAC * P = Q * R, where JAC is the final calculated jacobian, Q is
!    orthogonal (not stored), and R is upper triangular with diagonal
!    elements of nonincreasing magnitude.  Column J of P is column IPVT(J)
!    of the identity matrix.
!
!    Output, real ( kind = 8 ) QTF(N), the first N elements of Q'*FVEC.
!
  implicit none

  integer (I4B), INTENT(IN) :: ldfjac,m,n
  integer (I4B) :: i,iter,j,l
  integer (I2B) :: iflag
  integer (I4B), INTENT(OUT) :: nfev
  integer (I2B), INTENT(OUT) :: info
  integer (I4B),DIMENSION(n),INTENT(OUT) :: ipvt
  integer (I4B),INTENT(IN) :: maxfev,mode,nprint

  real (DP),INTENT(IN) :: epsfcn,factor,ftol,gtol,xtol
  real (DP) :: actred,delta,dirder,epsmch,fnorm,fnorm1,gnorm,par
  real (DP) :: pnorm,prered,ratio,sum2,temp,temp1,temp2,xnorm
  real (DP), DIMENSION(n),INTENT(INOUT) :: x,diag
  real (DP),DIMENSION(ldfjac,n),INTENT(OUT) :: fjac
  real (DP),DIMENSION(m),INTENT(OUT) :: fvec
  real (DP),DIMENSION(n),INTENT(OUT) :: qtf
  real (DP),DIMENSION(n) :: wa1,wa2,wa3
  real (DP),DIMENSION(m) :: wa4
  external fcn

  epsmch = epsilon ( epsmch )

  xnorm = 0
  info = 0
  iflag = 0
  nfev = 0


  if ( n <= 0 ) then
    go to 300
  else if ( m < n ) then
    go to 300
  else if ( ldfjac < m ) then
    go to 300
  else if ( ftol < 0.0D+00 ) then
    go to 300
  else if ( xtol < 0.0D+00 ) then
    go to 300
  else if ( gtol < 0.0D+00 ) then
    go to 300
  else if ( maxfev <= 0 ) then
    go to 300
  else if ( factor <= 0.0D+00 ) then
    go to 300
  end if

  if ( mode == 2 ) then
    do j = 1, n
      if ( diag(j) <= 0.0D+00 ) then
        go to 300
      end if
    end do
  end if
!
!  Evaluate the function at the starting point and calculate its norm.
!


  iflag = 1
  call fcn ( m, n, x, fvec)
  nfev = 1

  if ( iflag < 0 ) then
    go to 300
  end if

  fnorm = enorm ( m, fvec )
!
!  Initialize Levenberg-Marquardt parameter and iteration counter.
!
  par = 0.0D+00
  iter = 1
!
!  Beginning of the outer loop.
!
30 continue
!
!  Calculate the jacobian matrix.
!
  iflag = 2
  call fdjac2 ( fcn, m, n, x, fvec, fjac, ldfjac, iflag, epsfcn)
  nfev = nfev + n

  if ( iflag < 0 ) then
    go to 300
  end if
!
!  If requested, call FCN to enable printing of iterates.
!
     if ( 0 < nprint ) then
       iflag = 0
       if ( mod ( iter-1_I2B, nprint ) == 0 ) then
         call fcn ( m, n, x, fvec)
       end if
       if ( iflag < 0 ) then
         go to 300
       end if
     end if
!
!  Compute the QR factorization of the jacobian.
!
     call qrfac ( m, n, fjac, ldfjac, .true., ipvt, n, wa1, wa2 )
!
!  On the first iteration and if MODE is 1, scale according
!  to the norms of the columns of the initial jacobian.
!
     if ( iter == 1 ) then

       if ( mode /= 2 ) then
         diag(1:n) = wa2(1:n)
         do j = 1, n
           if ( wa2(j) == 0.0D+00 ) then
             diag(j) = 1.0D+00
           end if
         end do
       end if
!
!  On the first iteration, calculate the norm of the scaled X
!  and initialize the step bound DELTA.
!
       wa3(1:n) = diag(1:n) * x(1:n)
       xnorm = enorm ( n, wa3 )
       delta = factor * xnorm
       if ( delta == 0.0D+00 ) then
         delta = factor
       end if
     end if
!
!  Form Q' * FVEC and store the first N components in QTF.
!
     wa4(1:m) = fvec(1:m)

     do j = 1, n

       if ( fjac(j,j) /= 0.0D+00 ) then
         sum2 = dot_product ( wa4(j:m), fjac(j:m,j) )
         temp = - sum2 / fjac(j,j)
         wa4(j:m) = wa4(j:m) + fjac(j:m,j) * temp
       end if

       fjac(j,j) = wa1(j)
       qtf(j) = wa4(j)

     end do
!
!  Compute the norm of the scaled gradient.
!
     gnorm = 0.0D+00

     if ( fnorm /= 0.0D+00 ) then

       do j = 1, n

         l = ipvt(j)

         if ( wa2(l) /= 0.0D+00 ) then
           sum2 = 0.0D+00
           do i = 1, j
             sum2 = sum2 + fjac(i,j) * ( qtf(i) / fnorm )
           end do
           gnorm = max ( gnorm, abs ( sum2 / wa2(l) ) )
         end if

       end do

     end if
!
!  Test for convergence of the gradient norm.
!
     if ( gnorm <= gtol ) then
       info = 4
       go to 300
     end if
!
!  Rescale if necessary.
!
     if ( mode /= 2 ) then
       do j = 1, n
         diag(j) = max ( diag(j), wa2(j) )
       end do
     end if
!
!  Beginning of the inner loop.
!
200  continue
!
!  Determine the Levenberg-Marquardt parameter.
!
        call lmpar ( n, fjac, ldfjac, ipvt, diag, qtf, delta, par, wa1, wa2 )
!
!  Store the direction P and X + P.
!  Calculate the norm of P.
!
        wa1(1:n) = -wa1(1:n)
        wa2(1:n) = x(1:n) + wa1(1:n)
        wa3(1:n) = diag(1:n) * wa1(1:n)

        pnorm = enorm ( n, wa3 )
!
!  On the first iteration, adjust the initial step bound.
!
        if ( iter == 1 ) then
          delta = min ( delta, pnorm )
        end if
!
!  Evaluate the function at X + P and calculate its norm.
!
        iflag = 1
        call fcn ( m, n, wa2, wa4)
        nfev = nfev + 1_I2B
        if ( iflag < 0 ) then
          go to 300
        end if
        fnorm1 = enorm ( m, wa4 )
!
!  Compute the scaled actual reduction.
!
        if ( 0.1D+00 * fnorm1 < fnorm ) then
          actred = 1.0D+00 - ( fnorm1 / fnorm )**2
        else
          actred = -1.0D+00
        end if
!
!  Compute the scaled predicted reduction and the scaled directional derivative.
!
        do j = 1, n
          wa3(j) = 0.0D+00
          l = ipvt(j)
          temp = wa1(l)
          wa3(1:j) = wa3(1:j) + fjac(1:j,j) * temp
        end do

        temp1 = enorm ( n, wa3 ) / fnorm
        temp2 = ( sqrt ( par ) * pnorm ) / fnorm
        prered = temp1**2 + temp2**2 / 0.5D+00
        dirder = - ( temp1**2 + temp2**2 )
!
!  Compute the ratio of the actual to the predicted reduction.
!
        ratio = 0.0D+00
        if ( prered /= 0.0D+00 ) then
          ratio = actred / prered
        end if
!
!  Update the step bound.
!
        if ( ratio <= 0.25D+00 ) then

           if ( actred >= 0.0D+00 ) then
             temp = 0.5D+00
           endif

           if ( actred < 0.0D+00 ) then
             temp = 0.5D+00 * dirder / ( dirder + 0.5D+00 * actred )
           end if

           if ( 0.1D+00 * fnorm1 >= fnorm .or. temp < 0.1D+00 ) then
             temp = 0.1D+00
           end if

           delta = temp * min ( delta, pnorm / 0.1D+00  )
           par = par / temp

        else

           if ( par == 0.0D+00 .or. ratio >= 0.75D+00 ) then
             delta = 2.0D+00 * pnorm
             par = 0.5D+00 * par
           end if

        end if
!
!  Test for successful iteration.
!

!
!  Successful iteration. update X, FVEC, and their norms.
!
        if ( 0.0001D+00 <= ratio ) then
          x(1:n) = wa2(1:n)
          wa2(1:n) = diag(1:n) * x(1:n)
          fvec(1:m) = wa4(1:m)
          xnorm = enorm ( n, wa2 )
          fnorm = fnorm1
          iter = iter + 1_I2B
        end if
!
!  Tests for convergence.
!
        if ( abs ( actred) <= ftol .and. prered <= ftol &
          .and. 0.5D+00 * ratio <= 1.0D+00 ) then
          info = 1
        end if

        if ( delta <= xtol * xnorm ) then
          info = 2
        end if

        if ( abs ( actred) <= ftol .and. prered <= ftol &
          .and. 0.5D+00 * ratio <= 1.0D+00 .and. info == 2 ) info = 3

        if ( info /= 0 ) then
          go to 300
        end if
!
!  Tests for termination and stringent tolerances.
!
        if ( maxfev <= nfev ) then
          info = 5
        end if

        if ( abs ( actred) <= epsmch .and. prered <= epsmch &
          .and. 0.5D+00 * ratio <= 1.0D+00 ) then
          info = 6
        end if

        if ( delta <= epsmch * xnorm ) then
          info = 7
        end if

        if ( gnorm <= epsmch ) then
          info = 8
        end if

        if ( info /= 0 ) then
          go to 300
        end if
!
!  End of the inner loop.  Repeat if iteration unsuccessful.
!
        if ( ratio < 0.0001D+00 ) then
          go to 200
        end if
!
!  End of the outer loop.
!
     go to 30

300 continue
!
!  Termination, either normal or user imposed.
!
  if ( iflag < 0 ) then
    info = iflag
  end if

  iflag = 0

  if ( 0 < nprint ) then
    call fcn ( m, n, x, fvec)
  end if

  return
end subroutine lmdif
!#################################
subroutine lmdif1 ( fcn, m, n, x, fvec, tol, info)

!Corrado comment: I did some changes to this routine. I removed the iflag in
!the functions because useless. 

!*****************************************************************************80
!
!! LMDIF1 minimizes M functions in N variables using Levenberg-Marquardt method.
!
!  Discussion:
!
!    LMDIF1 minimizes the sum of the squares of M nonlinear functions in
!    N variables by a modification of the Levenberg-Marquardt algorithm.
!    This is done by using the more general least-squares solver LMDIF.
!    The user must provide a subroutine which calculates the functions.
!    The jacobian is then calculated by a forward-difference approximation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions.  The routine should have the form:
!
!      subroutine fcn ( m, n, x, fvec, iflag )
!      integer ( kind = 4 ) n
!      real fvec(m)
!      integer ( kind = 4 ) iflag
!      real x(n)
!
!
!    The value of IFLAG should not be changed by FCN unless
!    the user wants to terminate execution of the routine.
!    In this case set IFLAG to a negative integer.
!
!    Input, integer ( kind = 4 ) M, the number of functions.
!
!    Input, integer ( kind = 4 ) N, the number of variables.  
!    N must not exceed M.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real ( kind = 8 ) FVEC(M), the functions evaluated at the output X.
!
!    Input, real ( kind = 8 ) TOL.  Termination occurs when the algorithm
!    estimates either that the relative error in the sum of squares is at
!    most TOL or that the relative error between X and the solution is at
!    most TOL.  TOL should be nonnegative.
!
!    Output, integer ( kind = 4 ) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See description
!    of FCN.  Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, algorithm estimates that the relative error in the sum of squares
!       is at most TOL.
!    2, algorithm estimates that the relative error between X and the
!       solution is at most TOL.
!    3, conditions for INFO = 1 and INFO = 2 both hold.
!    4, FVEC is orthogonal to the columns of the jacobian to machine precision.
!    5, number of calls to FCN has reached or exceeded 200*(N+1).
!    6, TOL is too small.  No further reduction in the sum of squares
!       is possible.
!    7, TOL is too small.  No further improvement in the approximate
!       solution X is possible.
!
  implicit none

  integer (I4B), INTENT(IN) :: m,n
  integer (I2B), INTENT(OUT) :: info
  integer (I4B) :: ldfjac,mode,nfev,nprint
  integer (I4B) :: maxfev
  integer (I4B),DIMENSION(n) :: ipvt(n)

  real (DP),INTENT(IN) :: tol
  real (DP),DIMENSION(n),INTENT(INOUT) :: x
  real (DP),DIMENSION(m),INTENT(OUT) :: fvec
  real (DP),DIMENSION(n) :: diag,qtf
  real (DP) :: epsfcn,factor,ftol,gtol,xtol
  real (DP),DIMENSION(m,n) :: fjac

  external fcn

  info = 0

  if ( n <= 0 ) then
    return
  else if ( m < n ) then
    return
  else if ( tol < 0.0D+00 ) then
    return
  end if

  factor = 100.0D+00
  maxfev = 200_I4B * ( n + 1_I4B )
  ftol = tol
  xtol = tol
  gtol = 0.0D+00
  epsfcn = 0.0D+00
  mode = 1
  nprint = 0
  ldfjac = m



  call lmdif ( fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, &
    diag, mode, factor, nprint, info, nfev, fjac, ldfjac, ipvt, qtf)

  if ( info == 8 ) then
    info = 4
  end if

  return
end subroutine lmdif1

subroutine lmpar ( n, r, ldr, ipvt, diag, qtb, delta, par, x, sdiag )

!*****************************************************************************80
!
!! LMPAR computes a parameter for the Levenberg-Marquardt method.
!
!  Discussion:
!
!    Given an M by N matrix A, an N by N nonsingular diagonal
!    matrix D, an M-vector B, and a positive number DELTA,
!    the problem is to determine a value for the parameter
!    PAR such that if X solves the system
!
!      A*X = B,
!      sqrt ( PAR ) * D * X = 0,
!
!    in the least squares sense, and DXNORM is the euclidean
!    norm of D*X, then either PAR is zero and
!
!      ( DXNORM - DELTA ) <= 0.1 * DELTA,
!
!    or PAR is positive and
!
!      abs ( DXNORM - DELTA) <= 0.1 * DELTA.
!
!    This subroutine completes the solution of the problem
!    if it is provided with the necessary information from the
!    QR factorization, with column pivoting, of A.  That is, if
!    A*P = Q*R, where P is a permutation matrix, Q has orthogonal
!    columns, and R is an upper triangular matrix with diagonal
!    elements of nonincreasing magnitude, then LMPAR expects
!    the full upper triangle of R, the permutation matrix P,
!    and the first N components of Q'*B.  On output
!    LMPAR also provides an upper triangular matrix S such that
!
!      P' * ( A' * A + PAR * D * D ) * P = S'* S.
!
!    S is employed within LMPAR and may be of separate interest.
!
!    Only a few iterations are generally needed for convergence
!    of the algorithm.  If, however, the limit of 10 iterations
!    is reached, then the output PAR will contain the best
!    value obtained so far.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 2014
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of R.
!
!    Input/output, real ( kind = 8 ) R(LDR,N),the N by N matrix.  The full
!    upper triangle must contain the full upper triangle of the matrix R.
!    On output the full upper triangle is unaltered, and the strict lower
!    triangle contains the strict upper triangle (transposed) of the upper
!    triangular matrix S.
!
!    Input, integer ( kind = 4 ) LDR, the leading dimension of R.  LDR must be
!    no less than N.
!
!    Input, integer ( kind = 4 ) IPVT(N), defines the permutation matrix P 
!    such that A*P = Q*R.  Column J of P is column IPVT(J) of the 
!    identity matrix.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal elements of the matrix D.
!
!    Input, real ( kind = 8 ) QTB(N), the first N elements of the vector Q'*B.
!
!    Input, real ( kind = 8 ) DELTA, an upper bound on the euclidean norm
!    of D*X.  DELTA should be positive.
!
!    Input/output, real ( kind = 8 ) PAR.  On input an initial estimate of the
!    Levenberg-Marquardt parameter.  On output the final estimate.
!    PAR should be nonnegative.
!
!    Output, real ( kind = 8 ) X(N), the least squares solution of the system
!    A*X = B, sqrt(PAR)*D*X = 0, for the output value of PAR.
!
!    Output, real ( kind = 8 ) SDIAG(N), the diagonal elements of the upper
!    triangular matrix S.
!
  implicit none

  integer (I4B),INTENT(IN) :: ldr,n
  integer (I4B),DIMENSION(n),INTENT(IN) :: ipvt
  integer (I4B) :: iter,j,k,l,nsing

  real (DP),INTENT(IN) :: delta
  real (DP),DIMENSION(n),INTENT(IN) :: diag,qtb
  real (DP),INTENT(INOUT) :: par
  real (DP),DIMENSION(ldr,n),INTENT(INOUT) :: r
  real (DP),DIMENSION(n),INTENT(OUT) :: sdiag,x
  real (DP) :: dwarf,dxnorm,gnorm,fp,parc,parl,paru,sum2,temp
  real (DP),DIMENSION(n) :: wa1,wa2

!
!  DWARF is the smallest positive magnitude.
!
  dwarf = tiny ( dwarf )
!
!  Compute and store in X the Gauss-Newton direction.
!
!  If the jacobian is rank-deficient, obtain a least squares solution.
!
  nsing = n

  do j = 1, n
    wa1(j) = qtb(j)
    if ( r(j,j) == 0.0D+00 .and. nsing == n ) then
      nsing = j - 1_I2B
    end if
    if ( nsing < n ) then
      wa1(j) = 0.0D+00
    end if
  end do

  do k = 1, nsing
    j = nsing - k + 1_I2B
    wa1(j) = wa1(j) / r(j,j)
    temp = wa1(j)
    wa1(1:j-1) = wa1(1:j-1) - r(1:j-1,j) * temp
  end do

  do j = 1, n
    l = ipvt(j)
    x(l) = wa1(j)
  end do
!
!  Initialize the iteration counter.
!  Evaluate the function at the origin, and test
!  for acceptance of the Gauss-Newton direction.
!
  iter = 0
  wa2(1:n) = diag(1:n) * x(1:n)
  dxnorm = enorm ( n, wa2 )
  fp = dxnorm - delta

  if ( fp <= 0.1D+00 * delta ) then
    if ( iter == 0 ) then
      par = 0.0D+00
    end if
    return
  end if
!
!  If the jacobian is not rank deficient, the Newton
!  step provides a lower bound, PARL, for the zero of
!  the function.
!
!  Otherwise set this bound to zero.
!
  parl = 0.0D+00

  if ( n <= nsing ) then

    do j = 1, n
      l = ipvt(j)
      wa1(j) = diag(l) * ( wa2(l) / dxnorm )
    end do

    do j = 1, n
      sum2 = dot_product ( wa1(1:j-1), r(1:j-1,j) )
      wa1(j) = ( wa1(j) - sum2 ) / r(j,j)
    end do

    temp = enorm ( n, wa1 )
    parl = ( ( fp / delta ) / temp ) / temp

  end if
!
!  Calculate an upper bound, PARU, for the zero of the function.
!
  do j = 1, n
    sum2 = dot_product ( qtb(1:j), r(1:j,j) )
    l = ipvt(j)
    wa1(j) = sum2 / diag(l)
  end do

  gnorm = enorm ( n, wa1 )
  paru = gnorm / delta

  if ( paru == 0.0D+00 ) then
    paru = dwarf / min ( delta, 0.1D+00 )
  end if
!
!  If the input PAR lies outside of the interval (PARL, PARU),
!  set PAR to the closer endpoint.
!
  par = max ( par, parl )
  par = min ( par, paru )
  if ( par == 0.0D+00 ) then
    par = gnorm / dxnorm
  end if
!
!  Beginning of an iteration.
!
  do
 
    iter = iter + 1_I2B
!
!  Evaluate the function at the current value of PAR.
!
    if ( par == 0.0D+00 ) then
      par = max ( dwarf, 0.001D+00 * paru )
    end if

    wa1(1:n) = sqrt ( par ) * diag(1:n)

    call qrsolv ( n, r, ldr, ipvt, wa1, qtb, x, sdiag )

    wa2(1:n) = diag(1:n) * x(1:n)
    dxnorm = enorm ( n, wa2 )
    temp = fp
    fp = dxnorm - delta
!
!  If the function is small enough, accept the current value of PAR.
!
    if ( abs ( fp ) <= 0.1D+00 * delta ) then
      exit
    end if
!
!  Test for the exceptional cases where PARL
!  is zero or the number of iterations has reached 10.
!
    if ( parl == 0.0D+00 .and. fp <= temp .and. temp < 0.0D+00 ) then
      exit
    else if ( iter == 10 ) then
      exit
    end if
!
!  Compute the Newton correction.
!
    do j = 1, n
      l = ipvt(j)
      wa1(j) = diag(l) * ( wa2(l) / dxnorm )
    end do

    do j = 1, n
      wa1(j) = wa1(j) / sdiag(j)
      temp = wa1(j)
      wa1(j+1:n) = wa1(j+1:n) - r(j+1:n,j) * temp
    end do

    temp = enorm ( n, wa1 )
    parc = ( ( fp / delta ) / temp ) / temp
!
!  Depending on the sign of the function, update PARL or PARU.
!
    if ( 0.0D+00 < fp ) then
      parl = max ( parl, par )
    else if ( fp < 0.0D+00 ) then
      paru = min ( paru, par )
    end if
!
!  Compute an improved estimate for PAR.
!
    par = max ( parl, par + parc )
!
!  End of an iteration.
!
  end do
!
!  Termination.
!
  if ( iter == 0 ) then
    par = 0.0D+00
  end if

  return
end subroutine lmpar

END MODULE minimize
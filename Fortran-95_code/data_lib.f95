!    This module contains information on the grid employed by SP_Ace 
!    to move through the parameter space and some absorption lines that 
!    behave differently from the all the others contained in the line list.
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


MODULE data_lib
        USE num_type
        IMPLICIT NONE

        REAL(DP), DIMENSION(16), PARAMETER :: temp_grid=&
          &(/4000.,4200.,4400.,4600.,4800.,5000.,5200.,5400.,5600. &
          & ,5800.,6000.,6200.,6400.,6600.,6800.,7000./)
        REAL(DP), DIMENSION(10), PARAMETER :: logg_grid=&
          &(/1.0,1.4,1.8,2.2,2.6,3.0,3.4,3.8,4.2,4.6/)
        REAL(DP), DIMENSION(11), PARAMETER :: met_grid=&
          &(/-2.0,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,1e-4/)
        REAL(DP), DIMENSION(20), PARAMETER :: temp_gridL=&
          &(/3600.,3800.,4000.,4200.,4400.,4600.,4800.,5000.,5200.,5400.,5600. &
          & ,5800.,6000.,6200.,6400.,6600.,6800.,7000.,7200.,7400./)
        REAL(DP), DIMENSION(14), PARAMETER :: logg_gridL=&
          &(/0.2,0.6,1.0,1.4,1.8,2.2,2.6,3.0,3.4,3.8,4.2,4.6,5.0,5.4/)
        REAL(DP), DIMENSION(15), PARAMETER :: met_gridL=&
          &(/-2.4,-2.2,-2.0,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,1e-4,0.2,0.4/)

        !['1', 'x0', 'x1', 'x0^2', 'x0 x1', 'x1^2']
        !coeff computed with EWs expressed in A
!        REAL(DP), DIMENSION(6), PARAMETER :: gamG_coeff=&
!        &(/-0.07580377,0.29194086,0.75961454,0.01065768,-0.0416968,0.10474119/)

        !['1', 'x0', 'x1', 'x2', 'x0^2', 'x0 x1', 'x0 x2', 'x1^2', 'x1 x2', 'x2^2']
        !coeff computed with EWs expressed in A
!        REAL(DP), DIMENSION(10), PARAMETER :: gamL_coeff=&
!         &(/-0.07017293,0.6511216,0.01754621,-0.01243886,0.0078968,-0.07566144 &
!         &,0.03365214,-0.00207022,-0.00414377,0.0063777/)

        REAL(DP), DIMENSION(15), PARAMETER :: gamL_coeff=&
         &(/4.09178746e-02,9.71670522e-01,4.21578108e-04,7.89898167e-03 &
         &,-4.63726126e-05,1.32250692e-02,-7.40693971e-02,3.59427282e-02 &
         &,-9.09381428e-05,-2.56193818e-03,-4.54543975e-03,4.74718812e-06 &
         &,6.86968477e-03,-5.33724384e-06,4.92086560e-09/)

        REAL(DP), DIMENSION(10), PARAMETER :: gamL_coeff_H=&
        &(/-0.80416003,1.26632894,0.33366253,0.19390628,-0.03366942,-0.02230882 & !voigt
        &,0.10934166,0.02322909,0.0157008,-0.04058045/)
!         &(/-0.22626112,0.67650894,0.16448995,0.1579183,-0.0155744,-0.02943888 & !lorentzian only (test)
!         &,0.09563914,0.03604656,0.01045936,-0.04568341/)

        REAL(DP), DIMENSION(10), PARAMETER :: gamL_coeff_Na=&
         &(/-1.71216020e-01,6.88899191e-01,-8.13271006e-03,4.91015950e-02 &
         &,-2.37826204e-04,2.73305256e-03,7.33152275e-03,1.51602683e-03 &
         &,-6.04538703e-04,-1.04184217e-03/)


         REAL(DP), DIMENSION(2), PARAMETER :: w_rej_op=(/5227.17,5270.0/)
         REAL(DP), DIMENSION(2), PARAMETER :: r_rej_op=(/0.5,5.0/)

         REAL(DP), DIMENSION(10), PARAMETER :: w_rej_nlte=&
         &(/4861.3,6562.8,8413.322,8437.959,8467.258,8498.023,8542.091,8598.396,8662.141,8862.787/)!,8750.476/)
         REAL(DP), DIMENSION(10), PARAMETER :: r_rej_nlte=(/50.,50.,5.,5.,5.,10.,10.,10.,10.,10./)!,10./)

         REAL(DP), DIMENSION(15), PARAMETER :: w_rej_unknown=&
         &(/4849.4,5017.9,5053.6,5087.8,5215.5,5226.2,5371.5,5390.55,5654.5,5693.6,5899.5,6343.0,6449.1,6758.1,6765.45/)
         REAL(DP), DIMENSION(15), PARAMETER :: r_rej_unknown=&
          &(/0.2,0.2,0.4,0.2,0.2,0.2,0.4,0.2,0.2,0.2,0.2,1.0,0.5,0.2,0.2/)

         REAL(DP), DIMENSION(32), PARAMETER :: w_rej_bad=&
         &(/4920.52,4957.35,5175.,5727.027,5798.171,5816.372,5889.951,5895.924,5909.972,&
         &5918.536,5941.752,5944.660,&
         &5949.346,5958.27,5976.776,6029.900,6277.334,6278.234,6279.753,6280.617,6282.634,6292.824,6305.657,&
         &6309.920,6315.306,6355.028,6475.624,6483.944,6504.165,6516.077,6598.598,6743.122/)
         REAL(DP), DIMENSION(32), PARAMETER :: r_rej_bad=&
         &(/1.0,0.5,15.,0.2,0.2,0.2,10.,10.,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,&
         &0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2/)

        CHARACTER(LEN=3),DIMENSION(101), PARAMETER :: ELE_symb=&
        &(/'  H',' He',' Li',' Be','  B','  C',&
        &'  N','  O','  F',' Ne',' Na',' Mg',' Al',&
        &' Si','  P','  S',' Cl',' Ar','  K',' Ca',&
        &' Sc',' Ti','  V',' Cr',' Mn',' Fe',' Co',&
        &' Ni',' Cu',' Zn',' Ga',' Ge',' As',' Se',&
        &' Br',' Kr',' Rb',' Sr','  Y',' Zr',' Nb',&
        &' Mo',' Tc',' Ru',' Rh',' Pd',' Ag',' Cd',&
        &' In',' Sn',' Sb',' Te','  I',' Xe',' Cs',&
        &' Ba',' La',' Ce',' Pr',' Nd',' Pm',' Sm',&
        &' Eu',' Gd',' Tb',' Dy',' Ho',' Er',' Tm',&
        &' Yb',' Lu',' Hf',' Ta',' Wl',' Re',' Os',&
        &' Ir',' Pt',' Au',' Hg',' Tl',' Pb',' Bi',&
        &' Po',' At',' Rn',' Fr',' Ra',' Ac',' Th',&
        &' Pa','  U','  m','  a',' CH',' NH','MgH',&
        &'SiH',' CC',' CN','SiO'/)



END MODULE data_lib

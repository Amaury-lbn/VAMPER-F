

program test_fonctions

  use Principal, only : Vamper

  implicit none
  real, dimension(:), allocatable :: Temp, snw_totals
  real, dimension(:,:), allocatable :: Soil_temp
  integer :: unit_number
 ! allocate(Soil_temp(1:z_num,1:t_num))
 ! allocate(snw_totals(1:z_num))
 ! allocate(Temp(1:z_num)) 
  !write(*,*) "coucou"
  call Vamper(Temp, Soil_temp, snw_totals)

  open(newunit=unit_number,file="/home/users/alambin/VAMPER-F/Resultats/Temp_sol2.txt",status="old",action='write')

  !write(*,*) "coucou", dt

!  write(*,*) "coucou", Cp
  
!  write(*,*) "coucou", pori
  
!  write(*,*) "coucou", porf

!  do ll = 1,2

!     call AppHeatCapacity(z_num,T,T_freeze,n, organic_ind, Cp, porf, pori)
 
!     do kk=1,z_num-1
!
!     h_pori(kk) = (pori(kk) + pori(kk+1))/2
!     h_porf(kk) = (porf(kk) + porf(kk+1))/2
!     h_n(kk) = (n(kk) + n(kk+1))/2

!     call ThermalConductivity(kk,h_n(kk),h_pori(kk),h_porf(kk), organic_ind, T(kk), Kp(kk))

!     end do

!  end do

  write(unit_number,*) Temp

  write(*,*) Temp

  write(*,*) Soil_temp(:,5)

  !write(*,*) snw_totals

 ! write(*,*) "coucou", Soil_temp

  !write(*,*) "coucou", snw_totals

 end program test_fonctions



program test_fonctions

  use Principal, only : Vamper

  implicit none
  real :: T1,T2, t_temp
  real, dimension(:), allocatable :: Temp, snw_totals
  real, dimension(:,:), allocatable :: Soil_temp
  integer :: unit_number,kk, ll
  real, dimension(51) :: Temp_moy

  call Vamper(Temp, Soil_temp, snw_totals)

  !open(newunit=unit_number,file="/home/users/alambin/VAMPER-F/Resultats/Temp_sol2.txt",status="old",action='write')

  !write(unit_number,*) Temp

  !write(*,*) Temp

  !write(*,*) Soil_temp
     
  do ll=1,51
     
     do kk=1,12
     
        t_temp = t_temp + Soil_temp(ll,kk)
     
     end do

     Temp_moy(ll) = t_temp/12.0
      write(*,*) Temp_moy(ll)
     t_temp = 0.0

  end do

  write(*,*) Temp
  
 end program test_fonctions

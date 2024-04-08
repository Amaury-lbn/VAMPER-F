

program test_fonctions

  use Principal, only : Vamper_init,Lecture_forcing, Vamper_step
  use Fonction_init


  use Parametrisation, only : z_num,TotTime,Timestep,YearType,z_num,Depth,GridType,PorosityType,T_init,Bool_glacial 
  use Parametrisation, only : Bool_Organic,organic_depth,Gfx, T_freeze, EQ_Tr, EQ1_EQ2, Bool_delta,t_fin, alpha
  use Parametrisation, only : Bool_layer_temp,Forcage_Month_day,Bool_Swe_Snw,Bool_Model_Snow
  use Fonction_temp, only : AppHeatCapacity, ThermalConductivity
  use Fonction_init, only : Porosity_init, GeoHeatFlow, Glacial_index
  use Para_fonctions, only : t_disc, z_disc
  use Model_snow, only : snw_average_swe, snw_proc, snw_average_snw, snw_average_snw_tot
  use Fonction_implicit, only : Implicit_snow, Implicit

  implicit none
  real :: T1,T2, t_temp
  real, dimension(:), allocatable ::  snw_totals
  real, dimension(:,:), allocatable :: Soil_temp
  integer :: unit_number,kk, ll,organic_ind,spy, nb_lines,t_num,dim_temp,dim_swe,t_step
  real, dimension(51) :: Temp_moy
  real,dimension(:),allocatable :: time_gi, glacial_ind
  real :: Tb,dt
  real, dimension(:),allocatable:: T_air,Temp,Kp,n,dz,D,Cp,swe_f_t,pori,porf

  !call Glacial_index(time_gi,glacial_ind,nb_lines)

  !write(*,*) glacial_ind
  !write(*,*) time_gi

  !call Vamper(Temp, Soil_temp, snw_totals)

  !open(newunit=unit_number,file="/home/users/alambin/VAMPER-F/Resultats/Temp_sol2.txt",status="old",action='write')

  !write(unit_number,*) Temp

  !write(*,*) Temp

  !write(*,*) Soil_temp
     
  !do ll=1,51
     
     !do kk=1,12
     
!        t_temp = t_temp + Soil_temp(ll,kk)
        

 !    end do

  !   Temp_moy(ll) = t_temp/12.0
   !  t_temp = 0.0

  !end do
  !write(*,*) "Temp_moy = " , Temp_moy
  !do kk=1,12
     
  !write(*,*) "Nouveau mois",Soil_temp (:,kk)
     
 ! end do
  
 ! write(*,*) "Température finale =",Temp




  t_step = 100


  call t_disc(TotTime,Timestep,YearType,dt,spy,t_num)            
  call z_disc(z_num, Depth, GridType, dz, D)
  write(*,*)spy

  allocate(Kp(1:z_num-1))
  allocate(Cp(1:z_num))
  allocate(Temp(1:z_num))
  allocate(pori(1:z_num))
  allocate(porf(1:z_num))


  call Vamper_init(z_num,dz,D,Temp,time_gi,glacial_ind,nb_lines,Kp,Cp,n,organic_ind,Tb)

  call Lecture_forcing(z_num,T_air,swe_f_t,Temp,dim_temp,dim_swe)

  write(*,*) D

  write(*,*) Temp

  do kk = 1,99

     call Vamper_step(T_air,swe_f_t,Temp,Tb,Cp,Kp,n,organic_ind,glacial_ind,nb_lines,dim_temp,dim_swe,z_num,dz,dt,t_step, &
porf,pori,kk)

     write(*,*) Temp

  end do

  t_step = 20

  do kk = 1,4

     call Vamper_step(T_air,swe_f_t,Temp,Tb,Cp,Kp,n,organic_ind,glacial_ind,nb_lines,dim_temp,dim_swe,z_num,dz,dt,t_step, &
porf,pori,kk)

     write(*,*) Temp

  end do


  

 end program test_fonctions



program test_fonctions


!~   use Fonction_init, only: !dmr unused                                                                              [TBRMD]
!~   use Parametrisation, only : z_num,TotTime,Timestep,YearType,z_num,Depth,GridType,PorosityType,T_init,Bool_glacial [TBRMD] 
!~   use Parametrisation, only : Bool_Organic,organic_depth,Gfx, T_freeze, EQ_Tr, EQ1_EQ2, Bool_delta,t_fin, alpha     [TBRMD]
!~   use Parametrisation, only : Bool_layer_temp,Forcage_Month_day,Bool_Swe_Snw,Bool_Model_Snow                        [TBRMD]
!~   use Fonction_implicit, only : Implicit_snow, Implicit                                                             [TBRMD]

  use Parametrisation, only: z_num

  !dmr [2024-06-28] Functions used in the main
  use Principal, only : Vamper_init,Lecture_forcing, Vamper_step
  
  use Fonction_temp, only : AppHeatCapacity, ThermalConductivity
  use Fonction_init, only : Porosity_init, GeoHeatFlow, Glacial_index
  
  
  use Para_fonctions, only : t_disc, z_disc

  use Model_snow, only : snw_average_swe, snw_proc, snw_average_snw, snw_average_snw_tot

  implicit none
  
  integer :: kk, ll,organic_ind,spy, nb_lines,t_num,dim_temp,dim_swe,t_step,t_deb
  real,dimension(:),allocatable :: time_gi, glacial_ind
  real :: Tb,dt
  real, dimension(:),allocatable:: T_air,Temp,Kp,n,dz,D,Cp,swe_f_t,pori,porf,snow_dp_t,rho_snow_t,T_snw_t


  t_deb = 0
  kk=1
  ll=1

  !dmr [2024-06-28] [ADDING COMMENTS]
  !dmr Discretization routines

  !dmr Inputs to t_disc:
  !dmr
  !dmr intent(out)               t_num     the number of timesteps to perform from:     t_num = floor(model_secs/dt)
  !dmr intent(out)               spy       probably the number of steps per year from above
  !dmr intent(out)               dt        dt is the delta time step of the model, in seconds

  call t_disc(dt,spy,t_num)

!dmr [2024-06-28] Removed dependency to internal constants here. These intent(in) are parameters
  !dmr intent(in)                TotTime  defines the number of years (to run I presume)
  !dmr intent(in)                Timestep contains 1, 15 or 30 (from branching values)
  !dmr                                    30 seems to define monthly => spy = 12 and dt = real(Timestep)*60.0*60.0*24.0 in seconds
  !dmr                                    15 defines? [NOTA UNCLEAR] / spy = 24 that is two steps per months ???
  !dmr                                     1 defines a form of daily, with spy = 360 and dt as above. Why is there a Daily switch? [NOTA UNCLEAR] 
  !dmr intent(in)                YearType defines the number of days in years, 360 or 365
!~   call t_disc(TotTime,Timestep,YearType,dt,spy,t_num)
!dmr [2024-06-28] [TBRMD]


  !dmr Inputs to z_disc:
  !dmr 
  !dmr intent(out) (allocatable) dz thickness of the layer considered 
  !dmr intent(out) (allocatable) D = depth of the layer considered

  call z_disc(dz, D)
  
!dmr [2024-06-28] Removed dependency to internal constants here. These intent(in) are parameters

  !dmr intent(in)                Depth maximum depth, in meters, from the parametrisation file, now 1000 meters
  !dmr intent(in)                Gridtype if 2 : linspace else: logspace 
  !dmr intent(in)                z_num number of vertical layers, 51 or (now) 101 
!~   call z_disc(z_num, Depth, GridType, dz, D)
!dmr [2024-06-28] [TBRMD]

  write(*,*) "[MAIN] spy: ", spy

  allocate(Kp(1:z_num-1))
  allocate(Cp(1:z_num))
  allocate(Temp(1:z_num))
  allocate(pori(1:z_num))
  allocate(porf(1:z_num))

  !dmr [2024-06-28] [ADDING COMMENTS]
  !dmr No spatial dependence till here ...
  !dmr Inputs to Vamper_init:
  !dmr intent(in)                z_num == number of vertical slices
  !dmr intent(in)                dz = vertical stepping
  !dmr intent(in)                D is provided to porosity_init as depth_layer(z_num)
  !dmr intent(out) (allocatable) Temp = a form of temperature (which one?)
  !dmr intent(out) (allocatable) time_gi = years B.P. in the glacial index file (I think)
  !dmr intent(out) (allocatable) glacial_ind = array for glacial indexing, if not, then glacial_ind(1) = 0
  !dmr intent(out)               nb_lines = number of lines in the glacial index file, read in that file
  !dmr intent(out)               Kp = constant over the depth, current value is 2
  !dmr intent(out)               Cp 
  !dmr intent(out) (allocatable) n -> allocated in Porosity_init to z_num, contains porosity profile [NOTA: BAD_NAME] 
  !dmr intent(out)               organic_ind = depth of the organic layer? integer value in vertical index
  !dmr intent(out)               Tb = Temperature Bottom, lower boundary condition ... computed from GeoHeatFlow
  call Vamper_init(dz,D,Temp,time_gi,glacial_ind,nb_lines,Kp,Cp,n,organic_ind,Tb)

  !dmr [2024-06-28] [ADDING COMMENTS]
  !dmr This subroutine will be reading the forcing from external files. There is a suite of options as to how the forcing is done.
  !dmr Basically, it needs to fill in:
  !dmr
  !dmr 
  
  !dmr intent(inout) (z_num)     Temp         Temperature of the soil for each layer, read in external file unit_nb_3, one timestep (initial?)

  !dmr intent(out)               dim_swe      Length of the SWE forcing, calculated from the length of the input text file
  !dmr intent(out) (allocatable) swe_f_t      SWE forcing, allocated to (1:dim_swe) and values read in external text file (unit_nb_2)

  !dmr intent(out)               dim_temp     Length of the temp forcing, calculated from the length of the input text file      
  !dmr intent(out) (allocatable) T_air        temperature of the air forcing, allocated to (1:dim_temp) and read in external file (unit_nb_1)


  !dmr Those three are allocated to (1:dim_temp) and if BESSI, read from external files.

  !dmr intent(out) (allocatable) snw_dp_t     if BESSI, unit_nb_4 else set to zero
  !dmr intent(out) (allocatable) rho_snow_t   if BESSI, unit_nb_5 else set to zero
  !dmr intent(out) (allocatable) T_snw        if BESSI, unit_nb_6 else ... commented read, nothing done [NOTA: UNINITIALIZED]

  call Lecture_forcing(z_num,T_air,swe_f_t,snow_dp_t,rho_snow_t,T_snw_t,Temp,dim_temp,dim_swe)

  write(*,*) "[MAIN] D: ", D
  !write(*,*) T_air
  !write(*,*) dz

  write(*,*) "[MAIN] 1|Temp: ",Temp


  t_step = dim_temp

  !do kk = 1,800

  
  !dmr [2024-06-28] [ADDING COMMENTS]
  !dmr Main stepping of the VAMPER model
  !dmr
  
  !dmr intent(inout) (z_num)     dz
  !dmr intent(inout) (z_num)     n
  !dmr intent(inout) (z_num)     porf
  !dmr intent(inout) (z_num)     pori
  !dmr intent(inout) (z_num)     Kp
  !dmr intent(inout) (z_num)     Cp
  !dmr intent(inout) (z_num)     D
  !dmr intent(inout) (z_num)     Temp
  
  !dmr intent(in)    (dim_swe)   swe_f_t       / SWE forcing data
  !dmr intent(in)    (dim_temp)  T_air         / Air temperature forcing data
    
  !dmr intent(in)    (dim_temp)  snw_dp_t      / SNOW depth forcing
  !dmr intent(in)    (dim_temp)  rho_snow_t    / DENSITY of snow forcing
  !dmr intent(in)    (dim_temp)  T_snw_t       / TEMP of snow forcing
  !dmr intent(in)    (nb_lines)  glacial_ind   / glacial index modifier
  
  call Vamper_step(T_air,swe_f_t,Temp,Tb,Cp,Kp,n,organic_ind,glacial_ind,nb_lines,dim_temp,dim_swe,z_num,dz,dt,t_step, &
porf,pori,t_deb,rho_snow_t,snow_dp_t,T_snw_t,D)

  write(*,*) "[MAIN] 2|Temp: ",Temp

  write(*,*) "ok"
  

 end program test_fonctions

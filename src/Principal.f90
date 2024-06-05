module Principal


  use Parametrisation, only : z_num,TotTime,Timestep,YearType,Depth,GridType,PorosityType,T_init,Bool_glacial 
  use Parametrisation, only : Bool_Organic,organic_depth,Gfx, T_freeze, EQ_Tr, EQ1_EQ2, Bool_delta,t_fin, alpha
  use Parametrisation, only : Bool_layer_temp,Forcage_Month_day,Bool_Swe_Snw,Bool_Model_Snow,Bool_Bessi,s_l_max
  use Fonction_temp, only : AppHeatCapacity, ThermalConductivity, Permafrost_Depth
  use Fonction_init, only : Porosity_init, GeoHeatFlow, Glacial_index
  use Para_fonctions, only : t_disc, z_disc
  use Model_snow, only : snw_average_swe, snw_proc, snw_average_snw, snw_average_snw_tot
  use Fonction_implicit, only : Implicit_snow, Implicit
 
  Implicit none
  
  integer :: u_n_23,u_n_53,u_n_93,u_n_143,u_n_250,u_n_350,u_n_550,u_n_900
  integer :: layer_temp23,layer_temp53,layer_temp93,layer_temp143,layer_temp250,layer_temp350,layer_temp550,layer_temp900
  integer :: unit_nb_1,unit_nb_2,unit_nb_3,unit_nb_4,unit_nb_5,unit_nb_6

# include "constant.h"



contains

  subroutine Vamper(Temp, Soil_temp, snw_totals)
    
    integer :: t_num, kk, organic_ind, ll, nb_lines,spy, indice_tab,ind_days,snw_d,s_l_t
    real :: dt,Cp_snow,frac_snw,K_s,rho_snow,snw_tot,swe_f,Tb,snw_dp,swe_tot,T_soil, days,nb_days,snw_old,snw_f,G0_23,G0_18,itg
    real :: dz_snow


    integer, dimension(12) :: Day_per_month
    real, dimension(:),allocatable :: T_layer23,T_layer53,T_layer93,T_layer143,T_layer250,T_layer350,T_layer550,T_layer900
    real, dimension(z_num) :: T_old, Cp, porf, pori, Soil_temp_moy, delta_t_moy
    real, dimension(z_num-1) ::  h_n, h_pori, h_porf
    real, dimension(365*14) ::  snw_tot_t
    real, dimension(s_l_max+1) ::  Tsnw

    real, dimension(:),allocatable:: T_air,Timp,Kp,n,dz,D,Cp_t,swe_f_t,time_gi,glacial_ind,snw_f_t
    real, dimension(:,:),allocatable:: delta_T
    real, dimension(:),allocatable,intent(out) :: Temp
    real, dimension(:,:),allocatable,intent(out) :: Soil_temp
    real, dimension(:),allocatable,intent(out) :: snw_totals


    call t_disc(TotTime,Timestep,YearType,dt,spy,t_num)            ! Discretization of time -> this function find t_num(number of iteration) and dt(time step)
    call z_disc(z_num, Depth, GridType, dz, D)                     ! Discretization of space ->this function find D(depth of each layer) and dz(Distance between layer(x) and layer(x+1))
    
    allocate(Kp(1:z_num-1))
    allocate(Temp(1:z_num))

    ! ---- Allocation of different space in function of the case -----!

    if (Bool_delta==1)then

       allocate(delta_T(1:z_num,1:t_num-1))

    end if

    if (EQ_Tr == 1)then       
       
       allocate(Soil_temp(1:z_num,1:t_num))
       allocate(snw_totals(1:t_num))
       allocate(T_air(1:t_num))

    end if

    if (Forcage_Month_day == 1)then
       
       allocate(snw_f_t(1:365))
       allocate(T_air(1:365))
       allocate(Soil_temp(1:z_num,1:365))
       allocate(T_layer23(1:365*20))
       allocate(T_layer53(1:365*20))
       allocate(T_layer93(1:365*20))
       allocate(T_layer143(1:365*20))
       allocate(T_layer250(1:365*20))
       allocate(T_layer350(1:365*20))
       allocate(T_layer550(1:365*20))
       allocate(T_layer900(1:365*20))

    elseif(Forcage_Month_day == 0)then

       allocate(T_air(1:12))
       allocate(snw_f_t(1:12))
       allocate(Soil_temp(1:z_num,1:12))
       allocate(T_layer23(1:12*20))
       allocate(T_layer53(1:12*20))
       allocate(T_layer93(1:12*20))
       allocate(T_layer143(1:12*20))
       allocate(T_layer250(1:12*20))
       allocate(T_layer350(1:12*20))
       allocate(T_layer550(1:12*20))
       allocate(T_layer900(1:12*20))
    
    end if
    
    write(*,*) D
    !write(*,*) dz
    Day_per_month(1) = 31
    Day_per_month(2) = 28
    Day_per_month(3) = 31
    Day_per_month(4) = 30
    Day_per_month(5) = 31
    Day_per_month(6) = 30
    Day_per_month(7) = 31
    Day_per_month(8) = 31
    Day_per_month(9) = 30
    Day_per_month(10) = 31
    Day_per_month(11) = 30
    Day_per_month(12) = 31

    days = Day_per_month(1)
    ind_days = 1
    nb_days = 0
    Tsnw = -4.0
    frac_snw = 1.0
    rho_snow = 0.0
    snw_tot = 0.0
    swe_tot = 0.0

    !--------- OPENING OF FILE WHICH ARE READ OR WRITTEN ---------!

    if (Bool_layer_temp==1)then

       layer_temp23 = 14
       layer_temp53 = 18
       layer_temp93 = 20
       layer_temp143 = 22
       layer_temp250 = 24
       layer_temp350 = 26
       layer_temp550 = 28
       layer_temp900 = 30

       if (z_num == 101) then

          layer_temp23 = 28
          layer_temp53 = 35
          layer_temp93 = 40
          layer_temp143 = 44
          layer_temp250 = 49
          layer_temp350 = 51
          layer_temp550 = 55
          layer_temp900 = 60

       end if
       open(newunit=u_n_23,file="/home/users/alambin/VAMPER-F/Resultats/Temp_layer_23.txt",status="replace",action='write')
       open(newunit=u_n_53,file="/home/users/alambin/VAMPER-F/Resultats/Temp_layer_53.txt",status="replace",action='write')
       open(newunit=u_n_93,file="/home/users/alambin/VAMPER-F/Resultats/Temp_layer_93.txt",status="replace",action='write')
       open(newunit=u_n_143,file="/home/users/alambin/VAMPER-F/Resultats/Temp_layer_143.txt",status="replace",action='write')
       open(newunit=u_n_250,file="/home/users/alambin/VAMPER-F/Resultats/Temp_layer_250.txt",status="replace",action='write')
       open(newunit=u_n_350,file="/home/users/alambin/VAMPER-F/Resultats/Temp_layer_350.txt",status="replace",action='write')
       open(newunit=u_n_550,file="/home/users/alambin/VAMPER-F/Resultats/Temp_layer_550.txt",status="replace",action='write')
       open(newunit=u_n_900,file="/home/users/alambin/VAMPER-F/Resultats/Temp_layer_900.txt",status="replace",action='write') 
       open(newunit=snw_d,file="/home/users/alambin/VAMPER-F/Resultats/Snw_depth.txt",status="replace",action='write') 
    end if

    if (EQ_Tr == 0)then

      open(newunit=unit_nb_1,file="/home/users/alambin/VAMPER-F/Donnee/Temp_under_ice.txt",status="old",action='read')
      

      if(EQ1_EQ2==2 .and. z_num ==51)then
         open(newunit=unit_nb_3,file="/home/users/alambin/VAMPER-F/Donnee/Temp_soil_0.txt",status="old",action='read') 
      else
         open(newunit=unit_nb_3,file="/home/users/alambin/VAMPER-F/Donnee/Temp_soil_0_z101.txt",status="old",action='read') 
      end if
      
      if (Forcage_Month_day == 0)then
         open(newunit=unit_nb_2,file="/home/users/alambin/VAMPER-F/Donnee/Snow_EQ.txt",status="old",action='read')
         open(newunit=unit_nb_4,file="/home/users/alambin/VAMPER-F/Donnee/Temp_Bayevla.txt",status="old",action='read')
      else
         open(newunit=unit_nb_2,file="/home/users/alambin/VAMPER-F/Donnee/Snow_fresh_day.txt",status="old",action='read')
         open(newunit=unit_nb_4,file="/home/users/alambin/VAMPER-F/Donnee/Temperature_moyenne_jour.txt",status="old",action='read')
      end if

      if (Bool_Model_Snow==0)then
         open(newunit=unit_nb_5,file="/home/users/alambin/VAMPER-F/Donnee/Snow_tot.txt",status="old",action='read')
      end if

    elseif(EQ_Tr==1)then

       open(newunit=unit_nb_3,file="/home/users/alambin/VAMPER-F/Donnee/Temp_soil_1998.txt",status="old",action='read')
       open(newunit=unit_nb_1,file="/home/users/alambin/VAMPER-F/Donnee/Temp_18y_day.txt",status="old",action='read')
       open(newunit=unit_nb_2,file="/home/users/alambin/VAMPER-F/Donnee/Swe_18y_day.txt",status="old",action='read')

    end if
    

!----------INITIALISATION----------!   

    call Porosity_init(z_num, PorosityType, D, Bool_Organic, organic_depth, n, organic_ind )  !CALCULATION OF POROSITY
    write(*,*) n

    if (Bool_glacial == 1)then

       call Glacial_index(time_gi,glacial_ind,nb_lines)    ! GLACIAL INDEX FOR THE FORCING OF TEMPERATURE

    end if

    do kk=1,z_num

       Cp(kk)=1

    end do

   ! RECOVERY OF FORCING SNOW AND FORCING TEMPERATURE !

    if (EQ_Tr == 0)then                              !For an equilibrum run
       
       if (Forcage_Month_day == 0)then                 ! FOR A MONTHLY FORCING

          do ll =1,12

             read(unit_nb_2,*) snw_f_t(ll)            ! RECOVERY OF SWE

             if (Bool_glacial == 1)then               ! RECOVERY OF AIR TEMPERATURE
                read(unit_nb_4,*) T_air(ll)
             elseif(Bool_glacial == 0)then
                read(unit_nb_1,*) T_air(ll)
             end if

          end do

       elseif(Forcage_Month_day == 1)then             ! For a daily forcing 

          do ll =1,365
             if (Bool_Model_Snow ==1)then
                read(unit_nb_2,*) snw_f_t(ll)           ! Snow data
             end if
             read(unit_nb_4,*) T_air(ll)             ! Air temperature data
          end do

       end if
       
       if (EQ1_EQ2 ==1 ) then                                       ! If there is no initial condition for the ground temperature

          do kk=1,z_num-1

             Kp(kk)=2
          
          end do
          
          call GeoHeatFlow(Gfx, Kp, dz, T_init, z_num, Temp)            

       elseif(EQ1_EQ2==2)then                                        ! If there is initial condition for the ground temperature

          do ll =1,z_num
             read(unit_nb_3,*) Temp(ll)
          end do

       end if

    else                                             !For a transient run

       do ll =1,z_num

          read(unit_nb_3,*) Temp(ll)

       end do
       
       read(unit_nb_1,*) T_air(1)
       if(Bool_Model_Snow==1)then
          !read(unit_nb_2,*) swe_f
       end if
       do kk = 1,365

          read(unit_nb_2,*) snw_f_t(kk)

       end do

    end if

    if (Bool_Model_Snow == 0)then

       do kk=1,365*14
          read(unit_nb_5,*) snw_tot_t(kk) 
       end do

    end if

    Tb = Temp(z_num)                         ! Lower boundary condition 
    write(*,*) Temp

    do ll = 1,2
       
       call AppHeatCapacity(z_num,Temp,T_freeze,n, organic_ind, Cp, porf, pori)         !Calculation of heat capacity of soil
       
       do kk=1,z_num-1
          
          h_pori(kk) = (pori(kk) + pori(kk+1))/2
          h_porf(kk) = (porf(kk) + porf(kk+1))/2
          h_n(kk) = (n(kk) + n(kk+1))/2
          
          call ThermalConductivity(kk,h_n(kk),h_pori(kk),h_porf(kk), organic_ind, Temp(kk), Kp(kk))    !Calculation of thermal condutivity of soil
          
       end do
       
    end do

    if (Bool_delta==1)then
       delta_t(1:z_num,1) = Temp(1:z_num)
    end if
             

                   !  -------- MAIN LOOP --------- !
    do ll=1,t_num

       snw_old = snw_tot
       nb_days = nb_days + 1

       if (nb_days > days) then
          nb_days = 1
          ind_days = ind_days + 1

          if (ind_days == 13) then
             ind_days = 1
          end if

          days = Day_per_month(ind_days)

       end if

       ! ---------- Reading of forcing tamperature and snow/swe forcing --------------- !

       if (EQ_Tr == 0)then
          
          if( Forcage_Month_day == 0 .and. Bool_glacial == 0)then
             if (EQ1_EQ2 == 1)then
                T_soil = (T_air(mod(ll,12)+1))
             elseif (EQ1_EQ2==2)then
                T_soil = T_air(mod(ll,12)+1)
             end if
          end if

          if (Bool_Swe_Snw == 0)then
             if( Forcage_Month_day == 1)then
                swe_f = snw_f_t(mod(ll,365)+1)
             else
                swe_f =  snw_f_t(mod(ll,12)+1)
             end if
          elseif (Bool_Swe_Snw == 1)then
             if( Forcage_Month_day == 1)then
                snw_f = snw_f_t(mod(ll,365)+1)
             else
                snw_f = snw_f_t(mod(ll,12)+1)
             end if
          end if
          
          if (Bool_Model_Snow == 0)then
             snw_tot = snw_tot_t(mod(ll,365*14)+1) 
          end if
   
       elseif(EQ_Tr==1)then

          read(unit_nb_1,*) T_air(ll)
          !read(unit_nb_2,*) swe_f
          swe_f = snw_f_t(mod(ll,365)+1)

       end if

       !write(*,*) T_air(ll),swe_f,snw_tot

       if (Bool_glacial==1)then

          if( Forcage_Month_day == 0)then

             indice_tab = nb_lines-floor((-(ll/12.0)+t_fin+TotTime)/100.0)
             T_soil=alpha*(glacial_ind(indice_tab-1)+mod((ll/12.0),100.0)*(glacial_ind(indice_tab)-glacial_ind(indice_tab-1))/100.0)
             T_soil = (T_air(mod(ll,12)+1)+T_soil)
             write(*,*) T_soil

          elseif( Forcage_Month_day == 1)then

             indice_tab = nb_lines-floor((-(ll/(365.0))+t_fin+TotTime)/100.0)
             T_soil=alpha*(glacial_ind(indice_tab-1)+mod((ll/12.0),100.0)*(glacial_ind(indice_tab)-glacial_ind(indice_tab-1))/100.0)
             T_soil = T_air(mod(ll,365)+1)+T_soil

          end if

       end if

       T_old(1:z_num) = Temp(1:z_num)
       

       ! ----- Parametrisation of snow layers based on fresh/old snow depth --------!
       
       if (snw_tot > 0.000001 .or. swe_f > 0.000001 .or. snw_f>0.00001) then

          if (Bool_Swe_Snw == 0)then
             call snw_average_swe(swe_f, swe_tot, snw_tot, rho_snow)
          elseif(Bool_Swe_Snw == 1)then
             call snw_average_snw(snw_f, swe_tot, snw_tot, rho_snow,swe_f)
          end if
          
          if(Bool_Model_snow==0)then
             call  snw_average_snw_tot(snw_tot,snw_old, rho_snow,swe_tot,dt)
          end if
          
          if (abs(swe_tot - swe_f )< 0.000001 .and. Bool_Model_snow==1) then

             if (EQ_Tr == 0)then
                Tsnw = T_soil
             else
                Tsnw = T_air(ll) 
             end if
             K_s = 0.07
             frac_snw = 1

          end if

          if (abs(snw_old)< 0.00001 .and. Bool_Model_snow==0) then

             if (EQ_Tr == 0)then
                Tsnw = T_soil
             else
                Tsnw = T_air(ll) 
             end if
             K_s = 0.07
             frac_snw = 1

          end if
       
       end if

       !-------------- Numerical difference routine when there is snow or not --------!

       if (snw_tot > 0.00001) then

          if (EQ_Tr == 0)then

             call Implicit_snow(snw_tot,rho_snow,Tsnw,T_old,T_soil,dt,dz,n,organic_ind,Temp,Cp_t,Kp,Cp_snow,s_l_t)
          
          else

             call Implicit_snow(snw_tot,rho_snow,Tsnw,T_old,T_air(ll),dt,dz,n,organic_ind,Temp,Cp_t,Kp,Cp_snow,s_l_t)

          end if

          call snw_proc(Tsnw(1),Temp(1), snw_tot, swe_tot, frac_snw, Cp_snow, rho_snow, dt)
          !if(Tsnw > 0)then
           !  frac_snw = 0.0
            ! snw_tot = 0.0
            ! Tsnw = 0.0
          !end if
          
       else

          swe_tot = 0.0
          snw_tot = 0.0
          
          if (EQ_Tr == 0)then

             call Implicit(T_old,T_soil,Tb,dt,dz,n,organic_ind,Temp,Cp_t,Kp) 
          
          else

             call Implicit(T_old,T_air(ll),Tb,dt,dz,n,organic_ind,Temp,Cp_t,Kp)

          end if
    
       end if

       ! -------- Filing table to show soil temperature ----------!
          
       do kk=1,z_num
          
          if (EQ_Tr==0) then
             
             if( Forcage_Month_day == 1)then
                Soil_temp(kk,mod(ll,365)+1) = Temp(kk)
             elseif( Forcage_Month_day == 0)then
                Soil_temp(kk,mod(ll,12)+1) = Temp(kk)
             end if

             if (Bool_delta==1)then
                delta_t(kk,ll) = Temp(kk) - T_old(kk)
             end if

          else
             
             Soil_temp(kk,ll) = Temp(kk)
             
             if (Bool_delta==1)then
                delta_t(kk,ll) = Temp(kk) - T_old(kk)
             end if

          end if

       end do
                 
       if (mod(ll,120000) == 0) then

          write(*,*) indice_tab, nb_lines

          do kk=1,z_num
             
             if (Bool_delta==1)then

                delta_t_moy(kk)=sum(delta_t(kk,(ll-11):ll))
                
             end if
                
          end do
          write(*,*) ll , Temp
          
       end if
       
       !write(*,*) snw_tot, Cp_Snow, T_soil, Temp(1), rho_snow
       
       if (Forcage_month_day == 1 .and. Bool_layer_temp==1) then

          if (mod(ll,1000*365) == 0) then

             do kk=1,z_num
                
                Soil_temp_moy(kk) = sum(Soil_temp(kk,1:365))/365.0
                
             end do
             write(*,*) ll
             write(*,*) "Soil_temp_moy =" , Soil_temp_moy
             !write(*,*) "delta_t_moy =" , delta_t_moy 

          end if
          
          if (ll>t_num-7300)then
             T_layer23(ll+7300-t_num) = Temp(layer_temp23)
             T_layer53(ll+7300-t_num) = Temp(layer_temp53)
             T_layer93(ll+7300-t_num) = Temp(layer_temp93)
             T_layer143(ll+7300-t_num) = Temp(layer_temp143)
             T_layer250(ll+7300-t_num) = (Temp(layer_temp250)+Temp(layer_temp250 - 1))/2.0
             T_layer350(ll+7300-t_num) = Temp(layer_temp350)
             T_layer550(ll+7300-t_num) = Temp(layer_temp550)
             T_layer900(ll+7300-t_num) = Temp(layer_temp900)
             write(snw_d,*) snw_tot
          end if
          
       elseif (Forcage_month_day == 0 .and. Bool_layer_temp==1)then

          if (ll>t_num-240)then
             T_layer23(240+ll-t_num) = Temp(layer_temp23)
             T_layer53(240+ll-t_num) = Temp(layer_temp53)
             T_layer93(240+ll-t_num) = Temp(layer_temp93)
             T_layer143(240+ll-t_num) = Temp(layer_temp143)
             T_layer250(240+ll-t_num) = (Temp(layer_temp250)+Temp(layer_temp250+1))/2.0
             T_layer350(240+ll-t_num) = Temp(layer_temp350)
             T_layer550(240+ll-t_num) = Temp(layer_temp550)
             T_layer900(240+ll-t_num) = Temp(layer_temp900)
             write(snw_d,*) snw_tot
          end if

       end if

       itg=0

       do kk = 1,26
          itg = itg + Cp(kk)*(Temp(26)-T_old(26))*dz(kk)/dt
       end do
       G0_18 = Kp(26)*(Temp(26)-Temp(1))/D(26) + itg

       G0_23 = Kp(layer_temp23)*(Temp(layer_temp23)-Temp(1))/D(layer_temp23)
       
       write(*,*) G0_23, G0_18,Cp_snow, T_soil,snw_tot

       
       
    end do
    write(*,*) D

    
    if (Bool_layer_temp==1)then
       write(u_n_23,*) T_layer23
       write(u_n_53,*) T_layer53
       write(u_n_93,*) T_layer93
       write(u_n_143,*) T_layer143
       write(u_n_250,*) T_layer250
       write(u_n_350,*) T_layer350
       write(u_n_550,*) T_layer550
       write(u_n_900,*) T_layer900
       
    end if

  end subroutine Vamper













  subroutine Vamper_init(z_num,dz,D,Temp,time_gi,glacial_ind,nb_lines,Kp,Cp,n,organic_ind,Tb)

    integer, intent(in) :: z_num
    real, dimension(z_num),intent(in) :: dz,D
    real, dimension(z_num),intent(out):: Cp
    real, dimension(z_num-1),intent(out):: Kp
    real, dimension(:),allocatable,intent(out) :: time_gi,glacial_ind,n,Temp
    real, intent(out) :: Tb
    integer, intent(out) :: nb_lines
    integer, intent(out) :: organic_ind

    real, dimension(z_num-1) ::  h_n, h_pori, h_porf
    real, dimension(z_num) :: porf,pori
    integer :: kk,ll
    
    


    call Porosity_init(z_num, PorosityType, D, Bool_Organic, organic_depth, n, organic_ind )  !CALCULATION OF POROSITY
    !write(*,*) n

    !do kk=1,z_num

     !  if(D(kk)>1.4)then

      !    n(kk) = n(kk) + 0.5

       !end if

   ! end do

    if (Bool_glacial == 1)then

       call Glacial_index(time_gi,glacial_ind,nb_lines)    ! GLACIAL INDEX FOR THE FORCING OF TEMPERATURE

    else

       allocate(glacial_ind(1:1))
       glacial_ind(1) = 0

    end if
                                
    
    do kk=1,z_num-1
       
       Kp(kk)=2
       
    end do
          
    call GeoHeatFlow(Gfx, Kp, dz, T_init, z_num, Temp)            

    Tb = Temp(z_num)                         ! Lower boundary condition 
    
    !write(*,*) Temp

    do ll = 1,2
       
       call AppHeatCapacity(z_num,Temp,T_freeze,n, organic_ind, Cp, porf, pori)         !Calculation of heat capacity of soil

       do kk=1,z_num-1
          
          h_pori(kk) = (pori(kk) + pori(kk+1))/2
          h_porf(kk) = (porf(kk) + porf(kk+1))/2
          h_n(kk) = (n(kk) + n(kk+1))/2
          call ThermalConductivity(kk,h_n(kk),h_pori(kk),h_porf(kk), organic_ind, Temp(kk), Kp(kk))    !Calculation of thermal condutivity of soil
          
       end do
       
    end do


  end subroutine Vamper_init









  subroutine Lecture_forcing(z_num,T_air,swe_f_t,snw_dp_t,rho_snow_t,T_snw,Temp,dim_temp,dim_swe)
    
    
    integer, intent(in) :: z_num
    real, dimension(z_num),intent(inout) :: Temp
    integer, intent(out) :: dim_temp,dim_swe
    real, dimension(:), allocatable, intent(out):: T_air, swe_f_t, snw_dp_t,rho_snow_t,T_snw
    integer :: kk,ii,ll
    real :: ligne

    ll = 0
    dim_swe = 0
    dim_temp = 0

    if (EQ_Tr == 0)then


      if(EQ1_EQ2==2 .and. z_num ==51)then
         open(newunit=unit_nb_3,file="/home/users/alambin/VAMPER-F/Donnee/Temp_soil_0.txt",status="old",action='read') 
      else
         open(newunit=unit_nb_3,file="/home/users/alambin/VAMPER-F/Donnee/Temp_soil_0_z101.txt",status="old",action='read') 
      end if

      
# if Daily == 0
         open(newunit=unit_nb_2,file="/home/users/alambin/VAMPER-F/Donnee/Snow_EQ.txt",status="old",action='read')
         open(newunit=unit_nb_1,file="/home/users/alambin/VAMPER-F/Donnee/Temp_Bayevla.txt",status="old",action='read')
# else
         open(newunit=unit_nb_2,file="/home/users/alambin/VAMPER-F/Donnee/Snow_fresh_day.txt",status="old",action='read')
         open(newunit=unit_nb_1,file="/home/users/alambin/VAMPER-F/Donnee/Temperature_moyenne_jour.txt",status="old",action='read')
# endif


    elseif(EQ_Tr==1)then

       open(newunit=unit_nb_3,file="/home/users/alambin/VAMPER-F/Init_Svalbard/Ts_Sv_2.0_0.5_0.2.txt",&
status="old",action='read')
       open(newunit=unit_nb_1,file="/home/users/alambin/VAMPER-F/Donnee/T_snw_d.txt",status="old",action='read')
       open(newunit=unit_nb_2,file="/home/users/alambin/VAMPER-F/Donnee/Snow_tot.txt",status="old",action='read')

    end if

    if (Bool_Bessi==1)then

       open(newunit=unit_nb_4,file="/home/users/alambin/VAMPER-F/Donnee/Snow_dp_1998.txt",status="old",action='read')
       open(newunit=unit_nb_5,file="/home/users/alambin/VAMPER-F/Donnee/Rho_snow_1998.txt",status="old",action='read')

    end if

    !open(newunit=unit_nb_6,file="/home/users/alambin/VAMPER-F/Donnee/T_snw.txt",status="old",action='read')
    
    if(EQ1_EQ2==2)then
       
       do kk=1,z_num 
          read(unit_nb_3,*) Temp(kk)
       end do

       do kk=1,z_num 
          Temp(kk) = Temp(kk) 
       end do  

       close(unit_nb_3)

    end if

    
    do

       read(unit_nb_2,*,iostat=ii) ligne

       if(ii/=0)exit

       dim_swe = dim_swe + 1

    end do
    
    allocate(swe_f_t(1:dim_swe))
    rewind(unit_nb_2)
    
    do kk = 1,dim_swe
       
       read(unit_nb_2,*) swe_f_t(kk)

    end do


    do

       read(unit_nb_1,*,iostat=ii) ligne

       if(ii/=0)exit

       dim_temp = dim_temp + 1

    end do

    allocate(T_air(1:dim_temp))    
    rewind(unit_nb_1)

    do kk = 1,dim_temp
       
       read(unit_nb_1,*) T_air(kk)

    end do


    if (Bool_Bessi==1)then
       allocate(rho_snow_t(1:dim_temp)) 
       allocate(snw_dp_t(1:dim_temp))
       allocate(T_snw(1:dim_temp))
       
       do kk = 1,dim_temp
       
          read(unit_nb_5,*) rho_snow_t(kk)
          read(unit_nb_4,*) snw_dp_t(kk)
          read(unit_nb_6,*) T_snw(kk)

       end do

    else

       allocate(rho_snow_t(1:dim_temp)) 
       allocate(snw_dp_t(1:dim_temp)) 
       allocate(T_snw(1:dim_temp))
       
       do kk = 1,dim_temp
       
          rho_snow_t(kk) = 0
          snw_dp_t(kk) = 0
          !read(unit_nb_6,*) T_snw(kk)

       end do
 
    end if
    

    close(unit_nb_1)    
    close(unit_nb_2)
    

  end subroutine Lecture_forcing



  subroutine Vamper_step(T_air,swe_f_t,Temp,Tb,Cp,Kp,n,organic_ind,glacial_ind,nb_lines,dim_temp,dim_swe,z_num,dz,dt,t_step, &
porf,pori,t_deb,rho_snow_t,snw_dp_t,T_snw_t,D)

    integer, intent(inout) ::  organic_ind, nb_lines, dim_swe, dim_temp, t_step,t_deb
    real, intent(in) :: dt, Tb
    integer, intent(in) :: z_num

    real,dimension(z_num),intent(inout) :: dz,n,porf,pori
    real,dimension(z_num),intent(inout) :: Kp,Cp,D
    real,dimension(dim_swe),intent(in) :: swe_f_t
    real,dimension(dim_temp),intent(in) :: T_air,snw_dp_t,rho_snow_t,T_snw_t
    real,dimension(nb_lines),intent(in) :: glacial_ind
    real, dimension(z_num),intent(inout) :: Temp

    integer :: kk, ll, indice_tab, snw_d,t_num,s_l_t,ind_snw
    real :: T_soil, T_glacial, swe_tot,snw_tot,rho_snow,swe_f,frac_snw,k_s,Cp_snow,spy,snw_old,dz_snow,Per_depth
    real, dimension(s_l_max) :: Tsnw
    real, dimension(z_num-1) ::  h_n, h_pori, h_porf
    real, dimension(z_num) :: T_old
    real, dimension(:),allocatable :: T_layer23,T_layer53,T_layer93,T_layer143,T_layer250,T_layer350,T_layer550,T_layer900

    
     if (Bool_layer_temp==1)then

       open(newunit=u_n_23,file="/home/users/alambin/VAMPER-F/Resultats/Tl_23_2.0_0.1_0.0.txt",&
status="replace",action='write')
       open(newunit=u_n_53,file="/home/users/alambin/VAMPER-F/Resultats/Tl_53_2.0_0.1_0.0.txt",&
status="replace",action='write')
       open(newunit=u_n_93,file="/home/users/alambin/VAMPER-F/Resultats/Tl_93_2.0_0.1_0.0.txt",&
status="replace",action='write')
       open(newunit=u_n_143,file="/home/users/alambin/VAMPER-F/Resultats/Tl_143_2.0_0.1_0.0.txt",&
status="replace",action='write')
       open(newunit=u_n_250,file="/home/users/alambin/VAMPER-F/Resultats/Tl_250_2.0_0.1_0.0.txt",&
status="replace",action='write')
       open(newunit=u_n_350,file="/home/users/alambin/VAMPER-F/Resultats/Tl_350_2.0_0.1_0.0.txt",&
status="replace",action='write')
       open(newunit=u_n_550,file="/home/users/alambin/VAMPER-F/Resultats/Tl_550_2.0_0.1_0.0.txt",&
status="replace",action='write')
       open(newunit=u_n_900,file="/home/users/alambin/VAMPER-F/Resultats/Tl_900_2.0_0.1_0.0.txt",&
status="replace",action='write') 
       open(newunit=snw_d,file="/home/users/alambin/VAMPER-F/Resultats/Snw_depth.txt",status="replace",action='write') 

       layer_temp23 = 28
       layer_temp53 = 35
       layer_temp93 = 40
       layer_temp143 = 44
       layer_temp250 = 49
       layer_temp350 = 51
       layer_temp550 = 55
       layer_temp900 = 60
       
       !open(newunit=u_n_23,file="/home/users/alambin/VAMPER-F/Resultats/Tl_75_2.0_1.0_1.0.txt",status="replace",action='write')
       !open(newunit=u_n_53,file="/home/users/alambin/VAMPER-F/Resultats/Tl_175_2.0_1.0_1.0.txt",status="replace",action='write')
       !open(newunit=u_n_93,file="/home/users/alambin/VAMPER-F/Resultats/Tl_375_2.0_1.0_1.0.txt",status="replace",action='write')
       !open(newunit=u_n_143,file="/home/users/alambin/VAMPER-F/Resultats/Tl_675_2.0_1.0_1.0.txt",status="replace",action='write')
       !open(newunit=u_n_250,file="/home/users/alambin/VAMPER-F/Resultats/Tl_875_2.0_1.0_1.0.txt",status="replace",action='write')
       !open(newunit=u_n_350,file="/home/users/alambin/VAMPER-F/Resultats/Tl_1275_2.0_1.0_1.0.txt",status="replace",action='write')
       !open(newunit=u_n_550,file="/home/users/alambin/VAMPER-F/Resultats/Tl_1675_2.0_1.0_1.0.txt",status="replace",action='write')
       !open(newunit=u_n_900,file="/home/users/alambin/VAMPER-F/Resultats/Tl_2475_2.0_1.0_1.0.txt",status="replace",action='write') 
       !open(newunit=snw_d,file="/home/users/alambin/VAMPER-F/Resultats/Snw_depth.txt",status="replace",action='write')

       !layer_temp23 = 57
       !layer_temp53 = 68
       !layer_temp93 = 78
       !layer_temp143 = 85
       !layer_temp250 = 89
       !layer_temp350 = 94
       !layer_temp550 = 97
       !layer_temp900 = 101
       
       allocate(T_layer23(1:365*20))
       allocate(T_layer53(1:365*20))
       allocate(T_layer93(1:365*20))
       allocate(T_layer143(1:365*20))
       allocate(T_layer250(1:365*20))
       allocate(T_layer350(1:365*20))
       allocate(T_layer550(1:365*20))
       allocate(T_layer900(1:365*20))

    end if


    swe_tot = 0
    snw_tot = 0
    snw_old = 0
    ind_snw = 0
    s_l_t = 1
    write(*,*) organic_ind
    
    do kk =1,s_l_max

       Tsnw(kk) = -4

    end do

    dz_snow= 0.01
    D(1) = dz_snow*s_l_max
    do kk = 1,s_l_max-1
       D(kk+1) = D(kk) - dz_snow
    end do
    
    do kk=1,z_num

       Cp(kk) = 1

    end do

    do ll = 1,2
       
       call AppHeatCapacity(z_num,Temp,T_freeze,n, organic_ind, Cp, porf, pori)         !Calculation of heat capacity of soil
       
       do kk=1,z_num-1
          
          h_pori(kk) = (pori(kk) + pori(kk+1))/2
          h_porf(kk) = (porf(kk) + porf(kk+1))/2
          h_n(kk) = (n(kk) + n(kk+1))/2
          
          call ThermalConductivity(kk,h_n(kk),h_pori(kk),h_porf(kk), organic_ind, Temp(kk), Kp(kk))    !Calculation of thermal condutivity of soil
          
       end do
       
    end do

# if Daily == 0

       spy = 12
       t_num = t_step*12

#else 

       spy = 365
       t_num = t_step 
       !t_num = 2

#endif

    !write(*,*) n

    do ll=1,t_num

       !write(*,*) "ok"

       snw_old = snw_tot

       T_soil = T_air(mod(ll,dim_temp)+1)
       swe_f  = swe_f_t(mod(ll,dim_swe)+1)
       snw_tot = swe_f_t(mod(ll,dim_swe)+1)
       !snw_tot = 0

       if (Bool_Bessi==1)then
          rho_snow = rho_snow_t(mod(ll,dim_temp)+1)
          snw_tot = snw_dp_t(mod(ll,dim_temp)+1)
       end if
       !write(*,*) snw_tot,T_soil,rho_snow
       if (Bool_glacial==1)then

          indice_tab = nb_lines-floor((-(ll/spy)+t_deb)/100.0)
          T_glacial=alpha*(glacial_ind(indice_tab-1)+(100-mod((-(ll/spy)+t_deb),100.0))*(glacial_ind(indice_tab)- &
glacial_ind(indice_tab-1))/100.0)
          T_soil = (T_glacial+T_soil)
          !write(*,*) indice_tab,T_soil,glacial_ind(indice_tab-1),glacial_ind(indice_tab),T_glacial,mod((-(ll/spy)+t_deb),100.0)

       end if

       if (snw_tot > 0.000001 .or. swe_f > 0.000001  ) then

          if (Bool_Bessi==0) then
             call snw_average_swe(swe_f, swe_tot, snw_tot, rho_snow)
          end if
          
          if (abs(swe_tot - swe_f )< 0.000001 .and. Bool_Bessi==0) then

             Tsnw = T_soil
             K_s = 0.07
             frac_snw = 1

          end if
          
          if (abs(snw_old )< 0.000001) then

             Tsnw(1) = T_soil
             !Tsnw(2) = (
             K_s = 0.07
             frac_snw = 1

          end if
       end if
       
       if (T_soil>0.0) then

          snw_tot = 0

       end if

       T_old(1:z_num) = Temp(1:z_num)
       !write(*,*) Temp

       !-------------- Numerical difference routine when there is snow or not --------!

       if (snw_tot > 0.00001) then
          
          s_l_t = 1
          
          !call Implicit_snow(snw_tot,rho_snow,Tsnw,T_old,T_soil,dt,dz,n,organic_ind,Temp,Cp,Kp,Cp_snow,s_l_t)
       

          !T_soil = T_snw_t(mod(ll,dim_temp)+1)

          call Implicit(T_old,T_soil,Tb,dt,dz,n,organic_ind,Temp,Cp,Kp) 

          if (Bool_Bessi==0) then
             call snw_proc(Tsnw(1),Temp(1), snw_tot, swe_tot, frac_snw, Cp_snow, rho_snow, dt)
          end if
          
       else

          swe_tot = 0.0
          snw_tot = 0.0
          
          do kk =1,s_l_max

             Tsnw(kk) = 0

          end do

          !T_soil = T_snw_t(mod(ll,dim_temp)+1)

          call Implicit(T_old,T_soil,Tb,dt,dz,n,organic_ind,Temp,Cp,Kp) 
          

       end if
       
       !write(*,*) Kp(5)*dt/(Cp(5)*dz(5)*dz(5)) 

       if (Bool_layer_temp==1)then
         !write(*,*) "ok"
         if (ll>t_num-7300)then
             T_layer23(ll+7300-t_num) = Temp(layer_temp23)
             T_layer53(ll+7300-t_num) = Temp(layer_temp53)
             T_layer93(ll+7300-t_num) = Temp(layer_temp93)
             T_layer143(ll+7300-t_num) = (Temp(layer_temp143)+Temp(layer_temp143+1))/2.0
             !T_layer143(ll+7300-t_num) = Temp(layer_temp143)
             !T_layer250(ll+7300-t_num) = (Temp(layer_temp250)+Temp(layer_temp250 - 1))/2.0
             T_layer250(ll+7300-t_num) = Temp(layer_temp250)
             T_layer350(ll+7300-t_num) = Temp(layer_temp350)
             T_layer550(ll+7300-t_num) = Temp(layer_temp550)
             T_layer900(ll+7300-t_num) = Temp(layer_temp900)
            ! write(snw_d,*) snw_tot
          end if
          !write(*,*) "ok"
       end if

       call Permafrost_Depth(Temp,D,Per_depth)
       write(*,*) Per_depth
       
       !write(*,*) Kp

    end do

    t_deb = t_deb - floor(t_step/365.0)

    write(*,*) indice_tab,t_num,organic_ind
    !write(u_n_23,*) Kp

   
    if (Bool_layer_temp==1)then
       write(u_n_23,*) T_layer23
       write(u_n_53,*) T_layer53
       write(u_n_93,*) T_layer93
       write(u_n_143,*) T_layer143
       write(u_n_250,*) T_layer250
       write(u_n_350,*) T_layer350
       write(u_n_550,*) T_layer550
       write(u_n_900,*) T_layer900
    end if

    close(u_n_23)
   
    !write(*,*) "ok"

  end subroutine Vamper_step

  
  
end module Principal

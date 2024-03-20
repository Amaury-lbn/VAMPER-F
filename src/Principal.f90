module Principal


  use Parametrisation, only : z_num,TotTime,Timestep,YearType,z_num,Depth,GridType,PorosityType,T_init,Bool_glacial
  use Parametrisation, only : Bool_Organic,organic_depth,Gfx, T_freeze, EQ_Tr, EQ1_EQ2, Bool_delta,t_fin, alpha,Bool_1998
  use Parametrisation, only : Bool_layer_temp
  use Fonction_temp, only : AppHeatCapacity, ThermalConductivity
  use Fonction_init, only : Porosity_init, GeoHeatFlow, Glacial_index
  use Para_fonctions, only : t_disc, z_disc
  use Model_snow, only : snw_average, snw_proc
  use Fonction_implicit, only : Implicit_snow, Implicit
  
  Implicit none
  
contains

  subroutine Vamper(Temp, Soil_temp, snw_totals)
    integer :: spy, unit_number_1, unit_number_2,  unit_number_3,  unit_number_4, indice_tab
    real :: dt, Cp_snow, frac_snw, K_s, rho_snow, snw_tot, swe_f, Tb, snw_dp, swe_tot, Tsnw, T_soil, days, nb_days
    integer :: t_num, kk, organic_ind, ll, nb_lines,layer_temp23,layer_temp53,layer_temp93,layer_temp143
    integer :: layer_temp250,layer_temp350,layer_temp550,layer_temp900, ind_days
    real, dimension(z_num) :: T_old, Cp, porf, pori, Soil_temp_moy, delta_t_moy
    real, dimension(z_num-1) ::  h_n, h_pori, h_porf
    real, dimension(:),allocatable,intent(out) :: Temp
    real, dimension(:,:),allocatable,intent(out) :: Soil_temp
    real, dimension(:),allocatable,intent(out) :: snw_totals
    real, dimension(:),allocatable:: T_air, Timp, Kp, n, dz, D, Cp_t, swe_f_t, time_gi,glacial_ind
    real, dimension(:,:),allocatable:: delta_T
    character(len=20) :: ligne
    real, dimension(365*20) :: T_layer23,T_layer53,T_layer93,T_layer143,T_layer250,T_layer350,T_layer550,T_layer900
    integer, dimension(12) :: Day_per_month
    integer :: u_n_23,u_n_53,u_n_93,u_n_143,u_n_250,u_n_350,u_n_550,u_n_900

    call t_disc(TotTime,Timestep,YearType,dt,spy,t_num)
    call z_disc(z_num, Depth, GridType, dz, D)

    if (EQ_Tr == 0)then
       if (Bool_1998==0)then
          allocate(T_air(1:12))
       else
          allocate(T_air(1:365))
       end if
       allocate(Soil_temp(1:z_num,1:365))
       allocate(swe_f_t(1:12))
       
    else
       allocate(Soil_temp(1:z_num,1:t_num))
       allocate(snw_totals(1:t_num))
       allocate(T_air(1:t_num))
    end if
    allocate(Kp(1:z_num-1))
    allocate(Temp(1:z_num))

    if (Bool_delta==1)then

       allocate(delta_T(1:z_num,1:t_num-1))

    end if
    
    
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

    if (Bool_layer_temp==1)then

       layer_temp23 = 14
       layer_temp53 = 18
       layer_temp93 = 20
       layer_temp143 = 22
       layer_temp250 = 24
       layer_temp350 = 26
       layer_temp550 = 28
       layer_temp900 = 30
       open(newunit=u_n_23,file="/home/users/alambin/VAMPER-F/Resultats/Temp_layer_23.txt",status="replace",action='write')
       open(newunit=u_n_53,file="/home/users/alambin/VAMPER-F/Resultats/Temp_layer_53.txt",status="replace",action='write')
       open(newunit=u_n_93,file="/home/users/alambin/VAMPER-F/Resultats/Temp_layer_93.txt",status="replace",action='write')
       open(newunit=u_n_143,file="/home/users/alambin/VAMPER-F/Resultats/Temp_layer_143.txt",status="replace",action='write')
       open(newunit=u_n_250,file="/home/users/alambin/VAMPER-F/Resultats/Temp_layer_250.txt",status="replace",action='write')
       open(newunit=u_n_350,file="/home/users/alambin/VAMPER-F/Resultats/Temp_layer_350.txt",status="replace",action='write')
       open(newunit=u_n_550,file="/home/users/alambin/VAMPER-F/Resultats/Temp_layer_550.txt",status="replace",action='write')
       open(newunit=u_n_900,file="/home/users/alambin/VAMPER-F/Resultats/Temp_layer_900.txt",status="replace",action='write') 

    end if

    if (EQ_Tr == 0)then
       
       if(Bool_glacial==1)then
          open(newunit=unit_number_3,file="/home/users/alambin/VAMPER-F/Donnee/Temp_soil_0_x2.txt",status="old",action='read') 
       else
          open(newunit=unit_number_3,file="/home/users/alambin/VAMPER-F/Donnee/Temp_soil_0.txt",status="old",action='read') 
       end if


       open(newunit=unit_number_1,file="/home/users/alambin/VAMPER-F/Donnee/Temp_EQ.txt",status="old",action='read')
       open(newunit=unit_number_2,file="/home/users/alambin/VAMPER-F/Donnee/Snow_EQ.txt",status="old",action='read')
    open(newunit=unit_number_4,file="/home/users/alambin/VAMPER-F/Donnee/Temperature_moyenne_jour.txt",status="old",action='read')

    elseif(EQ_Tr==1)then

       open(newunit=unit_number_3,file="/home/users/alambin/VAMPER-F/Donnee/Temp_EQ_100k_3k.txt",status="old",action='read')
       open(newunit=unit_number_1,file="/home/users/alambin/VAMPER-F/Donnee/Temp_EXP1.txt",status="old",action='read')
       open(newunit=unit_number_2,file="/home/users/alambin/VAMPER-F/Donnee/Snow_EXP1.txt",status="old",action='read')

    end if
    

   
    call Porosity_init(z_num, PorosityType, D, Bool_Organic, organic_depth, n, organic_ind )
     
    
    do kk=1,z_num-1

       Kp(kk)=2
          
    end do

    if (Bool_glacial == 1)then

       call Glacial_index(time_gi,glacial_ind,nb_lines)
       do ll =1,z_num

          read(unit_number_3,*) Temp(ll)

       end do
       write(*,*) Temp

    end if

    if (EQ_Tr == 0)then
       
       if(Bool_1998==1)then
          call Glacial_index(time_gi,glacial_ind,nb_lines)
          do ll =1,365
             read(unit_number_4,*) T_air(ll)
             !write(*,*) T_air(ll)
          end do

          do ll = 1,z_num
             
             read(unit_number_3,*) Temp(ll)

          end do

       end if
       
       do ll =1,12
          
          if (Bool_glacial == 1)then
             read(unit_number_4,*) T_air(ll)
          elseif(Bool_glacial == 0)then
             read(unit_number_1,*) T_air(ll)
          end if

          read(unit_number_2,*) swe_f_t(ll)
          
       end do
       
       if (EQ1_EQ2 ==1 ) then
          
          call GeoHeatFlow(Gfx, Kp, dz, T_init, z_num, Temp)
          write(*,*) Temp

       elseif(EQ1_EQ2==2)then

          do ll =1,z_num

             read(unit_number_3,*) Temp(ll)

          end do
          

       end if
    else

       do ll =1,z_num

          read(unit_number_3,*) Temp(ll)

       end do
       
       read(unit_number_1,*) T_air(1)
       read(unit_number_2,*) swe_f
       
    end if

    !write(*,*) sum(T_air)/365.0
    !write(*,*) "ok"

    Tb = Temp(z_num)
    write(*,*) Tb

    do ll = 1,2
       
       call AppHeatCapacity(z_num,Temp,T_freeze,n, organic_ind, Cp, porf, pori)
       
       do kk=1,z_num-1
          
          h_pori(kk) = (pori(kk) + pori(kk+1))/2
          h_porf(kk) = (porf(kk) + porf(kk+1))/2
          h_n(kk) = (n(kk) + n(kk+1))/2
          
          call ThermalConductivity(kk,h_n(kk),h_pori(kk),h_porf(kk), organic_ind, Temp(kk), Kp(kk))
          
       end do
       
    end do

    if (Bool_delta==1)then
       delta_t(1:z_num,1) = Temp(1:z_num)
    end if
             

    do ll=1,t_num
       
       !write(*,*)  T_air(mod(ll,12)+1)+9.0
       
       nb_days = nb_days + 1

       if (nb_days > days) then
          nb_days = 1
          ind_days = ind_days + 1

          if (ind_days == 13) then
             ind_days = 1
          end if

          days = Day_per_month(ind_days)

       end if

       if (EQ_Tr == 0)then
          
          if (EQ1_EQ2 == 1)then
             T_soil = (T_air(mod(ll,12)+1)+7.71484184)
          elseif (EQ1_EQ2==2)then
             T_soil = T_air(mod(ll,12)+1)
          end if

          if(Bool_1998==1)then
             !write(*,*) "ok"
             indice_tab = nb_lines-floor((-(ll/(365.0))+t_fin+TotTime)/100.0)
             T_soil=alpha*(glacial_ind(indice_tab-1)+mod((ll/12.0),100.0)*(glacial_ind(indice_tab)-glacial_ind(indice_tab-1))/100.0)
             T_soil = T_air(mod(ll,365)+1)+T_soil
             !write(*,*) indice_tab

             !write(*,*) "ok"
          end if

          !write(*,*) int(days)

          swe_f = swe_f_t(ind_days)/(Day_per_month(ind_days))
          
       else

          read(unit_number_1,*) T_air(ll)
          read(unit_number_2,*) swe_f

       end if

       if (Bool_glacial==1)then

          !write(*,*) "ok", ll
          indice_tab = nb_lines-floor((-(ll/12.0)+t_fin+TotTime)/100.0)
          T_soil=alpha*(glacial_ind(indice_tab-1)+mod((ll/12.0),100.0)*(glacial_ind(indice_tab)-glacial_ind(indice_tab-1))/100.0)
          !write(*,*) "ok", T_soil
          T_soil = T_air(mod(ll,12)+1)+T_soil

          !T_soil = T_air(mod(ll,12)+1) + alpha*glacial_ind(indice_tab)

          !write(*,*) "ok", T_soil

       end if

       T_old(1:z_num) = Temp(1:z_num)
       
       
       if (snw_tot > 0.000001 .or. swe_f > 0.000001) then
          
          call snw_average(swe_f, swe_tot, snw_tot, rho_snow)

          !snw_totals(ll) = snw_tot
          
          if (abs(swe_tot - swe_f )< 0.00001) then

             if (EQ_Tr == 0)then

                Tsnw = T_soil

             else
                
                Tsnw = T_air(ll)
                
             end if

             K_s = 0.07
             frac_snw = 1
          
          end if
       
       end if

       if (snw_tot > 0.000001) then

          if (EQ_Tr == 0)then

             call Implicit_snow(snw_tot,rho_snow,Tsnw,T_old,T_soil,dt,dz,n,organic_ind,Temp,Cp_t,Kp,Cp_snow)
          
          else

             call Implicit_snow(snw_tot,rho_snow,Tsnw,T_old,T_air(ll),dt,dz,n,organic_ind,Temp,Cp_t,Kp,Cp_snow)

          end if

          !write(*,*) Cp_snow

          call snw_proc(Tsnw, snw_tot, swe_tot, frac_snw, Cp_snow, rho_snow, dt)
          
       else

          swe_tot = 0.0
          snw_tot = 0.0
          
          if (EQ_Tr == 0)then

            call Implicit(T_old,T_soil,Tb,dt,dz,n,organic_ind,Temp,Cp_t,Kp) 
          
         else

            call Implicit(T_old,T_air(ll),Tb,dt,dz,n,organic_ind,Temp,Cp_t,Kp)

         end if
    
       end if

       !write(*,*) ll, T_soil
          
       do kk=1,z_num
          
          if (EQ_Tr==0) then

             Soil_temp(kk,mod(ll,365)+1) = Temp(kk)

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

       !write(*,*) snw_tot

       if (Bool_glacial == 1) then
          
          if (mod(ll,1200*365) == 0) then

             do kk=1,z_num
                
                Soil_temp_moy(kk) = sum(Soil_temp(kk,1:12))/12.0

                if (Bool_delta==1)then

                   delta_t_moy(kk)=sum(delta_t(kk,(ll-11):ll))

                end if
                
             end do
             write(*,*) ll
             write(*,*) "Soil_temp_moy =" , Soil_temp_moy
             !write(*,*) "delta_t_moy =" , delta_t_moy 

          end if
          
          if (ll>t_num-240)then
             T_layer23(240+ll-t_num) = Temp(layer_temp23)
             T_layer53(240+ll-t_num) = Temp(layer_temp53)
             T_layer93(240+ll-t_num) = Temp(layer_temp93)
             T_layer143(240+ll-t_num) = Temp(layer_temp143)
             T_layer250(240+ll-t_num) = (Temp(layer_temp250)+Temp(layer_temp250+1))/2.0
             T_layer350(240+ll-t_num) = Temp(layer_temp350)
             T_layer550(240+ll-t_num) = Temp(layer_temp550)
             T_layer900(240+ll-t_num) = Temp(layer_temp900)
          end if

       end if

       if (Bool_1998 == 1) then

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
             T_layer250(ll+7300-t_num) = (Temp(layer_temp250)+Temp(layer_temp250 + 1))/2.0
             T_layer350(ll+7300-t_num) = Temp(layer_temp350)
             T_layer550(ll+7300-t_num) = Temp(layer_temp550)
             T_layer900(ll+7300-t_num) = Temp(layer_temp900)
          end if

       end if
       
    end do

    !write(*,*) Temp
    
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

    !write(*,*) "T_layer23 = ",T_layer23
    !write(*,*) "T_layer53 = ",T_layer53 
    !write(*,*) "T_layer93 = ",T_layer93
    !write(*,*) "T_layer143 = ",T_layer143
    !write(*,*) "T_layer250 = ",T_layer250
    !write(*,*) "T_layer350 = ",T_layer350
    !write(*,*) "T_layer550 = ",T_layer550
    !write(*,*) "T_layer900 = ",T_layer900 
    !write(*,*) D
    !write(*,*) t_num


  end subroutine Vamper
  
  
end module Principal

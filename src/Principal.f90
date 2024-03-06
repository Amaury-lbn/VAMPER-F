module Principal


  use Parametrisation, only : z_num,TotTime,Timestep,YearType,z_num,Depth,GridType,PorosityType,T_init
  use Parametrisation, only : Bool_Organic,organic_depth,Gfx, T_freeze, EQ_Tr
  use Fonction_temp, only : AppHeatCapacity, ThermalConductivity
  use Fonction_init, only : Porosity_init, GeoHeatFlow
  use Para_fonctions, only : t_disc, z_disc
  use Model_snow, only : snw_average, snw_proc
  use Fonction_implicit, only : Implicit_snow, Implicit
  
  Implicit none
  
contains

  subroutine Vamper(Temp, Soil_temp, snw_totals)
    integer :: spy, unit_number_1, unit_number_2,  unit_number_3
    real :: dt, Cp_snow, frac_snw, K_s, rho_snow, snw_tot, swe_f, Tb, snw_dp, swe_tot, Tsnw, T_soil
    integer :: t_num, kk, organic_ind, ll
    real, dimension(z_num) :: T_old, Cp, porf, pori
    real, dimension(z_num-1) ::  h_n, h_pori, h_porf
    real, dimension(:),allocatable,intent(out) :: Temp
    real, dimension(:,:),allocatable,intent(out) :: Soil_temp
    real, dimension(:),allocatable,intent(out) :: snw_totals
    real, dimension(:),allocatable:: T_air, Timp, Kp, n, dz, D, Cp_t, swe_f_t
    character(len=20) :: ligne
    
    call t_disc(TotTime,Timestep,YearType,dt,spy,t_num)
    call z_disc(z_num, Depth, GridType, dz, D)

    if (EQ_Tr == 0)then
       allocate(T_air(1:12))
       allocate(Soil_temp(1:z_num,1:12))
       allocate(swe_f_t(1:12))
    else
       allocate(Soil_temp(1:z_num,1:t_num))
       allocate(snw_totals(1:t_num))
       allocate(T_air(1:t_num))
    end if
    allocate(Kp(1:z_num-1))
    allocate(Temp(1:z_num))
    
    
    Tsnw = -4.0
    frac_snw = 1.0
    rho_snow = 0.0
    snw_tot = 0.0
    swe_tot = 0.0


    if (EQ_Tr == 0)then
    
       open(newunit=unit_number_3,file="/home/users/alambin/VAMPER-F/Donnee/Temp_EQ_NoPF_Porolin.txt",status="old",action='read') 
       open(newunit=unit_number_1,file="/home/users/alambin/VAMPER-F/Donnee/Temp_EQ.txt",status="old",action='read')
       open(newunit=unit_number_2,file="/home/users/alambin/VAMPER-F/Donnee/Snow_EQ.txt",status="old",action='read')

    else

       open(newunit=unit_number_3,file="/home/users/alambin/VAMPER-F/Donnee/Timp.txt",status="old",action='read')
       open(newunit=unit_number_1,file="/home/users/alambin/VAMPER-F/Donnee/Temp_EXP1.txt",status="old",action='read')
       open(newunit=unit_number_2,file="/home/users/alambin/VAMPER-F/Donnee/Snow_EXP1.txt",status="old",action='read')

    end if
    
   
    call Porosity_init(z_num, PorosityType, D, Bool_Organic, organic_depth, n, organic_ind )
     
    
    do kk=1,z_num-1

       Kp(kk)=2
          
    end do
       
    if (EQ_Tr == 0)then
       
       call GeoHeatFlow(Gfx, Kp, dz, T_init, z_num, Temp)
       
       do ll =1,12

          read(unit_number_1,*) T_air(ll)
          read(unit_number_2,*) swe_f_t(ll)
          
       end do

       !do ll =1,z_num

        !  read(unit_number_3,*) Temp(ll)

       !end do
    else

       do ll =1,z_num

          read(unit_number_3,*) Temp(ll)

       end do
       
       read(unit_number_1,*) T_air(1)
       read(unit_number_2,*) swe_f
       
    end if

    write(*,*) "ok"

    Tb = Temp(z_num)

    do ll = 1,2
       
       call AppHeatCapacity(z_num,Temp,T_freeze,n, organic_ind, Cp, porf, pori)
       
       do kk=1,z_num-1
          
          h_pori(kk) = (pori(kk) + pori(kk+1))/2
          h_porf(kk) = (porf(kk) + porf(kk+1))/2
          h_n(kk) = (n(kk) + n(kk+1))/2
          
          call ThermalConductivity(kk,h_n(kk),h_pori(kk),h_porf(kk), organic_ind, Temp(kk), Kp(kk))
          
       end do
       
    end do


    do ll=2,t_num+9
       
       write(*,*) ll
       

       if (EQ_Tr == 0)then

          T_soil = T_air(mod(ll,12)+1)+9.0
          swe_f = swe_f_t(mod(ll,12)+1)
          
       else

          read(unit_number_1,*) T_air(ll)
          read(unit_number_2,*) swe_f

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
          
       do kk=1,z_num

          Soil_temp(kk,mod(ll,12)+1) = Temp(kk)

       end do

       !snw_totals(ll) = snw_tot
       
    end do




  end subroutine Vamper
  
  
end module Principal

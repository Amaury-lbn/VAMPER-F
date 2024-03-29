

module Fonction_temp

  use Parametrisation, only : C_ice,C_dry_soil,C_organic,C_water,freezing_range,Latent_Heat,rho_ice,z_num
  use Parametrisation, only :rho_organic,rho_soil,rho_water,K_other_minerals,K_organic,K_quartz,q_quartz

  Implicit none

  contains



  subroutine AppHeatCapacity(z_num, T, Tf, n, org_ind, Cp, porf, pori)
    implicit none
    integer, intent(in) :: z_num, org_ind
    real, intent(in) :: Tf
    real, dimension(z_num), intent(in) :: T,n
    real, dimension(z_num), intent(inout) :: Cp
    real, dimension(z_num), intent(out) ::  porf, pori
    integer :: kk
    real :: dTheta, theta, a, Csoil
    
    
    if (org_ind > 1) then

       Csoil = (1.0 - n(1)) * rho_organic * C_organic
    else
       Csoil = (1.0 - n(1)) * rho_soil * C_dry_soil
    end if
    do kk = 1, z_num
       if (kk <= 0) then
          Csoil=((1 - n(kk))* rho_organic * C_organic)
       else
          Csoil= 1 * 10.0**6 
       end if
       
       if (T(kk) < Tf) then
          a = - (((T(kk) - Tf) / freezing_range) ** 2.0)
          theta = exp(a)
          a = -2.0 / (freezing_range * freezing_range)
          dTheta = a * (T(kk) - Tf) * theta
          porf(kk) = n(kk) * theta
          pori(kk) = n(kk) - porf(kk)
          Cp(kk) = Csoil + (pori(kk) * C_ice * rho_ice) + (porf(kk) * C_water *rho_water)+(n(kk)*rho_water*Latent_Heat*dTheta)
       else
          porf(kk) = n(kk)
          pori(kk) = 0.0
          Cp(kk) = Csoil + (porf(kk) * C_water * rho_water)
          
       end if

    end do

  end subroutine AppHeatCapacity


  subroutine ThermalConductivity(layer , h_n, h_pori, h_porf, org_ind, Temp, Ther_cond )
    
    integer, intent(in) :: layer, org_ind
    real, intent(in) :: h_n, h_pori, h_porf, Temp
    real, intent(out) :: Ther_cond
    real :: Ksoil, Kice, Kfluids
    
    if (org_ind > layer) then

       Ksoil = K_organic
       
    else

       Ksoil = (K_quartz ** q_quartz) * (K_other_minerals ** (1.0-q_quartz))


    end if

    Kfluids = 0.1145 + 0.0016318 * (273.15 + Temp)

    Kice = 0.4865 + 488.19/(273.15 +Temp)

    Ther_cond = 0.5 * (Ksoil**(1-h_n) * Kice**(h_pori) * Kfluids**(h_porf))

  end subroutine ThermalConductivity

  
end module Fonction_temp

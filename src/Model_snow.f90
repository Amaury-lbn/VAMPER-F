

module Model_snow

  use Parametrisation

  contains


    subroutine snw_average(swe_f, swe_tot, snw_tot, rho_snow)
      
      real, intent(in) :: swe_f
      real, intent(inout) :: swe_tot, snw_tot, rho_snow
      real :: snw_f
      
      
      if (swe_f > 0.000001) then

         snw_f = swe_f * (rho_water/rho_snow_freeze)
         snw_tot = snw_tot + snw_f
         swe_tot = swe_tot + swe_f
         
         if (snw_tot > 0.000001) then
            
            rho_snow = swe_tot*rho_water/snw_tot
            
         else
            
            rho_snow = rho_snow_freeze

         end if

      end if


    end subroutine snw_average


    subroutine snw_proc(Tsnw, snw_tot, swe_tot, frac_snw, Cp_snow, rho_snow, dt)

      real, intent(inout) :: Tsnw, snw_tot, swe_tot, frac_snw
      real, intent(in) :: Cp_snow, rho_snow, dt
      real :: frac_snw_old, rho_snow_old, swe_melt, H_1, H_2, N, var, rho_snow_new

      swe_melt = 0.0
      frac_snw_old = frac_snw
      rho_snow_old = rho_snow

      if ( Tsnw >= 0 ) then
         
         H_1 = (Tsnw + 273.15)*Cp_snow*snw_tot
         H_2 = rho_water * Latent_heat * swe_tot

         if (H_1 <= H_2) then
            
            frac_snw = 1-(H_1/H_2)
            snw_tot = frac_snw/frac_snw_old * snw_tot
            Tsnw = 0

         else
            
            frac_snw = 0
            snw_tot = 0
            swe_tot = 0
            Tsnw = 0
         end if
        
      else
         
         frac_snw = frac_snw_old
         
      end if
         
      N = 0.5*snw_tot*rho_snow*5.0*10**(-8)*rho_snow*gravity
      var = 14.643 - 4000/min(Tsnw+273.16,273.16)-0.02*rho_snow
      rho_snow_new = rho_snow + N*dt*exp(var)

      snw_tot = (rho_snow_old/rho_snow)*snw_tot
      
    end subroutine snw_proc

   
end module Model_snow

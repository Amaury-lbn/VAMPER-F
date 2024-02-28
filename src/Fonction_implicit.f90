module Fonction_implicit


  use Parametrisation, only : z_num, s_l, rho_ice, Gfx, T_freeze
  use Fonction_temp, only : AppHeatCapacity, ThermalConductivity

  Implicit none

  contains

    subroutine Implicit_snow(snw_dp,rho_snow,Tsnw,T_old,Tu,dt,dz,n,org_ind,Timp,Cp,Kp,Cp_snow)
      
      integer, intent(in) :: org_ind
      real, intent(in) ::  snw_dp, rho_snow, dt,Tu
      real, dimension(:), intent(in) :: T_old, n, dz
      real, intent(inout) :: Tsnw
      real, intent(out) :: Cp_snow
      real, dimension(:), allocatable, intent(out) :: Timp, Cp, Kp
    
      ! Variables locales !

      real, dimension(z_num+s_l) ::  Kp_s, T_last, n_s, dz_s, K_s
      real, dimension(:), allocatable :: Cp_s, porf, pori
      real, dimension(:), allocatable :: Cp_temp, porf_temp, pori_temp
      real, dimension(5,z_num+s_l) :: T_iter
      real, dimension(z_num+s_l,z_num+s_l) :: MM
      real, dimension(1,z_num+s_l) :: Solution
      real :: m_Gfx, A, B, C, Z1
      integer :: kk, ll, z_s

      allocate(Timp(1:z_num))
      allocate(Cp(1:z_num))
      allocate(Kp(1:z_num))
      allocate(porf(1:z_num+s_l))
      allocate(pori(1:z_num+s_l))
      allocate(Cp_s(1:z_num+s_l))
      allocate(Cp_temp(1:z_num))
      allocate(porf_temp(1:z_num))
      allocate(pori_temp(1:z_num))

      MM=0
      z_s = z_num+s_l
      m_Gfx = Gfx/1000.0

      T_last(s_l+1:z_s) = T_old(1:z_num)
      dz_s(1:s_l) = snw_dp

      dz_s(s_l+1:z_s) = dz
      n_s(s_l+1:z_s) = n
      
      Cp_s(1:s_l) = (1.9*10**6)*rho_snow/rho_ice
      K_s(1:s_l) = (2.9*10*(-6))*(rho_snow**2)

      T_last(1:s_l) = Tsnw
      
      
      do kk=1,5 
         
         if (kk==1) then
            
            T_iter(1, 1:z_s) = T_last(1:z_s)
            
         else

            T_iter(1,1:z_s)=0.5*(T_iter(1,1:z_s)+T_last(1:z_s))
         end if


         call AppHeatCapacity(z_num,T_iter(1,s_l+1:z_s),T_freeze,n_s(s_l+1:z_s),org_ind+s_l,Cp_temp,pori_temp,porf_temp)

         do ll=s_l+1,z_s
            
            call ThermalConductivity(ll,n_s(ll),pori_temp(ll-1),porf_temp(ll-1),org_ind+s_l,T_iter(1,ll),K_s(ll))
            
         end do

         Kp_s(1:z_s-1) = (K_s(1:z_s-1)+K_s(2:z_s))/2.0

         Solution(1,1) = Tu
         
         do ll=2,z_s-1
            
            Z1 = T_last(ll)
            
            A=(dt/((dz_s(ll-1)+dz_s(ll))*0.5*dz_s(ll)) * (Kp_s(ll-1)/Cp_temp(ll-1)))
            B=(dt/((dz_s(ll+1)+dz_s(ll))*0.5*dz_s(ll)) * (Kp_s(ll)/Cp_temp(ll-1)))
            C= 1+A+B

            MM(ll,ll-1) = -A
            MM(ll,ll) = C
            MM(ll,ll+1) = -B
            Solution(1,ll) = Z1

         end do
           
         A=(dt/((dz_s(z_num)+dz_s(z_num))*0.5*dz_s(z_num)) * (Kp_s(z_num)/Cp_temp(z_num)))
         C=1+A
         MM(z_s,z_s-1)=-A
         MM(z_s,z_s)=C

         Solution(1,z_s) = T_last(z_s) + ((dt/Cp_temp(z_num))*m_Gfx/dz_s(z_num))
         
      end do

    end subroutine Implicit_snow


    subroutine Implicit(T_old,Tu,Tb,dt,dz,n,org_ind,Timp,Cp,Kp)
      
      integer, intent(in) :: org_ind
      real, intent(in) :: dt,Tu,Tb
      real, dimension(:), intent(in) :: T_old, n, dz
      real, dimension(:), allocatable, intent(out) :: Timp, Cp, Kp
      
      real, dimension(:),allocatable :: pori, porf
      real, dimension(z_num,z_num) :: MM
      real, dimension(1,z_num) :: Solution
      real :: m_Gfx, A, B, C, Z1
      integer :: kk, ll
      real, dimension(z_num) :: T_last
      real, dimension(5,z_num) :: T_iter
      
      m_Gfx = gfx/1000.0

      allocate(Timp(1:z_num))
      allocate(Cp(1:z_num))
      allocate(Kp(1:z_num))
      allocate(pori(1:z_num))
      allocate(porf(1:z_num))

      T_last = T_old

      do kk=1,5 
         
         if (kk==1) then
            
            T_iter(1, 1:z_num) = T_last(1:z_num)
            
         else

            T_iter(1,1:z_num)=0.5*(T_iter(1,1:z_num)+T_last(1:z_num))

         end if


         call AppHeatCapacity(z_num,T_iter(1,1:z_num),T_freeze,n,org_ind,Cp,porf,pori)

         do ll=1,z_num
            
            call ThermalConductivity(ll,n(ll),pori(ll-1),porf(ll-1),org_ind,T_iter(1,ll),Kp(ll))
            
         end do
      

         do ll=2,z_num-1
            
            Z1 = T_last(ll)
            
            A=(dt/((dz(ll-1)+dz(ll))*0.5*dz(ll)) * (Kp(ll-1)/Cp(ll)))
            B=(dt/((dz(ll+1)+dz(ll))*0.5*dz(ll)) * (Kp(ll)/Cp(ll)))
            C= 1+A+B

            MM(ll,ll-1) = -A
            MM(ll,ll) = C
            MM(ll,ll+1) = -B
            Solution(1,ll) = Z1

         end do
           
         Solution(1,z_num)=Tb
         MM(z_num,z_num)=1
         MM(1,1)=1

      end do
   

    end subroutine Implicit


end module Fonction_implicit

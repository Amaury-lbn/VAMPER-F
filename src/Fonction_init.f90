

module Fonction_init

  use Parametrisation, only : n_organic, Porosity_soil
  use Para_fonctions

  Implicit none


  contains

  !Routine paramettrant la porosité du sol en fonction de la profondeur!

  subroutine Porosity_init(z_num, Porosity_Type, Depth_layer, Bool_Organic, organic_depth, n, organic_ind)

    ! Entrées et sorties !

    integer, intent(in) :: z_num, Porosity_Type, Bool_Organic
    real, intent(in) ::  organic_depth
    real, dimension(:), intent(in) :: Depth_layer
    integer, intent(out) :: organic_ind
    real, dimension(:), allocatable, intent(out) :: n
    
    ! Variables locales !

    real, dimension(:), allocatable :: DepthCalcOrg
    real :: pente
    real :: origine
    integer :: pas_z

    allocate(n(1:z_num))
    allocate(DepthCalcOrg(1:z_num))
     
    if (Bool_organic == 1) then
       

       do pas_z = 1,z_num

          DepthCalcOrg(pas_z) = abs(Depth_layer(pas_z) - organic_depth)
       
       end do

       call indice_minimum(DepthCalcOrg, z_num, organic_ind)
       
       do pas_z = 1,organic_ind
          
          n(pas_z) = n_organic

       end do

    else
       
       n(1) = Porosity_soil
       organic_ind = 1

    end if
    
    

    if (Porosity_Type == 1) then
       
       
       pente = (Porosity_soil -0.15)/(Depth_layer(organic_ind) - Depth_layer(z_num))
       origine = Porosity_soil - pente*Depth_layer(organic_ind)

       do pas_z = organic_ind, z_num
          
          n(pas_z) = pente * Depth_layer(pas_z) + origine

       end do

    else
       
       do pas_z = organic_ind, z_num
          
          n(pas_z) = Porosity_soil*exp(0.000395*Depth_layer(pas_z))

       end do

    end if

  end subroutine Porosity_init




  subroutine GeoHeatFlow(Gfx, Kp, dz, T0, z_num, T)
   
    integer, intent(in) :: z_num
    real, intent(in) :: Gfx, T0
    real, dimension(:), allocatable, intent(in) :: dz, Kp 
    real, dimension(:), allocatable, intent(out) :: T
    integer :: kk
    
    allocate(T(1:z_num))

    T(1) = T0
  
    do kk = 2, z_num
       T(kk) = T(kk-1) + (((Gfx/1000.0)/Kp(kk-1))*dz(kk-1))
    end do



  end subroutine GeoHeatFlow


end module Fonction_init
   



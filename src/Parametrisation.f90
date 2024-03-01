module Parametrisation

  
  Implicit none

  integer,parameter :: TotTime = 451         !temps total en année
  integer,parameter :: Timestep = 30          !nombre de jour entre chaque pas de temps 
  integer,parameter :: YearType = 360         !nombre de jour par an
  integer,parameter :: z_num = 100             !nombre de couche étudiée
  integer,parameter :: EXPE = 1                !de 1 à 4, quelle expérience va être réalisée
  integer,parameter :: GridType = 1           !(1) Log-generated, (2) Linear-generated 
  integer,parameter :: PorosityType = 1       !(1) linéaire, (2)constante, (3)exponentiellement décroissante en fonction de la profondeur
  integer,parameter :: Bool_Snow = 1          ! forçage en neige ou non (1 ou 0)
  integer,parameter :: Bool_Organic = 1       ! prise en compte de la couche organique ou non (1 ou 0)


  real,parameter :: Depth = 1000              !profondeur de la modélisation
  real,parameter :: T_init = 0                !température initiale a la surface
  real,parameter :: T_freeze = 0             !température où l'eau est considérée comme gelée
  real,parameter :: freezing_range = 1        
  real,parameter :: Gfx = 65                  ! flux géothermique de la terre (a modifier peut être)
  real,parameter :: Porosity_soil = 0.5       ! porosité du sol
  real,parameter :: organic_depth = 0.025     ! profondeur de la couche organique
  real,parameter :: n_organic = 0.5           ! porosité de la couche organique
  real,parameter :: n_soil_bot = 0.5          ! porosité en bas de la couche de sol


  !       DENSITÉ DE DIFFÉRENTES MATIÈRES (en kg/m³) ! 

  real,parameter :: rho_snow_freeze = 350            ! densité de la neige
  real,parameter :: rho_water = 1000          ! densité de l'eau
  real,parameter :: rho_ice = 917             ! densité de la glace
  real,parameter :: rho_organic = 1300        ! densité de la matière organique
  real,parameter :: rho_soil = 1600           ! densité du sol

  !      Capacité thermique massique  (en J/(K*kg))    !

  real,parameter :: C_water = 4180          ! Capacité thermique massique de l'eau
  real,parameter :: C_ice = 2100            ! Capacité thermique massique de la glace
  real,parameter :: C_organic = 1920        ! Capacité thermique massique de la matière organique
  real,parameter :: C_dry_soil = 850        ! Capacité thermique massique du sol


  !      Conductivité thermique (en W/(m*K))          !

  real,parameter :: K_other_minerals = 2    ! Conductivité thermique des autres minéraux
  real,parameter :: K_quartz = 6            ! conductivité thermique du quartz
  real,parameter :: K_organic = 0.25        ! conductivité thermique de la matière organique
  real,parameter :: K_ice = 2.24            ! conductivité thermique de la glace
  real,parameter :: K_fluids = 0.56         ! conductivité thermique des fluides

  real, parameter :: q_quartz=0.5            ! pourcentage de quartz dans le sol

  real, parameter :: gravity = 9.81          ! accéleration gravitationnelle 
  real, parameter :: Latent_heat = 333700    ! en J/kg

  integer, parameter :: s_l = 1      ! nombre de couche de neige (marche que avec 1)
  
end module Parametrisation























program test_fonctions

  use Parametrisation
  use Fonction_init
  use Para_fonctions
  use Fonction_temp

  implicit none
  integer spy
  real dt
  integer t_num, kk, organic_ind, ll
  real, dimension(:), allocatable :: dz,D,Kp,T,n
  real, dimension(:), allocatable :: Cp, porf, pori, h_n, h_pori, h_porf
  
  allocate(dz(1:z_num-1))
  allocate(D(1:z_num))
  allocate(Kp(1:z_num-1))
  allocate(T(1:z_num))
  allocate(n(1:z_num))
  allocate(Cp(1:z_num))
  allocate(porf(1:z_num))
  allocate(pori(1:z_num))
  allocate(h_porf(1:z_num-1))
  allocate(h_pori(1:z_num-1))
  allocate(h_n(1:z_num-1))

  do kk=1,z_num-1
     Kp(kk)=2
  end do

  call z_disc(z_num, Depth, GridType, dz, D)

  call GeoHeatFlow(Gfx, Kp, dz, T_init, z_num, T)

  call Porosity_init(z_num, PorosityType, D, Bool_Organic, organic_depth, n, organic_ind )



  !write(*,*) "coucou", dt

!  write(*,*) "coucou", Cp
  
!  write(*,*) "coucou", pori
  
!  write(*,*) "coucou", porf

  do ll = 1,2

     call AppHeatCapacity(z_num,T,T_freeze,n, organic_ind, Cp, porf, pori)
 
     do kk=1,z_num-1

     h_pori(kk) = (pori(kk) + pori(kk+1))/2
     h_porf(kk) = (porf(kk) + porf(kk+1))/2
     h_n(kk) = (n(kk) + n(kk+1))/2

     call ThermalConductivity(kk,h_n(kk),h_pori(kk),h_porf(kk), organic_ind, T(kk), Kp(kk))

     end do

  end do

  write(*,*) "coucou", Kp

  write(*,*) "coucou", Cp

  write(*,*) "coucou", T

 end program test_fonctions

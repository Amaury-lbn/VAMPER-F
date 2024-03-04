module Para_fonctions

Implicit none

contains




  !sous routine trouvant les parametre permettant la discrétisation en temps!   ! a changer plus tard !

  subroutine t_disc(TotTime,Timestep,YearType,dt,spy,t_num)   !fonctionnelle!
    
    Implicit none

    integer,intent(in) ::  TotTime,Timestep,YearType
    integer,intent(out) :: t_num, spy
    real, intent(out) :: dt
    real :: yrs2days
    real :: model_secs

    if (YearType == 360) then
       yrs2days = 360
    else
       yrs2days = 365
    end if

    
    model_secs = TotTime*60*60*24*yrs2days
    dt = Timestep*60*60*24
    t_num = floor(model_secs/dt)


    if (timestep == 15) then
       spy = 24
    elseif (timestep == 30) then
       spy = 12
    elseif (timestep == 1) then
       spy = 360
    else
       spy = 1
    end if

  end subroutine t_disc
    
  

!sous routine trouvant les parametre permettant la discrétisation selon z (la profondeur)!

  subroutine z_disc(z_num, Depth, GridType, dz, D)     !fonctionnelle!
    
    Implicit none 

    real, intent(in) :: Depth
    integer, intent(in) :: Gridtype, z_num
    real, dimension(:), allocatable, intent(out) :: dz,D
    real, dimension(:), allocatable :: D_temp
    integer :: kk


    allocate(dz(1:z_num))
    allocate(D(1:z_num))
    allocate(D_temp(1:z_num+1))


    if ( GridType == 2 ) then
        call linspace(real(0),Depth,z_num+1,D)
    else
         call logspace(real(-2),log10(Depth),z_num+1,D_temp)
    end if


    do kk = 1, z_num
        dz(kk)=D_temp(kk+1)-D_temp(kk)
    end do


    D(1)=dz(1)


    do kk = 2, z_num
        D(kk)=D(kk-1) + dz(kk)
    end do

    !write(*,*) "coucou", dz
    !write(*,*) "coucou", D
    
  end subroutine z_disc




! Routine permettant de créer un tableau de n valeur logarithmiquement espacé entre 10^deb et 10^fin !

  subroutine logspace(d_val, f_val, num, res)  !fonctionnelle!

    implicit none

    real, intent(in) :: d_val, f_val
    integer, intent(in) :: num
    real, dimension(:), allocatable, intent(out) :: res

    integer :: i

    allocate(res(1:num))

    do i = 1, num
        res(i) = 10.0 ** (d_val + (f_val - d_val) * real(i - 1) / real(num - 1))
    end do


  end subroutine logspace




! Routine permettant de créer un tableau de n valeur espacées linéairement entre deb et fin  !

  subroutine linspace(d_val, f_val, num, res)

    implicit none

    real, intent(in) :: d_val, f_val
    integer, intent(in) :: num
    real, dimension(:), allocatable, intent(out) :: res

    integer :: i

    allocate(res(1:num))


    do i = 1, num
        res(i) = (d_val + (f_val - d_val) * real(i - 1) / real(num - 1))
    end do


  end subroutine linspace
    

  

  subroutine indice_minimum(tab, taille, ind_min)
    
    implicit none
    
    integer, intent(in) :: taille
    real, dimension(:), intent(in) :: tab
    integer, intent(out) :: ind_min

    integer :: kk

    ind_min = 1

    do kk = 2,taille,1

       if(tab(kk)<tab(ind_min)) then

          ind_min = kk

       end if

    end do



  end subroutine indice_minimum
end module Para_fonctions
    

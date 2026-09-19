!   PROJET DE MODELISATION DES OPERATIONS UNITAIRES - ENSGTI 3 A
!   Réalisé par :                               Encadrant :
!   Aboubacar Sidiki BANGOURA                   Frédéric MARIAS
!
!   Parcours : Conception des Procédés Assistés par Ordinateur CPAO
!   Année Universitaire: 2022-2023
! ----------------------------------------------------------------------------
!   But: Calcul modelisation d'une colonne de distillation biphasique a plusieurs
!                        plateaux en regime permanent
! ----------------------------------------------------------------------------
!   Ficiher : Subroutine calcul de la temperature de rosée
! ----------------------------------------------------------------------------

subroutine CalculTrosee(Yi,Xi,Trosee)
    use dimensions
    Use Thermos
    use operatoires
    implicit none

    double precision :: dh, alpha, ERREUR
    double precision, dimension(NC),intent(in) :: Yi
    doubleprecision, dimension(NC) :: Gama, hold, Psat
    doubleprecision, dimension(NC) :: Gama_dh, Psat_dh, dh_array ! variables sur la variation au niveau de la jacobienne
    doubleprecision, dimension(NC+1) :: Xsol, F ! Variable et fonction objective
    doubleprecision, dimension(NC+1,NC+1) :: Jaco
    integer :: i, IER, INDIC, it, itmax
    doubleprecision :: Temp
    double precision, dimension(NC), intent(out) :: Xi
    doubleprecision, intent(out) :: Trosee

    ! Initialisation

    IER = 0
    INDIC = 1
    alpha = 0.1
    itmax = 1000
    it = 0
    ERREUR = 1d0
    dh = 1d-3
    Jaco = 0

    ! comment avoir le Temp

    if (Nb_alim==1) then
        Temp=Talim(Num_alim(1))
    else
        Temp=Talim(Num_alim(Nb_alim))
    endif

    do i = 1,NC
        Xsol(i) = Yi(i)
    end do

    Xsol(NC+1)= Temp

    call CalculGama(Xsol(NC+1),Xsol(1:NC),Gama)
    call Psaturation(Xsol(NC+1),Psat)

    do while (it < itmax .and. ERREUR.GT.1E-5)

        it = it + 1

        F = 0

        do i = 1,NC
            F(i) = Xsol(i)*Psat(i)*Gama(i) - Yi(i)*P
            F(NC+1) = F(NC+1) - Xsol(i)
        end do
        F(NC+1) = F(NC+1) + 1

        ERREUR = norm2(F)

        ! Jacobienne du système

        call CalculGama(Xsol(NC+1),Xsol(1:NC),Gama)

        do i = 1,NC
            dh_array = 1.0d0
            dh_array(i) = 1.0d0 + dh
            call CalculGama(Xsol(NC+1),Xsol(1:NC)*dh_array,hold)
            Gama_dh(i) = hold(i)
        end do

        call Psaturation(Xsol(NC+1),Psat)

        do i = 1,NC
            Jaco(i,i) = ( (Xsol(i)*(1+dh)*Psat(i)*Gama_dh(i)) - (Xsol(i)*Psat(i)*Gama(i)) )/(Xsol(i)*(dh))
            Jaco(NC+1,i) = -1
        end do

        call CalculGama(Xsol(NC+1)*(1+dh),Xsol(1:NC),Gama_dh)
        call Psaturation(Xsol(NC+1)*(1+dh),Psat_dh)

        do i = 1,NC
            Jaco(i,NC+1) = ( (Xsol(i)*Psat_dh(i)*Gama_dh(i)) - (Xsol(i)*Psat(i)*Gama(i)) )/(Xsol(NC+1)*(dh))
        end do

        call MRSL01(Jaco,F,NC+1,NC+1,IER,INDIC)

        Xsol = Xsol - alpha*F

    end do

  ! Pour finaliser le programme
   do i=1,NC
       Xi(i)=Xsol(i)
   enddo
   Trosee = Xsol(NC+1)


end subroutine


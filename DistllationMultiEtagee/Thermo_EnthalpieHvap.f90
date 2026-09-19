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
!   Ficiher : Subroutine calcul de l'enthalpie vapeur
! ----------------------------------------------------------------------------

subroutine EnthalpieHvap(Temp,Yi,Hvap)
    Use Dimensions
    Use Thermos
    implicit none

    integer :: i
    double precision , dimension(NC), intent(in) :: Yi
    double precision, intent(in) :: Temp
    double precision, intent(out) :: Hvap

    Hvap=0

    do i=1,NC
        Hvap=Hvap+Yi(i)*Cp_vap(i)*(Temp-Tref)
    end do

endsubroutine

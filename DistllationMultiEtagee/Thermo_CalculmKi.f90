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
!   Ficiher : Subroutine calcul de la constante d'équilibre
! ----------------------------------------------------------------------------

subroutine CalculmKi(gama,Psat,mKi)
    Use Dimensions
    Use Thermos
    implicit none

    integer :: i
    double precision , dimension(NC), intent(in) :: gama, Psat
    doubleprecision , dimension(NC), intent(out) :: mKi

    do i=1,NC
        mKi(i)=gama(i)*Psat(i)/P
    enddo

endsubroutine

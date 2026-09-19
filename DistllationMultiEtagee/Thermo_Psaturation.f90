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
!   Ficiher : Subroutine calcul de la pression de vapeur saturante
! ----------------------------------------------------------------------------

subroutine Psaturation(Temp,Psat)
    Use Dimensions
    Use Thermos
    implicit none

    integer :: i
    doubleprecision, intent(in) :: Temp
    double precision, dimension(NC), intent(out):: Psat

    do i=1,NC
        Psat(i)=exp(ANT(i,1)+(ANT(i,2)/Temp)+ANT(i,3)*log(Temp)+ANT(i,4)*Temp**ANT(i,5))
    enddo

end subroutine

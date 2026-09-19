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
!   Ficiher : Subroutine Affectation des dimensions
! ----------------------------------------------------------------------------

subroutine AffecteDimensions
    Use Dimensions
    implicit none

    open(unit=10,file=TRIM("DonneesDimensions.txt"),status="old")
        read(10,*)NP
        read(10,*)NC
    close (10)
end subroutine

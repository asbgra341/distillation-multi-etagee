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
!   Ficiher : Module des données thermodynamiques
! ----------------------------------------------------------------------------

module Thermos
    implicit none

    doubleprecision, allocatable, dimension(:,:):: ANT
    doubleprecision, allocatable, dimension(:,:):: Cij, Aij
    doubleprecision, allocatable, dimension(:):: Teb, Cp_vap, Cp_liq, deltaH

    doubleprecision, parameter :: R = 8.3145D0
    doubleprecision, parameter :: P=101325
    doubleprecision, parameter :: Tref = 298.15D0
    doubleprecision, parameter :: Tbullemax=1000D0    ! Valeur max de la temperature de bulle
    doubleprecision, parameter :: Preci= 1E-9        ! Epsillon pour le calcul de Tbulle



end module

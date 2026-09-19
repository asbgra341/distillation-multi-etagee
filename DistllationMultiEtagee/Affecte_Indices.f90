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
!   Ficiher : Subroutine Affectation des indices
! ----------------------------------------------------------------------------

subroutine AffecteIndices
    Use Dimensions
    Use Indices
    implicit none

    IV=1
    IYvap=IV+1
    ITemp=IYvap+NC
    IXliq=ITemp+1
    IL=IXliq+NC
    BMT=1
    BMP=BMT+1
    BE=BMP+NC
    EQ=BE+1
    SOM=EQ+NC

end subroutine

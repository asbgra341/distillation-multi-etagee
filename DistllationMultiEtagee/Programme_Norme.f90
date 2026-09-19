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
!   Ficiher : Subroutine calcul de la norme jacobienne
! ----------------------------------------------------------------------------

subroutine Norme(Res,ERREUR)
    Use Dimensions

    doubleprecision, dimension(NP*(2*NC+3)),intent(in)::Res
    doubleprecision,intent(out):: ERREUR

    ERREUR= norm2(Res)

endsubroutine

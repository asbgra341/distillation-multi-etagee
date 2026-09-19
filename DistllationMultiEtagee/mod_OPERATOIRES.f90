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
!   Ficiher : Module opératoires
! ----------------------------------------------------------------------------

module Operatoires
    implicit none

    double precision :: D, W
    double precision, allocatable, dimension(:):: Falim, Dvap, Dliq, Talim, Qchal, Halim
    double precision, allocatable,dimension(:,:):: Zalim
    double precision, allocatable, dimension(:):: Val_alim, Val_temp, Val_Sliq, Val_Svap
    doubleprecision, allocatable, dimension(:) :: Val_Zalim, Val_Tbulle, Val_gama
    integer, allocatable, dimension(:):: Num_alim, Num_Sliq, Num_Svap
    integer :: Nb_alim, Nb_Sliq, Nb_Svap

end module

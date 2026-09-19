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
!   Ficiher : Subroutine Jacobienne du solveur MRSL01
! ----------------------------------------------------------------------------

subroutine Jacobien(Res,Ysol,Jaco)
    Use Dimensions
    implicit none

    doubleprecision :: deltaX
    doubleprecision, dimension(NP*(2*NC+3)),intent(in):: Res, Ysol
    doubleprecision, dimension(NP*(2*NC+3),NP*(2*NC+3)), intent(out):: Jaco
    doubleprecision, allocatable, dimension(:):: ResMOD, YsolMOD
    integer :: i,j

    allocate(YsolMOD(NP*(2*NC+3)))
    allocate(ResMOD(NP*(2*NC+3)))

    YsolMOD=Ysol
    ResMOD=0D0
    deltaX=0.001D0

    do j=1,NP*(2*NC+3)
        YsolMOD(j)=Ysol(j)*(1+deltaX)
        call Residus(YsolMOD,ResMOD)
        do i=1,NP*(2*NC+3)
            Jaco(i,j)=(ResMOD(i)-Res(i))/(Ysol(j)*deltaX)
        enddo
        YsolMOD(j)=Ysol(j)
    enddo

endsubroutine

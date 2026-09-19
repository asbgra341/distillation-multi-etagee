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
!   Ficiher : Subroutine Jacobienne du solveur MRSL21
! ----------------------------------------------------------------------------

subroutine Jacobien21(Res,Ysol,A,B,C)
    use Dimensions
    use Operatoires
    Use Thermos
    use Indices
    implicit none

    integer ::i, j, l, k
    doubleprecision,dimension(2*NC+3,2*NC+3,NP):: A, B, C
    doubleprecision :: deltaX
    doubleprecision,dimension(NP*(2*NC+3)):: Ysol, Res
    doubleprecision,allocatable,dimension(:):: YsolMOD, ResMOD

    allocate(ResMOD(NP*(2*NC+3)))
    allocate(YsolMOD(NP*(2*NC+3)))

    A=0
    B=0
    C=0
    deltaX=1e-3
    YsolMOD=Ysol
    ResMOD=0D0
    deltaX=0.001D0
    k=0

    do j=1,(2*NC+3)

        k=k+1
        YsolMOD(k)=Ysol(k)*(1+deltaX)
        call Residus(YsolMOD,ResMOD)

        do i=1,(2*NC+3)
            B(i,j,1)=(ResMOD(i)-Res(i))/(Ysol(k)*deltaX)
            A(i,j,2)=(ResMOD(i+2*NC+3)-Res(i+2*NC+3))/(Ysol(k)*deltaX)
        end do
        YsolMOD(k)=Ysol(k)
    end do

    do l=2,NP-1
        do j=1,(2*NC+3)
            k=k+1
            YsolMOD(k)=Ysol(k)*(1+deltaX)
            call Residus(YsolMOD,ResMOD)
            do i=1,(2*NC+3)
                C(i,j,l-1)=(ResMOD((l-2)*(2*NC+3)+i)-Res((l-2)*(2*NC+3)+i))/(Ysol(k)*deltaX)
                B(i,j,l)=(ResMOD((l-1)*(2*NC+3)+i)-Res((l-1)*(2*NC+3)+i))/(Ysol(k)*deltaX)
                A(i,j,l+1)=(ResMOD(l*(2*NC+3)+i)-Res(l*(2*NC+3)+i))/(Ysol(k)*deltaX)
            end do
            YsolMOD(k)=Ysol(k)
        end do
    end do

    do j=1,(2*NC+3)
        k=k+1
        YsolMOD(k)=Ysol(k)*(1+deltaX)
        call Residus(YsolMOD,ResMOD)
        do i=1,(2*NC+3)
            C(i,j,NP-1)=(ResMOD((NP-2)*(2*NC+3)+i)-Res((NP-2)*(2*NC+3)+i))/(Ysol(k)*deltaX)
            B(i,j,NP)=(ResMOD((NP-1)*(2*NC+3)+i)-Res((NP-1)*(2*NC+3)+i))/(Ysol(k)*deltaX)
        end do
        YsolMOD(k)=Ysol(k)
    end do

 end subroutine

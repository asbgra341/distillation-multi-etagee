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
!   Ficiher : Subroutine Programme des résidus
! ----------------------------------------------------------------------------

subroutine Residus(Ysol,Res)
    Use Dimensions
    Use Thermos
    Use Operatoires
    Use Indices
    implicit none

    doubleprecision, dimension(NP*(2*NC+3)),intent(in)::Ysol
    doubleprecision, dimension(NP*(2*NC+3)),intent(out)::Res
    doubleprecision, dimension(NC):: Gama, mKi, Psat
    doubleprecision, dimension(NP):: Hliq, Hvap
    doubleprecision :: somX, somY, SomX1
    doubleprecision :: somY1, somYN, SomXN
    integer :: i,j


!                                       CALCUL DES ENTHALPIES
!===============================================================================================================================
    do i=1,NP
        call EnthalpieHliq(Ysol(ITemp+(i-1)*(2*NC+3)),Ysol(IXliq+(i-1)*(2*NC+3)),Hliq(i))
        call EnthalpieHvap(Ysol(ITemp+(i-1)*(2*NC+3)),Ysol(IYvap+(i-1)*(2*NC+3)),Hvap(i))
    enddo

!                                       POUR LE CONDENSEUR
!===============================================================================================================================
    Res(BMT)=-(Ysol(IL)+Dliq(1))-Ysol(IV)+Ysol(IV+(2*NC+3))+Falim(1)
    do i=1,NC
        Res(BMP+(i-1))=-(Ysol(IL)+Dliq(1))*Ysol(IXliq+(i-1))-Ysol(IV)*Ysol(IYvap+(i-1))+ Ysol(IV+(2*NC+3))&
        &*Ysol(IYvap+(i-1)+(2*NC+3))+Falim(1)*Zalim(1,i)
    enddo
    Res(BE)=-(Ysol(IL)+Dliq(1))*Hliq(1)-Ysol(IV)*Hvap(1)+Ysol(IV+(2*NC+3))*Hvap(2)+Falim(1)*Halim(1)-Qchal(1)

    call CalculGama(Ysol(ITemp),Ysol(IXliq),Gama)
    call Psaturation(Ysol(ITemp),Psat)
    call CalculmKi(Gama,Psat,mKi)
    do i=1,NC
        Res(EQ+(i-1))=Ysol(IYvap+(i-1))-mKi(i)*Ysol(IXliq+(i-1))
    enddo
    SomX1=0
    somY1=0
    do i=1,NC
        somX1=somX1+Ysol(IXliq+(i-1))
        somY1=somY1+Ysol(IYvap+(i-1))
    enddo
    Res(SOM)=somX1-somY1

!                                       POUR LE BOUILLEUR
!===============================================================================================================================
    Res(BMT+(NP-1)*(2*NC+3))=Ysol(IL+(NP-2)*(2*NC+3))-Ysol(IL+(NP-1)*(2*NC+3))-(Ysol(IV+(NP-1)*(2*NC+3))+Dvap(NP))+Falim(NP)
    do i=1,NC
        Res(BMP+(NP-1)*(2*NC+3)+(i-1))=Ysol(IL+(NP-2)*(2*NC+3))*Ysol(IXliq+(NP-2)*(2*NC+3)+(i-1))&
        &-Ysol(IL+(NP-1)*(2*NC+3))*Ysol(IXliq+(NP-1)*(2*NC+3)+(i-1))-(Ysol(IV+(NP-1)*(2*NC+3))+Dvap(NP))&
        &*Ysol(IYvap+(NP-1)*(2*NC+3)+(i-1))+Falim(NP)*Zalim(NP,i)
    enddo
    Res(BE+(NP-1)*(2*NC+3))=Ysol(IL+(NP-2)*(2*NC+3))*Hliq(NP-1)-Ysol(IL+(NP-1)*(2*NC+3))*Hliq(NP)&
    &-(Ysol(IV+(NP-1)*(2*NC+3))+Dvap(NP))*Hvap(NP)+Falim(NP)*Halim(NP)+Qchal(NP)

    call CalculGama(Ysol(ITemp+(NP-1)*(2*NC+3)),Ysol(IXliq+(NP-1)*(2*NC+3)),Gama)
    call Psaturation(Ysol(ITemp+(NP-1)*(2*NC+3)),Psat)
    call CalculmKi(Gama,Psat,mKi)

    do i=1,NC
        Res(EQ+(NP-1)*(2*NC+3)+(i-1))=Ysol(IYvap+(NP-1)*(2*NC+3)+(i-1))-mKi(i)*Ysol(IXliq+(NP-1)*(2*NC+3)+(i-1))
    enddo
    SomXN=0
    somYN=0
    do i=1,NC
        somXN=somXN+Ysol(IXliq+(NP-1)*(2*NC+3)+(i-1))
        somYN=somYN+Ysol(IYvap+(NP-1)*(2*NC+3)+(i-1))
    enddo
    Res(SOM+(NP-1)*(2*NC+3))=SomXN-somYN

!                                       POUR LES PLATEAUX INTERMEDIAIRES
!===============================================================================================================================

    do i=2,NP-1
        Res(BMT+(i-1)*(2*NC+3))=Ysol(IL+(i-2)*(2*NC+3))-(Ysol(IL+(i-1)*(2*NC+3))+Dliq(i))-(Ysol(IV+(i-1)*(2*Nc+3))+Dvap(i))+&
        &Ysol(IV+(i)*(2*NC+3))+Falim(i)         ! retour à la ligne

        do j=1,NC
            Res(BMP+(i-1)*(2*NC+3)+(j-1))=Ysol(IL+(i-2)*(2*NC+3))*Ysol(IXliq+(i-2)*(2*NC+3)+(j-1))-&
            &(Ysol(IL+(i-1)*(2*NC+3))+Dliq(i))*Ysol(IXliq+(i-1)*(2*NC+3)+(j-1))-(Ysol(IV+(i-1)*(2*NC+3))+Dvap(i))&
            &*Ysol(IYvap+(i-1)*(2*NC+3)+(j-1))+Ysol(IV+(i)*(2*NC+3))*Ysol(IYvap+(i)*(2*NC+3)+(j-1))+Falim(i)*Zalim(i,j)
        enddo

        Res(BE+(i-1)*(2*NC+3))=Ysol(IL+(i-2)*(2*NC+3))*Hliq(i-1)-(Ysol(IL+(i-1)*(2*NC+3))+Dliq(i))*Hliq(i)&
        &-(Ysol(IV+(i-1)*(2*Nc+3))+Dvap(i))*Hvap(i)+Ysol(IV+(i)*(2*NC+3))*Hvap(i+1)+ Falim(i)*Halim(i)-Qchal(i)

        call CalculGama(Ysol(ITemp+(i-1)*(2*NC+3)),Ysol(IXliq+(i-1)*(2*NC+3)),Gama)
        call Psaturation(Ysol(ITemp+(i-1)*(2*NC+3)),Psat)
        call CalculmKi(Gama,Psat,mKi)

        do j=1,NC
            Res(EQ+(i-1)*(2*NC+3)+(j-1))=Ysol(IYvap+(i-1)*(2*NC+3)+(j-1))-mKi(j)*Ysol(IXliq+(i-1)*(2*NC+3)+(j-1))
        enddo

        somX=0
        somY=0

        do j=1,NC
            somX=somX+Ysol(IXliq+(i-1)*(2*NC+3)+(j-1))
            somY=somY+Ysol(IYvap+(i-1)*(2*NC+3)+(j-1))
        enddo
        Res(SOM+(i-1)*(2*NC+3))=somX-somY

    enddo

endsubroutine

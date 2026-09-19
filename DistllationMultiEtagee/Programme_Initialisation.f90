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
!   Ficiher : Subroutine Initialisation du programme
! ----------------------------------------------------------------------------

subroutine Initialise(Ysol)
    Use Dimensions
    Use Thermos
    Use Operatoires
    Use Indices
    implicit none

    double precision, dimension (NP*(2*NC+3)),intent(out)::Ysol
    double precision, dimension (NC):: Xliqbulle, Yvapbulle,Xliqrosee,Yvaprosee
    double precision, dimension (NC):: Xliqbulle1, Yvapbulle1       ! Utiliser dans le cas ou il y'a plusieurs alimentations
    double precision, dimension (NC):: gama, mKi,Psat
    double precision :: Tbulle, Trosee,Temp
    double precision :: Hvap,Hliq,Sum_Falim
    integer::i,j


   if (Nb_alim==1) then
!                            POUR UN SEUL PLATEAU ALIMENTE
!=======================================================================================
        Tbulle=Talim(Num_alim(1))
        do i=1,NC
             Xliqbulle(i)=(Zalim(Num_alim(1),i))
        end do
        call CalculGama(Tbulle,Xliqbulle,gama)
        call Psaturation(Tbulle,Psat)
        call CalculmKi(gama,Psat,mKi)
        do i=1,NC
            Yvapbulle(i)=Xliqbulle(i)*mKi(i)
            Yvaprosee(i)=Xliqbulle(i)
        enddo
        call CalculTrosee(Yvaprosee,Xliqrosee,Trosee)

    else
!                            POUR PLUSIEURS PLATEAUX ALIMENTES
!=======================================================================================

! -------------------------- Pour le Condenseur -----------------------------------------
        Tbulle=Talim(Num_alim(1))
        do i=1,NC
             Xliqbulle(i)=(Zalim(Num_alim(1),i))
        end do
        call CalculGama(Tbulle,Xliqbulle,gama)
        call Psaturation(Tbulle,Psat)
        call CalculmKi(gama,Psat,mKi)
        do i=1,NC
            Yvapbulle(i)=Xliqbulle(i)*mKi(i)
        enddo
! ------------------------- Pour le Bouilleur -------------------------------------------
        Temp=Talim(Num_alim(Nb_alim))
        do i=1,NC
             Xliqbulle1(i)=(Zalim(Num_alim(Nb_alim),i))
        end do
        call CalculGama(Temp,Xliqbulle1,gama)
        call Psaturation(Temp,Psat)
        call CalculmKi(gama,Psat,mKi)
        do i=1,NC
            Yvapbulle1(i)=Xliqbulle1(i)*mKi(i)
            Yvaprosee(i)=Xliqbulle1(i)
        enddo
        call CalculTrosee(Yvaprosee,Xliqrosee,Trosee)
    endif

!                            INITILISATION DU PROGRAMME
!=======================================================================================

!--------------------------Pour les Températures ------------------------------
    Ysol(Itemp)=Tbulle
    Ysol(Itemp + (NP-1)*(2*NC+3))=Trosee
    do i=2,NP-1
        Ysol(Itemp + (i-1)*(2*NC+3))=Ysol(Itemp + (i-2)*(2*NC+3))+((Trosee-Tbulle)/(NP-1))
    end do
!----------------------------Pour les Compositions -----------------------

    do j=1,NC
        Ysol(IXliq+j-1)=Xliqbulle(j)
        Ysol(IYvap+j-1)=Yvapbulle(j)
        Ysol(IXliq+ (NP-1)*(2*NC+3)+j-1)=Xliqrosee(j)
        Ysol(IYvap+ (NP-1)*(2*NC+3)+j-1)=Yvaprosee(j)
    end do

    do i=2,NP-1
        do j=1,NC
            Ysol(IXliq+(i-1)*(2*NC+3)+(j-1))=Ysol(IXliq+(i-2)*(2*NC+3)+(j-1))+((Xliqrosee(j)-Xliqbulle(j))/(NP-1))
            Ysol(IYvap+(i-1)*(2*NC+3)+(j-1))=Ysol(IYvap+(i-2)*(2*NC+3)+(j-1))+((Yvaprosee(j)-Yvapbulle(j))/(NP-1))
        end do
    end do

!-------------------------- Pour les Débits ---------------------------------
    Sum_Falim=0
    do i=1,NP
        Sum_Falim=Sum_Falim+Falim(i)*Zalim(i,1)
    end do
    W=(Sum_Falim-sum(Falim)*Yvapbulle(1))/(Xliqrosee(1)-Yvapbulle(1))

    D=sum(Falim)-W

    call EnthalpieHliq(Tbulle,Xliqbulle,Hliq)
    call EnthalpieHvap(Tbulle,Yvapbulle,Hvap)

    Ysol(IL)=(Qchal(1)/(Hvap-Hliq))-D
    do i=2,NP
        Ysol(IL+(i-1)*(2*NC+3))=Ysol(IL+(i-2)*(2*NC+3))+Falim(i)
    end do
    Ysol(IV)=D
    do i=2,NP
        Ysol(IV+(i-1)*(2*NC+3))=Ysol(IL)+D
    end do

endsubroutine

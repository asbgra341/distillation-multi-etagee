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
!   Ficiher : Subroutine Affectation des données opératoires
! ----------------------------------------------------------------------------

subroutine AffecteOperatoire
    Use Dimensions
    Use Thermos
    Use Operatoires
    implicit none

    integer ::i
    integer ::j

    allocate(Talim(NP))
    allocate(Halim(NP))
    allocate(Falim(NP))
    allocate(Zalim(NP,NC))
    allocate(Qchal(NP))
    allocate(Dvap(NP))
    allocate(Dliq(NP))

    Zalim=0
    Dvap=0
    Dliq=0
    Falim=0
    Talim=0
    Halim=0

    open(11, file="DonneesOperatoires.txt", status="old")
        read(11,*)Nb_alim
        allocate (Num_alim(Nb_alim))
        allocate (Val_alim(Nb_alim))
        allocate (Val_temp(Nb_alim))
        allocate (Val_Tbulle(Nb_alim))
        allocate (Val_Zalim(NC))
        allocate (Val_gama(NC))

        Num_alim=0
        Val_alim=0
        val_temp=0

        read(11,*)Num_alim
        read(11,*)Val_alim
        read(11,*)Val_temp
        do i=1,NC
            read(11,*)(Zalim(Num_alim(j),i),j=1,Nb_alim)
        end do
        Qchal=0
        read(11,*)Qchal(NP)
        read(11,*)Qchal(1)
        read(11,*)Nb_Sliq
        allocate(Num_Sliq(Nb_Sliq))
        allocate(Val_Sliq(Nb_Sliq))
        Num_Sliq=0
        Val_Sliq=0
        read(11,*)Num_Sliq
        read(11,*)Val_Sliq
        read(11,*)Nb_Svap
        allocate(Num_Svap(Nb_Svap))
        allocate(Val_Svap(Nb_Svap))
        Num_Svap=0
        Val_Svap=0
        read(11,*)Num_Svap
        read(11,*)Val_Svap

    close(11)


    do i=1,Nb_alim
        Falim(Num_alim(i))=Val_alim(i)
    end do

    do i=1,Nb_Sliq
        Dliq(Num_Sliq(i))=Val_Sliq(i)
    end do
    do i=1,Nb_Svap
        Dvap(Num_Svap(i))=Val_Svap(i)
    end do

!                        POUR IMPOSER LES TEMPERATURES D'ALIMENTATIONS
!=======================================================================================

! ----------------------- Les plateaux sont alimentés à Tbulle --------------------------
    do i=1,Nb_alim
        do j=1,NC
            Val_Zalim(j)=Zalim(Num_alim(i),j)
        end do
        call CalculGama(Val_temp(i),Val_Zalim,Val_gama)
        call CalculTbulle(Val_gama,Val_Zalim,Val_Tbulle(i))
        Talim(Num_alim(i))=Val_Tbulle(i)
        call EnthalpieHliq(Talim(Num_alim(i)),Val_Zalim,Halim(Num_alim(i)))
    enddo




end subroutine

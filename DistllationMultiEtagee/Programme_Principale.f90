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
!   Ficiher : Subroutine Programme principale
! ----------------------------------------------------------------------------

program DISTILLATION

    Use Dimensions
    Use Thermos
    Use Operatoires
    Use Indices
    Use Messages
    implicit none

    integer :: i,j, k
    integer :: MODE,MRSL, err
    doubleprecision,allocatable, dimension(:) :: Ysol, Res
    doubleprecision,allocatable, dimension(:,:) :: Jaco
    doubleprecision,allocatable, dimension(:,:,:) :: A, B, C
    doubleprecision :: alpha, ERREUR
    doubleprecision :: IER, INDIC, IRAZ, IMET
    doubleprecision :: temps_debut , temps_fin
    common MODE

    call messagedacceuil

    print*," Etape 1. Choix du comportement decrit par le modele Thermodynamique NRTL "
    print*
     33 print*,"tapez 1 pour une solution ideale et tapez 2 pour une solution non ideale"
    print*
    read(*,*,iostat=err)MODE
    call MessageErreur(err)

    if ((MODE/=1).AND.(MODE/=2)) then
        goto 33
    endif
    print*

    call AffecteDimensions
    call AffecteTHERMO
    call AffecteOperatoire
    call AffecteIndices

    allocate(Ysol(NP*(2*NC+3)))
    allocate(Res(NP*(2*NC+3)))
    allocate(Jaco(NP*(2*NC+3),NP*(2*NC+3)))
    allocate(A(2*NC+3,2*NC+3,NP))
    allocate(B(2*NC+3,2*NC+3,NP))
    allocate(C(2*NC+3,2*NC+3,NP))

    Ysol=0
    Res=0
    Jaco=0
    k=0
    IER=0
    INDIC=0
    IRAZ=1
    IMET=0
    alpha=0.1
    ERREUR=10D0

    call Initialise(Ysol)

    print*
    print*," Etape 2. Choix du Solveur  "
    print*
     22 print*,"tapez 1 pour le Solveur MRSL01 , tapez 2 pour le Solveur MRSL21"
    print*
    read(*,*,iostat=err)MRSL
    call MessageErreur(err)

    if ((MRSL/=1).AND.(MRSL/=2)) then
        goto 22
    endif
    print*
    print*
    print*," Etape 3. Debut du iterations "
    print*
    if (MRSL==1) then
        call cpu_time(temps_debut)

        do while (ERREUR.GT.1E-6)
            k=k+1
            call Residus(Ysol,Res)
            call Norme(Res,ERREUR)
            call Jacobien(Res,Ysol,Jaco)
            call MRSL01(Jaco,Res,NP*(2*NC+3),NP*(2*NC+3),IER,INDIC)
            Ysol=Ysol-alpha*Res
            print*, "Iterations : ",k, "Norme : ",ERREUR,"IER : ",IER," SIMULATION EN COURS "
        enddo
        open (unit=20, file="ResultatsResMRSL01.csv", status="unknown")
        do i=1,NP
            write(20,*)i,';',(Res((i-1)*(2*NC+3)+j),j=1,2*NC+3)
        enddo
        close(20)

        open (unit=20, file="ResultatsYsolMRSL01.csv", status="unknown")
        do i=1,NP
            write(20,*)i,';',(Ysol((i-1)*(2*NC+3)+j),j=1,2*NC+3)
        enddo
        close(20)

        call cpu_time(temps_fin)
    else
        call cpu_time(temps_debut)

        do while (ERREUR.GT.1E-6)
            k=k+1
            call Residus(Ysol,Res)
            call Norme(Res,ERREUR)
            call Jacobien21(Res,Ysol,A,B,C)
            call MRSL21(A,B,C,Res,2*NC+3,NP,1,2*NC+3,2*NC+3,IER,IRAZ,IMET)
            Ysol=Ysol-alpha*Res
            print*, "Iterations : ",k, "Norme : ",ERREUR,"IER : ",IER," SIMULATION EN COURS "
        enddo
        open (unit=20, file="ResultatsResMRSL21.csv", status="unknown")
        do i=1,NP
            write(20,*)i,';',(Res((i-1)*(2*NC+3)+j),j=1,2*NC+3)
        enddo
        close(20)

        open (unit=20, file="ResultatsYsolMRSL21.csv", status="unknown")
        do i=1,NP
            write(20,*)i,';',(Ysol((i-1)*(2*NC+3)+j),j=1,2*NC+3)
        enddo
        close(20)

        call cpu_time(temps_fin)
    endif
    print*
    print*
    print*," Etape 4. Messages du programme "
    print*
    print*,"Temps de simulation du solveur en seconde : ", temps_fin-temps_debut
    call messageResultat
    call messagedefin

end program



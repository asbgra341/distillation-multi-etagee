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
!   Ficiher : Subroutine Affectation des données thermodynamiques
! ----------------------------------------------------------------------------

subroutine AffecteTHERMO
    Use Dimensions
    Use Thermos
    implicit none

    integer :: i
    integer :: j

    allocate(ANT(NC,5))
    allocate(Cij(NC,NC))
    allocate(Aij(NC,NC))
    allocate(deltaH(NC))
    allocate(Teb(NC))
    allocate(Cp_liq(NC))
    allocate(Cp_vap(NC))

    ! Lecture des donnees du coefficient
    open(unit=11,file=TRIM("DonneesCoefAntoine.txt"),status="old")
        do i=1,NC
            read(11,*)(ANT(i,j),j=1,5)
        end do
    close (11)
    ! Lecture des donneesThermo
    open(unit=12,file=TRIM("DonneesThermos.txt"),status="old")
        read(12,*)(Teb(i),i=1,NC)
        read(12,*)(deltaH(i),i=1,NC)
        read(12,*)(Cp_vap(i),i=1,NC)
        read(12,*)(Cp_liq(i),i=1,NC)
    close (12)

   ! Lecture des données d'interaction Binaires
    open(unit=14,file=TRIM("DonneesCoefInteraction.txt"),status="old")
        do i=1,NC
            read(14,*)(Aij(i,j), j=1,NC)
        end do
       do i=1,NC
           read(14,*)(Cij(i,j), j=1,NC)
       end do
    close (14)

end subroutine

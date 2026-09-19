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
!   Ficiher : Subroutine calcul de la temperature de bulle
! ----------------------------------------------------------------------------

subroutine CalculTbulle(gama, Xi,Tebullition)
    Use Dimensions
    Use Thermos
    implicit none

    integer :: i
    doubleprecision :: Tmin, Tmax,Tmoyen, somme_min, somme_max
    doubleprecision , intent(out) :: Tebullition
    double precision , dimension(NC):: Psat_min, Psat_max, Psat_moy
    double precision , dimension(NC),intent(in) :: gama, Xi

    somme_min =0.
    somme_max =0.
    Tmin=Tref
    Tmax=Tbullemax

    do while((Tmax - Tmin) > preci)

        Tmoyen = (Tmin + Tmax )/2.
        call Psaturation(Tmin,Psat_min)
        call Psaturation(Tmax,Psat_max)

        do i=1,NC
            somme_min = somme_min + Xi(i)*gama(i)*Psat_min(i)
            somme_max = somme_max + Xi(i)*gama(i)*Psat_max(i)
        enddo

        if ( (somme_min - P) * (somme_max - P) < 0) then
            somme_max = 0.
            call Psaturation(Tmoyen,Psat_moy)
            do i = 1, nc
                somme_max = somme_max + Xi(i)*gama(i)*Psat_moy(i)
            enddo

            if( (somme_min - P) * (somme_max - P) < 0 ) then
                Tmax = Tmoyen
            else
                Tmin = Tmoyen
            endif

        else

            write(*,*)"Pas de solution dans l'intervalle de T donnee"
            stop

        endif

        somme_min = 0.
        somme_max = 0.

    enddo
    Tebullition=Tmoyen

endsubroutine

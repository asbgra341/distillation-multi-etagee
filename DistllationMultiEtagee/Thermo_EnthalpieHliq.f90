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
!   Ficiher : Subroutine calcul de l'enthalpie liquide
! ----------------------------------------------------------------------------

subroutine EnthalpieHliq(Temp,Xi,Hliq)
    Use Dimensions
    Use Thermos
    implicit none

    integer :: i, MODE
    doubleprecision , dimension(NC), intent(in) :: Xi
    doubleprecision, intent(in) :: Temp
    doubleprecision, intent(out) :: Hliq
    doubleprecision :: Hliq_exces, deltaT, Temp2, deltaLog
    doubleprecision , dimension(NC) :: gama, gama2
    common MODE

    Hliq=0
    Hliq_exces=0
    deltaT=0.001
    deltaLog=0

    if (MODE.eq.1) then ! solution idéale
        do i=1,NC
            Hliq=Hliq+Xi(i)*(Cp_vap(i)*(Teb(i)-Tref)- deltaH(i) + Cp_liq(i)*(Temp-Teb(i)))
        end do

    elseif (MODE.eq.2) then ! solution non idéale

        Temp2=Temp+deltaT
        call CalculGama(Temp,Xi,gama)
        call CalculGama(Temp2,Xi,gama2)

        do i=1,NC
            Hliq=Hliq+Xi(i)*(Cp_vap(i)*(Teb(i)-Tref)- deltaH(i) + Cp_liq(i)*(Temp-Teb(i)))-&
                                    & R*(Temp**2)*Xi(i)*(log(gama2(i))-log(gama(i)))/deltaT
        end do

    endif

endsubroutine

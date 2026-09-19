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
!   Ficiher : Subroutine calcul de gama
! ----------------------------------------------------------------------------

subroutine CalculGama(Temp,Xi,gama)
    Use Dimensions
    Use Thermos
    implicit none

    integer :: MODE
    integer :: i, j , l
    doubleprecision, intent(in) :: Temp
    doubleprecision , dimension(NC,NC) :: Toij
    doubleprecision , dimension(NC,NC) :: Gij
    doubleprecision , dimension(NC), intent(in) :: Xi
    doubleprecision , dimension(NC),intent(out) :: gama
    double precision :: somme1, somme2, somme3, somme4, somme6
    common MODE

    if (MODE.eq.1) then
        gama=1
    elseif (MODE.eq.2) then
!=========================================================================
        do i=1,NC
            do j=1,NC
                if(i==j) then
                    Toij(i,j)=0
                else
                    Toij(i,j)=(1/(R*Temp))*Cij(i,j)
                endif
            enddo
        enddo

        do i=1,NC
            do j=1,NC
                if(i==j) then
                    Gij(i,j)=1
                else
                    Gij(i,j)=exp(-(Toij(i,j)*Aij(i,j)))
                endif
            enddo
        enddo
!=========================================================================
        gama=0
        do i=1,NC
            somme1=0
            somme2=0
            somme6=0
            do l=1,NC
                somme1=somme1 + Toij(l,i)*Gij(l,i)*Xi(l)
                somme2=somme2 + Gij(l,i)*Xi(l)
            end do

            do j=1,NC
                somme3=0
                somme4=0
                do l=1,NC
                    somme3=somme3 + Toij(l,j)*Gij(l,j)*Xi(l)
                    somme4=somme4 + Gij(l,j)*Xi(l)
                end do
                    somme6=somme6 + ((Xi(j)*Gij(i,j))/somme4)*(toij(i,j)-(somme3/somme4))
            end do
            gama(i)=exp((somme1/somme2)+somme6)
        enddo

    endif

endsubroutine


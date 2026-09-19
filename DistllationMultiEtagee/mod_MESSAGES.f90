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
!   Ficiher : Module messages
! ----------------------------------------------------------------------------

module Messages

contains

subroutine messagedacceuil()

        implicit none
        print*,'**************************************************************************************************************'
        print*,'**************************************************************************************************************'
        print*,'                                               DEBUT DU PROGRAMME'
        print*,'**************************************************************************************************************'
        print*
        print*,"                                Modelisation dune Colonne de Distillation Biphasique"
        print*,""
        print*
        print*,"Bienvenue dans votre programme de modelisation d'une colonne de distillation biphasique a plusieurs plateaux"
        print*,"theoriques en regime permanent avec un condenseur partiel."
        print*
        print*

    endsubroutine messagedacceuil

    subroutine messagedefin()

        implicit none
        print*
        print*,'**************************************************************************************************************'
        print*,'                                                   FIN DU PROGRAMME'
        print*,'**************************************************************************************************************'
        print*,'**************************************************************************************************************'
        print*

    endsubroutine messagedefin

    subroutine messageResultat()

        implicit none
        print*
        print*,'Fichiers ResultatsYsol et ResultatsRes creer avec succes !!!!!        '
        print*
        print*,'Pour lexploitation du fichier : ResultatsYsol (Identification des colonnes)'
        print*
        print*,'Numero Plateau ; Debit Vapeur ; Composition Vapeur (NC) ; Temperature ; Debit Liquide ; Composition Liquide (NC)'
        print*
        print*,'Pour lexploitation du fichier : ResultatsRes (Identification des colonnes)'
        print*
        print*,'Numero Plateau ; B.M Total ; B.M.P (NC) ; B. Energetique ; Equilibre L/V (NC) ; Sommation'


    endsubroutine messageResultat


    subroutine MessageErreur(err)

        implicit none
        integer :: err

        if(err /= 0) THEN
            print*
            print*," PROGRAM ERROR : Le nombre saisi est incoherent "
            print*
            print*,' '
            call messagedefin()
            stop
        endif
    endsubroutine


endmodule

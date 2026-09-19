subroutine mKi(gama,Psat,Ki)
    Use Dimensions
    implicit none
    integer :: i
    double precision , dimension(NC), intent(in) :: gama, Psat
    doubleprecision , dimension(NC), intent(out) :: Ki

    do i=1,NC
        Ki(i)=(gama(i)*Psat(i))/P
    enddo
    Print*,"Valeur de Ki"
    print*, Ki

endsubroutine

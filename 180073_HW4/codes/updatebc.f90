! AMAN PAREKH - 180073 - ME630 - Monsoon 2021

subroutine updatebc(n, phi)

        implicit none

        integer :: n, i, j
        real*16 :: phi(n,n)

        do j = 2,n-1
            phi(1,j) = phi(n-1,j)
            phi(n,j) = phi(2,j)
            phi(j,1) = phi(j,n-1)
            phi(j,n) = phi(j,2)
        end do

end subroutine updatebc
 

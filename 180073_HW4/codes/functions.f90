! AMAN PAREKH - 180073 - ME630 - Monsoon 2021

real*16 function residual(n, phi, phi_old)

        implicit none
        
        integer :: n, i, j
        real*16 :: phi(n,n), phi_old(n,n), s

        s = 0
        do j = 2,n-1
          do i = 2,n-1
            s = s + (phi(i,j) - phi_old(i,j))**2
          end do
        end do

        residual = sqrt(s)

end function residual



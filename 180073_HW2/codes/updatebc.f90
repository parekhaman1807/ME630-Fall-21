! AMAN PAREKH - 180073 - Fall 2021 - ME630

subroutine updatebc(n, phi, mesh)

        implicit none

        integer :: n, i, j
        real*8 :: x, y
        real*8 :: phi(n,n), mesh(4,n*n)

        ! Left Boundary
        do j = 2,n-1
            y = mesh(4,(j-1)*n + 1)
            phi(1,j) = 2.0*(0.25*sinh(-5.0) + ((y-0.5)**2)*sinh(10.0*(y-0.5)) + 1) - phi(2,j)
        end do

        ! Right Boundary
        do j = 2,n-1
            y = mesh(4,(j)*n - 1)
            phi(n,j) = 2.0*(0.25*sinh(5.0) + ((y-0.5)**2)*sinh(10.0*(y-0.5)) + exp(2.0*y)) - phi(n-1,j)
        end do
        
        ! Bottom Boundary
        do i = 2,n-1           
            x = mesh(3,i)
            phi(i,1) = 2.0*(0.25*sinh(-5.0) + ((x-0.5)**2)*sinh(10.0*(x-0.5)) + 1) - phi(i,2)
        end do

        ! Top Boundary
        do i = 2,n-1
            x = mesh(3,(n-1)*n + i)
            phi(i,n) = 2.0*(0.25*sinh(5.0) + ((x-0.5)**2)*sinh(10.0*(x-0.5)) + exp(2.0*x)) - phi(i,n-1)
        end do

end subroutine updatebc
 

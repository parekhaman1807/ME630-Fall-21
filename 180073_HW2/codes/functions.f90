! AMAN PAREKH - 180073 - Fall 2021 - ME630

real*8 function g(x,y)         ! Function for Source Term

        implicit none

        real*8 :: x,y

        g = 2.0*sinh(10.0*(x-0.5)) + 40.0*(x-0.5)*cosh(10.0*(x-0.5))  &
           + 100.0*((x-0.5)**2)*sinh(10.0*(x-0.5))                 & 
           + 2.0*sinh(10.0*(y-0.5)) + 40.0*(y-0.5)*cosh(10.0*(y-0.5)) &
           + 100.0*((y-0.5)**2)*sinh(10.0*(y-0.5))                 &
           + 4.0*(x**2 +y**2)*exp(2.0*x*y)

end function g

real*8 function analytical(x,y)        ! Calculates Analytical Solution

        implicit none

        real*8 :: x,y

        analytical = ((x-0.5)**2)*sinh(10.0*(x-0.5)) + ((y-0.5)**2)*sinh(10.0*(y-0.5)) + exp(2*x*y)

end function analytical         

real*8 function residual(n, phi, phi_old)

        implicit none
        
        integer :: n, i, j
        real*8 :: phi(n,n), phi_old(n,n), s

        s = 0
        do j = 2,n-1
          do i = 2,n-1
            s = s + (phi(i,j) - phi_old(i,j))**2
          end do
        end do

        residual = sqrt(s)

end function residual



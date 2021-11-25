real*16 function g(x,y)         ! Function for Source Term

        implicit none

        real :: x,y

        g = 2.0*sinh(10.0*(x-0.5)) + 40.0*(x-0.5)*cosh(10.0*(x-0.5))  &
           + 100.0*((x-0.5)**2)*sinh(10.0*(x-0.5))                 & 
           + 2.0*sinh(10.0*(y-0.5)) + 40.0*(y-0.5)*cosh(10.0*(y-0.5)) &
           + 100.0*((y-0.5)**2)*sinh(10.0*(y-0.5))                 &
           + 4.0*(x**2 +y**2)*exp(2.0*x*y)

end function g

real*16 function analytical(x,y)        ! Calculates Analytical Solution

        implicit none

        real :: x,y

        analytical = ((x-0.5)**2)*sinh(10.0*(x-0.5)) + ((y-0.5)**2)*sinh(10.0*(y-0.5)) + exp(2*x*y)

end function analytical         

real*16 function residual(n, phi, phi_old)

        implicit none
        
        integer :: n, i, j
        real :: phi(n,n), phi_old(n,n), s

        s = 0
        do j = 1,n
          do i = 1,n
            s = s + (phi(i,j) - phi_old(i,j))**2
          end do
        end do

        residual = sqrt(s)

end function residual

subroutine meshgen(n, mesh)     ! Subroutine for Mesh Generation

      implicit none
      integer :: n, i, j, c  
      real :: mesh(4,n*n)
      real :: edge1, edge2, total_length, del

      edge1 = -0.5/(n-2.0)
      edge2 = 1.0 + (0.5/(n-2.0))
      total_length = edge2 - edge1
      del = (total_length*1.0)/(n-1.0)

      ! mesh = [x_index y_index x_loc y_loc]

      c = 1
      do j = 1,n
        do i = 1,n
          mesh(1,c) = i
          mesh(2,c) = j
          mesh(3,c) = edge1 + ((i-1)*1.0)*del
          mesh(4,c) = edge1 + ((j-1)*1.0)*del 
          c = c + 1
        end do  
      end do

end subroutine meshgen    

subroutine updatebc(n, phi, mesh)

        implicit none

        integer :: n, i, j
        real :: x, y
        real :: phi(n,n), mesh(4,n*n)

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
 
subroutine jacobi(mesh, phi, n, tol)

        implicit none
        real :: mesh(4,n*n), delx, dely, edge1, edge2, total_length, old_residual
        real :: phi(n,n), phi_old(n,n)
        integer :: n, niter, i, j
        real*16 :: g, residual, tol
        open(unit=11, file = 'convergence.dat', position = 'append')

        edge1 = -0.5/(n-2.0)
        edge2 = 1.0 + (0.5/(n-2.0))
        total_length = edge2 - edge1
        delx = (total_length*1.0)/(n-1.0)
        dely = (total_length*1.0)/(n-1.0)


        phi = 0.0
        phi_old = 1.0
        niter = 0

        call updatebc(n, phi, mesh)

        do while(abs(residual(n,phi,phi_old) - old_residual).gt.tol)
            old_residual = residual(n,phi,phi_old)
            phi_old = phi
            do j = 2,n-1
              do i = 2,n-1
                phi(i,j) = (g(mesh(3,((j-1)*n)+i),mesh(4,((j-1)*n)+i)) &
                            - (phi_old(i,j-1)/(dely**2)) &
                            - (phi_old(i-1,j)/(delx**2)) & 
                            - (phi_old(i+1,j)/(delx**2)) & 
                            - (phi_old(i,j+1)/(dely**2)))/(-2.0*((1.0/(delx**2)) + (1.0/(dely**2))))
              end do
            end do

            niter = niter + 1

            print *, "Iteration Count = ", niter
            print *, "Residual = ", old_residual
            write(11,*) niter, old_residual 

            call updatebc(n, phi, mesh)
        end do

        print *,"Number of Iterations= ", niter
        close(11)

end subroutine jacobi        

program main

      implicit none

      real*16 , external :: g, residual, analytical
      real*16 :: tol
      real :: mesh(4,43*43), phi(43,43), phi_a(43,43), error(43,43)
      integer :: n,i,j,c
      open(unit=10, file = 'numerical.dat')
      open(unit=12, file = 'error.dat')
      open(unit=13, file='analytical.dat')
      n = 43
      tol = 1e-6

      call meshgen(n,mesh)

      do j = 1,n
         do i = 1,n
           phi_a(i,j) = analytical(mesh(3,(j-1)*n+i),mesh(4,(j-1)*n+i))
         end do
      end do

      ! Jacobi      
      call jacobi(mesh, phi, n, tol)
      error = phi - phi_a

      do j = n,1,-1
          write (12,*) ((error(i,j)), i=1,n)
      end do

      do j = n,1,-1
          write (10,*)  ((phi(i,j)), i=1,n)
      end do

      do j = n,1,-1
          write (13,*) ((phi_a(i,j)), i=1,n)
      end do

      close(10)
      close(12)
      close(13)

end program main      

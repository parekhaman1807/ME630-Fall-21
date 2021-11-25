! AMAN PAREKH - 180073 - Fall 2021 - ME630

subroutine jacobi(mesh, phi, n, tol)

        implicit none
        real*8 :: mesh(4,n*n), delx, dely, edge1, edge2, total_length, old_residual
        real*8 :: phi(n,n), phi_old(n,n)
        integer :: n, niter, i, j
        real*8 :: g, residual, tol
        open(unit=11, file = 'convergence.dat', position = 'append')

        edge1 = -0.5/(n-2.0)
        edge2 = 1.0 + (0.5/(n-2.0))
        total_length = edge2 - edge1
        delx = (total_length*1.0)/(n-1.0)
        dely = (total_length*1.0)/(n-1.0)


        phi = 0.0
        phi_old = 1.0
        niter = 0
 
        do while(abs(residual(n,phi,phi_old)).gt.tol)
            old_residual = residual(n,phi,phi_old)
            phi_old = phi

            call updatebc(n, phi, mesh)

            ! Jacobi Loop starts
            do j = 2,n-1
              do i = 2,n-1
                phi(i,j) = (g(mesh(3,((j-1)*n)+i),mesh(4,((j-1)*n)+i)) &
                            - (phi_old(i,j-1)/(dely**2)) &
                            - (phi_old(i-1,j)/(delx**2)) & 
                            - (phi_old(i+1,j)/(delx**2)) & 
                            - (phi_old(i,j+1)/(dely**2)))/(-2.0*((1.0/(delx**2)) + (1.0/(dely**2))))
              end do
            end do
            ! Jacobi Loop ends

            niter = niter + 1

            print *, "Iteration Count = ", niter
            print *, "Residual = ", old_residual
            write(11,*) niter, old_residual 

        end do

        print *,"Number of Iterations= ", niter
        close(11)

end subroutine jacobi        

program main

      implicit none

      real*8 , external :: g, residual, analytical
      real*8 :: tol
      real*8 :: mesh(4,43*43), phi(43,43), phi_a(43,43), error(43,43)
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

      do j = n-1,2,-1
          write (12,*) ((error(i,j)), i=2,n-1)
      end do

      do j = n-1,2,-1
          write (10,*)  ((phi(i,j)), i=2,n-1)
      end do

      do j = n-1,2,-1
          write (13,*) ((phi_a(i,j)), i=2,n-1)
      end do


      close(10)
      close(12)
      close(13)

end program main      

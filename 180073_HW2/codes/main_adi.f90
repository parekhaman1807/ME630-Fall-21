! AMAN PAREKH - 180073 - Fall 2021 - ME630

subroutine adi(mesh, phi, n, tol)

        implicit none
        real*8 :: mesh(4,n*n), delx, dely, edge1, edge2, total_length, old_residual
        real*8 :: phi(n,n), phi_old(n,n), a, b, c, rhs(2:n-1)
        integer :: n, niter, i, j, k
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

            a = -(1.0/delx**2)
            b = 2.0*((1.0/delx**2) + (1.0/dely**2))
            c = -(1.0/delx**2)
                
            ! Row Wise Sweep Loop starts
            do j = 2,n-1

            ! Calculating RHS Vector
              rhs(2)   = ((1.0/delx**2)*phi(1,j)) &
                        + ((1.0/dely**2)*phi(2,j-1)) &
                        + ((1.0/dely**2)*phi(2,j+1)) &
                         - g(mesh(3,(j-1)*n + 2),mesh(4,(j-1)*n + 2))
              rhs(n-1) = ((1.0/delx**2)*phi(n,j)) &
                         + ((1.0/dely**2)*phi(n-1,j-1)) &
                         + ((1.0/dely**2)*phi(n-1,j+1)) &
                         - g(mesh(3,(j-1)*n + n-1),mesh(4,(j-1)*n + n-1))


              do k = 3,n-2,1
                rhs(k) = ((1/dely**2)*phi(k,j-1)) &
                         + ((1/dely**2)*phi(k,j+1)) &
                         - g(mesh(3,(j-1)*n + k),mesh(4,(j-1)*n + k))
              end do
            
            ! Solving Implicit System using Thomas Algorithm
              call TDMA(a,b,c,rhs,phi,n,j,2,delx,dely)  

            end do
            ! Row Wise Sweep Loop ends

            call updatebc(n,phi,mesh)

            a = -(1.0/dely**2)
            b = 2.0*((1.0/delx**2) + (1.0/dely**2))
            c = -(1.0/dely**2)

            ! Column Wise Sweep Loop Starts
            do i = 2,n-1

            ! Calculating RHS Vector
              rhs(2)   = ((1.0/dely**2)*phi(i,1)) &
                         + ((1.0/delx**2)*phi(i-1,2)) &
                         + ((1.0/delx**2)*phi(i+1,2)) &
                         - g(mesh(3,n+i),mesh(4,n+i))
              rhs(n-1) = ((1.0/dely**2)*phi(i,n)) &
                         + ((1.0/delx**2)*phi(i-1,n-1)) &
                         + ((1.0/delx**2)*phi(i+1,n-1)) &
                         - g(mesh(3,(n-2)*n + i),mesh(4,(n-2)*n + i))

              do k = 3,n-2,1
                 rhs(k) = ((1/delx**2)*phi(i-1,k)) &
                          + ((1/delx**2)*phi(i+1,k)) &
                          - g(mesh(3,(k-1)*n + i),mesh(4,(k-1)*n + i))                  
              end do
            
            ! Solving Implicit System using Thomas Algorithm
              call TDMA(a,b,c,rhs,phi,n,i,1,delx,dely)

            end do
            ! Column Wise Sweep Loop ends

            niter = niter + 1

            print *, "Iteration Count = ", niter
            print *, "Residual = ", old_residual
            write(11,*) niter, old_residual 

        end do

        print *,"Number of Iterations= ", niter
        close(11)

end subroutine adi

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

      ! ADI      
      call adi(mesh, phi, n, tol)
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

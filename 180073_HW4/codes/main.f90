! AMAN PAREKH - 180073 - ME630 - Monsoo 2021

program main

      implicit none
      integer :: n, re, i, j, k
      real*16 :: dt, delx, dely, edge1, edge2, total_length, tol
      real*16 :: mesh(4,64*64), psi(64,64), omega(64,64), u(64,64), v(64,64)
      open(unit=10, file = 're100_64.dat', position='append')

      n = 64
      re = 100
      dt = 1e-4
      tol = 1E-5
      call meshgen(n,mesh)

      edge1 = -0.5/(n-2.0)
      edge2 = 1.0 + (0.5/(n-2.0))
      total_length = edge2 - edge1
      delx = (total_length*1.0)/(n-1.0)
      dely = (total_length*1.0)/(n-1.0)

      ! Initializing Psi using random values
      do j=1,n
        do k=1,n
          call random_number(psi(j,k))
          psi(j,k) = (psi(j,k)-0.5)*2.0
        end do
      end do

      ! Calculating u,v from Psi
      do j=2,n-1
        do k=2,n-1
          u(j,k) = (psi(j,k+1) - psi(j,k-1))/(2*dely)
          v(j,k) = (psi(j-1,k) - psi(j+1,k))/(2*delx)
        end do
      end do

      ! Calculating omega from u,v
      do j=2,n-1
        do k=2,n-1
          omega(j,k) = -((u(j,k+1)-u(j,k-1))/dely) &
                        +((v(j+1,k)-v(j-1,k))/delx)
        end do
      end do

      do i=1,1000       ! Time Loop
        
        call rk3(mesh, omega, u, v, re, n, dt, delx, dely)

        call seidel(mesh, psi, omega, n, tol)
        
        ! Updating u,v using Psi
        do j=2,n-1
          do k=2,n-1
            u(j,k) = (psi(j,k+1) - psi(j,k-1))/(2.0*dely)
            v(j,k) = (psi(j-1,k) - psi(j+1,k))/(2.0*delx)
          end do
        end do

        ! Updating BC
        call updatebc(n, omega)
        call updatebc(n, u)
        call updatebc(n, v)

        ! Printing Time Step 
        print *, 'Time Iteration: ', i
        
        ! Soring Data
        do k=2,n-1
          write (10,*) ((omega(j,k)), j=2,n-1)
        end do

      end do

      close(10)

end program main

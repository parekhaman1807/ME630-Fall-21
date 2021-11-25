! AMAN PAREKH - 180073 - ME630 - Monsoon 2021

subroutine TDMA(a,b,c,d,u,nx)

      implicit none

      integer :: nx, i
      real*16 :: a(nx), b(nx), c(nx), d(nx), u(nx)

      c(1) = c(1)/b(1)
      d(1) = d(1)/b(1)

      ! Forward Elimination
      do i = 2,nx
        c(i) = (c(i))/(b(i) - a(i)*c(i-1))
        d(i) = (d(i) - a(i)*d(i-1))/(b(i) - a(i)*c(i-1))
      end do

      ! Backward Substitution
      u(nx) = d(nx)

      do i = nx-1,1,-1
        u(i) = d(i) - u(i+1)*c(i)
      end do

end subroutine TDMA

program main

      implicit none

      integer :: re, i, j, n, nx
      real*16 :: x(201), u(201), a(201), b(201), c(201), d(201)
      real :: dx, vis, t, dt, s
      open(unit=10, file = 'b.dat', position = 'append')

      re = 10
      vis = 1.0/re

      nx = 201
      dx = 2.0/(nx-1)

      t = 0.47
      dt = 0.01
      n = t/dt

      s = (vis*dt)/(dx**2)

      do i = 1,nx
        x(i) = (i-1)*dx - 1.0
      end do
      
      write(10,*) x  

      ! Initializing u
      do i = 1,nx
        if(x(i).gt.0.0) then
          u(i) = 0.0
        else
          u(i) = 1.0
        end if
      end do

      do i = 1,n

        ! Setting Array Entries to implement Boundary Conditions in Time      
        b(1) = 1.0
        b(n) = 1.0

        a(n) = 0.0
        c(1) = 0.0

        d(1) = 1.0
        d(n) = 0.0
        
       ! Calculating a array
        do j = 2,nx-1
          a(j) = -(0.25*(dt/dx)*u(j-1)) - (0.5*s)
        end do

       ! Calculating b array
        do j = 2,nx-1
          b(j) = 1 + s
        end do

       ! Calculating c array
        do j = 2,nx-1
          c(j) = (0.25*(dt/dx)*u(j+1)) - (0.5*s)
        end do

       ! Calculating d array
        do j = 2,nx-1
          d(j) = (0.5*s*u(j-1)) + ((1-s)*u(j)) +  (0.5*s*u(j+1))
        end do

        call TDMA(a,b,c,d,u,nx)

        write(10,*) u
      end do

      close(10)

end program main


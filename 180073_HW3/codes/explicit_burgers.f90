! AMAN PAREKH  - 180073 - ME630 - Monsoon 2021

program main
      
      implicit none

      integer :: re, i, j, n, nx
      real*16 :: x(201), u(201), u_old(201)
      real :: dx, vis, t, dt, c, s
      open(unit=10, file = 'e.dat', position = 'append')

      re = 10
      vis = 1.0/re

      nx = 201
      dx = 2.0/(nx-1)

      t = 0.47
      dt = 0.0001
      n = t/dt

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
        u_old = u
        do j = 2,nx-1
          u(j) = u_old(j) &             ! Explicit Stencil (BS)
                 - (dt/(dx))*u_old(j)*(u_old(j) - u_old(j-1)) &
                 + ((vis*dt)/dx**2)*(u_old(j+1) - 2*u_old(j) + u_old(j-1))

!           u(j) = u_old(j) &             ! Explicit Stencil (CS)
!                 - (dt/(2*dx))*u_old(j)*(u_old(j+1) - u_old(j-1)) &
!                 + ((vis*dt)/dx**2)*(u_old(j+1) - 2*u_old(j) + u_old(j-1))

!           u(j) = u_old(j) &             ! Explicit Stencil (FS)
!                 - (dt/(dx))*u_old(j)*(u_old(j+1) - u_old(j)) &
!                 + ((vis*dt)/dx**2)*(u_old(j+1) - 2*u_old(j) + u_old(j-1))
        end do
        write(10,*) u
      end do

      close(10)

end program main

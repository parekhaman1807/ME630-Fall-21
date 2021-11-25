! AMAN PAREKH - 180073 - ME630 - Monsoon 2021

real*16 function rhs_func(uij, vij, oij, oim1j, oip1j, oijm1, oijp1, delx, dely, re)    ! RHS of RK3 equation

      implicit none
      integer :: re
      real*16 :: uij, vij, oij, oim1j, oip1j, oijm1, oijp1, delx, dely

      rhs_func = - (uij*((oip1j-oim1j)/(2.0*delx))) &
                 - (vij*((oijp1-oijm1)/(2.0*dely))) &
                 + (((oip1j - 2.0*oij + oim1j)/(delx**2)) + ((oijp1 - 2.0*oij + oijm1)/(dely**2)))/(1.0*re)

end function rhs_func

subroutine rk3(mesh, omega, u, v, re, n, dt, delx, dely)

      implicit none
      integer :: re, n, j, k
      real*16 :: mesh(4,n*n), omega(n,n), u(n,n), v(n,n), dt, delx, dely
      real*16 :: omegap(n,n), k1(n,n)
      real*16, external :: rhs_func

      ! RK-1
      omegap = omega

      do j=2,n-1
        do k=2,n-1
          k1(j,k) = rhs_func(u(j,k), v(j,k), omegap(j,k), omegap(j-1,k), &
                        omegap(j+1,k), omegap(j,k-1), omegap(j,k+1), delx, dely, re)
        end do
      end do

      ! RK-2
      do j=2,n-1
        do k=2,n-1
          omegap(j,k) = omegap(j,k) + ((dt/3.0)*k1(j,k))
        end do
      end do

      do j=2,n-1
        do k=2,n-1
          k1(j,k) = ((-5.0/9.0)*k1(j,k)) &
                    + rhs_func(u(j,k), v(j,k), omegap(j,k), omegap(j-1,k),&
                         omegap(j+1,k), omegap(j,k-1), omegap(j,k+1), delx, dely, re)
        end do
      end do

      ! RK-3
      do j=2,n-1
        do k=2,n-1
          omegap(j,k) = omegap(j,k) + (((15.0*dt)/16.0)*k1(j,k))
        end do
      end do

      do j=2,n-1
        do k=2,n-1
          k1(j,k) = ((-153.0/128.0)*k1(j,k)) &
                    + rhs_func(u(j,k), v(j,k), omegap(j,k), omegap(j-1,k), &
                        omegap(j+1,k), omegap(j,k-1), omegap(j,k+1), delx, dely, re)
        end do
      end do

      ! Calculating Final omega
      do j=2,n-1
        do k=2,n-1
          omega(j,k) = omegap(j,k) + (((8.0*dt)/15.0)*k1(j,k))
        end do
      end do

end subroutine rk3

! AMAN PAREKH - 180073 - ME630 - Monsoon 2021

subroutine seidel(mesh, psi, omega, n, tol)

        implicit none
        real*16 :: mesh(4,n*n), delx, dely, edge1, edge2, total_length, old_residual
        real*16 :: psi(n,n), psi_old(n,n), omega(n,n)
        integer :: n, niter, i, j
        real*16 :: residual, tol

        edge1 = -0.5/(n-2.0)
        edge2 = 1.0 + (0.5/(n-2.0))
        total_length = edge2 - edge1
        delx = (total_length*1.0)/(n-1.0)
        dely = (total_length*1.0)/(n-1.0)

        psi = 0.0
        psi_old = 1.0
        niter = 0

        do while(abs((residual(n,psi,psi_old))-abs(old_residual)).gt.tol)
            old_residual = residual(n,psi,psi_old)
            psi_old = psi
            call updatebc(n, psi)

            ! Siedel Loop starts
            do j = 2,n-1
              do i = 2,n-1
                psi(i,j) = 0.25*(psi(i+1,j) + psi(i-1,j) + psi(i,j+1) + psi(i,j-1) &
                                        + (omega(i,j)*delx*delx))
              end do
            end do
            ! Siedel Loop ends

            niter = niter + 1

        end do

        print *,"Number of Iterations= ", niter, "| Converged Residual= ", old_residual

end subroutine seidel

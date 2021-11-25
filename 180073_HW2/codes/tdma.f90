! AMAN PAREKH - 180073 - Fall 2021 - ME630

subroutine TDMA(a,b,c,rhs,phi,n,loc,ind,delx,dely)

      implicit none
      real*8 :: a, b, c, phi(n,n), delx, dely
      real*8 :: rhs(2:n-1), At(2:n-1), Bt(2:n-1), Ct(2:n-1)
      integer :: i, loc, ind, n

      At(2:n-1) = a
      Bt(2:n-1) = b
      Ct(2:n-1) = c

      ! Forward Elimination 
         
      Ct(2) = Ct(2)/Bt(2)
      rhs(2) = rhs(2)/Bt(2)

      do i= 3,n-1
          Ct(i) = Ct(i)/(Bt(i) - At(i)*Ct(i-1))
          rhs(i) = (rhs(i) - At(i)*rhs(i-1))/(Bt(i) - At(i)*Ct(i-1))
      end do

 
      if(ind.eq.1) then ! This case is for Row Wise Sweep

        ! Backward Substitution

        phi(loc,n-1) = rhs(n-1)

        do i= n-2,2,-1
            phi(loc,i) = rhs(i) - phi(loc,i+1)*Ct(i)
        end do

      else  ! This case is for Column Wise Sweep
        
        ! Backward Substitution

        phi(n-1,loc) = rhs(n-1)

        do i = n-2,2,-1
            phi(i,loc) = rhs(i) - phi(i+1,loc)*Ct(i)
        end do

      end if

end subroutine TDMA

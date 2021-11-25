! AMAN PAREKH - 180073 - ME630 - Monsoon 2021

subroutine meshgen(n, mesh)     ! Subroutine for Mesh Generation

      implicit none
      integer :: n, i, j, c  
      real*16 :: mesh(4,n*n)
      real*16 :: edge1, edge2, total_length, del

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



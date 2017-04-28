program main
  implicit none

  integer, parameter   :: n1 = 1024, n2 = 1024, niter = 1
  integer              :: i, j, n
  real(8), parameter   :: epsilon = 0.1
  real(8), allocatable :: a(:,:), b(:,:)
  real(8), allocatable :: x(:)  , y(:)

  allocate(a(n1,n2))
  allocate(b(n1,n2))
  allocate(x(n1))
  allocate(y(n1))
  
  do i = 1, n1
     x(i) = 1./dble(n1-1)*dble(i-1)
  end do

  do j = 1, n2
     y(j) = 1./dble(n2-1)*dble(j-1)
  end do

  do j = 1, n2
     do i = 1, n1
        if (x(i) .lt. 0.5) then
           a(i,j) = cos(x(i)+y(j))
        else
           a(i,j) = sin(x(i)+y(j))
        end if
        b(i,j) = a(i,j)
     end do
  end do

  do n = 1, niter
     do j = 2, n2-1
        do i = 2, n1-1
           a(i,j) = b(i,j)+ epsilon*( &
                b(i-1,j+1)+         b(i,j+1)+b(i+1,j+1) + &
                b(i-1,j  )-dble(8.)*b(i,j  )+b(i+1,j  ) + &
                b(i-1,j-1)+         b(i,j-1)+b(i+1,j-1))
        end do
     end do
     b = a
  end do

  print *, "norm(a) = ", sqrt(sum(a**2))
  print *, 1,a(512,512)
  print *, 2,a(513,512)
  print *, 3,a(512,513)
  print *, 4,a(513,513)
  
  deallocate(a,b,x,y)

end program main

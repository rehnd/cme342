program main
  implicit none
  include 'mpif.h'

  integer, parameter   :: n1 = 1024, n2 = 1024, niter = 100
  integer              :: i, j, n
  integer              :: my_id, size, ierr
  integer              :: np1, np2
  real(8), parameter   :: epsilon = 0.1
  real(8), allocatable :: a(:,:), b(:,:)
  real(8), allocatable :: x(:)  , y(:)
  integer              :: nid(4) ! neighbor id's

  np1 = 4
  np2 = 4

  allocate(a(n1,n2))
  allocate(b(n1,n2))
  allocate(x(n1))
  allocate(y(n1))

  call mpi_init(ierr)
  
  call mpi_comm_size(mpi_comm_world, size, ierr)
  call mpi_comm_rank(mpi_comm_world, my_id, ierr)
  my_id = my_id + 1

  if ((.not. np1*np2  == size)) then
     if (my_id == 1) print *, 'Error: Domain decomp (np1*np2) not equal to MPI_RANK'
     stop('')
  end if

  call initialize_arrays(a,b,x,y)

  call get_neighbor_ids(my_id, nid, np1, np2)

  deallocate(a,b,x,y)

  call mpi_finalize(ierr)

contains

  
  subroutine get_neighbor_ids(my_id, nid, np1, np2)
    integer, intent(in)    :: my_id, np1, np2
    integer, intent(inout) :: nid(4)
    integer                :: ip1, ip2 ! indices of my_id in the global array

    ip2 = (my_id-1)/np1 + 1
    ip1 =  my_id - (ip2-1)*np2
    
    !print *, 'my_id, ip1, ip2 ', my_id, ip1, ip2

    nid(:) = -1

    if (my_id - np1  >  0)           nid(1) = my_id - np1    
    if (my_id + 1    <= ip2*np1)     nid(2) = my_id + 1
    if (my_id + np1  <= np1*np2)     nid(3) = my_id + np1
    if (my_id - 1    > (ip2-1)*np1)  nid(4) = my_id - 1
    
    !print *, 'my_id, np', my_id, nid(:)

  end subroutine get_neighbor_ids

  
  subroutine initialize_arrays(a,b,x,y)
    real(8), intent(inout) :: a(:,:), b(:,:), x(:), y(:)

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
             a(i,j) = b(i,j) + epsilon*( &
                  b(i-1,j+1)+         b(i,j+1)+b(i+1,j+1) + &
                  b(i-1,j  )-dble(8.)*b(i,j  )+b(i+1,j  ) + &
                  b(i-1,j-1)+         b(i,j-1)+b(i+1,j-1))
          end do
       end do
       b = a
    end do

  end subroutine initialize_arrays

end program main

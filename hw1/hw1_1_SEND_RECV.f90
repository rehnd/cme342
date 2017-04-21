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
  integer              :: nid(8) ! neighbor id's

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


  deallocate(a,b,x,y)

  call mpi_finalize(ierr)

contains

  
  subroutine get_neighbor_ids(my_id, nid, np1, np2)
    integer, intent(in)    :: my_id, np1, np2
    integer, intent(inout) :: nid(8)
    integer                :: ip1, ip2   ! indices of my_id in the global array
    character(len=10)      :: fmt = "(1I,8I)"

    ! Fill the nid array based on nearest neighbors. The grid for 16 processors is:
    !
    !       -------->  j (ip2, np2)
    !        ________________________
    !   |   |     |     |     |      |
    !   |   |  1  |  5  |  9  |  13  |
    !   |   |_____|_____|_____|______|
    !   v   |     |     |     |      |
    !       |  2  |  6  |  10 |  14  |
    !   i   |_____|_____|_____|______|
    ! (ip1) |     |     |     |      | 
    ! (np1) |  3  |  7  |  11 |  15  |
    !       |_____|_____|_____|______|
    !       |     |     |     |      |
    !       |  4  |  8  |  12 |  16  |
    !       |_____|_____|_____|______|
    !
    ! nid filling:
    !    Use nid(1,...,8) = (l, d, r, u, lu, ld, rd, ru)
    !    where:  l = left, d = down, r = right, u = up,  lu = left+up, etc.
    !
    ! Examples:
    !    nid(1,...,8) for my_id = 6 should be: (2, 7,10,5,  1, 3,11, 9)
    !    nid(1,...,8) for my_id = 8 should be: (4,-1,12,7,  3,-1,-1,11)

    ! First determine index in i,j of my_id
    ip2 = (my_id-1)/np1 + 1
    ip1 =  my_id - (ip2-1)*np2
    
    nid(:) = -1

    ! Set the first 4 elements of nid (edges)
    if (my_id - np1  >  0)           nid(1) = my_id - np1    
    if (my_id + 1    <= ip2*np1)     nid(2) = my_id + 1
    if (my_id + np1  <= np1*np2)     nid(3) = my_id + np1
    if (my_id - 1    > (ip2-1)*np1)  nid(4) = my_id - 1

    ! Set the last 4 elements of nid (corners)
    if (nid(4) > 0 .and. nid(1) > 0) nid(5) = my_id-np1-1
    if (nid(1) > 0 .and. nid(2) > 0) nid(6) = my_id-np1+1
    if (nid(2) > 0 .and. nid(3) > 0) nid(7) = my_id+np1+1
    if (nid(3) > 0 .and. nid(4) > 0) nid(8) = my_id+np1-1

    ! write (*,fmt)  my_id, nid(:)

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

  end subroutine initialize_arrays

end program main

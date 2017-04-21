program main
  implicit none
  include 'mpif.h'

  integer, parameter   :: n1 = 1024, n2 = 1024, niter = 100
  integer              :: i, j, n
  integer              :: my_id, np, ierr
  real(8), parameter   :: epsilon = 0.1
  real(8), allocatable :: a(:,:), b(:,:)
  real(8), allocatable :: x(:)  , y(:)
  integer              :: nid(8) ! neighbor id's
  integer              :: i_f, i_l, j_f, j_l ! => i_first, i_last, j_first, j_last
  integer              :: np1, np2  
  np1 = 4
  np2 = 4

  !  MPI Setup, get rank, size, and run tests
  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world, np, ierr)
  call mpi_comm_rank(mpi_comm_world, my_id, ierr)
  my_id = my_id + 1
  call check_num_processes(my_id, np, np1, np2)

  ! Set up arrays on each processor, get the indices this process will work on
  call allocate_arrays(a,b,x,y,n1,n2)
  call initialize_arrays(a,b,x,y)
  call get_my_ij(i_f, i_l, j_f, j_l, np1, np2, my_id, n1, n2)
  !call get_neighbor_ids(my_id, nid, np1, np2)

  ! Main loop. Inside, update the subgrid on each process and send data
  do n = 1, niter
     call update_subgrid(a,b,epsilon, i_f, i_l, j_f, j_l)

     !call send_recv_data()
  end do

  call deallocate_arrays(a,b,x,y)

  call mpi_finalize(ierr)

contains
  
  subroutine send_recv_data(np1, np2, my_id)
    integer, intent(in) :: np1, np2, my_id
    integer  :: ip, jp, curr_id, p

     do jp = 1, np2
        do ip = 1, np1
           curr_id = np1*(jp-1) + i
           do p = 1, np
              if (my_id == curr_id) then
                 !call mpi_send()
              else
                 !call mpi_recv()
              end if
           end do
        end do
     end do


  end subroutine send_recv_data


  subroutine update_subgrid(a,b,epsilon, i_f, i_l, j_f, j_l)
    real(8), intent(inout) :: a(:,:), b(:,:)
    integer, intent(in)    :: i_f, i_l, j_f, j_l
    real(8), intent(in)    :: epsilon
    integer                :: i, j

    do j = j_f+1, j_l-1
       do i = i_f+1, i_f-1
          a(i,j) = b(i,j) + epsilon*( &
               b(i-1,j+1)+         b(i,j+1)+b(i+1,j+1) + &
               b(i-1,j  )-dble(8.)*b(i,j  )+b(i+1,j  ) + &
               b(i-1,j-1)+         b(i,j-1)+b(i+1,j-1))
       end do
    end do
    b = a
    ! or b(i_f:i_l, j_f:j_l) = a(i_f:i_l, j_f:j_l)

  end subroutine update_subgrid

  subroutine allocate_arrays(a,b,x,y,n1,n2)
    real(8), allocatable, intent(inout) :: a(:,:), b(:,:), x(:), y(:)
    integer, intent(in)                 :: n1, n2

    allocate(a(n1,n2))
    allocate(b(n1,n2))
    allocate(x(n1))
    allocate(y(n2))
  end subroutine allocate_arrays

  subroutine deallocate_arrays(a,b,x,y)
    real(8), allocatable, intent(inout) :: a(:,:), b(:,:), x(:), y(:)
    deallocate(a,b,x,y)
  end subroutine deallocate_arrays
  subroutine get_my_ij(i_f, i_l, j_f, j_l, np1, np2, my_id, n1, n2)
    integer, intent(inout) :: i_f, i_l, j_f, j_l
    integer, intent(in)    :: np1, np2, my_id, n1, n2
    integer                :: ip, jp
    character(len=10)      :: fmt = '(5I)'

    jp = (my_id-1)/np1 + 1
    ip =  my_id - (jp-1)*np1
    
    i_f = (n1/np1)*(ip-1) + 1
    i_l = i_f + (n1/np1)  - 1
    j_f = (n2/np2)*(jp-1) + 1
    j_l = j_f + (n2/np2)  - 1

    !if (my_id == 1) write(*,*) '       my_id     i_f      i_l      j_f      j_l'
    !write(*,fmt) my_id, i_f, i_l, j_f, j_l

  end subroutine get_my_ij
  
  subroutine check_num_processes(my_id, np, np1, np2)
    integer, intent(in) :: my_id, np, np1, np2

    if ((.not. np1*np2  == np)) then
       if (my_id == 1) print *, 'Error: Domain decomp (np1*np2) not equal to MPI_RANK'
       stop('')
    end if
  end subroutine check_num_processes


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

    !write (*,fmt)  my_id, nid(:)

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

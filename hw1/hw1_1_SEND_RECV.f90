program main
  implicit none
  include 'mpif.h'
  integer, parameter    :: n1 = 1024, n2 = 1024, niter = 100
  integer               :: i, j, n
  integer               :: my_id, np, ierr
  real(8), parameter    :: epsilon = 0.1
  real(8), allocatable  :: a(:,:), b(:,:)
  real(8), allocatable  :: x(:)  , y(:)
  integer               :: nid(4) ! neighbor id's
  integer               :: i_f, i_l, j_f, j_l ! => i_first, i_last, j_first, j_last
  integer               :: np1, np2  
  np1 = 4
  np2 = 4

  !  MPI Setup, get rank, size, and run tests
  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world,    np, ierr)
  call mpi_comm_rank(mpi_comm_world, my_id, ierr)
  my_id = my_id + 1
  call check_num_processes()

  ! Set up arrays on each processor, get the indices this process will work on
  call allocate_arrays(n1/np1,n2/np2)
  call get_my_ij()
  call initialize_arrays()
  nid = get_neighbor_ids(my_id)

  ! Main loop. Inside, update the interior points on each process and send data
  do n = 1, niter
     !call send_receive()
     call update_interior()
     call update_edges()
     b=a
  end do

  call deallocate_arrays()

  call mpi_finalize(ierr)

contains


  subroutine send_receive()
    integer  :: ip, jp, p
    integer  :: cid         ! current id
    integer  :: cid_nid(4)  ! current id's neighbors

     do jp = 1, np2
        do ip = 1, np1
           cid     = np1*(jp-1) + i
           cid_nid = get_neighbor_ids(cid)
           do p = 1, np
              if (my_id == cid) then
                 !call mpi_send()
              else
                 !call mpi_recv()
              end if
           end do
        end do
     end do

  end subroutine send_receive


  subroutine update_interior()
    integer                :: i, j, n1me, n2me
    n1me = n1/np1
    n2me = n2/np2

    do j = 2, n2me-1
       do i = 2, n1me-1
          a(i,j) = b(i,j) + epsilon*( &
               b(i-1,j+1)+         b(i,j+1)+b(i+1,j+1) + &
               b(i-1,j  )-dble(8.)*b(i,j  )+b(i+1,j  ) + &
               b(i-1,j-1)+         b(i,j-1)+b(i+1,j-1))
       end do
    end do

  end subroutine update_interior


  subroutine update_edges()
    integer :: n1me, n2me

    n1me = n1/np1
    n2me = n2/np2

    if (nid(1) > 0)  then
       a(:,1)        = a(:,1)        + b(:,2)        - dble(8.)*b(:,1)
       a(1:n1me-1,1) = a(1:n1me-1,1) + b(2:n1me  ,2)
       a(2:n1me  ,1) = a(2:n1me  ,1) + b(1:n1me-1,2)
    end if

    if (nid(2) > 0) then
       a(n1me,:)        = a(n1me,:)        + b(n1me-1,:)     - dble(8.)*b(n1me,:)
       a(n1me,1:n2me-1) = a(n1me,1:n2me-1) + b(n1me-1,2:n2me)
       a(n1me,2:n2me)   = a(n1me,2:n2me)   + b(n1me-1,1:n2me-1)
    end if

    if (nid(3) > 0) then
       a(:,n2me)        = a(:,n2me)        + b(:,n2me-1)     - dble(8.)*b(:,n2me)
       a(1:n1me-1,n2me) = a(1:n1me-1,n2me) + b(2:n1me  ,n2me-1)
       a(2:n1me  ,n2me) = a(2:n1me  ,n2me) + b(1:n1me-1,n2me-1)

    end if
    
    if (nid(4) > 0) then
       a(1,:)        = a(1,:)        + b(2,:)     - dble(8.)*b(1,:)
       a(1,1:n2me-1) = a(1,1:n2me-1) + b(2,2:n2me)
       a(1,2:n2me)   = a(1,2:n2me)   + b(2,1:n2me-1)       
    end if

  end subroutine update_edges


  subroutine allocate_arrays(n1me,n2me)
    integer, intent(in)                 :: n1me, n2me

    allocate(a(n1me,n2me))
    allocate(b(n1me,n2me))
    allocate(x(n1me))
    allocate(y(n2me))
  end subroutine allocate_arrays


  subroutine deallocate_arrays()
    deallocate(a,b,x,y)
  end subroutine deallocate_arrays


  subroutine get_my_ij()
    integer                :: ip, jp
    character(len=10)      :: fmt = '(5I)'
    logical                :: printtest = .false.

    jp = (my_id-1)/np1 + 1
    ip =  my_id - (jp-1)*np1
    
    i_f = (n1/np1)*(ip-1) + 1
    i_l = i_f + (n1/np1)  - 1
    j_f = (n2/np2)*(jp-1) + 1
    j_l = j_f + (n2/np2)  - 1

    if (printtest) then
       if (my_id == 1) write(*,'(5A12)') 'my_id', 'i_f', 'i_l', 'j_f', 'j_l'
       write(*,fmt) my_id, i_f, i_l, j_f, j_l
    end if
  end subroutine get_my_ij
  

  subroutine check_num_processes()
    if ((.not. np1*np2  == np)) then
       if (my_id == 1) print *, 'Error: Domain decomp (np1*np2) not equal to MPI_RANK'
       stop('')
    end if
  end subroutine check_num_processes


  pure function get_neighbor_ids(id) 
    integer, intent(in)  :: id
    integer              :: ip1, ip2   ! indices of my_id in the global array
    integer              :: get_neighbor_ids(4)
    !character(len=10)    :: fmt = "(1I,4I)"

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
    !    Use nid(1,...,8) = (l, d, r, u)
    !    where:  l = left, d = down, r = right, u = up,  lu = left+up, etc.
    !
    ! Examples:
    !    nid(1,...,8) for id = 6 should be: (2, 7,10,5)
    !    nid(1,...,8) for id = 8 should be: (4,-1,12,7)

    ! First determine index in i,j of id
    ip2 = (id-1)/np1 + 1
    ip1 =  id - (ip2-1)*np2
    
    get_neighbor_ids(:) = -1

    ! Set the first 4 elements of get_neighbor_ids (edges)
    if (id - np1  >  0)           get_neighbor_ids(1) = id - np1    
    if (id + 1    <= ip2*np1)     get_neighbor_ids(2) = id + 1
    if (id + np1  <= np1*np2)     get_neighbor_ids(3) = id + np1
    if (id - 1    > (ip2-1)*np1)  get_neighbor_ids(4) = id - 1

    !write (*,fmt)  id, get_neighbor_ids(:)

  end function get_neighbor_ids

  
  subroutine initialize_arrays()
    logical           :: testvals = .false.
    character(len=20) :: fmt = "(1I,2F12.2)"

    do i = 1, n1/np1
       x(i) = 1./dble(n1-1)*dble(i+i_f-2)
    end do
    
    do j = 1, n2/np2
       y(j) = 1./dble(n2-1)*dble(j+j_f-2)
    end do

    do j = 1, n2/np2
       do i = 1, n1/np1
          if (x(i) .lt. 0.5) then
             a(i,j) = cos(x(i+i_f-1)+y(j+j_f-1))
          else
             a(i,j) = sin(x(i+i_f-1)+y(j+j_f-1))
          end if
          b(i,j) = a(i,j)
       end do
    end do

    if (testvals) then
       if (my_id == 1) write(*,'(3A12)') 'my_id', 'x(1)', 'y(1)'
       write(*,fmt) my_id, x(1), y(1)
    end if

  end subroutine initialize_arrays

end program main

program main
  implicit none
  include 'mpif.h'

  integer, parameter            :: n1 = 1024, n2 = 1024, niter = 1
  integer                       :: i, j, n
  integer                       :: my_id, np, ierr, status
  double precision              :: norm
  double precision, parameter   :: epsilon = 0.1
  double precision, allocatable :: a(:,:), b(:,:), b1(:),b2(:),b3(:),b4(:)
  double precision, allocatable :: x(:)  , y(:)
  integer                       :: nid(4) ! neighbor id's
  integer                       :: i_f, i_l, j_f, j_l ! => i_first, i_last, j_first, j_last
  integer                       :: np1, np2
  integer                       :: n1me, n2me
  integer,allocatable           :: top(:), left(:)

  np1 = 1
  np2 = 2
  n1me = n1/np1
  n2me = n2/np2

  !  MPI Setup, get rank, size, and run tests
  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world,    np, ierr)
  call mpi_comm_rank(mpi_comm_world, my_id, ierr)
  my_id = my_id + 1
  call check_num_processes()

  ! Set up arrays on each processor, get the indices this process will work on
  call allocate_arrays()
  call get_my_ij()
  call initialize_arrays()
  nid = get_neighbor_ids(my_id)

  ! Main loop. Inside, update the interior points on each process and send data
  do n = 1, niter
     call send_receive()
     call update_interior()
     call update_edges()
     b=a
  end do

  call get_total_norm()
  if (my_id == 1) print *, "norm(a) = ", norm

  call deallocate_arrays()
  call mpi_finalize(ierr)


contains


  subroutine send_receive()
    integer  :: ip, jp, p, tag=2001
    integer  :: cid         ! current id
    integer  :: cid_nid(4)  ! current id's neighbors

    ! Send and receive N to S (nid(4) = N, nid(2) = S)
    if (any(top == my_id) .and. nid(2) > 0) then
       b2 = b(n1me,:)
       call mpi_send(b2,n2me,mpi_double_precision,nid(2)-1,tag,mpi_comm_world, ierr)
       call mpi_recv(b2,n2me,mpi_double_precision,nid(2)-1,tag,mpi_comm_world,status, ierr)
       a(n1me,:) = a(n1me,:) + epsilon*b2

    elseif (nid(4) > 0 .and. nid(2) > 0) then
       call mpi_recv(b4,n2me,mpi_double_precision,nid(4)-1,tag,mpi_comm_world, status, ierr)
       a(1,:) = a(1,:) + epsilon*b4
       b4     = b(1,:)
       call mpi_send(b4,n2me,mpi_double_precision,nid(4)-1,tag,mpi_comm_world, ierr)
       
       b2 = b(n1me,:)
       call mpi_send(b2,n2me,mpi_double_precision,nid(2)-1,tag,mpi_comm_world, ierr)
       call mpi_recv(b2,n2me,mpi_double_precision,nid(2)-1,tag,mpi_comm_world,status,ierr)
       a(n1me,:) = a(n1me,:) + epsilon*b2

    elseif (nid(4) > 0) then
       call mpi_recv(b4,n2me,mpi_double_precision,nid(4)-1,tag,mpi_comm_world, status, ierr)
       a(1,:) = a(1,:) + epsilon*b4
       b4     = b(1,:)
       call mpi_send(b4,n2me,mpi_double_precision,nid(4)-1,tag,mpi_comm_world, ierr)

    end if

    ! Send and receive W to E (nid(1) = W, nid(3) = E)
    if (any(left == my_id) .and. nid(3) > 0) then
       b3 = b(:,n2me)
       call mpi_send(b3,n1me,mpi_double_precision,nid(3)-1,tag,mpi_comm_world, ierr)
       call mpi_recv(b3,n1me,mpi_double_precision,nid(3)-1,tag,mpi_comm_world,status, ierr)
       a(2:n1me-1,n2me) = a(2:n1me-1,n2me) + epsilon*b3(2:n1me-1)


    elseif (nid(1) > 0 .and. nid(3) > 0) then
       call mpi_recv(b1,n1me,mpi_double_precision,nid(1)-1,tag,mpi_comm_world, status, ierr)
       a(:,1) = a(:,1) + epsilon*b1
       b1     = b(:,1)
       call mpi_send(b1,n1me,mpi_double_precision,nid(1)-1,tag,mpi_comm_world, ierr)
       
       b3 = b(:,n2me)
       call mpi_send(b3,n1me,mpi_double_precision,nid(3)-1,tag,mpi_comm_world, ierr)
       call mpi_recv(b3,n1me,mpi_double_precision,nid(3)-1,tag,mpi_comm_world,status,ierr)
       a(:,n2me) = a(:,n2me) + epsilon*b3

    elseif (nid(1) > 0) then
       call mpi_recv(b1,n1me,mpi_double_precision,nid(1)-1,tag,mpi_comm_world, status, ierr)
       a(:,1) = a(:,1) + epsilon*b1
       b1     = b(:,1)
       call mpi_send(b1,n1me,mpi_double_precision,nid(1)-1,tag,mpi_comm_world, ierr)

    end if

  end subroutine send_receive


  subroutine update_interior()
    integer :: i, j

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

    if (nid(1) > 0)  then
       a(:,1)        = a(:,1)        + epsilon*( b(:,2)        - dble(8.)*b(:,1))
       a(1:n1me-1,1) = a(1:n1me-1,1) + epsilon*b(2:n1me  ,2)
       a(2:n1me  ,1) = a(2:n1me  ,1) + epsilon*b(1:n1me-1,2)
    end if

    if (nid(2) > 0) then
       a(n1me,:)        = a(n1me,:)        + epsilon*(b(n1me-1,:)     - dble(8.)*b(n1me,:))
       a(n1me,1:n2me-1) = a(n1me,1:n2me-1) + epsilon*b(n1me-1,2:n2me)
       a(n1me,2:n2me)   = a(n1me,2:n2me)   + epsilon*b(n1me-1,1:n2me-1)
    end if

    if (nid(3) > 0) then
       a(:,n2me)        = a(:,n2me)        + epsilon*(b(:,n2me-1)     - dble(8.)*b(:,n2me))
       a(1:n1me-1,n2me) = a(1:n1me-1,n2me) + epsilon*b(2:n1me  ,n2me-1)
       a(2:n1me  ,n2me) = a(2:n1me  ,n2me) + epsilon*b(1:n1me-1,n2me-1)

    end if
    
    if (nid(4) > 0) then
       a(1,:)        = a(1,:)        + epsilon*(b(2,:)     - dble(8.)*b(1,:))
       a(1,1:n2me-1) = a(1,1:n2me-1) + epsilon*b(2,2:n2me)
       a(1,2:n2me)   = a(1,2:n2me)   + epsilon*b(2,1:n2me-1)       
    end if

  end subroutine update_edges


  subroutine get_total_norm()
    integer :: ip, tag = 3001
    double precision :: norm_ip
    ! First sum squares of entries
    norm = sum(a**2)
    
    if (my_id == 1) then
       do ip = 2,np
          call mpi_recv(norm_ip,1,mpi_double_precision,ip-1,tag,mpi_comm_world,status,ierr)
          norm = norm + norm_ip
       end do
       norm = sqrt(norm)
    else
       call mpi_send(norm,1,mpi_double_precision,0,tag,mpi_comm_world,ierr)
    end if
    
  end subroutine get_total_norm


  subroutine allocate_arrays()
    allocate(a(n1me,n2me))
    allocate(b(n1me,n2me))
    allocate(x(n1me))
    allocate(y(n2me))

    allocate(b1(n1me))
    allocate(b2(n2me))
    allocate(b3(n1me))
    allocate(b4(n2me))

    allocate( top(np2))
    allocate(left(np1))

  end subroutine allocate_arrays


  subroutine deallocate_arrays()
    deallocate(a,b,x,y,b1,b2,b3,b4)
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

    do i=1,np2
       top(i) = np1*(i-1)+1
    end do

    do i=1,np1
       left(i) = i
    end do

    if (testvals) then
       if (my_id == 1) write(*,'(3A12)') 'my_id', 'x(1)', 'y(1)'
       write(*,fmt) my_id, x(1), y(1)
    end if

  end subroutine initialize_arrays

end program main

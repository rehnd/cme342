program main
  implicit none
  include 'mpif.h'
  integer                       :: n1 = 1024, n2 = 1024, np1=1, np2=1, niter = 1
  integer                       :: n1me, n2me
  integer                       :: nid(4)             ! neighbor id's
  integer                       :: my_id, np, ierr, status
  integer                       :: i, j, n
  integer                       :: i_f, i_l, j_f, j_l ! i_first, i_last, j_first, j_last
  double precision              :: norm
  double precision, parameter   :: epsilon = 0.1
  double precision, allocatable :: a(:,:),b(:,:)
  double precision, allocatable :: b1w(:),b1me(:),b2s(:),b2me(:),b3e(:),b3me(:),b4n(:),b4me(:)
  double precision, allocatable :: x(:), y(:)
  integer,          allocatable :: top(:), left(:)    ! nodes at the top, left of processor grid

  !  MPI Setup, get rank, size. Increment my_id by 1 for Fortran 1-based indexing
  call mpi_init(ierr)
  call mpi_comm_size(mpi_comm_world,    np, ierr)
  call mpi_comm_rank(mpi_comm_world, my_id, ierr)
  my_id = my_id + 1

  ! Read in input values and run tests
  call read_input()
  call check_num_processes()

  ! Set up arrays on each processor, get the indices this process will work on
  call allocate_arrays()
  call get_my_ij()
  call initialize_arrays()
  nid = get_neighbor_ids(my_id)

  ! Main loop. Inside, update the interior & edge points on each processor
  do n = 1, niter
     call update_edges()
     call update_interior()
     b=a
  end do

  ! Compute the norm as a quick means of testing correctness
  call get_total_norm()
  if (my_id == 1) print *, "norm(a) = ", norm

  call deallocate_arrays()
  call mpi_finalize(ierr)


contains


  subroutine update_edges()
    integer  :: ip, jp, p, tag=2001
    integer  :: cid         ! current id
    integer  :: cid_nid(4)  ! current id's neighbors

    ! Send and receive N to S (nid(4) = N, nid(2) = S)
    if (any(top == my_id) .and. nid(2) > 0) then

       b2me = b(n1me,:)
       call mpi_sendrecv( &
            b2me,n2me,mpi_double_precision,nid(2)-1,tag, &
            b2s, n2me,mpi_double_precision,nid(2)-1,tag, &
            mpi_comm_world,status,ierr)
       call update_S_edge()

    elseif (nid(4) > 0 .and. nid(2) > 0) then

       b4me = b(1,:)
       call mpi_sendrecv( &
            b4me,n2me,mpi_double_precision,nid(4)-1,tag,&
            b4n ,n2me,mpi_double_precision,nid(4)-1,tag,&
            mpi_comm_world,status,ierr)
       call update_N_edge()
       
       b2me = b(n1me,:)
       call mpi_sendrecv( &
            b2me,n2me,mpi_double_precision,nid(2)-1,tag, &
            b2s ,n2me,mpi_double_precision,nid(2)-1,tag, &
            mpi_comm_world,status,ierr)
       call update_S_edge()

    elseif (nid(4) > 0) then

       b4me = b(1,:)
       call mpi_sendrecv( &
            b4me,n2me,mpi_double_precision,nid(4)-1,tag, &
            b4n, n2me,mpi_double_precision,nid(4)-1,tag, &
            mpi_comm_world,status,ierr)
       call update_N_edge()
       
    end if

    ! Send and receive W to E (nid(1) = W, nid(3) = E)
    if (any(left == my_id) .and. nid(3) > 0) then

       b3me(2:n1me+1) = b(:,n2me)
       call mpi_sendrecv( &
            b3me,n1me+2,mpi_double_precision,nid(3)-1,tag, &
            b3e, n1me+2,mpi_double_precision,nid(3)-1,tag, &
            mpi_comm_world,status,ierr)
       call update_E_edge()

    elseif (nid(1) > 0 .and. nid(3) > 0) then

       b1me(2:n1me+1) = b(:,1)
       call mpi_sendrecv( &
            b1me,n1me+2,mpi_double_precision,nid(1)-1,tag, &
            b1w, n1me+2,mpi_double_precision,nid(1)-1,tag, &
            mpi_comm_world,status,ierr)
       call update_W_edge()

       b3me(2:n1me+1) = b(:,n2me)
       call mpi_sendrecv( &
            b3me,n1me+2,mpi_double_precision,nid(3)-1,tag, &
            b3e, n1me+2,mpi_double_precision,nid(3)-1,tag, &
            mpi_comm_world,status,ierr)
       call update_E_edge()

    elseif (nid(1) > 0) then

       b1me(2:n1me+1)  = b(:,1)
       call mpi_sendrecv( &
            b1me,n1me+2,mpi_double_precision,nid(1)-1,tag, &
            b1w, n1me+2,mpi_double_precision,nid(1)-1,tag, &
            mpi_comm_world,status,ierr)
       call update_W_edge()

    end if

  end subroutine update_edges


  subroutine update_N_edge()
    ! Immediately fill b1(1), b3(1)
    b1me(1) = b4n(1)
    b3me(1) = b4n(n2me)

    ! Perform complete update of the interior edge elements
    a(1,2:n2me-1) = b(1,2:n2me-1) + epsilon*( &
         + b4n(2:n2me-1)          & ! above
         + b4n(1:n2me-2)          & ! above-left
         + b4n(3:n2me)            & ! above-right
         + b(1,1:n2me-2)          & ! left
         - b(1,2:n2me-1)*dble(8.) & ! yourself 
         + b(1,3:n2me)            & ! right
         + b(2,1:n2me-2)          & ! below-left
         + b(2,2:n2me-1)          & ! below
         + b(2,3:n2me)            & ! below-right
         )

    ! Treat node corners as special case
    if (nid(1) > 0) then
       a(1,1) = b(1,1) + epsilon*(      &
            + b4n(1) + b4n(2)           & ! up, up-right
            + b(1,2) - dble(8.)*b(1,1)  & ! yourself, right
            + b(2,1) + b(2,2) )           ! below, below-right
    end if

    if (nid(3) > 0) then
       a(1,n2me) = b(1,n2me) + epsilon*(      &
            + b4n(n2me-1) + b4n(n2me)         & ! up, up-left
            + b(1,n2me-1) -dble(8.)*b(1,n2me) & ! yourself, left
            + b(2,n2me-1) + b(2,n2me))          ! below, below-left
    end if

  end subroutine update_N_edge
  

  subroutine update_S_edge()
    ! Immediately fill b1(end), b3(end) to account for corner values
    b1me(n1me+2) = b2s(1)
    b3me(n1me+2) = b2s(n2me)

    ! Perform a complete update of the interior edge elements |_|x|x|...|x|x|_|
    a(n1me,2:n2me-1) = b(n1me,2:n2me-1) + epsilon*( & 
         + b2s(2:n2me-1)                            & ! below
         + b2s(1:n2me-2)                            & ! below-left
         + b2s(3:n2me)                              & ! below-right
         + b(n1me,1:n2me-2)                         & ! left
         - b(n1me,2:n2me-1)*dble(8.)                & ! yourself
         + b(n1me,3:n2me)                           & ! right
         + b(n1me-1,2:n2me-1)                       & ! above
         + b(n1me-1,1:n2me-2)                       & ! above-left
         + b(n1me-1,3:n2me)                         & ! above-right
         )

    ! Treat corners as a special case, doing a *partial* update
    if (nid(1) > 0) then
       a(n1me,1) = b(n1me,1) + epsilon*(       &
            b(n1me-1,1) + b(n1me-1,2) +        & ! up, up-right
            b(n1me,2)  - dble(8.)*b(n1me,1) +  & ! yourself, right
            b2s(2) + b2s(1))                     ! below, below-right
    end if

    if (nid(3) > 0) then
       a(n1me,n2me) = b(n1me,n2me) + epsilon*(        &
            + b(n1me-1,n2me) + b(n1me-1,n2me-1)       & ! up, up-left
            + b(n1me,n2me-1) - dble(8.)*b(n1me,n2me)  & ! me, left
            + b2s(n2me-1)+b2s(n2me))                    ! below, below-left
    end if

  end subroutine update_S_edge


  subroutine update_W_edge()
    ! Update interior points
    a(2:n1me-1,1) = b(2:n1me-1,1) + epsilon*( &
         + b1w(2:n1me-1)             & ! left-up
         + b1w(3:n1me)               & ! left
         + b1w(4:n1me+1)             & ! left-down
         + b(1:n1me-2,1)             & ! up
         - b(2:n1me-1,1)*dble(8.)    & ! yourself
         + b(3:n1me,  1)             & ! down
         + b(1:n1me-2,2)             & ! right-up
         + b(2:n1me-1,2)             & ! right
         + b(3:n1me,  2)             & ! right-down
         )

    ! Update top corner
    if (nid(4) > 0) a(1,1) = a(1,1) + epsilon*(b1w(1) + b1w(2) + b1w(3))

    ! Update bottom corner
    if (nid(2) > 0) a(n1me,1) = a(n1me,1) + epsilon*(b1w(n1me+2)+b1w(n1me+1)+b1w(n1me))

  end subroutine update_W_edge


  subroutine update_E_edge()
    ! Update interior points
    a(2:n1me-1,n2me) = b(2:n1me-1,n2me) + epsilon*( &
         + b3e(2:n1me-1)             & ! right-up
         + b3e(3:n1me)               & ! right
         + b3e(4:n1me+1)             & ! right-down
         + b(1:n1me-2,n2me)          & ! up
         - b(2:n1me-1,n2me)*dble(8.) & ! yourself
         + b(3:n1me,  n2me)          & ! down
         + b(1:n1me-2,n2me-1)        & ! left-up
         + b(2:n1me-1,n2me-1)        & ! left
         + b(3:n1me,  n2me-1)        & ! left-down
         )

    ! Update top corner
    if (nid(4) > 0) a(1,n2me) = a(1,n2me) + epsilon*(b3e(1) + b3e(2) + b3e(3))

    ! Update bottom corner
    if (nid(2) > 0) a(n1me,n2me) = a(n1me,n2me) + epsilon*(b3e(n1me+2)+b3e(n1me+1)+b3e(n1me))

  end subroutine update_E_edge


  subroutine update_interior()
    integer :: i, j

    do j = 2, n2me-1
       do i = 2, n1me-1
          a(i,j) = b(i,j) + epsilon*(                        &
               b(i-1,j+1) +          b(i,j+1) + b(i+1,j+1) + &
               b(i-1,j  ) - dble(8.)*b(i,j  ) + b(i+1,j  ) + &
               b(i-1,j-1) +          b(i,j-1) + b(i+1,j-1))
       end do
    end do

  end subroutine update_interior


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

    allocate(b1w (n1me+2))  ! +2 is for ghost cells for the 'corners'
    allocate(b2s (n2me))
    allocate(b3e (n1me+2))  ! +2 is for ghost cells for the 'corners'

    allocate(b4n(n2me))
    allocate(b1me(n1me+2))  ! +2 is for ghost cells for the 'corners'
    allocate(b2me(n2me))
    allocate(b3me(n1me+2))  ! +2 is for ghost cells for the 'corners'
    allocate(b4me(n2me))


    allocate( top(np2))
    allocate(left(np1))

  end subroutine allocate_arrays


  subroutine deallocate_arrays()
    deallocate(a,b,x,y,b1me,b2me,b3me,b4me,b1w,b2s,b3e,b4n)
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

    do i = 1, n1me
       x(i) = 1./dble(n1-1)*dble(i+i_f-2)
    end do
    
    do j = 1, n2me
       y(j) = 1./dble(n2-1)*dble(j+j_f-2)
    end do

    do j = 1, n2me
       do i = 1, n1me
          if (x(i) .lt. 0.5) then
             a(i,j) = cos(x(i)+y(j))
          else
             a(i,j) = sin(x(i)+y(j))
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


  subroutine read_input()
    character(20)  :: filename
    character(20)  :: arg
    integer :: i
    namelist /input_parameters/ n1, n2, np1, np2, niter

    if (iargc() == 0) then
       if (my_id == 1) then
          print *,"Error: Please specify a file name"
          print *, "Usage: "
          print *, "    mpirun -np <X> ./hw1_<i/ii/iii> -i FILENAME"
       end if
       stop
    else 
       do i=1,iargc()
          call getarg(i,arg)
          if (arg == '-i') then
             call getarg(i+1,filename)
             exit
          else
             filename = arg
          end if

       end do
    end if    
    
    open( unit=222, file=filename)
    read( unit=222, nml = input_parameters)
    close(unit=222)

    n1me = n1/np1
    n2me = n2/np2
  end subroutine read_input


end program main

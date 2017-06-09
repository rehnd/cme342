program main
  use iso_fortran_env
  implicit none
  integer                       :: n1 = 1024, n2 = 1024, np1=1, np2=1, niter = 1
  integer                       :: n1me, n2me
  integer                       :: nid(4)             ! neighbor id's
  integer                       :: np, ierr, status
  integer                       :: i, j, n
  integer                       :: i_f, i_l, j_f, j_l ! i_first, i_last, j_first, j_last
  real                          :: time
  double precision              :: norm
  double precision, parameter   :: epsilon = 0.1
  double precision, allocatable :: a(:,:)[:],b(:,:)[:],b1(:)[:],b2(:)[:],b3(:)[:],b4(:)[:]
  double precision, allocatable :: x(:)[:], y(:)[:]
  integer,          allocatable :: top(:)[:], left(:)[:] ! nodes at top & left of processor grid

  np = num_images()

  ! Read in input values and run tests
  call read_input()
  call check_num_processes()

  ! Start timer
  sync all
  
  ! Set up arrays on each processor, get the indices this process will work on
  call allocate_arrays()
  call get_my_ij()
  call initialize_arrays()
  nid = get_neighbor_ids(this_image())
  sync all
  
  ! Main loop. Inside, update the interior & edge points on each processor
  if (this_image() == 1) time = secnds(0.0)
  do n = 1, niter
     call update_edges()
     call update_interior()
     b=a
  end do

  sync all
  if (this_image() == 1) time = secnds(time)

  ! Compute the norm as a quick means of testing correctness
  call get_total_norm()

  if (this_image() == 1) then
     print *, "norm(a)   = ", norm
     write(*,"(A,F10.5,A)") " Wall time = ", time, " sec"
  end if

  call deallocate_arrays()


contains


  subroutine update_edges()

    ! Get N and S data
    if (nid(2) > 0) then
       b2 = b(1,:)[nid(2)] 
       call update_S_edge()
    end if
    if (nid(4) > 0) then 
       b4 = b(n1me,:)[nid(4)]
       call update_N_edge()
    end if

    sync all

    ! Get W and E data
    if (nid(3) > 0) then
       b3(2:n1me+1) = b(:,1)[nid(3)]
       call update_E_edge()
    end if
    if (nid(1) > 0) then
       b1(2:n1me+1) = b(:,n2me)[nid(1)]
       call update_W_edge()
    end if

  end subroutine update_edges

  
  subroutine update_N_edge()
    ! Immediately fill b1(1), b3(1)
    b1(1) = b4(1)
    b3(1) = b4(n2me)

    ! Perform complete update of the interior edge elements
    a(1,2:n2me-1) = b(1,2:n2me-1) + epsilon*( &
         + b4(2:n2me-1)           & ! above
         + b4(1:n2me-2)           & ! above-left
         + b4(3:n2me)             & ! above-right
         + b(1,1:n2me-2)           & ! left
         - b(1,2:n2me-1)*dble(8.)  & ! yourself 
         + b(1,3:n2me)             & ! right
         + b(2,1:n2me-2)           & ! below-left
         + b(2,2:n2me-1)           & ! below
         + b(2,3:n2me)             & ! below-right
         )

    ! Treat node corners as special case
    if (nid(1) > 0) then
       a(1,1) = b(1,1) + epsilon*(      &
            + b4(1) + b4(2)           & ! up, up-right
            + b(1,2) - dble(8.)*b(1,1)  & ! yourself, right
            + b(2,1) + b(2,2) )           ! below, below-right
    end if

    if (nid(3) > 0) then
       a(1,n2me) = b(1,n2me) + epsilon*(      &
            + b4(n2me-1) + b4(n2me)         & ! up, up-left
            + b(1,n2me-1) -dble(8.)*b(1,n2me) & ! yourself, left
            + b(2,n2me-1) + b(2,n2me))          ! below, below-left
    end if

  end subroutine update_N_edge
  

  subroutine update_S_edge()
    ! Immediately fill b1(end), b3(end) to account for corner values
    b1(n1me+2) = b2(1)
    b3(n1me+2) = b2(n2me)

    ! Perform a complete update of the interior edge elements |_|x|x|...|x|x|_|
    a(n1me,2:n2me-1) = b(n1me,2:n2me-1) + epsilon*( & 
         + b2(2:n2me-1)                            & ! below
         + b2(1:n2me-2)                            & ! below-left
         + b2(3:n2me)                              & ! below-right
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
            b2(2) + b2(1))                       ! below, below-right
    end if

    if (nid(3) > 0) then
       a(n1me,n2me) = b(n1me,n2me) + epsilon*(        &
            + b(n1me-1,n2me) + b(n1me-1,n2me-1)       & ! up, up-left
            + b(n1me,n2me-1) - dble(8.)*b(n1me,n2me)  & ! me, left
            + b2(n2me-1)+b2(n2me))                      ! below, below-left
    end if

  end subroutine update_S_edge


  subroutine update_W_edge()
    ! Update interior points
    a(2:n1me-1,1) = b(2:n1me-1,1) + epsilon*( &
         + b1(2:n1me-1)             & ! left-up
         + b1(3:n1me)               & ! left
         + b1(4:n1me+1)             & ! left-down
         + b(1:n1me-2,1)             & ! up
         - b(2:n1me-1,1)*dble(8.)    & ! yourself
         + b(3:n1me,  1)             & ! down
         + b(1:n1me-2,2)             & ! right-up
         + b(2:n1me-1,2)             & ! right
         + b(3:n1me,  2)             & ! right-down
         )

    ! Update top corner
    if (nid(4) > 0) a(1,1) = a(1,1) + epsilon*(b1(1) + b1(2) + b1(3))

    ! Update bottom corner
    if (nid(2) > 0) a(n1me,1) = a(n1me,1) + epsilon*(b1(n1me+2)+b1(n1me+1)+b1(n1me))

  end subroutine update_W_edge


  subroutine update_E_edge()
    ! Update interior points
    a(2:n1me-1,n2me) = b(2:n1me-1,n2me) + epsilon*( &
         + b3(2:n1me-1)             & ! right-up
         + b3(3:n1me)               & ! right
         + b3(4:n1me+1)             & ! right-down
         + b(1:n1me-2,n2me)          & ! up
         - b(2:n1me-1,n2me)*dble(8.) & ! yourself
         + b(3:n1me,  n2me)          & ! down
         + b(1:n1me-2,n2me-1)        & ! left-up
         + b(2:n1me-1,n2me-1)        & ! left
         + b(3:n1me,  n2me-1)        & ! left-down
         )

    ! Update top corner
    if (nid(4) > 0) a(1,n2me) = a(1,n2me) + epsilon*(b3(1) + b3(2) + b3(3))

    ! Update bottom corner
    if (nid(2) > 0) a(n1me,n2me) = a(n1me,n2me) + epsilon*(b3(n1me+2)+b3(n1me+1)+b3(n1me))

  end subroutine update_E_edge


  subroutine update_interior()
    integer :: i, j

    ! do concurrent(j = 2:n2me-1,i = 2:n1me-1)
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

    call co_sum(norm, result_image = 1)

    if (this_image() == 1)  norm = sqrt(norm)
    
  end subroutine get_total_norm


  subroutine allocate_arrays()
    allocate(a(n1me,n2me)[*])
    allocate(b(n1me,n2me)[*])
    allocate(x(n1me)[*])
    allocate(y(n2me)[*])

    allocate(b1(n1me+2)[*])  ! +2 is for ghost cells for the 'corners'
    allocate(b2(n2me)[*])
    allocate(b3(n1me+2)[*])  ! +2 is for ghost cells for the 'corners'
    allocate(b4(n2me)[*])

    allocate( top(np2)[*])
    allocate(left(np1)[*])

  end subroutine allocate_arrays


  subroutine deallocate_arrays()
    deallocate(a,b,x,y,b1,b2,b3,b4)
  end subroutine deallocate_arrays


  subroutine get_my_ij()
    integer                :: ip, jp
    character(len=10)      :: fmt = '(5I)'
    logical                :: printtest = .false.

    jp = (this_image()-1)/np1 + 1
    ip =  this_image() - (jp-1)*np1
    
    i_f = (n1/np1)*(ip-1) + 1
    i_l = i_f + (n1/np1)  - 1
    j_f = (n2/np2)*(jp-1) + 1
    j_l = j_f + (n2/np2)  - 1

    if (printtest) then
       if (this_image() == 1) write(*,'(5A12)') 'this_image()', 'i_f', 'i_l', 'j_f', 'j_l'
       write(*,fmt) this_image(), i_f, i_l, j_f, j_l
    end if
  end subroutine get_my_ij
  

  subroutine check_num_processes()
    if ((.not. np1*np2  == np)) then
       if (this_image() == 1) print *, 'Error: Domain decomp (np1*np2) not equal to MPI_RANK'
       stop
    end if
  end subroutine check_num_processes


  pure function get_neighbor_ids(id) 
    integer, intent(in)  :: id
    integer              :: ip1, ip2   ! indices of this_image() in the global array
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
       if (this_image() == 1) write(*,'(3A12)') 'this_image()', 'x(1)', 'y(1)'
       write(*,fmt) this_image(), x(1), y(1)
    end if

  end subroutine initialize_arrays


  subroutine read_input()
    character(20)  :: filename
    character(20)  :: arg
    integer :: i
    namelist /input_parameters/ n1, n2, np1, np2, niter

    if (iargc() == 0) then
       if (this_image() == 1) then
          print *,"Error: Please specify a file name"
          print *, "Usage: "
          print *, "    cafrun -np <X> ./coarray -i FILENAME"
       end if
       error stop
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

program main
  implicit none
  integer                       :: n1 = 1024, n2 = 1024, np1=1, np2=1, niter = 1
  integer                       :: n1me, n2me
  integer                       :: i, j, n, np
  integer                       :: nid(4)             ! neighbor id's
  integer                       :: i_f, i_l, j_f, j_l,  if_, il_, jf_, jl_
  real                          :: time
  double precision              :: norm
  double precision, parameter   :: epsilon = 0.1
  double precision, allocatable :: a(:,:)[:], b(:,:)[:], x(:)[:], y(:)[:]

  ! Read input values and run tests
  call read_input()
  call check_num_processes()

  ! Start timer
  sync all

  
  ! Set up arrays on each processor, get the indices this process will work on
  call allocate_arrays()
  call get_my_ij()
  call initialize_arrays()
  nid = get_neighbor_ids(this_image())
  call get_ifl_jfl()

  sync all
  if (this_image() == 1) time = secnds(0.0)

  ! Main loop. Inside, update the interior & edge points on each processor
  do n = 1, niter
     call send_edges()
     call update_interior()
     b=a
     sync all
  end do

  if (this_image() == 1) time = secnds(time)

  ! Compute the norm as a quick means of testing correctness
  call get_total_norm()

  if (this_image() == 1) then
     print *, "norm(a)   = ", norm
     write(*,"(A,F10.5,A)") " Wall time = ", time, " sec"
  end if

  call deallocate_arrays()


contains


  subroutine send_edges()

    ! Get N and S data
    if (nid(2) > 0)  b(n1me+2,:) = b(2,:)[nid(2)]
    if (nid(4) > 0)  b(1,:)      = b(n1me+1,:)[nid(4)]

    sync all

    ! Get W and E data
    if (nid(3) > 0)  b(:,n2me+2) = b(:,2)[nid(3)]
    if (nid(1) > 0)  b(:,1)      = b(:,n2me+1)[nid(1)]

    sync all
  end subroutine send_edges

  subroutine update_interior()
    integer :: i, j

    do concurrent(j = jf_:jl_)
       do concurrent(i=if_:il_)
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
    norm = sum(a(2:n1me+1,2:n2me+1)**2)

    call co_sum(norm, result_image = 1)

    if (this_image() == 1)  norm = sqrt(norm)
    
  end subroutine get_total_norm


  subroutine allocate_arrays()
    allocate(a(n1me+2,n2me+2)[*])
    allocate(b(n1me+2,n2me+2)[*])
    allocate(x(n1me)[*])
    allocate(y(n2me)[*])
  end subroutine allocate_arrays


  subroutine deallocate_arrays()
    deallocate(a,b,x,y)
  end subroutine deallocate_arrays

  subroutine get_ifl_jfl()
    logical :: test = .false.
    if_ = 2
    il_ = n1me+1
    jf_ = 2
    jl_ = n2me+1
    if (nid(1) < 0) jf_ = jf_ + 1
    if (nid(2) < 0) il_ = il_ - 1
    if (nid(3) < 0) jl_ = jl_ - 1
    if (nid(4) < 0) if_ = if_ + 1

    if(test) then
       print *, this_image(), if_, il_, jf_,jl_
    end if

  end subroutine get_ifl_jfl
  

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
    if ((.not. np1*np2  == num_images())) then
       if (this_image() == 1) print *, 'Error: Domain decomp (np1*np2) not equal to MPI_RANK'
       error stop
    end if
  end subroutine check_num_processes


  pure function get_neighbor_ids(id) 
    integer, intent(in)  :: id
    integer              :: ip1, ip2   ! indices of this_image() in the global array
    integer              :: get_neighbor_ids(4)

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

  end function get_neighbor_ids

  
  subroutine initialize_arrays()

    ! do concurrent(i = 1:n1me)
    do concurrent(i = 1:n1me)
       x(i) = 1./dble(n1-1)*dble(i+i_f-2)
    end do
    
    do concurrent(j=1:n2me)
       y(j) = 1./dble(n2-1)*dble(j+j_f-2)
    end do

    do concurrent(j=1:n2me)
       do concurrent(i=1:n1me)
          if (x(i) .lt. 0.5) then
             a(i+1,j+1) = cos(x(i)+y(j))
          else
             a(i+1,j+1) = sin(x(i)+y(j))
          end if
       end do
    end do
    
    b = a

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

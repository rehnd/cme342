program main
  use omp_lib
  implicit none
  integer                       :: n1 = 1024, n2 = 1024, np1=1, np2=1, niter = 1
  integer                       :: n1me, n2me, nthreads, my_id, np
  integer                       :: nid(4)             ! neighbor id's
  integer                       :: i, j, n
  integer                       :: i_f, i_l, j_f, j_l, i_ff, i_ll, j_ff,j_ll
  double precision              :: start, end
  double precision, parameter   :: epsilon = 0.1
  double precision, allocatable :: a(:,:),b(:,:)

  double precision, allocatable :: x(:), y(:)

  ! Read in input values and run tests
  call read_input()
  call omp_set_num_threads(nthreads)

  np = nthreads
  call check_num_processes()
  call allocate_arrays()

  !start = omp_get_wtime()
  call cpu_time(start)

  !$omp parallel shared(a,b,x,y) private(i,j, i_f,i_l,j_f,j_l, i_ff,i_ll,j_ff,j_ll, my_id) 
  my_id = omp_get_thread_num() + 1

  call get_my_ij()
  call set_ij_interior()

  call initialize_arrays()
  !$omp barrier

  ! Main loop. Inside, update the interior & edge points on each processor
  do n = 1, niter
     call update_interior()
     b(i_f:i_l, j_f:j_l) = a(i_f:i_l, j_f:j_l)
     !$omp barrier
  end do
  
  !$omp end parallel
  !  endt = omp_get_wtime()
  call cpu_time(end)
  print *, "norm(a)    = ", norm2(a)
  print *, "a(512,512) = ", a(512,512)
  print *, "time       = ", (end - start)


  call deallocate_arrays()


contains


  subroutine set_ij_interior()
    logical :: printtest = .false.
    character(len=10)      :: fmt = '(5I)'
    i_ff = i_f
    i_ll = i_l
    j_ff = j_f
    j_ll = j_l
    if (i_f == 1)  i_ff = 2
    if (i_l == n1) i_ll = n1-1

    if (j_f == 1)  j_ff = 2
    if (j_l == n2) j_ll = n2-1

    if (printtest) then
       if (my_id == 1) write(*,'(5A12)') 'my_id', 'i_ff', 'i_ll', 'j_ff', 'j_ll'
       write(*,fmt) my_id, i_ff, i_ll, j_ff, j_ll
    end if

  end subroutine set_ij_interior

  subroutine update_interior()
    integer :: i,j
    
    do j = j_ff, j_ll
       do i = i_ff, i_ll
          a(i,j) = b(i,j) + epsilon*(                        &
               b(i-1,j+1) +          b(i,j+1) + b(i+1,j+1) + &
               b(i-1,j  ) - dble(8.)*b(i,j  ) + b(i+1,j  ) + &
               b(i-1,j-1) +          b(i,j-1) + b(i+1,j-1))
       end do
    end do

  end subroutine update_interior


  subroutine allocate_arrays()
    allocate(a(n1,n2))
    allocate(b(n1,n2))
    allocate(x(n1))
    allocate(y(n2))
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
    if ((.not. np1*np2  == nthreads)) then
       print *, 'Error: Domain decomp (np1*np2) not equal to MPI_RANK'
       stop
    end if
  end subroutine check_num_processes

  
  subroutine initialize_arrays()

    do i = i_f,i_l
       x(i) = 1./dble(n1-1)*dble(i-1)
    end do
    
    do j = j_f, j_l
       y(j) = 1./dble(n2-1)*dble(j-1)
    end do
    
    do j = j_f, j_l
       do i = i_f, i_l
          if (x(i) .lt. 0.5) then
             a(i,j) = cos(x(i)+y(j))
          else
             a(i,j) = sin(x(i)+y(j))
          end if
          b(i,j) = a(i,j)
       end do
    end do

  end subroutine initialize_arrays


  subroutine read_input()
    character(20)  :: filename
    character(20)  :: arg
    integer :: i
    namelist /input_parameters/ n1, n2, np1, np2, niter, nthreads

    if (iargc() == 0) then
       print *,"Error: Please specify a file name"
       print *, "Usage: "
       print *, "   ./hw3_2 -i FILENAME"
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

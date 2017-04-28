program main
  implicit none

  integer              :: n1 = 1024, n2 = 1024, niter = 1, np1, np2
  integer              :: i, j, n
  real(8), parameter   :: epsilon = 0.1
  real(8), allocatable :: a(:,:), b(:,:)
  real(8), allocatable :: x(:)  , y(:)

  call read_input()

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
  ! print *, 1,a(512,512)
  ! print *, 2,a(513,512)
  ! print *, 3,a(512,513)
  ! print *, 4,a(513,513)
  
  deallocate(a,b,x,y)

contains

  subroutine read_input()
    character(20)  :: filename
    character(20)  :: arg
    integer :: i
    namelist /input_parameters/ n1, n2, np1, np2, niter

    np1 = 1
    np2 = 1
    ! Reset these for serial case

    if (iargc() == 0) then
       print *,"Error: Please specify a file name"
       print *, "Usage: "
       print *, "    ./serial -i FILENAME"
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

  end subroutine read_input

end program main

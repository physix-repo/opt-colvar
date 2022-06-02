!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      program debug
      !
      implicit none
      !!!!!!!!!!!!!!Variable Declaration!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision, allocatable :: t(:),x_init(:),cv(:,:)
      double precision, allocatable :: x_new(:),cv_temp(:)
      double precision, allocatable :: F(:),G(:),D(:),Q(:)
      integer n,i,io,j,k,o,clock,label,m,error,col
      double precision t_tmp
      double precision r,factor
      double precision col1,col2,col3,col4,col5
      double precision tau_new,tau_old
      double precision, allocatable :: Integral(:),z1(:),z2(:)
      double precision a,b,q0,dq,dtau,Fmin,tfinal
      double precision, allocatable ::s(:),w(:),dw(:)
      character str
      character(len=80) :: string
      real :: kBT,val
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !read colvar file and get number of columns (m) and lines (n)  
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      n=0
      o=0
      open(unit=3,file='colvar_x1')
        do
                read(3,*,iostat=io) !read the columns in colvar file
                if (io/=0) exit
                n=n+1
        end do
      print *, 'Number of lines in colvar is',n
      !
      close(3)
      open(unit=3,file='colvar_x1')
      read (3,'(a)') string
        do i =1,40   ! The very maximum that the string can contain
                read( string, *, iostat=error ) ( val, k=1,i )
                if ( error .ne. 0 ) then
                m = i - 1
                exit
                endif
        enddo
      close(3)
      print *, 'Number of columns in colvar file is',m
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !!!!!!!!!!!store and normalize colvar file!!!!!!!!!!!!!!!!! 
      ! 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      allocate(cv(n,m),t(n))
      open(unit=3,file='colvar_x1')
        do i = 1,n
                read(3,*) (cv(i,col),col=1,m) !read columns in colvar
        enddo
      close(3)
      t=cv(:,1)!store time in separate array
      tfinal=t(n)
      !
        do i = 2,m !Let's scale the colvar file between 0 and 1
                cv(:,i)= (cv(:,i) - minval(cv(:,i))) / (maxval(cv(:,i)-minval(cv(:,i))))
        end do
      !Store normalized cv in file
      open(unit=88,file='normalized_cv')
      do i=1,n
        if (t(i)==tfinal) then
                if (cv(i,2) .gt. cv(1,2)) then
                        write(88,*) (cv(i,col),col=1,m),2
                else
                        write(88,*) (cv(i,col),col=1,m),1
                endif
        else
                write(88,*) (cv(i,col),col=1,m)

        end if 
      end do
      print *, cv(1,2)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !INITIALIZATION
      CALL init_random_seed()  !generate random seed
      allocate(w(m-1))
      allocate(dw(m-1))
      allocate(s(m-1))
        do i=1,m-1
                call random_number(w(i)) !generate random weights
        end do
      tau_old=0
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !Suggest a move
      allocate(F(1000),G(1000),D(1000),Q(1000),Integral(1000))
        do j=1,1000 !MAIN LOOP!!!!!!!!!!!!!!!!!!
              CALL create_colvar(n,m,w,cv,t,dw,a,b,q0,s)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !run the opt LE code to get F and D of colvar
                call execute_command_line ('./optle')
                open(unit=11,file='PROFILES')
                open(unit=66,file='PROFILES-TEMP')
                read(11,*)
                do i=1,1000
                        read(11,*) col1,col2,col3,col4,col5
                        Q(i)=col1
                        F(i)=col3
                        G(i)=col4
                        D(i)=1/col4
                        write(66,*) Q(i),F(i),D(i)
                end do
                print*,F(1000),Q(1000),G(1000)
                close(11)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !calculate mfpt from integral
                Integral(1)=0
                do k=2,1000 !COMPUTING INNER INTEGRAL
                        if (Q(k)>a) then
                                Integral(k)=Integral(k-1)+(Q(k)-Q(k-1))*exp(-F(k))
                        else
                                Integral(k)=0
                        endif
                enddo
                tau_new=0
                do k=2,1000 !COMPUTING OUTER INTEGRAL
                        if (Q(k)>q0 .and. Q(k)<b) then 
                                dq=Q(k)-Q(k-1)
                                dtau=dq*(exp(F(k))/D(k))*Integral(k)
                                tau_new=tau_new+dtau
                        endif
                end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Accept tau with certain probability METROPOLIS
                !kBT=0.1
                call random_number(r)
                open(unit=13,file='tau')
                if (tau_old<=tau_new) then
                        write(13,*) j,s,tau_new,'A',sqrt(sum(s**2))
                        tau_old=tau_new
                        w=s !update weights
                elseif ((tau_new/tau_old)>r) then
                        write(13,*)j,s,tau_new,'B',sqrt(sum(s**2))
                        tau_old=tau_new
                        w=s !update weights
                endif
                print*,tau_old
        enddo
                close(13)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        contains
        !
        SUBROUTINE create_colvar(n,m,w,cv,t,dw,a,b,q0,s)
        INTEGER :: i,n,m,k,col
        DOUBLE PRECISION factor1,a,b,q0
        DOUBLE PRECISION, ALLOCATABLE :: w(:),cv(:,:),dw(:),s(:),x_new(:)
        DOUBLE PRECISION, ALLOCATABLE :: t(:),z1(:),z2(:)
        CALL INIT_RANDOM_SEED()
        allocate(x_new(n))
        allocate(z1(0),z2(0))
        do i=1,m-1
                CALL random_number(dw(i))
        enddo
        dw=-0.01+(0.01-(-0.01))*dw !random increment between -0.01 and 0.01
        s=w+dw
        factor1=1/sqrt(dot_product(s,s)) ! XXX check
        open(unit=10,file='colvar')
        do i=1,n
                x_new(i)=factor1*sum(s*cv(i,2:m))  
                if (t(i) .eq. tfinal ) then !label colvar file
                        if (cv(i,2) .gt. cv(1,2)) then !shootings start from barrier top
                                write(10,*) t(i),x_new(i),2
                                z2=[z2,x_new(i)]
                        else
                                write(10,*) t(i),x_new(i), 1
                                z1=[z1,x_new(i)]
                                end if
                        else
                                write(10,*) t(i), x_new(i)
                        end if
        enddo
        close(10)
        s=factor1*s
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Find integral bounds and Free energy minimum a,b,q0 to calculate
        !mfpt from integral 
                a=minval(x_new) ! reflecting bound is lowest q value
                b=sum(z2)/size(z2) !average value of second minima absorbing bound
                q0=sum(z1)/size(z1) !average value of first minima
                print*, 'integral bounds are=',a,b,q0
        RETURN
        END SUBROUTINE
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!
       SUBROUTINE init_random_seed() !Generate random seed
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
            DOUBLE PRECISION dw
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))

            CALL SYSTEM_CLOCK(COUNT=clock)

            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)
            DEALLOCATE(seed)
       END SUBROUTINE
      !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!
      end program debug
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

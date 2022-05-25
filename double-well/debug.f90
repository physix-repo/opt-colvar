!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      program debug
      !
      implicit none
      !!!!!!!!!!!!!!Variable Declaration!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision, allocatable :: t(:),x(:),y(:)
      double precision, allocatable :: x_new(:)
      double precision, allocatable :: F(:),G(:),D(:),Q(:)
      integer n,i,io,j,k,o,clock,label
      double precision t_tmp,x_tmp,y_tmp
      double precision w,dw,r,w1,w2
      double precision col1,col2,col3,col4,col5
      double precision tau_new,tau_old
      double precision, allocatable :: Integral(:),z1(:),z2(:)
      double precision a,b,q0,dq,dtau,Fmin
      integer ub,lb,i_min,temp(1)
      character str
      real kBT
      !
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call random_seed()
      !
      !read and store input and get number of lines in file n
      n=0
      o=0
      open(unit=3,file='colvar-short')
        do
                read(3,*,iostat=io) !read the columns in colvar file
                if (io/=0) exit
                n=n+1
        end do
      close(3)
      print *, 'Number of lines in colvar is',n
      allocate(t(n),x(n),y(n))
      open(unit=3,file='colvar-short')
        do i=1,n
                read(3,*) t_tmp,x_tmp,y_tmp !read the columns in colvar file
                t(i)=t_tmp
                x(i)=x_tmp
                y(i)=y_tmp
        end do
      print *, t(n),x(n),y(n)
      close(3)
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !INITIALIZATION
        call init_random_seed()  !generate random number
        call random_number(w)
        tau_old=1 !initialization
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !Suggest a move
      allocate(x_new(n))
      allocate(F(1000),G(1000),D(1000),Q(1000),Integral(1000))
      allocate(z1(0),z2(0)) !we have 100 shootings
      Print *, n
        do j=1,1000 !MAIN LOOP!!!!!!!!!!!!!!!!!!
                CALL RANDOM_NUMBER(dw)
                dw=-0.01+(0.01-(-0.01))*dw
                w1=w+dw
                w2=1-w1 !normalization
                !keeps weights between 0 and 1
                if (w1.lt.0) then
                        w1=w1*(-1)**cmplx(0,1)
                endif
                if (w1.gt.1) then
                        w1=1-(w1-1)
                endif
                if (w2.lt.0) then
                        w2=w2*(-1)**cmplx(0,1)
                endif
                if (w2.gt.1) then
                        w2=1-(w2-1)
                endif

                open(unit=10,file='colvar')
                do i=1,n
                        x_new(i)=(w1)*x(i)+(w2)*y(i)
                        if (t(i) .eq. 3 ) then !label colvar file
                                if (x_new(i) .gt. 1.7) then 
                                       write(10,*) t(i),x_new(i),y(i), 2
                                       z2=[z2,x_new(i)]
                                else
                                       write(10,*) t(i),x_new(i),y(i), 1
                                       z1=[z1,x_new(i)]
                                end if
                        else
                                write(10,*) t(i), x_new(i), y(i)
                        end if
                end do
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
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Find integral bounds and Free energy minimum a,b,q0 to calculate
      !mfpt from integral 
                a=minval(x_new) ! reflecting bound is lowest q value
                b=sum(z2)/size(z2) !average value of second minima absorbing bound
                q0=sum(z1)/size(z1) !average value of first minima
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
                !kBT=0.3
                call random_number(r)
                open(unit=13,file='tau')
                if (tau_old<=tau_new) then
                        write(13,*) j,w+dw,tau_new,'A'
                        tau_old=tau_new
                        w=w+dw
                elseif ((tau_new/tau_old)>r) then
                        write(13,*) j,w+dw,tau_new,'B'
                        tau_old=tau_new
                        w=w+dw
                endif
                print*,tau_old
                close(10)
                if (w .gt. 1) then
                        exit
                endif
        enddo
                close(13)
      end program debug
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

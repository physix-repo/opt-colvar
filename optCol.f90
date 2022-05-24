!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      program debug
      !
      implicit none
      !!!!!!!!!!!!!!Variable Declaration!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision, allocatable :: t(:),x_init(:)
      double precision, allocatable :: cv1(:),cv2(:),cv3(:)
      double precision, allocatable :: x_new(:)
      double precision, allocatable :: F(:),G(:),D(:),Q(:)
      integer n,i,io,j,k,o,clock,label
      double precision t_tmp,cv1_tmp,cv2_tmp,cv3_tmp
      double precision w1,w2,w3,dw1,dw2,dw3,r,factor
      double precision col1,col2,col3,col4,col5
      double precision tau_new,tau_old
      double precision, allocatable :: Integral(:),z1(:),z2(:)
      double precision a,b,q0,dq,dtau,Fmin
      character str
      real kBT
      !
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !read and store input and get number of lines in file n
      n=0
      o=0
      open(unit=3,file='colvar_x1')
        do
                read(3,*,iostat=io) !read the columns in colvar file
                if (io/=0) exit
                n=n+1
        end do
      close(3)
      print *, 'Number of lines in colvar is',n
      !
      allocate(t(n),cv1(n),cv2(n),cv3(n),x_init(n))
      open(unit=3,file='colvar_x1')
        do i=1,n
                read(3,*) t_tmp,cv1_tmp,cv2_tmp,cv3_tmp !read the columns in colvar file
                t(i)=t_tmp
                cv1(i)=cv1_tmp
                cv2(i)=cv2_tmp
                cv3(i)=cv3_tmp
        end do
      print *, t(n),cv1(n),cv2(n),cv3(n)
      close(3)
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Let's scale the colvar file between 0 and 1
      cv1= (cv1 - minval(cv1)) / (maxval(cv1)- minval(cv1))
      cv2= (cv2 - minval(cv2)) / (maxval(cv2)- minval(cv2))
      cv3= (cv3 - minval(cv3)) / (maxval(cv3)- minval(cv3))
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Store normalized cv in file
      open(unit=88,file='normalized_cv')
      do i=1,n
        if (t(i)==100) then
                if (cv1(i) .gt. 0.5) then
                        write(88,*) t(i),cv1(i),cv2(i),cv3(i),2
                else
                        write(88,*) t(i),cv1(i),cv2(i),cv3(i),1
                endif
        else
                write(88,*) t(i),cv1(i),cv2(i),cv3(i)

        end if 
      end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !INITIALIZATION
      CALL init_random_seed()  !generate random seed
      CALL random_number(w1)
      CALL random_number(w2)
      CALL random_number(w3)
      tau_old=0
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      !Suggest a move
      allocate(F(1000),G(1000),D(1000),Q(1000),Integral(1000))
        do j=1,1000 !MAIN LOOP!!!!!!!!!!!!!!!!!!
              CALL create_colvar(n,w1,w2,w3,cv1,cv2,cv3,t,dw1,dw2,dw3,a,b,q0)
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
                        write(13,*) j,w1+dw1,w2+dw2,w3+dw3,tau_new,'A',w1+dw1+w2+dw2+w3+dw3
                        tau_old=tau_new
                        w1=w1+dw1
                        w2=w2+dw2
                        w3=w3+dw3
                elseif ((tau_new/tau_old)<r) then
                        write(13,*)j,w1+dw1,w2+dw2,w3+dw3,tau_new,'B',w1+dw1+w2+dw2+w3+dw3
                        tau_old=tau_new
                        w1=w1+dw1
                        w2=w2+dw2
                        w3=w3+dw3
                endif
                print*,tau_old
        enddo
                close(13)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        contains
        !
        SUBROUTINE create_colvar(n,w1,w2,w3,cv1,cv2,cv3,t,dw1,dw2,dw3,a,b,q0)
        INTEGER :: i,n
        DOUBLE PRECISION dw1,dw2,dw3,w1,w2,w3,factor1,a,b,q0
        DOUBLE PRECISION, ALLOCATABLE :: cv1(:),cv2(:),cv3(:),x_new(:)
        DOUBLE PRECISION, ALLOCATABLE ::t(:),z1(:),z2(:)
        CALL INIT_RANDOM_SEED()
        CALL RANDOM_NUMBER(dw1)
        CALL RANDOM_NUMBER(dw2)
        CALL RANDOM_NUMBER(dw3)
        allocate(x_new(n))
        allocate(z1(0),z2(0))
        dw1=-0.01+(0.01-(-0.01))*dw1 !random increment between -0.01 and 0.01
        dw2=-0.01+(0.01-(-0.01))*dw2
        dw3=-0.01+(0.01-(-0.01))*dw3
        !factor1=1/sqrt((w1+dw1)**2+(w2+dw2)**2+(w3+dw3)**2)
        factor1=1/(w1+w2+w3+dw1+dw2+dw3)
        open(unit=10,file='colvar')
        do i=1,n
                x_new(i)=factor1*((w1+dw1)*cv1(i)+(w2+dw2)*cv2(i)+(w3+dw3)*cv3(i))
                if (t(i) .eq. 100 ) then !label colvar file
                        if (cv1(i) .gt. 0.5) then
                                write(10,*) t(i),x_new(i),2
                                z2=[z2,x_new(i)]
                        else
                                write(10,*) t(i),x_new(i), 1
                                z1=[z1,x_new(i)]
                                end if
                        else
                                write(10,*) t(i), x_new(i)
                        end if
                end do
        close(10)
        w1=factor1*w1
        w2=factor1*w2
        w3=factor1*w3
        dw1=factor1*dw1
        dw2=factor1*dw2
        dw3=factor1*dw3
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


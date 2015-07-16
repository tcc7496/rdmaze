
          implicit real*8 (a-h,o-z)
          parameter (nx=210)
          parameter (ny=210)

         real eps,dt,dx,Diff,delta,a,b,gamma,dlap
         real u,v
         real ut
         real c

         real time
         integer  ntime
c--       dimensions for the 2 variables
         dimension u(0:nx+1, 0:ny+1),v(0:nx+1, 0:ny+1)
         dimension ut(0:nx+1, 0:ny+1)
         dimension phi(0:nx+1, 0:ny+1)
         dimension c(0:nx+1, 0:ny+1)

c--        fixed parameters
          nf=19
           eps=0.05
c--------- eps = 0.08
           t=0
           gamma=0.8
           beta=.7
           dt=0.005
           dx=0.25
           Diff=1.
           dlap=Diff*dt/(dx*dx)
           ntime=15000


c---------------- Initial Conditions (rest state) ---------------------------------
        do i=0,nx+1
        do j=0,ny+1
        u(i,j)=-1.199
        ut(i,j)=u(i,j)
        phi(i,j)=1
        c(i,j)=4
        enddo
        enddo

c------ blocks/obstacles
c
c       do i=20,30
c       do j=20,40
c       phi(i,j)=0
c       ut(i,j) = 3
c       enddo
c       enddo
c
c       do i=25,35
c       do j=50,70
c       phi(i,j)=0
c       ut(i,j) = 3
c       enddo
c       enddo
c
c       do i=10,60
c       do j=50,60
c       phi(i,j)=0
c       enddo
c       enddo

c---------------Initialise perturbation----------------

        do i=1,15
        do j=1,15
         u(i,j)=1
        enddo
        enddo

        do i=1,nx+1
        do j=1,ny+1
         v(i,j)=-0.6242
        enddo
        enddo



c--      Iextt=1.90   ! external current to produce AP
c--      Opening up the wall points to get boundary

c       open(unit=10, file="wallPoints.txt")
c       do k=1,10000
c       read(10,*,end=10)i,j
c       phi(i,j)=0
c       enddo
c10      continue

c-------colour assignment of

c----------------------------------------------------
c           time integration
c---------------------------------------------------
         do nt=0,ntime
c--        Iext=0
c-----  this injects current
c--      if(mod(nt,40000).eq.1)Iext=Iextt

c--------- updating boundary conditions for zero flux
         do i=1,nx
         do j=1,ny
         u(i,0)=u(i,2)
         u(i,ny+1)=u(i,ny-1)
         u(0,j)=u(2,j)
         u(nx+1,j)=u(nx-1,j)
         enddo
         enddo
c-------- updating boundary conditions to rest state
c           do i=1,nx
c           do j=1,ny
c           u(i,0)=-1.199
c           u(i,ny+1)=-1.199
c           u(0,j)=-1.199
c           u(nx+1,j)=-1.199
c           enddo
c           enddo

c---updating boundary conditions for first derivative 0
c          do i=1,nx
c          do j=1,ny
c          u(i,0)=u(i,1)
c          u(i,ny+1)=u(i,ny)
c          u(0,j)=u(1,j)
c          u(nx+1,j)=u(nx,j)
c          enddo
c          enddo

c------- obstacle boundary conditions
c       do i=20,30
c       do j=20,40
c       u(i,20) = u(i,18)
c       u(i,40)=u(i,42)
c       u(20,j)=u(18,j)
c       u(30,j)=u(32,j)
c       enddo
c       enddo
c
c       do i=25,35
c       do j=50,70
c       u(i,50) = u(i,48)
c       u(i,70)=u(i,72)
c       u(25,j)=u(23,j)
c       u(35,j)=u(37,j)
c       enddo
c       enddo

c---------integration in space
        do i=1,nx
        do j=1,ny

        if(phi(i,j).eq.1)then
         v(i,j)=v(i,j)+eps*(u(i,j)-gamma*v(i,j)+beta)*dt
c--------updating the voltage to ut
         xlap=u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4*u(i,j)

         ut(i,j)=u(i,j)+(u(i,j)-u(i,j)**3/3.-v(i,j))*dt+xlap*dlap

        endif
        enddo
        enddo

c------- update u value
         do i=1,nx
         do j=1,ny
         u(i,j)=ut(i,j)
         if(mod(nt,100).eq.1)then
c         if(u(i,j).gt.-0.5)then
          write(nf,*)i,j,u(i,j),v(i,j),c(i,j)
c         endif
         endif
         enddo
         enddo
         if(mod(nt,100).eq.1)then
           close(nf)
c           write(6,*)nf
           nf=nf+1
         endif

         t=t+dt
        if(mod(nt,20).eq.1)write(16,*)t,u(5,5),u(50,50),u(90,90)
        if(mod(nt,20).eq.1)write(17,*)t,v(5,5),v(50,50),v(90,90)

         enddo
         end


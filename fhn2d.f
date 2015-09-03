
          implicit real*8 (a-h,o-z)
          parameter (nx = 100)
          parameter (ny = 100)

         real eps,dt,dx,Diff,delta,a,b,gamma,dlap
         real u,v
         real ut

         real time
         integer  ntime
c        dimensions for the 2 variables:
c        indices inclusive, default in 1 if m omitted in m:n
         dimension u(0:nx+1, 0:ny+1),v(nx,ny)
         dimension ut(0:nx+1, 0:ny+1)

c--        fixed parameters
          nf=19
           eps=0.0005
c--------- eps = 0.08
           t=0
           gamma=0.8
           beta=.7
           dt=0.005
           dx=0.25
           Diff=1.
           dlap=Diff*dt/(dx*dx)
           ntime=75000
c           ntime=10000

c---------------- Initial Conditions (rest state) ---------------------------------
        do i=0,nx+1
          do j=0,ny+1
          u(i,j)=-1.199
          ut(i,j)=u(i,j)
          enddo
        enddo
        
        do i=0,10
        do j=0,10
         u(i,j)=1
        enddo
        enddo

        do i=1,nx+1
        do j=1,ny+1
         v(i,j)=-0.6242
        enddo
        enddo
c       Iextt=1.90   ! external current to produce AP


c----------------------------------------------------
c           time integration
c---------------------------------------------------
         do nt=0,ntime
c         Iext=0
c-----  this injects current
c        if(mod(nt,40000).eq.1)Iext=Iextt

c--------- updating boundary conditions for zero flux
         do i=1,ny
            u(i,0) = u(i,2)
            u(i,nx+1)=u(i,nx-1)
         enddo
         do j = 1,nx
            u(0,j) = u(2,j)
            u(ny+1,j) = u(ny-1,j)
         enddo
c---------integration in space
        do i=1,nx
        do j=1,ny

        
         v(i,j)=v(i,j)+eps*(u(i,j)-gamma*v(i,j)+beta)*dt
c--------updating the voltage to ut
         xlap=u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4*u(i,j)

         ut(i,j)=u(i,j)+(u(i,j)-u(i,j)**3/3.-v(i,j))*dt+xlap*dlap

       
        enddo
        enddo         
cc---------integration in space
c        do i=1,nx
c         do j=1,ny
c         v(i,j)=v(i,j)+eps*(u(i,j)-gamma*v(i,j)+beta)*dt
cc--------updating the voltage to ut
c         xlap=u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4.*u(i,j)
c
c         ut(i,j)=u(i,j)+(u(i,j)-u(i,j)**3/3.-v(i,j))*dt+xlap*dlap
c         enddo
c        enddo

c------- update u value
         do i=1,nx
           do j=1,ny
         u(i,j)=ut(i,j)
         if(mod(nt,500).eq.1)then
           if(u(i,j).gt.0)then
             write(nf,*)i,j,u(i,j),v(i,j)
           endif
         endif
           enddo
         enddo
         if(mod(nt,500).eq.1)then
           close(nf)
           write(6,*)nf
           nf=nf+1
         endif

         t=t+dt
        if(mod(nt,20).eq.1)write(16,*)t,u(5,5),u(50,50),u(90,90)
        if(mod(nt,20).eq.1)write(17,*)t,v(5,5),v(50,50),v(90,90)

         enddo
         end


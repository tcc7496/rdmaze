
          implicit real*8 (a-h,o-z)
          parameter (nx=209)
          parameter (ny=209)

         real eps,dt,dx,Diff,delta,a,b,gamma,dlap
         real u,v
         real ut
         real pinvalue

         real time
         integer  ntime, n, m
c--       dimensions for the 2 variables
         dimension u(0:nx+1, 0:ny+1),v(0:nx+1, 0:ny+1)
         dimension ut(0:nx+1, 0:ny+1)
         dimension phi(0:nx+1, 0:ny+1)

c--------------------------------------------------------------------------------
c             Diffusion
c--------------------------------------------------------------------------------
c--        fixed parameters
          nf=19
           eps=0.0005
c--------- eps = 0.08, 0.0005
           t=0
           gamma=1.8
           beta=-0.2
           dt=0.005
           dx=0.25
           Diff=1.
           dlap=Diff*dt/(dx*dx)
           ntime=75000
c           ntime=10000


c---------------- Initial Conditions (rest state) ---------------------------------
        do i=0,nx+1
        do j=0,ny+1
        u(i,j)=-0.99999
        ut(i,j)=u(i,j)
        phi(i,j)=1
        enddo
        enddo
c -1.199
c------ blocks/obstacles

c       do i=20,30
c       do j=20,40
c       phi(i,j)=0
c       enddo
c       enddo

c       do i=25,35
c       do j=50,70
c       phi(i,j)=0
c       enddo
c       enddo

c       do i=10,60
c       do j=50,60
c       phi(i,j)=0
c       enddo
c       enddo

        do i=1,nx+1
        do j=1,ny+1
         v(i,j)=-0.66667
        enddo
        enddo
c -0.6242

c------ initial excitation
c------- 47/62
c       do i=0,10
c       do j=0,10
c        u(i,j)=1.26376
c        v(i,j)=0.590979
c       enddo
c       enddo

c------- exciting the end of the maze        
c       do i=nx-15,nx
c       do j=ny-15,ny
c        u(i,j)=1.26376
c        v(i,j)=0.590979
c       enddo
c       enddo        

c------- exciting the top left corner of the maze        
       do i=10,20
       do j=ny-20,ny-10
        u(i,j)=1.26376
        v(i,j)=0.590979
       enddo
       enddo 



c--      Iextt=1.90   ! external current to produce AP
c--      Opening up the wall points to get boundary

       open(unit=10, file="largeMaze2.txt")
       do k=1,(nx+1)*(ny+1)
       read(10,*,end=10)i,j
       phi(i,j)=0
       enddo
10       continue

       do i=0,15
       do j=0,15
        phi(i,j)=1
       enddo
       enddo

       do i=nx-15,nx
       do j=ny-15,ny
        phi(i,j)=1
       enddo
       enddo

c top left corner of the maze
       do i=10,20
       do j=ny-20,ny-10
        phi(i,j)=1
       enddo
       enddo

c       do i=0,nx+1
c       do j=0, ny+1
c       if(phi(i,j).eq.0)then
c       ut(i,j) = 3
c       endif
c       enddo
c       enddo



c----------------------------------------------------
c           time integration
c---------------------------------------------------
         do nt=0,ntime
c--        Iext=0
c-----  this injects current
c--      if(mod(nt,40000).eq.1)Iext=Iextt

c--------- updating boundary conditions for zero flux
c         do i=1,nx
c         do j=1,ny
c         u(i,0)=u(i,2)
c         u(i,ny+1)=u(i,ny-1)
c         u(0,j)=u(2,j)
c         u(nx+1,j)=u(nx-1,j)
c         enddo
c         enddo
c---updating boundary conditions for first derivative 0
          do i=1,nx
          do j=1,ny
          u(i,0)=u(i,1)
          u(i,ny+1)=u(i,ny)
          u(0,j)=u(1,j)
          u(nx+1,j)=u(nx,j)
          enddo
          enddo
c-------boundary conditions for maze
      do i=0, nx
      do j=0, ny
      n = phi(i,j)
      m = phi(i,j+1)
      if(n.gt.m)then
      u(i,j+1) = u(i,j)
      endif
      if(n.lt.m)then
      u(i,j) = u(i,j+1)
      endif
      n = phi(i,j)
      m = phi(i+1,j)
      if(n.gt.m)then
      u(i+1,j) = u(i,j)
      endif
      if(n.lt.m)then
      u(i,j) = u(i+1,j)
      endif
      enddo
      enddo
c------- keeping the beginning/end fixed at unexcited values
       do i=0,10
       do j=0,10
        u(i,j)=-0.99999
        v(i,j)=-0.66667
       enddo
       enddo
   
       do i=nx-15,nx
       do j=ny-15,ny
        u(i,j)=-0.99999
        v(i,j)=-0.66667
       enddo
       enddo  


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
c         if(u(i,j).gt.0)then
          write(nf,*)i,j,u(i,j),v(i,j)
c         endif
         endif
         enddo
         enddo
         if(mod(nt,100).eq.1)then
           close(nf)
           if(mod(nf,100).eq.0)then
             write(6,*)nf
           endif
           nf=nf+1
         endif

         t=t+dt
        if(mod(nt,20).eq.1)write(16,*)t,u(5,5),u(50,50),u(90,90)
        if(mod(nt,20).eq.1)write(17,*)t,v(5,5),v(50,50),v(90,90)

         enddo
         end
c---------------------------------------------------------------------------              Contraction
c--------------------------------------------------------------------------------

c Pin ends
c change beta to -ve beta
c
c      pinvalue = ut(nx-15,ny-15)
c      beta = -0.7
c
c----------------------------------------------------
c           time integration
c---------------------------------------------------
c
c         do nt=0,ntime
c--        Iext=0
c-----  this injects current
c--      if(mod(nt,40000).eq.1)Iext=Iextt
c
c--------- updating boundary conditions for zero flux
c         do i=1,nx
c         do j=1,ny
c         u(i,0)=u(i,2)
c         u(i,ny+1)=u(i,ny-1)
c         u(0,j)=u(2,j)
c         u(nx+1,j)=u(nx-1,j)
c         enddo
c         enddo
c
c-------- pin end values
c      do i = 0, 15
c      do j = 0, 15
c      u(i,j) = pinvalue
c      enddo
c      enddo
c
c      do i = nx-8, nx-4
c      do j = ny-8, ny-4
c      u(i,j) = pinvalue
c      enddo
c      enddo
c-------boundary conditions for maze
c      do i=0, nx
c      do j=0, ny
c      n = phi(i,j)
c      m = phi(i,j+1)
c      if(n.gt.m)then
c      u(i,j+1) = u(i,j)
c      endif
c      if(n.lt.m)then
c      u(i,j) = u(i,j+1)
c      endif
c      n = phi(i,j)
c      m = phi(i+1,j)
c      if(n.gt.m)then
c      u(i+1,j) = u(i,j)
c      endif
c      if(n.lt.m)then
c      u(i,j) = u(i+1,j)
c      endif
c      enddo
c      enddo

c---------integration in space
c        do i=1,nx
c        do j=1,ny

c        if(phi(i,j).eq.1)then
c         v(i,j)=v(i,j)+eps*(u(i,j)-gamma*v(i,j)+beta)*dt
c--------updating the voltage to ut
c         xlap=u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4*u(i,j)

c         ut(i,j)=u(i,j)+(u(i,j)-u(i,j)**3/3.-v(i,j))*dt+xlap*dlap

c        endif
c        enddo
c        enddo

c------- update u value
c         do i=1,nx
c         do j=1,ny
c         u(i,j)=ut(i,j)
c         if(mod(nt,100).eq.1)then
c         if(u(i,j).gt.-0.5)then
c          write(nf,*)i,j,u(i,j),v(i,j)
c         endif
c         endif
c         enddo
c         enddo
c         if(mod(nt,100).eq.1)then
c           close(nf)
c           if(mod(nf,100).eq.0)then
c             write(6,*)nf
c           endif
c           nf=nf+1
c         endif

c         t=t+dt
c        if(mod(nt,20).eq.1)write(16,*)t,u(5,5),u(50,50),u(90,90)
c        if(mod(nt,20).eq.1)write(17,*)t,v(5,5),v(50,50),v(90,90)

c         enddo
c         end





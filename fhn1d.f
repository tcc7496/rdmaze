
          implicit real*8 (a-h,o-z)
          parameter (nx=100)

         real eps,dt,dx,Diff,delta,a,b,gamma,dlap
         real u,v
         real ut   

         real time
         integer  ntime
c        dimensions for the 2 variables
         dimension u(0:nx+1),v(nx)
         dimension ut(0:nx+1)

c         fixed parameters
          nf=19
           eps=.08
           t=0
           gamma=0.8
           beta=.7
           dt=0.005
           dx=0.25
           Diff=1.
           dlap=Diff*dt/(dx*dx)
           ntime=10000
     

c---------------- Initial Conditions (rest state) ---------------------------------
        do i=0,nx+1
        u(i)=-1.199
        ut(i)=u(i)
        enddo
        do i=1,10
         u(i)=1
        enddo
        do i=1,nx
         v(i)=-0.6242
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
         u(0)=u(2)
         u(nx+1)=u(nx-1)
c---------integration in space
        do i=1,nx

         v(i)=v(i)+eps*(u(i)-gamma*v(i)+beta)*dt   
c--------updating the voltage to ut
         xlap=u(i+1)+u(i-1)-2.*u(i)
          
         ut(i)=u(i)+(u(i)-u(i)**3/3.-v(i))*dt+xlap*dlap

         enddo
          
c------- update u value
         do i=1,nx
         u(i)=ut(i)
         if(mod(nt,500).eq.1)then
          write(nf,*)i,u(i)
         endif
         enddo
         if(mod(nt,500).eq.1)then
           close(nf)
           write(6,*)nf
           nf=nf+1
         endif

         t=t+dt
        if(mod(nt,20).eq.1)write(16,*)t,u(5),u(50),u(90)
        if(mod(nt,20).eq.1)write(17,*)t,v(5),v(50),v(90)

         enddo
         end


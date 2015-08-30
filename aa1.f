       implicit real*8(a-h,o-z)
        dimension u(0:100,0:100),v(0:100,0:100),ut(0:100,0:100)
        dimension ph(0:100,0:100)
       a=2
       eps=10.
       b=0.4
       dt=0.025
       dx=.5
      
       do i=0,100
       do j=0,100
      u(i,j)= -0.8
      v(i,j)=-0.25
      ph(i,j)=0
        enddo
        enddo


        do i=1,20
        do j=1,20
        u(i,j)=0.1
        v(i,j)=0
       enddo
       enddo

       do n=1,200

        do i=0,100
         u(i,0)=u(i,2)
         u(i,100)=u(i,98)
         u(0,i)=u(2,i)
         u(100,i)=u(98,i)
       enddo

c      u(1,1)=1
c      v(1,1)=.2
c       write(6,*)n,u(1,1),u(1,3)
       do i=1,99
       do j=1,99
c      if(ph(i,j).eq.0)then

       xlap=u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4*u(i,j)
c         xlap=0
       ut(i,j)=u(i,j)+5.*eps*dt*(u(i,j)-u(i,j)**3-v(i,j))+
     %  dt*xlap/(dx*dx)
c      use ut instead of u because of time
       v(i,j)=v(i,j)+dt*(u(i,j)-a*v(i,j)+b)

c      endif
       enddo 
       enddo
 
      
        do i=1,99
        do j=1,99
        u(i,j)=ut(i,j)
        if(u(i,j).gt.0)then
          if(mod(n,2).eq.1)then
        write(n+19,*)i,j
        endif
        endif
        enddo
        enddo
        


      write(11,*)n,u(5,5),v(5,5)
       enddo




      end

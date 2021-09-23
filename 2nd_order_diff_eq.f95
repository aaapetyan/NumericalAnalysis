program thomasalgorithm
implicit none
integer n,k
real(8) ax,bx,h,alfx,betx,aa,bb,pi
real(8), allocatable, dimension(:) :: x,y,a,b,c,d,f,p,q,alf,bet
open(unit=1, file='results.dat')
n=10
pi=4*atan(1.0)
ax=0.0; bx=1.0; h=(bx-ax)/n; 
alfx=0; betx=-0.3; aa=1.0; bb=2.0   ! the boundary conditions data           
allocate(x(0:n),y(0:n), a(0:n), b(0:n), c(0:n), &
d(0:n), alf(0:n),bet(0:n), f(0:n),q(0:n),p(0:n))
  do k=0,n
    x(k) = ax + k*h
    f(k) = x(k)*(1-x(k))            ! defining 
    q(k) = x(k)/(1+sin(pi/2*x(k)))  ! the
    p(k) = 1+x(k)                   ! equation
      a(k) = 1/h**2-p(k)/(2*h)
      b(k) = 2/h**2-q(k)            ! solving
      c(k) = 1/h**2+p(k)/(2*h)      ! the
      d(k) = f(k)                   ! system
  enddo
      c(0) = 1/h+p(0)/2
      b(0) = 1/h+p(0)/2-h/2*q(0)+alfx
      d(0) = h/2*f(0) + aa
        a(n) = -1/h+p(n)/2
        b(n) = -1/h+p(n)/2+h/2*q(n)+betx
        d(n) = -h/2*f(n) + bb
          alf(0) = c(0)/b(0); bet(0) = - d(0)/b(0)
            do k=1,n 
	            alf(k) = c(k)/(b(k)-a(k)*alf(k-1)) 
	            bet(k) = (a(k)*bet(k-1)-d(k))/(b(k)-a(k)*alf(k-1))
            enddo
              y(n) = bet (n)
                do k=n-1,0,-1
	                y(k) = alf(k)*y(k+1)+bet(k)
                enddo                
                  write(1,*) "a, b: "
                    do k=0,n
	                    write(1,*) a(k),b(k)
                    enddo
                  write(1,*) "c, d: "
                    do k=0,n
	                    write(1,*) c(k),d(k)		! output
                    enddo
                  write(1,*) "alf, bet: "
                    do k=0,n
                    	write(1,*) alf(k), bet(k)
                    enddo
                  write(1,*) "x, y: "
                    do k=0,n
	                    write(1,*) x(k),y(k)
                    enddo
deallocate(x,y,a,b,c,d,f,p,q,alf,bet)
close(1)
end

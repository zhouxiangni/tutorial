
module interpolation 

! the result is not the same as MATLAB, but not wrong ,check it later 
! 03/21//2008! 03/21//2008! 03/21//2008

		implicit none 
		private 
	    integer , parameter :: maxm = 1e5
		double precision :: work(3*maxm)

		public :: intpol_10 
		public :: intpol_12
		public :: spline
		public :: seval 
		public :: linear_interpolate
contains


  subroutine intpol_10(m, s, x, mm, ss, xx)
    integer, intent(in) :: m , mm 
	double precision , intent(in) :: s(m),x(m), ss(mm)
	double precision  :: xx(mm) , tmp(4)
	integer :: i, m2, iflag 
!C
!C Input: nd - dimension of the config space
!C        m - number of images along the string
!C        s(m) - normalized arclength of old images
!C        x(m) - coordinates of old images
!C        mm - number of new images on the string
!C        ss(mm) - equal arclength
!C        work - working array of length at least 3*m
!C Output: 
!C        xx(mm) - coordinates of new images after interpolation
!C 
      !print *, 'interpolation'  , m, 1E6 

      m2=m*2
      call spline(m,s,x,work(1),work(m+1),work(m2+1))
!C 
      iflag = 1 
	  do i=1,mm
        call seval(m,ss(i),s,x,work(1),work(m+1),work(m2+1),iflag,tmp)
        xx(i)=tmp(1)
      enddo
!C

      end subroutine intpol_10


 subroutine intpol_12(m, s, x, mm, ss, xx, fstDer,snder )
      implicit none
!C
!C Input: nd - dimension of the config space
!C        m - number of images along the string
!C        s(m) - normalized arclength of old images
!C        x(nd,m) - coordinates of old images
!C        mm - number of new images on the string
!C        ss(mm) - equal arclength
!C        work - working array of length at least 3*m
!C Output: 
!C        xx(nd,mm) - coordinates of new images after interpolation
!C 
      integer :: m,   mm  
	  double precision :: s(m), x(m), ss(mm), xx(mm), fstDer(mm), snder(mm)   
      double precision :: tmp(4)
	  integer :: i, iflag 

	
      call spline(m,s,x,work(1:m),work(m+1:2*m),work(m*2+1:m*3))

	  iflag = 2 
 
	 
      do i=1,mm
        call seval(m,ss(i),s,x,work(1:m),work(m+1:2*m),work(m*2+1:3*m),iflag,tmp)
        xx(i)=tmp(1)
		fstDer(i) = tmp(2)
		snder(i) = tmp(3)
      enddo		
      return
      end subroutine 

!
!---------------------------------------------------------------
!---------------------------------------------------------------
!---------------------------------------------------------------
      subroutine spline (n, x, y, b, c, d)
      integer n
      double precision x(n), y(n), b(n), c(n), d(n)
!c
!! the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
!! for a cubi!interpolating spline
!c
!c    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!c
!c    for  x(i) .le. x .le. x(i+1)
!c
!c  input..
!c
!c    n = the number of data points or knots (n.ge.2)
!c    x = the abscissas of the knots in strictly increasing order
!c    y = the ordinates of the knots
!c
!c  output..
!c
!c    b, c, d  = arrays of spline coefficients as defined above.
!c
!c  using  p  to denote differentiation,
!c
!c    y(i) = s(x(i))
!c    b(i) = sp(x(i))
!c    c(i) = spp(x(i))/2
!c    d(i) = sppp(x(i))/6  (derivative from the right)
!c
!c  the accompanying function subprogram  seval  can be used
!c  to evaluate the spline.
!c
!c
      integer nm1, ib, i
      double precision t
 
      nm1 = n-1
      if ( n .lt. 2 ) return
      if ( n .lt. 3 ) go to 50
!c
!c  set up tridiagonal system
!c
!c  b = diagonal, d = offdiagonal, c = right hand side.
!c
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do 10 i = 2, nm1
         d(i) = x(i+1) - x(i)
         b(i) = 2.*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
   10 continue
!c
!c  end conditions.  third derivatives at  x(1)  and  x(n)
!c  obtained from divided differences
!c
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.
      c(n) = 0.
      if ( n .eq. 3 ) go to 15
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
!c
!c  forward elimination
!c
   15 do 20 i = 2, n
         t = d(i-1)/b(i-1)
         b(i) = b(i) - t*d(i-1)
         c(i) = c(i) - t*c(i-1)
   20 continue
!c
!c  back substitution
!c
      c(n) = c(n)/b(n)
      do 30 ib = 1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
   30 continue
!c
!c  c(i) is now the sigma(i) of the text
!c
!c  compute polynomial coefficients
!c
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
      do 40 i = 1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.*c(i)
   40 continue
      c(n) = 3.*c(n)
      d(n) = d(n-1)
      return
!c
   50 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.
      d(1) = 0.
      b(2) = b(1)
      c(2) = 0.
      d(2) = 0.
      return
      end subroutine spline 

!---------------------------------------------------------------
    
	
	
	subroutine seval(n, u, x, y, b, c, d, iflag, f)
      integer n, iflag
      double precision  u, x(n), y(n), b(n), c(n), d(n), f(4)
!c
!c  this subroutine evaluates the cubic spline function
!c
!c    f(1) = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
!C  and/or its derivatives
!c
!c    where  x(i) .lt. u .lt. x(i+1), using horner's rule
!c
!c  if  u .lt. x(1) then  i = 1  is used.
!c  if  u .ge. x(n) then  i = n  is used.
!c
!c  input..
!c
!c    n = the number of data points
!c    u = the abscissa at which the spline is to be evaluated
!c    x,y = the arrays of data abscissas and ordinates
!c    b,c,d = arrays of spline coefficients computed by spline
!C    iflag =1: function value
!C           2: 1st + 2nd derivative
!c
!c  if  u  is not in the same interval as the previous call, then a
!c  binary search is performed to determine the proper interval.
!c
      integer i, j, k
      double precision dx, tmp1, tmp2
      data i/1/
      if ( i .ge. n ) i = 1
      if ( u .lt. x(i) ) go to 10
      if ( u .le. x(i+1) ) go to 30
!c
!c  binary search
!c
   10 i = 1
      j = n+1
   20 k = (i+j)/2
      if ( u .lt. x(k) ) j = k
      if ( u .ge. x(k) ) i = k
      if ( j .gt. i+1 ) go to 20
!c
!c  evaluate spline
!c
   30 dx = u - x(i)
!C
      if(iflag.eq.1) then
        f(1) = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      elseif(iflag.eq.2) then
		f(1) = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
        tmp1=3.d0*dx*d(i)
        tmp2=2.d0*c(i) + tmp1
        f(2) = b(i) + dx*tmp2
        f(3) = tmp2 + tmp1
      endif
!C
      return
      end subroutine seval
!---------------------------------------------------------------

  subroutine linear_interpolate(m,t0,u0,n,t1,u1)
!C    linear interpolation for u: from mesh t0 to t1 
!C    the range of mesh t0 and t1 should with [0,1], maybe not contain 0 and 1. but 
 
!C    t0, u0 : old images, t0(:) increase 
!C    t1, u0 : new images, t1(:) increase 

        integer, intent(in)::   m, n 
        double precision, intent(in) :: t0(m), u0(m), t1(n) 
		double precision  :: u1(n)
        integer i0,i1, j, k, loc0, loc1 
 
		u1 = 0 
 ! find the region in old mesh where t1(1)  is 
 ! t0(loc0) <= t1(1) < t0(loc1) where loc1 - loc0 = 1
 ! t1(1) > t0(m) :   loc0 = m, loc1 = m+1 
       loc1 = m +1 
       do k = 1, m 
           if ( t0(k) > t1(1) ) then
              loc1 = k 
              exit 
           endif 
       enddo 
       loc0 = loc1 - 1 
	   	

		
		         
        do k = 1,n

          if ( loc0 == 0 ) then 
             i0 = 1 
             i1 = 2  
			 goto 200   ! use  exponential outerpolation 
          elseif (loc0 == m) then
             i0 = m -1  ! use exponential outerpolation 
             i1 = m 
			 goto 300
          else 
             i0 = loc0 
             i1 = loc1     
          endif 
		  
! use t0(i0) , t0(i1) to do linear interpolation
        
         u1(k) =  (t0(i1)- t1(k) )/(t0(i1)-t0(i0))*u0(i0) 
         u1(k) =  u1(k)+ ( t1(k) -t0(i0))/(t0(i1)-t0(i0))*(u0(i1))		
		

		 goto 500

		  
	

200       if(t0(1)==0.d0)  then 
			u1(k) = 0.d0 
		  else 
			u1(k) =    u0(1) +  log(t1(k)/t0(1))
		  endif 
		  go to 500

300       if ( t0(m) == 1.d0 ) then 
			u1(k) = 1.0d0
		  else
			u1(k) =   u0(m) - log( (t1(k) -1 )/ ( t0(m) -1 ) ) 
		  endif
       		   

500       continue 
			
		  
		  if (k==n)  then 
				return 
		  endif
! seek the regin for next t1(k+1)
          i1 =  loc1
          loc1 = m +1 
          do j = i1, m 
           if ( t0(j) >= t1(k+1) ) then
              loc1 = j 
              exit 
           endif 
          enddo 
          loc0 = loc1 - 1		 


        enddo        
		end subroutine 


! 
!---------------------------------------------------------------


end module 

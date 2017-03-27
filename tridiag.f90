
	module tri_sys_solver

contains
!---------------------------------------------------------------
    SUBROUTINE tridag_solver(a,b,c,r,u,n)
! NUMERICAL RECIPES FORTRAN 77
    INTEGER :: n,NMAX
    double precision ::  a(n),b(n),c(n),r(n),u(n)
    PARAMETER(NMAX=50000)
    INTEGER :: j
    double precision ::  bet, gam(NMAX)
! solve the linear tri-diagonal equation:   a(k)*u(k-1) + b(k)*u(k) + c(k)*u(k+1) = r(k) for k = 1..n   with  small modification at ends for u
! olves for a vector u(1:n) of length n the tridiagonal linear set given by equation (2.4.1).
!a(1:n), b(1:n), c(1:n), and r(1:n) are input vectors and are not modied.
!a(2:n)  is the low band ; c(1:n-1) is the upper band, b(i:n) is the diagonal
!	Parameter: NMAX is the maximum expected value of n.
!	One vector of workspace, gam is needed.
    if(b(1) == 0.)	pause 'tridag: rewrite equations'
!	If this happens then you should rewrite your equations as a set of order N . 1, with u2
! trivially eliminated.
    bet=b(1)
    u(1)=r(1)/bet
    do   j=2,n  !Decomposition and forward substitution.
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet == 0.) pause 'tridag failed' !Algorithm fails; see below.
        u(j)=(r(j)-a(j)*u(j-1))/bet
    enddo
    do   j=n-1,1,-1 ! Backsubstitution.
        u(j)=u(j)-gam(j+1)*u(j+1)
    enddo
    return
    END SUBROUTINE tridag_solver 


	function tridag_check(a,b,c,r,u,n) result(err)
		integer , intent(in) :: n
		double precision , intent(in) :: a(n),b(n),c(n),u(n),r(n)
		integer :: k, err, tmp
		err = 0  
		do k =2 , n-1
			tmp = a(k)*u(k-1)+ b(k)*u(k)+ c(k)*u(k+1)
			tmp = abs(tmp-r(k))
			err = max(tmp,err)
		enddo 
			k=1
			tmp =  b(k)*u(k)+ c(k)*u(k+1)
			tmp = abs(tmp-r(k))
			err = max(tmp,err)
			k=n
			tmp = a(k)*u(k-1)+ b(k)*u(k)
			tmp = abs(tmp-r(k))
			err = max(tmp,err)
			return 
	end function 


	SUBROUTINE tridag_inv(a,b,c,M,n)
!   directly seek the inverse matrix, M, of the tri-diagonal matrix defind by arrays: a,b,c
		INTEGER, INTENT(IN) :: N 
		DOUBLE PRECISION , INTENT(IN) :: A(N),B(N),C(N)
		DOUBLE PRECISION :: E(N), M(N,N)
		integer :: k 		
		do k = 1, n 
			E = 0.d0 ;	E(k) = 1.d0 
			call tridag_solver(a,b,c,E,M(:,k),n)
		enddo 
	END SUBROUTINE TRIDAG_INV



   !SUBROUTINE CYC_tridiag_precessor 
   ! reference :
   ! www-dinma.univ.trieste.it/nirftc/research/cyctrid/cyclic_trid.pdf

   !END SUBROUTINE CYC_tridiag_precessor

   !SUBROUTINE CYC_tridiag_solver 
   !END SUBROUTINE CYC_tridiag_solver 




end module tri_sys_solver

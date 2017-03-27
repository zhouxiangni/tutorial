module LU_module 

 ! Reference 
 ! http://www.math.buffalo.edu/~pitman/courses/mth437/na2/node6.html
     implicit none 

    !interface  LU_solver 
          !subroutine LU_solver(n,A,b)
                 !integer, intent(in) :: N
         !double precision :: A(n,n),b(n)
         !integer :: indx(n), D
      !end subroutine LU_solver 
    !end interface LU_solver  

        double precision, parameter :: tiny = 1e-20
contains 
        
!*********************************************************************72
    subroutine ludcmp(A,N,NP,INDX,D)
!*********************************************************************72
! GIVEN THE NxN MATRIX  A  WITH PHYSICAL DIMENSIONS NPxNP,
! FIND THE LU DECOMPOSITION OF  A
! USE WITH LUBKSB TO SOLVE LINEAR SYSTEM
! CROUT's METHOD WITH PARTIAL PIVOTING
! indx is output vector which records row permutations
! d is output scalar equal to +/- 1, depending on even or odd
! number of permutations
! in  "NUMERICAL RECIPES"  p 35 ff.

    implicit none
    
        integer , intent(in) :: n, np 
        double precision ::  A(np,np)
        integer :: indx(n), d 
        double precision :: vv(n), aamax,dum,sum 
        
        
        integer :: i, j, k,imax
        
    d=1
    do 12 i=1,n
        aamax=0.d0
        do 11 j=1,n
            if (dabs(A(i,j)) > aamax) aamax=dabs(A(i,j))
        11 ENDDO
        if (aamax == 0.d0) then
            write(6,666)
            666 format('','LUDCMP BOMBED')
             return 
        endif
        vv(i)=1.d0/aamax
    12 ENDDO

! CROUT's METHOD
    do 19 j=1,n
        do 14 i=1,j-1
            sum=A(i,j)
            do 13 k=1,i-1
                sum=sum-A(i,k)*A(k,j)
            13 ENDDO
            A(i,j)=sum
        14 ENDDO
    ! SEARCH FOR PIVOT ELEMENT
        aamax=0.d0
        do 16 i=j,n
            sum=A(i,j)
            do 15 k=1,j-1
                sum=sum-A(i,k)*A(k,j)
            15 ENDDO
            A(i,j)=sum
            dum=vv(i)*dabs(sum)
            if (dum >= aamax) then
                imax=i
                aamax=dum
            endif
        16 ENDDO
    ! INTERCHANGE ROWS???
        if (j /= imax) then
            do 17 k=1,n
                dum=A(imax,k)
                A(imax,k)=A(j,k)
                A(j,k)=dum
            17 ENDDO
            d=-d
            vv(imax)=vv(j)
        endif
        indx(j)=imax
        if (A(j,j) == 0.d0) A(j,j)=tiny
        if (j /= n) then
            dum=1.d0/A(j,j)
            do 18 i=j+1,n
                A(i,j)=A(i,j)*dum
            18 ENDDO
        endif
    19 ENDDO


    999 return
    end subroutine ludcmp

!*********************************************************************72
    subroutine lubksb(A,N,NP,INDX,B)
!*********************************************************************72
! GIVEN THE NxN MATRIX  A  WITH PHYSICAL DIMENSIONS NPxNP,
! WHICH HAS BEEN LU-DECOMPOSED IN ludcmp,
! SOLVE         AX = B
! indx is output vector which records row permutations
! B is input r.h.s. vector and on output contains the solution X
! in  "NUMERICAL RECIPES"  p 36 ff.

    implicit none
        integer, intent(in) :: n, np 
        double precision, intent(in) :: A(np,np)
    integer , intent(in) ::  indx(n)
    double precision :: B(N), sum 
        integer :: ii, ll, i,j 
        
    ii=0
    do 12 i=1,n
        ll=indx(i); sum=B(ll); B(ll)=B(i)
        if (ii /= 0) then
            do  j=ii,i-1
                sum=sum-A(i,j)*B(j)
             ENDDO
        else if (sum /= 0.d0) then
            ii=i
        endif
        B(i)=sum
    12 ENDDO

! BACKSUBSTITUTION
    do 14 i=n,1,-1
        sum=B(i)
        if (i < n) then
            do 13 j=i+1,n
                sum=sum-A(i,j)*B(j)
            13 ENDDO
        endif
        B(i)=sum/A(i,i)
    14 ENDDO


    return
    end subroutine lubksb
        
        
        function LU_solver2(n,A,b) result(err)

         integer, intent(in) :: N
         double precision :: A(n,n),b(n), A0(n,n), b0(n), err
         integer :: indx(n), D


         A0 =A; b0=b; 
         call ludcmp(A,N,N,INDX,D)
         call lubksb(A,N,N,INDX,B)
         b0 = matmul(A0,b)-b0;
         err = sqrt(dot_product(b0,b0)) 
         
    end function 
        
        subroutine LU_solver(n,A,b)

         integer, intent(in) :: N
         double precision :: A(n,n),b(n)
         integer :: indx(n), D
 
         call ludcmp(A,N,N,INDX,D)
         call lubksb(A,N,N,INDX,B)

         !call refine(A,ALU,N,N,INDX,B,x)
         
    end subroutine 
    
        
! *********************************************************************72
         subroutine refine(A,ALU,N,NP,INDX,B,X)
!*********************************************************************72
!        GIVEN THE NxN MATRIX  A  WITH PHYSICAL DIMENSIONS NPxNP,
!        AND LU DECOMPOSITION ALU, X THE APPROX SOLUTION
!        OF   AX = b  
!        USE ITERATIVE REFINEMENT TO IMPROVE THE SOLN
!        indx is output vector which records row permutations
!        in  "NUMERICAL RECIPES"  p 42 ff.

            implicit none
                    integer, intent(in) :: N,NP
                integer, intent(in) :: indx(n)
            double precision,intent(in) :: A(np,np), ALU(np,np),b(n)
                        double  precision :: x(n),r(n), spd
                    integer :: i , j 
            

      do  i=1,n
                        spd =-b(i)
                        do  j=1,n
             spd=spd+A(i,j)*x(j)
         enddo 
         r(i) = spd 
         enddo 


         call lubksb(ALU,n,np,indx,r)
                 
         x(1:n)=x(1:n)-r(1:n)
      
         return
         end subroutine refine 
end module LU_module

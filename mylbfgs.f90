!---------------------------------------------------------------
! in future : V1.21   
!  # 1      
!         get rid of tri_sys_solver module to let this modoule independently 
!  # 2    rewrite some 77 subrouinte in concise 90 intrinsic subrouinte 

! V1.2 MODIFIED BY XIANG ZHOU, 03-25-2008
!---------------------------------------------
! # 1
! IFLAG = 4 is introduced to make solution X is accessible outside
!              after one iteration of L-BFGS 
!              the origianal version gives trial X in line search
!              immediately after last iteration of L-BFGS    
! #2
! STP :    initial guess of STP for line search changed from 1 to 
!                  stepsize in the previous iteration 


! V1.1 MODIFY BY XIANG ZHOU , 03-25-2007 
! ---------------------------------------
!  # 1
! change the diagonal precondition matrix to tri diagonal matrix (the matrix b)
! DIAGCO -> PRECCO
! DIAG  -> INVERSE OF TRI_DIAG [SEARCH DIRECTION = DIAG*GRADIENT]
!         NOTE :  THE CHECK OF  POSITIVE-DIFINITE OF TRI_DIAG IS SKIPPED
!  # 2
! IFLAG =3  is introduced to distinguish normal termination 
!                               and initial state IFLAG=0
!  # 3
! work memory is provided by two subroutines : 
!                   ALLOCATE_LBFGS AND DEALLOCATE_LBFGS METHOD , 
!              instead of explicitly supply work memory for l-bfgs subroutine 
! WP          is introduced to deal with precondition case for solving TRI_DIAG 
!
! # 4
! NPT   is initialized to be zero 
!---------------------------------------------------------------


	MODULE MY_LBFGS_MODULE  

	IMPLICIT NONE 
	INTEGER, PARAMETER  :: MP = 6 , LP =6


! PARAMETERS FOR LINE SEARCH ROUTINE
      DOUBLE PRECISION :: FTOL= 1.0D-4
      INTEGER :: MAXFEV= 40
      !  default 20
	  DOUBLE PRECISION  :: GTOL=0.9D0,STPMIN=1.0D-20,STPMAX=1.0D+20
	  !GTOL default = 0.9

!------------
!WORK MEMRORY 
!------------
	DOUBLE PRECISION, ALLOCATABLE :: W(:), WA(:)  , WP(:)

	
! RESTRICT DATA ACCESS FROM OUTSIDE 	
	PRIVATE :: W, WA , WP
	PRIVATE :: MP, LP, GTOL, STPMIN,STPMAX

CONTAINS  
!-----------------------------------------------
! ALLOCATE AND DEALLOCATE WORK MEMORY FOR LBFGS 
!-----------------------------------------------
	SUBROUTINE ALLOCATE_LBFGS(N,M)
        ! M is the level of LBFGS
        ! N is the total no. of variables to minimize  over 
		INTEGER , INTENT(IN) :: N , M 
		INTEGER :: ALLOC_ERR
		ALLOCATE( W(1:N*(2*M+1)+2*M), WA(1:N), WP(1:N),STAT=ALLOC_ERR)
		IF (ALLOC_ERR >0) THEN 
			WRITE (*,*) 'ALLOCATE WORK MEMORY FAILS IN LBFGS';	PAUSE 
		ENDIF 
		W = 0.D0; 	WA = 0.D0 
	END SUBROUTINE ALLOCATE_LBFGS
	
	SUBROUTINE DEALLOCATE_LBFGS()
		DEALLOCATE(W); DEALLOCATE(WA,WP); 
	END SUBROUTINE DEALLOCATE_LBFGS

! ----------------------------------------------------------------------
! This file contains the LBFGS algorithm and supporting routines

! ****************
! LBFGS SUBROUTINE
! ****************

    SUBROUTINE LBFGS(N,M, ITER, X,F,G,PRECO,TRI_DIAG,IPRINT,EPS,XTOL,IFLAG)

    INTEGER , INTENT(IN)   :: N,M,IPRINT(2)
	INTEGER , INTENT(INOUT) :: IFLAG, ITER 
    double precision, INTENT(INOUT) :: X(N),F, G(N),TRI_DIAG(N,3) 
	DOUBLE PRECISION, INTENT(IN) :: EPS,XTOL
    LOGICAL , INTENT(IN) :: PRECO

! LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
! JORGE NOCEDAL
! *** July 1990 ***


! This subroutine solves the unconstrained minimization problem

! min F(x),    x= (x1,x2,...,xN),

! using the limited memory BFGS method. The routine is especially
! effective on problems involving a large number of variables. In
! a typical iteration of this method an approximation Hk to the
! inverse of the Hessian is obtained by applying M BFGS updates to
! a diagonal matrix Hk0, using information from the previous M steps.
! The user specifies the number M, which determines the amount of
! storage required by the routine. The user may also provide the
! diagonal matrices Hk0 if not satisfied with the default choice.
! The algorithm is described in "On the limited memory BFGS method
! for large scale optimization", by D. Liu and J. Nocedal,
! Mathematical Programming B 45 (1989) 503-528.

! The user is required to calculate the function value F and its
! gradient G. In order to allow the user complete control over
! these computations, reverse  communication is used. The routine
! must be called repeatedly under the control of the parameter
! IFLAG.

! The steplength is determined at each iteration by means of the
! line search routine MCVSRCH, which is a slight modification of
! the routine CSRCH written by More' and Thuente.

! The calling statement is

! CALL LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG)

! where

! N       is an INTEGER variable that must be set by the user to the
! number of variables. It is not altered by the routine.
! Restriction: N>0.

! M       is an INTEGER variable that must be set by the user to
! the number of corrections used in the BFGS update. It
! is not altered by the routine. Values of M less than 3 are
! not recommended; large values of M will result in excessive
! computing time. 3<= M <=7 is recommended. Restriction: M>0.

! X       is a DOUBLE PRECISION array of length N. On initial entry
! it must be set by the user to the values of the initial
! estimate of the solution vector. On exit with IFLAG=0, it
! contains the values of the variables at the best point
! found (usually a solution).

! F       is a DOUBLE PRECISION variable. Before initial entry and on
! a re-entry with IFLAG=1, it must be set by the user to
! contain the value of the function F at the point X.

! G       is a DOUBLE PRECISION array of length N. Before initial
! entry and on a re-entry with IFLAG=1, it must be set by
! the user to contain the components of the gradient G at
! the point X.

! DIAGCO  is a LOGICAL variable that must be set to .TRUE. if the
! user  wishes to provide the diagonal matrix Hk0 at each
! iteration. Otherwise it should be set to .FALSE., in which
! case  LBFGS will use a default value described below. If
! DIAGCO is set to .TRUE. the routine will return at each
! iteration of the algorithm with IFLAG=2, and the diagonal
! matrix Hk0  must be provided in the array DIAG.


! DIAG    is a DOUBLE PRECISION array of length N. If DIAGCO=.TRUE.,
! then on initial entry or on re-entry with IFLAG=2, DIAG
! it must be set by the user to contain the values of the
! diagonal matrix Hk0.  Restriction: all elements of DIAG
! must be positive.

! IPRINT  is an INTEGER array of length two which must be set by the
! user.

! IPRINT(1) specifies the frequency of the output:
! IPRINT(1) < 0 : no output is generated,
! IPRINT(1) = 0 : output only at first and last iteration,
! IPRINT(1) > 0 : output every IPRINT(1) iterations.

! IPRINT(2) specifies the type of output generated:
! IPRINT(2) = 0 : iteration count, number of function
! evaluations, function value, norm of the
! gradient, and steplength,
! IPRINT(2) = 1 : same as IPRINT(2)=0, plus vector of
! variables and  gradient vector at the
! initial point,
! IPRINT(2) = 2 : same as IPRINT(2)=1, plus vector of
! variables,
! IPRINT(2) = 3 : same as IPRINT(2)=2, plus gradient vector.


! EPS     is a positive DOUBLE PRECISION variable that must be set by
! the user, and determines the accuracy with which the solution
! is to be found. The subroutine terminates when

! ||G|| < EPS max(1,||X||),

! where ||.|| denotes the Euclidean norm.

! XTOL    is a  positive DOUBLE PRECISION variable that must be set by
! the user to an estimate of the machine precision (e.g.
! 10**(-16) on a SUN station 3/60). The line search routine will
! terminate if the relative width of the interval of uncertainty
! is less than XTOL.

! W       is a DOUBLE PRECISION array of length N(2M+1)+2M used as
! workspace for LBFGS. This array must not be altered by the
! user.

! ITER    is the number of iterations in LBFGS

! NFUN    is the number of evaluations of objective function and its gradient 
! NFEV    is the number of evaluations of obj and grad in one line search 

! IFLAG   is an INTEGER variable that must be set to 0 on initial entry
! to the subroutine. A return with IFLAG<0 indicates an error,
! and IFLAG=0 indicates that the routine has terminated without
! detecting errors. On a return with IFLAG=1, the user must
! evaluate the function F and gradient G. On a return with
! IFLAG=2, the user must provide the diagonal matrix Hk0.

!---------------------------------------------------------------
! IFLAG =3  SUCCESSFUL RETURN WITHOU ERROR 
!---------------------------------------------------------------

! The following negative values of IFLAG, detecting an error,
! are possible:

! IFLAG=-1  The line search routine MCSRCH failed. The
! parameter INFO provides more detailed information
! (see also the documentation of MCSRCH):

! INFO = 0  IMPROPER INPUT PARAMETERS.

! INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF
! UNCERTAINTY IS AT MOST XTOL.

! INFO = 3  MORE THAN (MAXFEV)(default=20) FUNCTION EVALUATIONS WERE
! REQUIRED AT THE PRESENT ITERATION.

! INFO = 4  THE STEP IS TOO SMALL.

! INFO = 5  THE STEP IS TOO LARGE.

! INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
! THERE MAY NOT BE A STEP WHICH SATISFIES
! THE SUFFICIENT DECREASE AND CURVATURE
! CONDITIONS. TOLERANCES MAY BE TOO SMALL.


! IFLAG=-2  The i-th diagonal element of the diagonal inverse
! Hessian approximation, given in DIAG, is not
! positive.

! IFLAG=-3  Improper input parameters for LBFGS (N or M are
! not positive).

! ON THE DRIVER:

! The program that calls LBFGS must contain the declaration:

! EXTERNAL LB2

! LB2 is a BLOCK DATA that defines the default values of several
! parameters described in the COMMON section.



! COMMON:

! The subroutine contains one common area, which the user may wish to
! reference:


! MP  is an INTEGER variable with default value 6. It is used as the
! unit number for the printing of the monitoring information
! controlled by IPRINT.

! LP  is an INTEGER variable with default value 6. It is used as the
! unit number for the printing of error messages. This printing
! may be suppressed by setting LP to a non-positive value.

! GTOL is a DOUBLE PRECISION variable with default value 0.9, which
! controls the accuracy of the line search routine MCSRCH. If the
! function and gradient evaluations are inexpensive with respect
! to the cost of the iteration (which is sometimes the case when
! solving very large problems) it may be advantageous to set GTOL
! to a small value. A typical small value is 0.1.  Restriction:
! GTOL should be greater than 1.D-04.

! STPMIN and STPMAX are non-negative DOUBLE PRECISION variables which
! specify lower and uper bounds for the step in the line search.
! Their default values are 1.D-20 and 1.D+20, respectively. These
! values need not be modified unless the exponents are too large
! for the machine being used, or unless the problem is extremely
! badly scaled (in which case the exponents should be increased).


! MACHINE DEPENDENCIES

! The only variables that are machine-dependent are XTOL,
! STPMIN and STPMAX.


! GENERAL INFORMATION

! Other routines called directly:  DAXPY, DDOT, LB1, MCSRCH

! Input/Output  :  No input; diagnostic messages on unit MP and
! error messages on unit LP.


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    DOUBLE PRECISION ::  ONE,ZERO,GNORM,DDOT,STP1,&
         STP,YS,YY,SQ,YR,BETA,XNORM
!WP's declaration is removed. V1.2
    INTEGER :: NFUN,POINT,ISPT,IYPT,INFO, &
        BOUND,NPT,CP,I,NFEV,INMC,IYCN,ISCN,J
    LOGICAL :: FINISH

    SAVE
    DATA ONE,ZERO/1.0D+0,0.0D+0/
	
	
! INITIALIZE
! ----------
        
!  NORMAL TERMINATION 
!     print *,'Start :iflag=', iflag, 'iter=', iter ; pause

      SELECT CASE (IFLAG)
        CASE (0)  
            GOTO 10      ! FIRST ENTRY 
        CASE (1)  
            GOTO 172     ! UPADTED F AND G IS PROVIDED 
        CASE (2)  
            GOTO 100     ! PRECOND IS PROVIDED 
        CASE (3)  
            RETURN        ! NORMAL TERMINATION 
!V1.2
       CASE (4)  
            GOTO 80      ! NEW ITERATION OF L-BFGS, LINE SEARCH WAS DONE 
      END SELECT 

! FIRST ENTRY 
    10 ITER= 0
    IF(N <= 0 .OR. M < 0) GO TO 196
	

    IF(GTOL <= 1.D-04) THEN
        IF(LP > 0) WRITE(LP,245)
        GTOL=9.D-01
    ENDIF
    NFUN= 1
    POINT= 0

!V1.1#4
!   NPT = N*POINT ! NPT NEED TO BE INITIALIZED BY POINT

    FINISH= .FALSE.
	
    IF(PRECO) THEN
! check positive definite 
!        DO 30 I=1,N
!            IF (DIAG(I) <= ZERO) GO TO 195
!        30 END DO
    ELSE			
!  default  for first iteration: steepest descent 
		tri_diag = 0.D0 		
        DO 40 I=1,N        	
			tri_diag(I,2)=1.D0				
        40 END DO
    ENDIF
	
	

! THE WORK VECTOR W IS DIVIDED AS FOLLOWS:
! ---------------------------------------
! THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND
! OTHER TEMPORARY INFORMATION.
! LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO.
! LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED
! IN THE FORMULA THAT COMPUTES H*G.
! LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH
! STEPS.
! LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M
! GRADIENT DIFFERENCES.

! THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A
! CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT.

	IF (M==0) THEN 
		ISPT = 0 
	ELSE 
		ISPT= N+2*M
	ENDIF 

    IYPT= ISPT+N*M 
	

 
! FIRST SEARCH DIRECTION AND STEP SIEZE
! --------------------------------------
!    DO 50 I=1,N
!        W(ISPT+I)= -G(I)*diag(I)
!    50 END DO

	CALL TRIDAG_solver(tri_diag(:,1),tri_diag(:,2),tri_diag(:,3),&
                            -G(1:N), W(ISPT+1:ISPT+N), N)		 
	
!V1.21
    GNORM= DSQRT(DOT_PRODUCT(G,G))
    !GNORM= DSQRT(DDOTP(N,G,1,G,1))
    STP1= ONE/GNORM    ! STPE SIZE : X1-X0 = STP1* SEARCH DIRECTION

    IF(IPRINT(1) >= 0) CALL LB1(IPRINT,ITER,NFUN,GNORM,N,M,X,F,G,STP,FINISH)
    

! --------------------
! MAIN ITERATION LOOP
! --------------------

    80  continue 

! V1.1
!      ITER = ITER +1 

!V1.0  
!     IF (ITER == 1) GO TO 165
!      BOUND=ITER-1  !initilzed for BFGS
!      IF (ITER > M) BOUND=M 
              
!V1.2  !first entry changed, iter++ after finished iteration
       IF(ITER==0) GO TO 165 
       BOUND=ITER !initilzed for BFGS
       IF (ITER+1 > M) BOUND=M 


      INFO=0  ! initialized for line search 



!-----------------------------------
! STEEPEST DESCENT METHOD FOR M == 0 
!-----------------------------------
    IF(M == 0) THEN 
        POINT = 0 ;
		IF ( .NOT.PRECO) THEN 
		    TRI_DIAG = 0.D0
		    TRI_DIAG(:,2) = 1.D0
		    GO TO 100
		ELSE 
		    IFLAG = 2 ;RETURN 
		ENDIF 
	ENDIF 

!V1.21
    YS= DOT_PRODUCT(W(IYPT+NPT+1:IYPT+NPT+N),W(ISPT+NPT+1:ISPT+NPT+N))

    !YS= DDOTP(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)

    IF( .NOT. PRECO) THEN
		! DEFAULTE DIAGONAL PRECONDITION MATRIX 
!V1.21 
        YY= DOT_PRODUCT(W(IYPT+NPT+1:IYPT+NPT+N),W(IYPT+NPT+1:IYPT+NPT+N))
        !YY= DDOTP(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
        tri_diag = 0.D0; tri_diag(:,2) =  YY/YS
    ELSE
        IFLAG=2; RETURN
		! KICK OUT OF THIS SUBROUTINE TO OBTAIN USER-DEFINED 
        !    PRE-CONDIDTION MATRIX AND JUMP TO 100 DIRECTLY 
    ENDIF
!	PRINT *, "PRE-CONDITION MATRIX IS PREPARED, NOW CHECK IT AND GO ON "
!----------------------------------------------------------


    100 CONTINUE

      IF(M==0) THEN 
          CALL TRIDAG_solver(tri_diag(:,1),tri_diag(:,2),tri_diag(:,3),&
                  -G(1:N), W(ISPT+1:ISPT+N),N); 
          GOTO 165
      ENDIF
!    IF(DIAGCO) THEN
!        DO 110 I=1,N
!            IF (DIAG(I) <= ZERO) GO TO 195
!        110 END DO
!    ENDIF

! COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
! "Updating quasi-Newton matrices with limited storage",
! Mathematics of Computation, Vol.24, No.151, pp. 773-782.
! ---------------------------------------------------------

    CP= POINT
    IF (POINT == 0) CP=M
    W(N+CP)= ONE/YS
    W(1:N)= -G(1:N)
    CP= POINT
    DO 125 I= 1,BOUND
        CP=CP-1
        IF (CP == -1)CP=M-1
!        SQ= DDOTP(N,W(ISPT+CP*N+1),1,W,1)
!V1.21
        SQ= DOT_PRODUCT(W(ISPT+CP*N+1:ISPT+CP*N+N),W(1:N))

        INMC=N+M+CP+1
        IYCN=IYPT+CP*N
        W(INMC)= W(N+CP+1)*SQ
!V1.21   
        W(1:N) = W(1:N) - W(INMC) * W(IYCN+1:IYCN+N);
        !CALL DAXPY(N,-W(INMC),W(IYCN+1),1,W,1)
    125 ENDDO

! SOLVE THE EUQATION R <- HK0*Q. Q IS STORED IN W(1:n), THEN W(1:N) IS UPDATED BY R    
!    DO 130 I=1,N
!        W(I)=diag(I)*W(I)
!    130 END DO
	if (PRECO)	 THEN 
		call TRIDAG_solver(tri_diag(:,1),tri_diag(:,2),tri_diag(:,3), &
                    W(1:N), WP(1:N), N) 
		W(1:N) = WP(1:N)
 	ELSE 	
      	W(1:n)=   W(1:n)  * (YS/ YY) 
	ENDIF 



    DO 145 I=1,BOUND
!        YR= DDOTP(N,W(IYPT+CP*N+1),1,W,1)
!V1.21
        YR= DOT_PRODUCT( W(IYPT+CP*N+1:IYPT+CP*N+N),W(1:N))
        
        BETA= W(N+CP+1)*YR
        INMC=N+M+CP+1
        BETA= W(INMC)-BETA
        ISCN=ISPT+CP*N
!V1.21  
        W(1:N) = W(1:N) + BETA*W(ISCN+1:ISCN+N)
        !CALL DAXPY(N,BETA,W(ISCN+1),1,W,1)
        CP=CP+1
        IF (CP == M)CP=0
    145 ENDDO

! STORE THE NEW SEARCH DIRECTION
! ------------------------------

    DO 160 I=1,N
        W(ISPT+POINT*N+I)= W(I)
    160 END DO
    

! OBTAIN THE ONE-DIMENSIONAL MINIMIZER OF THE FUNCTION
! BY USING THE LINE SEARCH ROUTINE MCSRCH
! ----------------------------------------------------
    165 continue 
    NFEV=0

!V1.1
!     STP=ONE
!    IF (ITER == 1) STP=STP1

!V1.2 ! change due to iflag 4 
     STP = ONE
     IF (ITER==0) STP = STP1

    IF(M>0) W(1:N)=G(1:N) 

    172 CONTINUE
	
    
    CALL MCSRCH(N,X,F,G,W(ISPT+POINT*N+1),STP, &
                FTOL, XTOL,MAXFEV,INFO,NFEV) 
 
	!---------------------------------------------------------------

    IF (INFO == -1) THEN
        IFLAG=1 
         !       print *, "Line Search unfinished, need evaluate f g"
         RETURN
		! LINE SEARCH UNFINISHED, NEED EVALUATIONS AT NEW TESTING STEPS
     ENDIF
     IF (INFO /= 1) GO TO 190 !   line search fail 
  
     NFUN= NFUN + NFEV
     !   print *, 'line search done, # of f=', nfev
 

! COMPUTE THE NEW STEP AND GRADIENT CHANGE
! -----------------------------------------
    NPT=POINT*N
    DO 175 I=1,N
        W(ISPT+NPT+I)= STP*W(ISPT+NPT+I)
        W(IYPT+NPT+I)= G(I)-W(I)
    175 END DO
    POINT=POINT+1
    IF (POINT == M)POINT=0

! TERMINATION TEST
! ----------------

!V0.0
    !GNORM= DSQRT(DDOTP(N,G,1,G,1))
    !XNORM= DSQRT(DDOTP(N,X,1,X,1))
!V1.21
    GNORM= DSQRT(DOT_PRODUCT(G,G))
    XNORM= DSQRT(DOT_PRODUCT(X,X))


    XNORM= DMax1(1.0D0,XNORM)    
    IF (GNORM/XNORM <= EPS) FINISH= .TRUE. 
!V1.2  ITER ++ only when new point updated by successful a few trial of line
!search 
     IFLAG =4 ; ITER= ITER+1
        

    IF(IPRINT(1) >= 0) CALL LB1(IPRINT,ITER,NFUN, &
                        GNORM,N,M,X,F,G,STP,FINISH)
    IF (FINISH) THEN
!        IFLAG=0
         IFLAG=3 ;  
         RETURN
        ! MODIFY TO DISTINGUISH INPUT IFLAG=0 AND 
        !    SUCCESSFUL RETURN IFALG=3
    ENDIF

!V1.2
     RETURN 
!V1.1    
!     GO TO 80

! ------------------------------------------------------------
! END OF MAIN ITERATION LOOP. ERROR EXITS.
! ------------------------------------------------------------

    190 IFLAG=-1
    IF(LP > 0) WRITE(LP,200) INFO
    RETURN
    195 IFLAG=-2
    IF(LP > 0) WRITE(LP,235) I
    RETURN
    196 IFLAG= -3
    IF(LP > 0) WRITE(LP,240)

! FORMATS
! -------

    200 FORMAT(/' IFLAG= -1 ',/' LINE SEARCH FAILED. SEE' &
    ' DOCUMENTATION OF ROUTINE MCSRCH',/' ERROR RETURN' &
    ' OF LINE SEARCH: INFO= ',I2,/ &
    ' POSSIBLE CAUSES: FUNCTION OR GRADIENT ARE INCORRECT',/, &
    ' OR INCORRECT TOLERANCES')
    235 FORMAT(/' IFLAG= -2',/' THE',I5,'-TH DIAGONAL ELEMENT OF THE',/, &
    ' INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE')
    240 FORMAT(/' IFLAG= -3',/' IMPROPER INPUT PARAMETERS (N OR M', &
    ' ARE NOT POSITIVE)')
    245 FORMAT(/'  GTOL IS LESS THAN OR EQUAL TO 1.D-04', &
    / ' IT HAS BEEN RESET TO 9.D-01')
    RETURN
    end SUBROUTINE LBFGS

! LAST LINE OF SUBROUTINE LBFGS


    SUBROUTINE LB1(IPRINT,ITER,NFUN, &
    GNORM,N,M,X,F,G,STP,FINISH)

! -------------------------------------------------------------
! THIS ROUTINE PRINTS MONITORING INFORMATION. THE FREQUENCY AND
! AMOUNT OF OUTPUT ARE CONTROLLED BY IPRINT.
! -------------------------------------------------------------

    INTEGER :: IPRINT(2),ITER,NFUN,N,M,i
    double precision :: X(N),G(N),F,GNORM,STP
    LOGICAL :: FINISH
  

    IF (ITER == 0)THEN
        WRITE(MP,10)
        IF (M>0) THEN 
        WRITE(MP,20) N,M
        ELSE
        WRITE(MP,21) N
        ENDIF 

        WRITE(MP,30)F,GNORM
        IF (IPRINT(2) >= 1)THEN
            WRITE(MP,40)
            WRITE(MP,50) (X(I),I=1,N)
            WRITE(MP,60)
            WRITE(MP,50) (G(I),I=1,N)
        ENDIF
        WRITE(MP,10)
        WRITE(MP,70)
    ELSE
        IF ((IPRINT(1) == 0) .AND. (ITER /= 1 .AND. .NOT. FINISH))RETURN
        IF (IPRINT(1) /= 0)THEN
            IF(MOD(ITER-1,IPRINT(1)) == 0 .OR. FINISH)THEN
                IF(IPRINT(2) > 1 .AND. ITER > 1) WRITE(MP,70)
                WRITE(MP,80)ITER,NFUN,F,GNORM,STP
            ELSE
                RETURN
            ENDIF
        ELSE
            IF( IPRINT(2) > 1 .AND. FINISH) WRITE(MP,70)
            WRITE(MP,80)ITER,NFUN,F,GNORM,STP
        ENDIF
        IF (IPRINT(2) == 2 .OR. IPRINT(2) == 3)THEN
            IF (FINISH)THEN
                WRITE(MP,90)
            ELSE
                WRITE(MP,40)
            ENDIF
            WRITE(MP,50)(X(I),I=1,N)
            IF (IPRINT(2) == 3)THEN
                WRITE(MP,60)
                WRITE(MP,50)(G(I),I=1,N)
            ENDIF
        ENDIF
        IF (FINISH) WRITE(MP,100)
    ENDIF

    10 FORMAT('*************************************************')
    20 FORMAT('  N=',I5,'   NUMBER OF CORRECTIONS=',I2, &
    /,  '       INITIAL VALUES')
    21 FORMAT('  N=',I5,'  M=0:|| STEEPEST DESCENT METHOD ||')
    30 FORMAT(' F= ',1PD12.5,'   GNORM= ',1PD12.5)
    40 FORMAT(' VECTOR X= ')
    50 FORMAT(6(2X,1PD12.5))
    60 FORMAT(' GRADIENT VECTOR G= ')
    70 FORMAT(/'   I   NFN',4X,'FUNC',8X,'GNORM',7X,'STEPLENGTH'/)
    80 FORMAT(2(I4,1X),3X,3(1PD12.5,2X))
    90 FORMAT(' FINAL POINT X= ')
    100 FORMAT(/' THE MINIMIZATION TERMINATED WITHOUT DETECTING ERRORS.', &
    /' IFLAG = 3')
!    /' IFLAG = 0')

    RETURN
    END SUBROUTINE LB1



! **************************
! LINE SEARCH ROUTINE MCSRCH
! **************************

    SUBROUTINE MCSRCH(N,X,F,G,S,STP,FTOL,XTOL,MAXFEV,INFO,NFEV)
    INTEGER :: N,MAXFEV,INFO,NFEV
    DOUBLE PRECISION :: F,STP,FTOL,XTOL
    DOUBLE PRECISION :: X(N),G(N),S(N)
   
    SAVE

! SUBROUTINE MCSRCH

! A slight modification of the subroutine CSRCH of More' and Thuente.
! The changes are to allow reverse communication, and do not affect
! the performance of the routine.

! THE PURPOSE OF MCSRCH IS TO FIND A STEP WHICH SATISFIES
! A SUFFICIENT DECREASE CONDITION AND A CURVATURE CONDITION.

! AT EACH STAGE THE SUBROUTINE UPDATES AN INTERVAL OF
! UNCERTAINTY WITH ENDPOINTS STX AND STY. THE INTERVAL OF
! UNCERTAINTY IS INITIALLY CHOSEN SO THAT IT CONTAINS A
! MINIMIZER OF THE MODIFIED FUNCTION

! F(X+STP*S) - F(X) - FTOL*STP*(GRADF(X)'S).

! IF A STEP IS OBTAINED FOR WHICH THE MODIFIED FUNCTION
! HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE DERIVATIVE,
! THEN THE INTERVAL OF UNCERTAINTY IS CHOSEN SO THAT IT
! CONTAINS A MINIMIZER OF F(X+STP*S).

! THE ALGORITHM IS DESIGNED TO FIND A STEP WHICH SATISFIES
! THE SUFFICIENT DECREASE CONDITION

! F(X+STP*S) .LE. F(X) + FTOL*STP*(GRADF(X)'S),

! AND THE CURVATURE CONDITION

! ABS(GRADF(X+STP*S)'S)) .LE. GTOL*ABS(GRADF(X)'S).

! IF FTOL IS LESS THAN GTOL AND IF, FOR EXAMPLE, THE FUNCTION
! IS BOUNDED BELOW, THEN THERE IS ALWAYS A STEP WHICH SATISFIES
! BOTH CONDITIONS. IF NO STEP CAN BE FOUND WHICH SATISFIES BOTH
! CONDITIONS, THEN THE ALGORITHM USUALLY STOPS WHEN ROUNDING
! ERRORS PREVENT FURTHER PROGRESS. IN THIS CASE STP ONLY
! SATISFIES THE SUFFICIENT DECREASE CONDITION.

! THE SUBROUTINE STATEMENT IS

! SUBROUTINE MCSRCH(N,X,F,G,S,STP,FTOL,XTOL, MAXFEV,INFO,NFEV,WA)
! WHERE

! N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
! OF VARIABLES.

! X IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
! BASE POINT FOR THE LINE SEARCH. ON OUTPUT IT CONTAINS
! X + STP*S.

! F IS A VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F
! AT X. ON OUTPUT IT CONTAINS THE VALUE OF F AT X + STP*S.

! G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
! GRADIENT OF F AT X. ON OUTPUT IT CONTAINS THE GRADIENT
! OF F AT X + STP*S.

! S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE
! SEARCH DIRECTION.

! STP IS A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN
! INITIAL ESTIMATE OF A SATISFACTORY STEP. ON OUTPUT
! STP CONTAINS THE FINAL ESTIMATE.

! FTOL AND GTOL ARE NONNEGATIVE INPUT VARIABLES. (In this reverse
! communication implementation GTOL is defined in a COMMON
! statement.) TERMINATION OCCURS WHEN THE SUFFICIENT DECREASE
! CONDITION AND THE DIRECTIONAL DERIVATIVE CONDITION ARE
! SATISFIED.

! XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS
! WHEN THE RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
! IS AT MOST XTOL.

! STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH
! SPECIFY LOWER AND UPPER BOUNDS FOR THE STEP. (In this reverse
! communication implementatin they are defined in a COMMON
! statement).

! MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION
! OCCURS WHEN THE NUMBER OF CALLS TO FCN IS AT LEAST
! MAXFEV BY THE END OF AN ITERATION.

! INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:

! INFO = 0  IMPROPER INPUT PARAMETERS.

! INFO =-1  A RETURN IS MADE TO COMPUTE THE FUNCTION AND GRADIENT.

! INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE
! DIRECTIONAL DERIVATIVE CONDITION HOLD.

! INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
! IS AT MOST XTOL.

! INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.

! INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN.

! INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX.

! INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
! THERE MAY NOT BE A STEP WHICH SATISFIES THE
! SUFFICIENT DECREASE AND CURVATURE CONDITIONS.
! TOLERANCES MAY BE TOO SMALL.

! NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF
! CALLS TO FCN.

! WA IS A WORK ARRAY OF LENGTH N.

! SUBPROGRAMS CALLED

! MCSTEP

! FORTRAN-SUPPLIED...ABS,MAX,MIN

! ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
! JORGE J. MORE', DAVID J. THUENTE

! **********
    INTEGER :: INFOC,J
    LOGICAL :: BRACKT,STAGE1
    double precision :: DG,DGM,DGINIT,DGTEST,DGX,DGXM,DGY,DGYM, &
    FINIT,FTEST1,FM,FX,FXM,FY,FYM,P5,P66,STX,STY, &
    STMIN,STMAX,WIDTH,WIDTH1,XTRAPF,ZERO
    DATA P5,P66,XTRAPF,ZERO /0.5D0,0.66D0,4.0D0,0.0D0/
    IF(INFO == -1) GO TO 45
    INFOC = 1

! CHECK THE INPUT PARAMETERS FOR ERRORS.

    IF (N <= 0 .OR. STP <= ZERO .OR. FTOL < ZERO .OR. &
    GTOL < ZERO .OR. XTOL < ZERO .OR. STPMIN < ZERO &
     .OR. STPMAX < STPMIN .OR. MAXFEV <= 0) RETURN

! COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION
! AND CHECK THAT S IS A DESCENT DIRECTION.

    DGINIT = ZERO

    DO 10 J = 1, N
        DGINIT = DGINIT + G(J)*S(J)
    10 ENDDO
    IF (DGINIT >= ZERO) then
        write(LP,15)
        15 FORMAT(/'  THE SEARCH DIRECTION IS NOT A DESCENT DIRECTION')
        RETURN
    ENDIF

! INITIALIZE LOCAL VARIABLES.

    BRACKT = .FALSE.
    STAGE1 = .TRUE.
    NFEV = 0
    FINIT = F
    DGTEST = FTOL*DGINIT
    WIDTH = STPMAX - STPMIN
    WIDTH1 = WIDTH/P5
    DO 20 J = 1, N
        WA(J) = X(J)
    20 ENDDO

! THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,
! FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
! THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,
! FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF
! THE INTERVAL OF UNCERTAINTY.
! THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP,
! FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.

    STX = ZERO
    FX = FINIT
    DGX = DGINIT
    STY = ZERO
    FY = FINIT
    DGY = DGINIT

! START OF ITERATION.

    30 CONTINUE

! SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND
! TO THE PRESENT INTERVAL OF UNCERTAINTY.

    IF (BRACKT) THEN
        STMIN = MIN(STX,STY)
        STMAX = MAX(STX,STY)
    ELSE
        STMIN = STX
        STMAX = STP + XTRAPF*(STP - STX)
    END IF

! FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.

    STP = MAX(STP,STPMIN)
    STP = MIN(STP,STPMAX)

! IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET
! STP BE THE LOWEST POINT OBTAINED SO FAR.

    IF ((BRACKT .AND. (STP <= STMIN .OR. STP >= STMAX)) &
     .OR. NFEV >= MAXFEV-1 .OR. INFOC == 0 &
     .OR. (BRACKT .AND. STMAX-STMIN <= XTOL*STMAX)) STP = STX

! EVALUATE THE FUNCTION AND GRADIENT AT STP
! AND COMPUTE THE DIRECTIONAL DERIVATIVE.
! We return to main program to obtain F and G.

    DO 40 J = 1, N
        X(J) = WA(J) + STP*S(J)
    40 ENDDO
    INFO=-1
    RETURN

! initial call info should not be -1 
! first call terminated, for second call to get function value of f(x)

    45 INFO=0
!  second call start from here after evaulation of vaule of funcion 
    NFEV = NFEV + 1

    DG = ZERO
    DO 50 J = 1, N
        DG = DG + G(J)*S(J)
    50 ENDDO
    FTEST1 = FINIT + STP*DGTEST
!  calcuate current DG and function value 

! TEST FOR CONVERGENCE.

    IF ((BRACKT .AND. (STP <= STMIN .OR. STP >= STMAX)) &
     .OR. INFOC == 0) INFO = 6
    IF (STP == STPMAX .AND. &
    F <= FTEST1 .AND. DG <= DGTEST) INFO = 5
    IF (STP == STPMIN .AND. &
    (F > FTEST1 .OR. DG >= DGTEST)) INFO = 4
    IF (NFEV >= MAXFEV) INFO = 3
    IF (BRACKT .AND. STMAX-STMIN <= XTOL*STMAX) INFO = 2
! the effect of XTOL in 1-d line search 
    IF (F <= FTEST1 .AND. ABS(DG) <= GTOL*(-DGINIT)) INFO = 1
!  desired result 

! CHECK FOR TERMINATION.

    IF (INFO /= 0) RETURN

! IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED
! FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.

    IF (STAGE1 .AND. F <= FTEST1 .AND. &
    DG >= MIN(FTOL,GTOL)*DGINIT) STAGE1 = .FALSE. 

! A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF
! WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED
! FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE
! DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN
! OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.

    IF (STAGE1 .AND. F <= FX .AND. F > FTEST1) THEN
    
    ! DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.
    
        FM = F - STP*DGTEST
        FXM = FX - STX*DGTEST
        FYM = FY - STY*DGTEST
        DGM = DG - DGTEST
        DGXM = DGX - DGTEST
        DGYM = DGY - DGTEST
    
    ! CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
    ! AND TO COMPUTE THE NEW STEP.
    
        CALL MCSTEP(STX,FXM,DGXM,STY,FYM,DGYM,STP,FM,DGM, &
        BRACKT,STMIN,STMAX,INFOC)
    
    ! RESET THE FUNCTION AND GRADIENT VALUES FOR F.
    
        FX = FXM + STX*DGTEST
        FY = FYM + STY*DGTEST
        DGX = DGXM + DGTEST
        DGY = DGYM + DGTEST
    ELSE
    
    ! CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
    ! AND TO COMPUTE THE NEW STEP.
    
        CALL MCSTEP(STX,FX,DGX,STY,FY,DGY,STP,F,DG, &
        BRACKT,STMIN,STMAX,INFOC)
    END IF

! FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE
! INTERVAL OF UNCERTAINTY.

    IF (BRACKT) THEN
        IF (ABS(STY-STX) >= P66*WIDTH1) &
        STP = STX + P5*(STY - STX)
        WIDTH1 = WIDTH
        WIDTH = ABS(STY-STX)
    END IF

! END OF ITERATION.

    GO TO 30

! LAST LINE OF SUBROUTINE MCSRCH.

    end SUBROUTINE MCSRCH



    SUBROUTINE MCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT, &
    STPMIN,STPMAX,INFO)
    INTEGER :: INFO
    double precision :: STX,FX,DX,STY,FY,DY,STP,FP,DP,STPMIN,STPMAX
    LOGICAL :: BRACKT,BOUND

! SUBROUTINE MCSTEP

! THE PURPOSE OF MCSTEP IS TO COMPUTE A SAFEGUARDED STEP FOR
! A LINESEARCH AND TO UPDATE AN INTERVAL OF UNCERTAINTY FOR
! A MINIMIZER OF THE FUNCTION.

! THE PARAMETER STX CONTAINS THE STEP WITH THE LEAST FUNCTION
! VALUE. THE PARAMETER STP CONTAINS THE CURRENT STEP. IT IS
! ASSUMED THAT THE DERIVATIVE AT STX IS NEGATIVE IN THE
! DIRECTION OF THE STEP. IF BRACKT IS SET TRUE THEN A
! MINIMIZER HAS BEEN BRACKETED IN AN INTERVAL OF UNCERTAINTY
! WITH ENDPOINTS STX AND STY.

! THE SUBROUTINE STATEMENT IS

! SUBROUTINE MCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT,
! STPMIN,STPMAX,INFO)

! WHERE

! STX, FX, AND DX ARE VARIABLES WHICH SPECIFY THE STEP,
! THE FUNCTION, AND THE DERIVATIVE AT THE BEST STEP OBTAINED
! SO FAR. THE DERIVATIVE MUST BE NEGATIVE IN THE DIRECTION
! OF THE STEP, THAT IS, DX AND STP-STX MUST HAVE OPPOSITE
! SIGNS. ON OUTPUT THESE PARAMETERS ARE UPDATED APPROPRIATELY.

! STY, FY, AND DY ARE VARIABLES WHICH SPECIFY THE STEP,
! THE FUNCTION, AND THE DERIVATIVE AT THE OTHER ENDPOINT OF
! THE INTERVAL OF UNCERTAINTY. ON OUTPUT THESE PARAMETERS ARE
! UPDATED APPROPRIATELY.

! STP, FP, AND DP ARE VARIABLES WHICH SPECIFY THE STEP,
! THE FUNCTION, AND THE DERIVATIVE AT THE CURRENT STEP.
! IF BRACKT IS SET TRUE THEN ON INPUT STP MUST BE
! BETWEEN STX AND STY. ON OUTPUT STP IS SET TO THE NEW STEP.

! BRACKT IS A LOGICAL VARIABLE WHICH SPECIFIES IF A MINIMIZER
! HAS BEEN BRACKETED. IF THE MINIMIZER HAS NOT BEEN BRACKETED
! THEN ON INPUT BRACKT MUST BE SET FALSE. IF THE MINIMIZER
! IS BRACKETED THEN ON OUTPUT BRACKT IS SET TRUE.

! STPMIN AND STPMAX ARE INPUT VARIABLES WHICH SPECIFY LOWER
! AND UPPER BOUNDS FOR THE STEP.

! INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
! IF INFO = 1,2,3,4,5, THEN THE STEP HAS BEEN COMPUTED
! ACCORDING TO ONE OF THE FIVE CASES BELOW. OTHERWISE
! INFO = 0, AND THIS INDICATES IMPROPER INPUT PARAMETERS.

! SUBPROGRAMS CALLED

! FORTRAN-SUPPLIED ... ABS,MAX,MIN,SQRT

! ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
! JORGE J. MORE', DAVID J. THUENTE

    double precision :: GAMMA,P,Q,R,S,SGND,STPC,STPF,STPQ,THETA
    INFO = 0

! CHECK THE INPUT PARAMETERS FOR ERRORS.

    IF ((BRACKT .AND. (STP <= MIN(STX,STY) .OR. &
    STP >= MAX(STX,STY))) .OR. &
    DX*(STP-STX) >= 0.0 .OR. STPMAX < STPMIN) RETURN

! DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN.

    SGND = DP*(DX/ABS(DX))

! FIRST CASE. A HIGHER FUNCTION VALUE.
! THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER
! TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN,
! ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN.

    IF (FP > FX) THEN
        INFO = 1
        BOUND = .TRUE.
        THETA = 3*(FX - FP)/(STP - STX) + DX + DP
        S = MAX(ABS(THETA),ABS(DX),ABS(DP))
        GAMMA = S*SQRT((THETA/S)**2 - (DX/S)*(DP/S))
        IF (STP < STX) GAMMA = -GAMMA
        P = (GAMMA - DX) + THETA
        Q = ((GAMMA - DX) + GAMMA) + DP
        R = P/Q
        STPC = STX + R*(STP - STX)
        STPQ = STX + ((DX/((FX-FP)/(STP-STX)+DX))/2)*(STP - STX)
        IF (ABS(STPC-STX) < ABS(STPQ-STX)) THEN
            STPF = STPC
        ELSE
            STPF = STPC + (STPQ - STPC)/2
        END IF
        BRACKT = .TRUE.
    
    ! SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF
    ! OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC
    ! STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP,
    ! THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN.
    
    ELSE IF (SGND < 0.0) THEN
        INFO = 2
        BOUND = .FALSE.
        THETA = 3*(FX - FP)/(STP - STX) + DX + DP
        S = MAX(ABS(THETA),ABS(DX),ABS(DP))
        GAMMA = S*SQRT((THETA/S)**2 - (DX/S)*(DP/S))
        IF (STP > STX) GAMMA = -GAMMA
        P = (GAMMA - DP) + THETA
        Q = ((GAMMA - DP) + GAMMA) + DX
        R = P/Q
        STPC = STP + R*(STX - STP)
        STPQ = STP + (DP/(DP-DX))*(STX - STP)
        IF (ABS(STPC-STP) > ABS(STPQ-STP)) THEN
            STPF = STPC
        ELSE
            STPF = STPQ
        END IF
        BRACKT = .TRUE.
    
    ! THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
    ! SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES.
    ! THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY
    ! IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC
    ! IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE
    ! EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO
    ! COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP
    ! CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN.
    
    ELSE IF (ABS(DP) < ABS(DX)) THEN
        INFO = 3
        BOUND = .TRUE.
        THETA = 3*(FX - FP)/(STP - STX) + DX + DP
        S = MAX(ABS(THETA),ABS(DX),ABS(DP))
    
    ! THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND
    ! TO INFINITY IN THE DIRECTION OF THE STEP.
    
        GAMMA = S*SQRT(MAX(0.0D0,(THETA/S)**2 - (DX/S)*(DP/S)))
        IF (STP > STX) GAMMA = -GAMMA
        P = (GAMMA - DP) + THETA
        Q = (GAMMA + (DX - DP)) + GAMMA
        R = P/Q
        IF (R < 0.0 .AND. GAMMA /= 0.0) THEN
            STPC = STP + R*(STX - STP)
        ELSE IF (STP > STX) THEN
            STPC = STPMAX
        ELSE
            STPC = STPMIN
        END IF
        STPQ = STP + (DP/(DP-DX))*(STX - STP)
        IF (BRACKT) THEN
            IF (ABS(STP-STPC) < ABS(STP-STPQ)) THEN
                STPF = STPC
            ELSE
                STPF = STPQ
            END IF
        ELSE
            IF (ABS(STP-STPC) > ABS(STP-STPQ)) THEN
                STPF = STPC
            ELSE
                STPF = STPQ
            END IF
        END IF
    
    ! FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
    ! SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES
    ! NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP
    ! IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN.
    
    ELSE
        INFO = 4
        BOUND = .FALSE.
        IF (BRACKT) THEN
            THETA = 3*(FP - FY)/(STY - STP) + DY + DP
            S = MAX(ABS(THETA),ABS(DY),ABS(DP))
            GAMMA = S*SQRT((THETA/S)**2 - (DY/S)*(DP/S))
            IF (STP > STY) GAMMA = -GAMMA
            P = (GAMMA - DP) + THETA
            Q = ((GAMMA - DP) + GAMMA) + DY
            R = P/Q
            STPC = STP + R*(STY - STP)
            STPF = STPC
        ELSE IF (STP > STX) THEN
            STPF = STPMAX
        ELSE
            STPF = STPMIN
        END IF
    END IF

! UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT
! DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE.

    IF (FP > FX) THEN
        STY = STP
        FY = FP
        DY = DP
    ELSE
        IF (SGND < 0.0) THEN
            STY = STX
            FY = FX
            DY = DX
        END IF
        STX = STP
        FX = FP
        DX = DP
    END IF

! COMPUTE THE NEW STEP AND SAFEGUARD IT.

    STPF = MIN(STPMAX,STPF)
    STPF = MAX(STPMIN,STPF)
    STP = STPF
    IF (BRACKT .AND. BOUND) THEN
        IF (STY > STX) THEN
            STP = MIN(STX+0.66*(STY-STX),STP)
        ELSE
            STP = MAX(STX+0.66*(STY-STX),STP)
        END IF
    END IF
    RETURN

! LAST LINE OF SUBROUTINE MCSTEP.

    end SUBROUTINE MCSTEP




!---------------------------------------------------------------
    SUBROUTINE tridag_solver(a,b,c,r,u,n)
! NUMERICAL RECIPES FORTRAN 77
    INTEGER,INTENT(IN) :: n
    DOUBLE PRECISION, INTENT(IN) ::   a(n), b(n), c(n), r(n)
    DOUBLE PRECISION, INTENT(OUT) ::  u(n)
    INTEGER :: j
    DOUBLE PRECISION ::  bet, gam(n)
! solve the linear tri-diagonal equation:
!        a(k)*u(k-1) + b(k)*u(k) + c(k)*u(k+1) = r(k) for k = 1..n   
! with  small modification at ends for u
! solves for a vector u(1:n) of length n the tridiagonal linear set given by equation (2.4.1).
! a(2:n)  is the low band ; c(1:n-1) is the upper band, b(i:n) is the diagonal
! If this happens then you should rewrite your equations as a set of order N . 1, with u2
! trivially eliminated.
    bet=b(1); u(1)=r(1)/bet
    DO   j=2,n  !Decomposition and forward substitution.
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet == 0.) pause 'tridag failed'   !Algorithm fails; see below.
        u(j)=(r(j)-a(j)*u(j-1))/bet
    ENDDO
    DO   j=n-1,1,-1 ! Backsubstitution.
        u(j)=u(j)-gam(j+1)*u(j+1)
    ENDDO
    RETURN
    END SUBROUTINE tridag_solver 


END MODULE MY_LBFGS_MODULE 

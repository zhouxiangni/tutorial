!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     sparse_matrix.f90 : implementation of the module sparse matrix
!

MODULE SPARSE_MATRIX_MODULE

	IMPLICIT NONE

	! public access for all subroutines
	PUBLIC::SPARSITY_PATTERN_INITIALIZE, &
		SPARSITY_PATTERN_ALLOCATE, &
		SPARSITY_PATTERN_DEALLOCATE, &
		SPARSITY_PATTERN_ADD_ENTRY, &
		SPARSITY_PATTERN_ADD_ENTRIES, &
		SPARSITY_PATTERN_COMPRESS, &
		SPARSITY_PATTERN_ASSIGNMENT, &
		SPARSITY_PATTERN_MULTIPLE, &
		SPARSITY_PATTERN_TRANSPOSE, &
		SPARSE_MATRIX_INITIALIZE, &
		SPARSE_MATRIX_ALLOCATE, &
		SPARSE_MATRIX_DEALLOCATE, &
		SPARSE_MATRIX_ADD_ENTRY, &
		SPARSE_MATRIX_ADD_ENTRIES, &
		SPARSE_MATRIX_ASSIGNMENT, &
		SPARSE_MATRIX_TRANSPOSE, &
		SPARSE_MATRIX_FULL, &
		SPARSE_MATRIX_GAUSS_SIDEL, &
		SPARSE_MATRIX_VECTOR_MULTIPLE, &
		SPARSE_MATRIX_VECTOR_MULT_P, &
		SPARSE_MATRIX_VECTOR_MULT_PS, &
		SPARSE_MATRIX_VECTOR_MULT_PSS, &
		SPARSE_MATRIX_VECTOR_T_MULTIPLE, &
		SPARSE_MATRIX_MULTIPLE

	! declaration of struct SPARSITY_PATTERN
	! The struct describe the non-zero structure for sparse matrix. The same
	! data structure as MATLAB is adopted here. The structure can be initialize
	! by given row, column and the maximum row length. After initialization,
	! corresponding memory are allocated for the sparsity pattern. The pattern
	! MUST be compressed before used to initialize a sparse matrix. The special
	! dealt is occured when the pattern is square, that the first entry for enery
	! row is diagonal entry. 
	TYPE, PUBLIC::SPARSITY_PATTERN
		INTEGER::N_ROW
		INTEGER::N_COLUMN
		INTEGER::N_NONZERO
		INTEGER::N_MAX_ROW_LENGTH
		INTEGER::N_MAX_VEC_LENGTH
		LOGICAL::IS_COMPRESSED
		INTEGER, DIMENSION(:), POINTER::ROW_START=>NULL()
		INTEGER, DIMENSION(:), POINTER::COLUMN_INDEX=>NULL()
	END TYPE SPARSITY_PATTERN

	! declaration of struct SPARSE_MATRIX
	! The struct describe a sparse matrix. The sparse matrix is associated with
	! a sparsity pattern. Only the values of the entries is stored in the structure.
	! The memory is allocated when the structure is initialized or the assignment
	! operator is called. DESTROY subroutine should be called to free its memory
	! after the matrix will be not used any more.
	TYPE, PUBLIC::SPARSE_MATRIX
		TYPE(SPARSITY_PATTERN), POINTER::PATTERN=>NULL()
		REAL, DIMENSION(:), POINTER::MATRIX_ENTRY=>NULL()
	END TYPE SPARSE_MATRIX

	CONTAINS

	! Initialize sparsity pattern, set the row, column and allocate memory
	! according maxinum row length
	SUBROUTINE SPARSITY_PATTERN_INITIALIZE(PATTERN, ROW, COL, MAXROWLEN)
		TYPE(SPARSITY_PATTERN), INTENT(INOUT)::PATTERN
		INTEGER, INTENT(IN)::ROW, COL
		INTEGER, INTENT(INOUT)::MAXROWLEN

		INTEGER::I, J
		
		CALL SPARSITY_PATTERN_DESTROY(PATTERN)

		IF (ROW == COL) THEN
			MAXROWLEN = MAX(1, MAXROWLEN)
		END IF

		PATTERN%N_ROW = ROW
		PATTERN%N_COLUMN = COL
		PATTERN%N_NONZERO = 0
		PATTERN%N_MAX_ROW_LENGTH = MAXROWLEN
		PATTERN%N_MAX_VEC_LENGTH = ROW*MAXROWLEN
		PATTERN%IS_COMPRESSED = .FALSE.
		ALLOCATE(PATTERN%ROW_START(1:ROW+1))
		ALLOCATE(PATTERN%COLUMN_INDEX(1:(ROW*MAXROWLEN)))
		
		DO I=1,ROW
			J = (I-1)*MAXROWLEN + 1
			PATTERN%ROW_START(I) = J
		END DO
		PATTERN%ROW_START(ROW+1) = ROW*MAXROWLEN + 1
		PATTERN%COLUMN_INDEX = -1

		IF (ROW == COL) THEN
			PATTERN%N_NONZERO = ROW
			DO I=1,ROW
				J = PATTERN%ROW_START(I)
				PATTERN%COLUMN_INDEX(J) = I
			END DO
		END IF
		
	END SUBROUTINE SPARSITY_PATTERN_INITIALIZE

	SUBROUTINE SPARSITY_PATTERN_DESTROY(PATTERN)
		TYPE(SPARSITY_PATTERN), INTENT(INOUT)::PATTERN
		
		IF (ASSOCIATED(PATTERN%ROW_START)) THEN
			DEALLOCATE(PATTERN%ROW_START)
		END IF
		IF (ASSOCIATED(PATTERN%COLUMN_INDEX)) THEN
			DEALLOCATE(PATTERN%COLUMN_INDEX)
		END IF
	END SUBROUTINE SPARSITY_PATTERN_DESTROY

	SUBROUTINE SPARSITY_PATTERN_ALLOCATE(PATTERN)
		TYPE(SPARSITY_PATTERN), POINTER::PATTERN
		
		ALLOCATE(PATTERN)
		NULLIFY(PATTERN%ROW_START)
		NULLIFY(PATTERN%COLUMN_INDEX)
	END SUBROUTINE SPARSITY_PATTERN_ALLOCATE
		
	SUBROUTINE SPARSITY_PATTERN_DEALLOCATE(PATTERN)
		TYPE(SPARSITY_PATTERN), POINTER::PATTERN
		
		CALL SPARSITY_PATTERN_DESTROY(PATTERN)
		DEALLOCATE(PATTERN)
	END SUBROUTINE SPARSITY_PATTERN_DEALLOCATE

	! Add a entry at (ROW, COL) for the sparsity pattern. Then operation MUST
	! be taken before the sparsity pattern is compressed
	SUBROUTINE SPARSITY_PATTERN_ADD_ENTRY(PATTERN, ROW, COL)
		TYPE(SPARSITY_PATTERN), INTENT(INOUT)::PATTERN
		INTEGER, INTENT(IN)::ROW, COL

		INTEGER::I, J
		IF (PATTERN%IS_COMPRESSED) THEN
			PRINT *, "Error: in SPARSITY_PATTERN_ADD_ENTRY, the sparsity pattern is compressed"
			READ *
		END IF 
			
		DO I=PATTERN%ROW_START(ROW),PATTERN%ROW_START(ROW+1)
			J = PATTERN%COLUMN_INDEX(I)
			IF (J == -1) THEN
				PATTERN%COLUMN_INDEX(I) = COL
				EXIT
			ELSE IF (J == COL) THEN
				EXIT
			END IF
		END DO
		IF (I == PATTERN%ROW_START(ROW+1)) THEN
			PRINT *, "Error: in SPARSITY_PATTERN_ADD_ENTRY, too much entry for the row"
			READ *
		END IF
	END SUBROUTINE SPARSITY_PATTERN_ADD_ENTRY

	! Add a series of entry
	SUBROUTINE SPARSITY_PATTERN_ADD_ENTRIES(PATTERN, ROW, COL)
		TYPE(SPARSITY_PATTERN), INTENT(INOUT)::PATTERN
		INTEGER, DIMENSION(:), INTENT(IN)::ROW, COL

		INTEGER::I, J, K, N
		
		N = SIZE(ROW)
		IF (PATTERN%IS_COMPRESSED) THEN
			PRINT *, "Error: in SPARSITY_PATTERN_ADD_ENTRY, the sparsity pattern is compressed"
			READ *
		END IF 

		DO K=1,N			
			DO I=PATTERN%ROW_START(ROW(K)),PATTERN%ROW_START(ROW(K)+1)
				J = PATTERN%COLUMN_INDEX(I)
				IF (J == -1) THEN
					PATTERN%COLUMN_INDEX(I) = COL(K)
					EXIT
				ELSE IF (J == COL(K)) THEN
					EXIT
				END IF
			END DO
			IF (I == PATTERN%ROW_START(ROW(K)+1)) THEN
				PRINT *, "Error: in SPARSITY_PATTERN_ADD_ENTRIES, too much entry for the row"
				READ *
			END IF
		END DO
	END SUBROUTINE SPARSITY_PATTERN_ADD_ENTRIES

	! Compress the sparsity pattern. Romove non-used memory in the sparsity
	! pattern to make the pattern in a compact formation. The pattern MUST
	! be taken before it is used to initialize a sparse matrix.
	SUBROUTINE SPARSITY_PATTERN_COMPRESS(PATTERN)
		TYPE(SPARSITY_PATTERN), INTENT(INOUT)::PATTERN

		INTEGER::I, J, K, M, N
		INTEGER, DIMENSION(:), POINTER::ROW_START_BUFFER=>NULL()
		INTEGER, DIMENSION(:), POINTER::COLUMN_INDEX_BUFFER=>NULL()
		
		I = 0
		M = 0
		DO J=1,PATTERN%N_ROW
			N = 0
			DO K=PATTERN%ROW_START(J), PATTERN%ROW_START(J+1)-1
				IF (PATTERN%COLUMN_INDEX(K) == -1) THEN
					EXIT
				ELSE
					N = N+1
				END IF
			END DO
			M = MAX(M, N)
			I = I+N
		END DO
		ROW_START_BUFFER=>PATTERN%ROW_START
		ALLOCATE(PATTERN%ROW_START(1:PATTERN%N_ROW+1))
		COLUMN_INDEX_BUFFER=>PATTERN%COLUMN_INDEX
		ALLOCATE(PATTERN%COLUMN_INDEX(1:I))
		PATTERN%N_NONZERO = I
		PATTERN%N_MAX_ROW_LENGTH = M
		PATTERN%N_MAX_VEC_LENGTH = I
		PATTERN%IS_COMPRESSED = .TRUE.
		PATTERN%ROW_START(1) = 1
		I = 0
		DO J=1,PATTERN%N_ROW
			DO K=ROW_START_BUFFER(J), ROW_START_BUFFER(J+1)-1
				IF (COLUMN_INDEX_BUFFER(K) == -1) THEN
					EXIT
				ELSE
					I = I+1
					PATTERN%COLUMN_INDEX(I) = COLUMN_INDEX_BUFFER(K)
				END IF
			END DO
			PATTERN%ROW_START(J+1) = I+1
		END DO
		DEALLOCATE(ROW_START_BUFFER)
		DEALLOCATE(COLUMN_INDEX_BUFFER)
	END SUBROUTINE SPARSITY_PATTERN_COMPRESS

	! Copy the sparsity pattern
	SUBROUTINE SPARSITY_PATTERN_ASSIGNMENT(P1, P2)
		TYPE(SPARSITY_PATTERN), INTENT(INOUT)::P1
		TYPE(SPARSITY_PATTERN), INTENT(IN)::P2

		CALL SPARSITY_PATTERN_DESTROY(P1)
		P1%N_ROW = P2%N_ROW
		P1%N_COLUMN = P2%N_COLUMN
		P1%N_NONZERO = P2%N_NONZERO
		P1%N_MAX_ROW_LENGTH = P2%N_MAX_ROW_LENGTH
		P1%N_MAX_VEC_LENGTH = P2%N_MAX_VEC_LENGTH
		P1%IS_COMPRESSED = P2%IS_COMPRESSED
		ALLOCATE(P1%ROW_START(1:P2%N_ROW+1))
		P1%ROW_START = P2%ROW_START
		ALLOCATE(P1%COLUMN_INDEX(1:P2%N_MAX_VEC_LENGTH))
		P1%COLUMN_INDEX = P2%COLUMN_INDEX
	END SUBROUTINE SPARSITY_PATTERN_ASSIGNMENT

	! Construct the sparsity pattern of the multiplication of two sparsity pattern
	SUBROUTINE SPARSITY_PATTERN_MULTIPLE(P1, P2, P3)
		TYPE(SPARSITY_PATTERN), INTENT(IN)::P1, P2
		TYPE(SPARSITY_PATTERN), INTENT(INOUT)::P3

		INTEGER::I, J, K, L, M
		INTEGER, DIMENSION(1:P2%N_COLUMN)::FLAG
		INTEGER, DIMENSION(1:P1%N_ROW)::N_ROW_NONZERO
		INTEGER::N_NONZERO, N_MAX_ROW_LENGTH

		CALL SPARSITY_PATTERN_DESTROY(P3)
		IF (P1%N_COLUMN /= P2%N_ROW) THEN
			PRINT *, "Error: in SPARSITY_PATTERN_MULTIPLE, sparsity pattern dimension not match"
			READ *
		END IF
		IF (.NOT. P1%IS_COMPRESSED) THEN
			PRINT *, "Error: in SPARSITY_PATTERN_MULTIPLE, sparsity pattern must be compressed"
			READ *
		END IF
		IF (.NOT. P2%IS_COMPRESSED) THEN
			PRINT *, "Error: in SPARSITY_PATTERN_MULTIPLE, sparsity pattern must be compressed"
			READ *
		END IF
		FLAG = 0
		N_ROW_NONZERO = 0
		N_NONZERO = 0
		DO I=1,P1%N_ROW
			DO J=P1%ROW_START(I),P1%ROW_START(I+1)-1
				L = P1%COLUMN_INDEX(J)
				DO K=P2%ROW_START(L),P2%ROW_START(L+1)-1
					M = P2%COLUMN_INDEX(K)
					IF (FLAG(M) /= I) THEN
						N_ROW_NONZERO(I) = N_ROW_NONZERO(I) + 1
						N_NONZERO = N_NONZERO + 1
						FLAG(M) = I
					END IF
				END DO
			END DO
		END DO
		ALLOCATE(P3%ROW_START(1:P1%N_ROW+1))
		ALLOCATE(P3%COLUMN_INDEX(1:N_NONZERO))
		P3%ROW_START(1) = 1
		N_MAX_ROW_LENGTH = 0
		N_NONZERO = 0
		FLAG = 0
		DO I=1,P1%N_ROW
			DO J=P1%ROW_START(I),P1%ROW_START(I+1)-1
				L = P1%COLUMN_INDEX(J)
				DO K=P2%ROW_START(L),P2%ROW_START(L+1)-1
					M = P2%COLUMN_INDEX(K)
					IF (FLAG(M) /= I) THEN
						N_NONZERO = N_NONZERO + 1
						P3%COLUMN_INDEX(N_NONZERO) = M
						FLAG(M) = I
					END IF
				END DO
			END DO
			N_MAX_ROW_LENGTH = MAX(N_MAX_ROW_LENGTH, N_ROW_NONZERO(I))
			P3%ROW_START(I+1) = P3%ROW_START(I) + N_ROW_NONZERO(I)
		END DO

		P3%N_ROW = P1%N_ROW
		P3%N_COLUMN = P2%N_COLUMN
		P3%N_NONZERO = N_NONZERO
		P3%N_MAX_ROW_LENGTH = N_MAX_ROW_LENGTH
		P3%N_MAX_VEC_LENGTH = N_NONZERO
		P3%IS_COMPRESSED = .TRUE.

		IF (P3%N_ROW == P3%N_COLUMN) THEN
			DO I=1,P3%N_ROW
				IF (P3%COLUMN_INDEX(P3%ROW_START(I)) == I) THEN
					CYCLE
				END IF
				DO J=P3%ROW_START(I)+1,P3%ROW_START(I+1)-1
					IF (P3%COLUMN_INDEX(J) == I) THEN
						P3%COLUMN_INDEX(J) = P3%COLUMN_INDEX(P3%ROW_START(I))
						P3%COLUMN_INDEX(P3%ROW_START(I)) = I
						EXIT
					END IF
				END DO
			END DO
		END IF
	END SUBROUTINE SPARSITY_PATTERN_MULTIPLE

	! Transpose the sparsity pattern
	SUBROUTINE SPARSITY_PATTERN_TRANSPOSE(P1, P2)
		TYPE(SPARSITY_PATTERN), INTENT(IN)::P1
		TYPE(SPARSITY_PATTERN), INTENT(INOUT)::P2

		INTEGER::I, J, K
		INTEGER, DIMENSION(1:P1%N_COLUMN)::N_ROW_NONZERO
		INTEGER::N_MAX_ROW_LENGTH

		CALL SPARSITY_PATTERN_DESTROY(P2)
		IF (.NOT. P1%IS_COMPRESSED) THEN
			PRINT *, "Error: in SPARSE_PATTERN_TRANSPOSE, sparsity pattern is not compressed"
			READ *
		END IF

		N_ROW_NONZERO = 0
		DO I=1,P1%N_ROW
			DO J=P1%ROW_START(I),P1%ROW_START(I+1)-1
				K = P1%COLUMN_INDEX(J)
				N_ROW_NONZERO(K) = N_ROW_NONZERO(K) + 1
			END DO
		END DO
		ALLOCATE(P2%ROW_START(1:P1%N_COLUMN+1))
		ALLOCATE(P2%COLUMN_INDEX(1:P1%N_MAX_VEC_LENGTH))
		P2%ROW_START(1) = 1
		N_MAX_ROW_LENGTH = 0
		DO I=1,P1%N_COLUMN
			N_MAX_ROW_LENGTH = MAX(N_MAX_ROW_LENGTH, N_ROW_NONZERO(I))
			P2%ROW_START(I+1) = P2%ROW_START(I) + N_ROW_NONZERO(I)
		END DO
		N_ROW_NONZERO = 0
		DO I=1,P1%N_ROW
			DO J=P1%ROW_START(I),P1%ROW_START(I+1)-1
				K = P1%COLUMN_INDEX(J)
				P2%COLUMN_INDEX(P2%ROW_START(K) + N_ROW_NONZERO(K)) = I
				N_ROW_NONZERO(K) = N_ROW_NONZERO(K) + 1
			END DO
		END DO

		P2%N_ROW = P1%N_COLUMN
		P2%N_COLUMN = P1%N_ROW
		P2%N_NONZERO = P1%N_NONZERO
		P2%N_MAX_ROW_LENGTH = N_MAX_ROW_LENGTH
		P2%N_MAX_VEC_LENGTH = P1%N_MAX_VEC_LENGTH
		P2%IS_COMPRESSED = .TRUE.

		IF (P2%N_ROW == P2%N_COLUMN) THEN
			DO I=1,P2%N_ROW
				IF (P2%COLUMN_INDEX(P2%ROW_START(I)) == I) THEN
					CYCLE
				END IF
				DO J=P2%ROW_START(I)+1,P2%ROW_START(I+1)-1
					IF (P2%COLUMN_INDEX(J) == I) THEN
						P2%COLUMN_INDEX(J) = P2%COLUMN_INDEX(P2%ROW_START(I))
						P2%COLUMN_INDEX(P2%ROW_START(I)) = I
						EXIT
					END IF
				END DO
			END DO
		END IF
	END SUBROUTINE SPARSITY_PATTERN_TRANSPOSE

	! Initialize a sparse matrix according the given sparsity pattern. It will allocate memory
	! to store the entry of the sparse matrix. The sparsity pattern MUST be a compressed one.
	SUBROUTINE SPARSE_MATRIX_INITIALIZE(MATRIX, PATTERN)
		TYPE(SPARSE_MATRIX), INTENT(INOUT)::MATRIX
		TYPE(SPARSITY_PATTERN), TARGET, INTENT(IN)::PATTERN

		IF (.NOT. PATTERN%IS_COMPRESSED) THEN
			PRINT *, "Error: in SPARSE_MATRIX_INITIALIZE, the sparsity pattern must be compressed " &
				//"before used to initialize sparse matrix"
			READ *
		END IF
		CALL SPARSE_MATRIX_DESTROY(MATRIX)
		MATRIX%PATTERN=>PATTERN
		ALLOCATE(MATRIX%MATRIX_ENTRY(1:PATTERN%N_MAX_VEC_LENGTH))
		MATRIX%MATRIX_ENTRY = 0.0E0
	END SUBROUTINE SPARSE_MATRIX_INITIALIZE

	SUBROUTINE SPARSE_MATRIX_DESTROY(MATRIX)
		TYPE(SPARSE_MATRIX), INTENT(INOUT)::MATRIX
		
		IF (ASSOCIATED(MATRIX%MATRIX_ENTRY)) THEN
			DEALLOCATE(MATRIX%MATRIX_ENTRY)
		END IF
	END SUBROUTINE SPARSE_MATRIX_DESTROY

	SUBROUTINE SPARSE_MATRIX_ALLOCATE(MATRIX)
		TYPE(SPARSE_MATRIX), POINTER::MATRIX
		
		ALLOCATE(MATRIX)
		NULLIFY(MATRIX%MATRIX_ENTRY)
	END SUBROUTINE SPARSE_MATRIX_ALLOCATE

	SUBROUTINE SPARSE_MATRIX_DEALLOCATE(MATRIX)
		TYPE(SPARSE_MATRIX), POINTER::MATRIX
		
		CALL SPARSE_MATRIX_DESTROY(MATRIX)
		DEALLOCATE(MATRIX)
	END SUBROUTINE SPARSE_MATRIX_DEALLOCATE

	! Add the value VAL on the entry at (ROW, COL) 
	SUBROUTINE SPARSE_MATRIX_ADD_ENTRY(MATRIX, ROW, COL, VAL)
		TYPE(SPARSE_MATRIX), INTENT(INOUT)::MATRIX
		INTEGER, INTENT(IN)::ROW, COL
		REAL, INTENT(IN)::VAL

		INTEGER::I
		DO I=MATRIX%PATTERN%ROW_START(ROW), MATRIX%PATTERN%ROW_START(ROW+1)
			IF (MATRIX%PATTERN%COLUMN_INDEX(I) == COL) THEN
				MATRIX%MATRIX_ENTRY(I) = MATRIX%MATRIX_ENTRY(I) + VAL
				EXIT
			END IF
		END DO
		IF (I == MATRIX%PATTERN%ROW_START(ROW+1)) THEN
			PRINT *, "Error: in SPARSE_MATRIX_ADD_ENTRY, no such entry in the matrix"
			READ *
		END IF
	END SUBROUTINE SPARSE_MATRIX_ADD_ENTRY

	! Add a series of value
	SUBROUTINE SPARSE_MATRIX_ADD_ENTRIES(MATRIX, ROW, COL, VAL, N)
		TYPE(SPARSE_MATRIX), INTENT(INOUT)::MATRIX
		INTEGER, INTENT(IN)::N
		INTEGER, DIMENSION(:), INTENT(IN)::ROW, COL
		REAL, DIMENSION(:), INTENT(IN)::VAL

		INTEGER::I, J
		DO I=1,N
			DO J=MATRIX%PATTERN%ROW_START(ROW(I)), MATRIX%PATTERN%ROW_START(ROW(I)+1)
				IF (MATRIX%PATTERN%COLUMN_INDEX(J) == COL(I)) THEN
					MATRIX%MATRIX_ENTRY(J) = MATRIX%MATRIX_ENTRY(J) + VAL(I)
					EXIT
				END IF
			END DO
			IF (J == MATRIX%PATTERN%ROW_START(ROW(I)+1)) THEN
				PRINT *, "Error: in SPARSE_MATRIX_ADD_ENTRIES, no such entry in the matrix"
				READ *
			END IF
		END DO
	END SUBROUTINE SPARSE_MATRIX_ADD_ENTRIES

	! Copy a sparse matrix
	SUBROUTINE SPARSE_MATRIX_ASSIGNMENT(M1, M2)
		TYPE(SPARSE_MATRIX), INTENT(INOUT)::M1
		TYPE(SPARSE_MATRIX), INTENT(IN)::M2

		CALL SPARSE_MATRIX_DESTROY(M1)
		M1%PATTERN=>M2%PATTERN
		ALLOCATE(M1%MATRIX_ENTRY(1:M2%PATTERN%N_MAX_VEC_LENGTH))
		M1%MATRIX_ENTRY = M2%MATRIX_ENTRY
	END SUBROUTINE SPARSE_MATRIX_ASSIGNMENT

	! Transpose a sparse matrix and construct its transposed sparsity pattern
	! at the same time
	SUBROUTINE SPARSE_MATRIX_TRANSPOSE(M1, M2, P2)
		TYPE(SPARSE_MATRIX), INTENT(IN)::M1
		TYPE(SPARSE_MATRIX), INTENT(INOUT)::M2
		TYPE(SPARSITY_PATTERN), TARGET, INTENT(INOUT)::P2

		INTEGER::I, J, K
		INTEGER, DIMENSION(1:M1%PATTERN%N_COLUMN)::N_ROW_NONZERO
		INTEGER::N_MAX_ROW_LENGTH
		REAL::D

		CALL SPARSE_MATRIX_DESTROY(M2)
		CALL SPARSITY_PATTERN_DESTROY(P2)

		N_ROW_NONZERO = 0
		DO I=1,M1%PATTERN%N_ROW
			DO J=M1%PATTERN%ROW_START(I),M1%PATTERN%ROW_START(I+1)-1
				K = M1%PATTERN%COLUMN_INDEX(J)
				N_ROW_NONZERO(K) = N_ROW_NONZERO(K) + 1
			END DO
		END DO
		ALLOCATE(P2%ROW_START(1:M1%PATTERN%N_COLUMN+1))
		ALLOCATE(P2%COLUMN_INDEX(1:M1%PATTERN%N_MAX_VEC_LENGTH))
		ALLOCATE(M2%MATRIX_ENTRY(1:M1%PATTERN%N_MAX_VEC_LENGTH))
		P2%ROW_START(1) = 1
		N_MAX_ROW_LENGTH = 0
		DO I=1,M1%PATTERN%N_COLUMN
			N_MAX_ROW_LENGTH = MAX(N_MAX_ROW_LENGTH, N_ROW_NONZERO(I))
			P2%ROW_START(I+1) = P2%ROW_START(I) + N_ROW_NONZERO(I)
		END DO
		N_ROW_NONZERO = 0
		DO I=1,M1%PATTERN%N_ROW
			DO J=M1%PATTERN%ROW_START(I),M1%PATTERN%ROW_START(I+1)-1
				K = M1%PATTERN%COLUMN_INDEX(J)
				P2%COLUMN_INDEX(P2%ROW_START(K) + N_ROW_NONZERO(K)) = I
				M2%MATRIX_ENTRY(P2%ROW_START(K) + N_ROW_NONZERO(K)) = M1%MATRIX_ENTRY(J)
				N_ROW_NONZERO(K) = N_ROW_NONZERO(K) + 1
			END DO
		END DO

		P2%N_ROW = M1%PATTERN%N_COLUMN
		P2%N_COLUMN = M1%PATTERN%N_ROW
		P2%N_NONZERO = M1%PATTERN%N_NONZERO
		P2%N_MAX_ROW_LENGTH = N_MAX_ROW_LENGTH
		P2%N_MAX_VEC_LENGTH = M1%PATTERN%N_MAX_VEC_LENGTH
		P2%IS_COMPRESSED = .TRUE.

		IF (P2%N_ROW == P2%N_COLUMN) THEN
			DO I=1,P2%N_ROW
				IF (P2%COLUMN_INDEX(P2%ROW_START(I)) == I) THEN
					CYCLE
				END IF
				DO J=P2%ROW_START(I)+1,P2%ROW_START(I+1)-1
					IF (P2%COLUMN_INDEX(J) == I) THEN
						D = M2%MATRIX_ENTRY(J)
						M2%MATRIX_ENTRY(J) = M2%MATRIX_ENTRY(P2%ROW_START(I))
						M2%MATRIX_ENTRY(P2%ROW_START(I)) = D
						P2%COLUMN_INDEX(J) = P2%COLUMN_INDEX(P2%ROW_START(I))
						P2%COLUMN_INDEX(P2%ROW_START(I)) = I
						EXIT
					END IF
				END DO
			END DO
		END IF

		M2%PATTERN=>P2
	END SUBROUTINE SPARSE_MATRIX_TRANSPOSE

	! Transform the sparse matrix to a full matrix
	SUBROUTINE SPARSE_MATRIX_FULL(M1, M2)
		TYPE(SPARSE_MATRIX), INTENT(IN)::M1
		REAL, DIMENSION(:,:), INTENT(INOUT)::M2

		INTEGER::I, J, K
		M2 = 0.0E0
		DO I=1,M1%PATTERN%N_ROW
			DO J=M1%PATTERN%ROW_START(I),M1%PATTERN%ROW_START(I+1)-1
				K = M1%PATTERN%COLUMN_INDEX(J)
				M2(I,K) = M1%MATRIX_ENTRY(J)
			END DO
		END DO
	END SUBROUTINE SPARSE_MATRIX_FULL

	! Gauss-Sidel iteration for certain times. The matrix is assumed to be square.
	SUBROUTINE SPARSE_MATRIX_GAUSS_SIDEL(M, X, R, S)
		TYPE(SPARSE_MATRIX), INTENT(IN)::M
		REAL, DIMENSION(:), INTENT(INOUT)::X
		REAL, DIMENSION(:), INTENT(IN)::R
		INTEGER, INTENT(IN)::S

		INTEGER::I, J, K, L
		REAL::TEMP
		DO I=1,S
			DO J=1,M%PATTERN%N_ROW
				TEMP = R(J)
				DO K=M%PATTERN%ROW_START(J)+1,M%PATTERN%ROW_START(J+1)-1
					L = M%PATTERN%COLUMN_INDEX(K)
					TEMP = TEMP - M%MATRIX_ENTRY(K) * X(L)
				END DO
				X(J) = TEMP / M%MATRIX_ENTRY(M%PATTERN%ROW_START(J))
			END DO
		END DO
	END SUBROUTINE SPARSE_MATRIX_GAUSS_SIDEL

	! Matrix vector multiplication: R = M * V
	SUBROUTINE SPARSE_MATRIX_VECTOR_MULTIPLE(M, V, R)
		TYPE(SPARSE_MATRIX), INTENT(IN)::M
		REAL, DIMENSION(:), INTENT(IN)::V
		REAL, DIMENSION(:), INTENT(INOUT)::R

		INTEGER::I, J
		R = 0.0E0
		DO I=1,M%PATTERN%N_ROW
			DO J=M%PATTERN%ROW_START(I),M%PATTERN%ROW_START(I+1)-1
				R(I) = R(I) + M%MATRIX_ENTRY(J) * V(M%PATTERN%COLUMN_INDEX(J))
			END DO
		END DO
	END SUBROUTINE SPARSE_MATRIX_VECTOR_MULTIPLE

	! R = R + M * V
	SUBROUTINE SPARSE_MATRIX_VECTOR_MULT_P(M, V, R)
		TYPE(SPARSE_MATRIX), INTENT(IN)::M
		REAL, DIMENSION(:), INTENT(IN)::V
		REAL, DIMENSION(:), INTENT(INOUT)::R

		INTEGER::I, J
		DO I=1,M%PATTERN%N_ROW
			DO J=M%PATTERN%ROW_START(I),M%PATTERN%ROW_START(I+1)-1
				R(I) = R(I) + M%MATRIX_ENTRY(J) * V(M%PATTERN%COLUMN_INDEX(J))
			END DO
		END DO
	END SUBROUTINE SPARSE_MATRIX_VECTOR_MULT_P

	! R = R + S * M * V
	SUBROUTINE SPARSE_MATRIX_VECTOR_MULT_PS(M, V, S, R)
		TYPE(SPARSE_MATRIX), INTENT(IN)::M
		REAL, DIMENSION(:), INTENT(IN)::V
		REAL, INTENT(IN)::S
		REAL, DIMENSION(:), INTENT(INOUT)::R

		INTEGER::I, J
		DO I=1,M%PATTERN%N_ROW
			DO J=M%PATTERN%ROW_START(I),M%PATTERN%ROW_START(I+1)-1
				R(I) = R(I) + S * M%MATRIX_ENTRY(J) * V(M%PATTERN%COLUMN_INDEX(J))
			END DO
		END DO
	END SUBROUTINE SPARSE_MATRIX_VECTOR_MULT_PS

	! R = T * R + S * M * V
	SUBROUTINE SPARSE_MATRIX_VECTOR_MULT_PSS(M, V, S, R, T)
		TYPE(SPARSE_MATRIX), INTENT(IN)::M
		REAL, DIMENSION(:), INTENT(IN)::V
		REAL, INTENT(IN)::S, T
		REAL, DIMENSION(:), INTENT(INOUT)::R

		INTEGER::I, J
		R = T*R
		DO I=1,M%PATTERN%N_ROW
			DO J=M%PATTERN%ROW_START(I),M%PATTERN%ROW_START(I+1)-1
				R(I) = R(I) + S * M%MATRIX_ENTRY(J) * V(M%PATTERN%COLUMN_INDEX(J))
			END DO
		END DO
	END SUBROUTINE SPARSE_MATRIX_VECTOR_MULT_PSS

	! Transpose matrix vector multiplication: R = (M^T) * V
	SUBROUTINE SPARSE_MATRIX_VECTOR_T_MULTIPLE(M, V, R)
		TYPE(SPARSE_MATRIX), INTENT(IN)::M
		REAL, DIMENSION(:), INTENT(IN)::V
		REAL, DIMENSION(:), INTENT(INOUT)::R

		INTEGER::I, J
		R = 0.0E0
		DO I=1,M%PATTERN%N_ROW
			DO J=M%PATTERN%ROW_START(I),M%PATTERN%ROW_START(I+1)-1
				R(M%PATTERN%COLUMN_INDEX(J)) = R(M%PATTERN%COLUMN_INDEX(J)) + M%MATRIX_ENTRY(J) * V(I)
			END DO
		END DO
	END SUBROUTINE SPARSE_MATRIX_VECTOR_T_MULTIPLE

	! R = R + (M^T) * V
	SUBROUTINE SPARSE_MATRIX_VECTOR_T_MULT_P(M, V, R)
		TYPE(SPARSE_MATRIX), INTENT(IN)::M
		REAL, DIMENSION(:), INTENT(IN)::V
		REAL, DIMENSION(:), INTENT(INOUT)::R

		INTEGER::I, J
		DO I=1,M%PATTERN%N_ROW
			DO J=M%PATTERN%ROW_START(I),M%PATTERN%ROW_START(I+1)-1
				R(M%PATTERN%COLUMN_INDEX(J)) = R(M%PATTERN%COLUMN_INDEX(J)) + M%MATRIX_ENTRY(J) * V(I)
			END DO
		END DO
	END SUBROUTINE SPARSE_MATRIX_VECTOR_T_MULT_P

	! R = R + S * (M^T) * V
	SUBROUTINE SPARSE_MATRIX_VECTOR_T_MULT_PS(M, V, S, R)
		TYPE(SPARSE_MATRIX), INTENT(IN)::M
		REAL, DIMENSION(:), INTENT(IN)::V
		REAL, INTENT(IN)::S
		REAL, DIMENSION(:), INTENT(INOUT)::R

		INTEGER::I, J
		DO I=1,M%PATTERN%N_ROW
			DO J=M%PATTERN%ROW_START(I),M%PATTERN%ROW_START(I+1)-1
				R(M%PATTERN%COLUMN_INDEX(J)) = R(M%PATTERN%COLUMN_INDEX(J)) + S*M%MATRIX_ENTRY(J) * V(I)
			END DO
		END DO
	END SUBROUTINE SPARSE_MATRIX_VECTOR_T_MULT_PS

	! R = T * R + S * (M^T) * V
	SUBROUTINE SPARSE_MATRIX_VECTOR_T_MULT_PSS(M, V, S, R, T)
		TYPE(SPARSE_MATRIX), INTENT(IN)::M
		REAL, DIMENSION(:), INTENT(IN)::V
		REAL, INTENT(IN)::S
		REAL, DIMENSION(:), INTENT(INOUT)::R
		REAL, INTENT(IN)::T

		INTEGER::I, J
		R = T*R
		DO I=1,M%PATTERN%N_ROW
			DO J=M%PATTERN%ROW_START(I),M%PATTERN%ROW_START(I+1)-1
				R(M%PATTERN%COLUMN_INDEX(J)) = R(M%PATTERN%COLUMN_INDEX(J)) + S*M%MATRIX_ENTRY(J) * V(I)
			END DO
		END DO
	END SUBROUTINE SPARSE_MATRIX_VECTOR_T_MULT_PSS

	! Sparse matrix multiplication. The sparsity pattern of M3 MUST be the multiplication
	! of the sparsity pattern of M1 and M2
	SUBROUTINE SPARSE_MATRIX_MULTIPLE(M1, M2, M3)
		TYPE(SPARSE_MATRIX), INTENT(IN)::M1, M2
		TYPE(SPARSE_MATRIX), INTENT(INOUT)::M3

		INTEGER::I, J, K, L, M
		REAL, DIMENSION(1:M2%PATTERN%N_COLUMN)::A

		A = 0.0E0
		DO I=1,M1%PATTERN%N_ROW
			DO J=M1%PATTERN%ROW_START(I),M1%PATTERN%ROW_START(I+1)-1
				L = M1%PATTERN%COLUMN_INDEX(J)
				DO K=M2%PATTERN%ROW_START(L),M2%PATTERN%ROW_START(L+1)-1
					M = M2%PATTERN%COLUMN_INDEX(K)
					A(M) = A(M) + M1%MATRIX_ENTRY(J) * M2%MATRIX_ENTRY(K)
				END DO
			END DO
			DO J=M3%PATTERN%ROW_START(I),M3%PATTERN%ROW_START(I+1)-1
				L = M3%PATTERN%COLUMN_INDEX(J)
				M3%MATRIX_ENTRY(J) = A(L)
				A(L) = 0.0E0
			END DO
		END DO
	END SUBROUTINE SPARSE_MATRIX_MULTIPLE

! the following code is added later for further utility

	subroutine sparsity_pattern_h_concat(p1, p2, p)
		type(sparsity_pattern), intent(in)::p1
		type(sparsity_pattern), intent(in)::p2
		type(sparsity_pattern), intent(inout)::p

		integer::i, j, m, n
		
		if (p1%n_row /= p2%n_row) then
			print *, "in sparsity_pattern_h_concat: number of row not match"
			read *
		end if
		if ((.not. p1%is_compressed) .or. (.not. p2%is_compressed)) then
			print *, "in sparsity_pattern_h_concat: sparsity pattern needed to be compressed"
			read *
		end if
		
		call sparsity_pattern_destroy(p)
		p%n_row = p1%n_row
		p%n_column = p1%n_column + p2%n_column
		p%n_nonzero = p1%n_nonzero + p2%n_nonzero
		p%n_max_vec_length = p1%n_max_vec_length + p2%n_max_vec_length
		p%is_compressed = .true.
		allocate(p%row_start(1:p%n_row+1))
		allocate(p%column_index(1:p%n_nonzero))
		p%n_max_row_length = 0
		p%row_start(1) = 1
		do i=1,p%n_row
			m = p1%row_start(i+1) - p1%row_start(i)
			n = p2%row_start(i+1) - p2%row_start(i)
			p%row_start(i+1) = p%row_start(i)+m+n
			p%n_max_row_length = max(p%n_max_row_length,m+n)
			p%column_index(p%row_start(i):p%row_start(i+1)-1) = (/p1%column_index(p1%row_start(i):p1%row_start(i+1)-1), &
				p1%n_column + p2%column_index(p2%row_start(i):p2%row_start(i+1)-1)/)
		end do
		p%row_start(p%n_row+1) = p%n_nonzero+1

		if (p%n_row == p%n_column) then
			do i=1,p%n_row
				if (p%column_index(p%row_start(i)) == i) then
					cycle
				end if
				do j=p%row_start(i)+1,p%row_start(i+1)-1
					if (p%column_index(j) == i) then
						p%column_index(j) = p%column_index(p%row_start(i))
						p%column_index(p%row_start(i)) = i
						exit
					end if
				end do
			end do
		end if
	end subroutine sparsity_pattern_h_concat
			
	subroutine sparsity_pattern_v_concat(p1, p2, p)
		type(sparsity_pattern), intent(in)::p1
		type(sparsity_pattern), intent(in)::p2
		type(sparsity_pattern), intent(inout)::p

		integer::i, j
		
		if (p1%n_column /= p2%n_column) then
			print *, "in sparsity_pattern_v_concat: number of row not match"
			read *
		end if
		if ((.not. p1%is_compressed) .or. (.not. p2%is_compressed)) then
			print *, "in sparsity_pattern_v_concat: sparsity pattern needed to be compressed"
			read *
		end if
		
		call sparsity_pattern_destroy(p)
		p%n_row = p1%n_row + p2%n_row
		p%n_column = p1%n_column
		p%n_nonzero = p1%n_nonzero + p2%n_nonzero
		p%n_max_row_length = max(p1%n_max_row_length, p2%n_max_row_length)
		p%n_max_vec_length = p1%n_max_vec_length + p2%n_max_vec_length
		p%is_compressed = .true.
		p%row_start(p%n_row+1) = p%n_nonzero+1

		allocate(p%row_start(1:p%n_row+1))
		allocate(p%column_index(1:p%n_nonzero))
		p%row_start(1:p1%n_row+1) = p1%row_start
		p%row_start(p1%n_row+1:p%n_row+1) = p1%n_nonzero + p2%row_start
		p%column_index(1:p1%n_nonzero) = p1%column_index
		p%column_index(p1%n_nonzero+1:p%n_nonzero) = p2%column_index

		if (p%n_row == p%n_column) then
			do i=1,p%n_row
				if (p%column_index(p%row_start(i)) == i) then
					cycle
				end if
				do j=p%row_start(i)+1,p%row_start(i+1)-1
					if (p%column_index(j) == i) then
						p%column_index(j) = p%column_index(p%row_start(i))
						p%column_index(p%row_start(i)) = i
						exit
					end if
				end do
			end do
		end if
	end subroutine sparsity_pattern_v_concat

	subroutine sparsity_pattern_d_concat(p1, p2, p)
		type(sparsity_pattern), intent(in)::p1
		type(sparsity_pattern), intent(in)::p2
		type(sparsity_pattern), intent(inout)::p

		integer::i, j
		
		if ((.not. p1%is_compressed) .or. (.not. p2%is_compressed)) then
			print *, "in sparsity_pattern_d_concat: sparsity pattern needed to be compressed"
			read *
		end if
		
		call sparsity_pattern_destroy(p)
		p%n_row = p1%n_row + p2%n_row
		p%n_column = p1%n_column + p2%n_column
		p%n_nonzero = p1%n_nonzero + p2%n_nonzero
		p%n_max_row_length = max(p1%n_max_row_length, p2%n_max_row_length)
		p%n_max_vec_length = p1%n_max_vec_length + p2%n_max_vec_length
		p%is_compressed = .true.

		allocate(p%row_start(1:p%n_row+1))
		allocate(p%column_index(1:p%n_nonzero))
		p%row_start(1:p1%n_row+1) = p1%row_start
		p%row_start(p1%n_row+1:p%n_row+1) = p1%n_nonzero + p2%row_start
		p%column_index(1:p1%n_nonzero) = p1%column_index
		p%column_index(p1%n_nonzero+1:p%n_nonzero) = p1%n_column + p2%column_index

		if (p%n_row == p%n_column) then
			do i=1,p%n_row
				if (p%column_index(p%row_start(i)) == i) then
					cycle
				end if
				do j=p%row_start(i)+1,p%row_start(i+1)-1
					if (p%column_index(j) == i) then
						p%column_index(j) = p%column_index(p%row_start(i))
						p%column_index(p%row_start(i)) = i
						exit
					end if
				end do
			end do
		end if
	end subroutine sparsity_pattern_d_concat
	
	subroutine sparse_matrix_h_concat(m1, m2, m, p)
		type(sparse_matrix), intent(in)::m1
		type(sparse_matrix), intent(in)::m2
		type(sparse_matrix), intent(inout)::m
		type(sparsity_pattern), target, intent(inout)::p
		
		integer::i, j, n1, n2
		real::d
		type(sparsity_pattern), pointer::p1=>NULL()
		type(sparsity_pattern), pointer::p2=>NULL()
		
		p1 => m1%pattern
		p2 => m2%pattern
		if (p1%n_row /= p2%n_row) then
			print *, "in sparse_matrix_h_concat: number of row not match"
			read *
		end if
		call sparsity_pattern_destroy(p)
		call sparse_matrix_destroy(m)
		
		p%n_row = p1%n_row
		p%n_column = p1%n_column + p2%n_column
		p%n_nonzero = p1%n_nonzero + p2%n_nonzero
		p%n_max_vec_length = p1%n_max_vec_length + p2%n_max_vec_length
		p%is_compressed = .true.
		m%pattern => p

		allocate(p%row_start(1:p%n_row+1))
		allocate(p%column_index(1:p%n_nonzero))
		allocate(m%matrix_entry(1:p%n_nonzero))
		p%n_max_row_length = 0
		p%row_start(1) = 1
		do i=1,p%n_row
			n1 = p1%row_start(i+1) - p1%row_start(i)
			n2 = p2%row_start(i+1) - p2%row_start(i)
			p%row_start(i+1) = p%row_start(i)+n1+n2
			p%n_max_row_length = max(p%n_max_row_length,n1+n2)
			p%column_index(p%row_start(i):p%row_start(i+1)-1) = (/p1%column_index(p1%row_start(i):p1%row_start(i+1)-1), &
				p1%n_column + p2%column_index(p2%row_start(i):p2%row_start(i+1)-1)/)
			m%matrix_entry(p%row_start(i):p%row_start(i+1)-1) = (/m1%matrix_entry(p1%row_start(i):p1%row_start(i+1)-1), &
				m2%matrix_entry(p2%row_start(i):p2%row_start(i+1)-1)/)
		end do

		if (p%n_row == p%n_column) then
			do i=1,p%n_row
				if (p%column_index(p%row_start(i)) == i) then
					cycle
				end if
				do j=p%row_start(i)+1,p%row_start(i+1)-1
					if (p%column_index(j) == i) then
						d = m%matrix_entry(j)
						m%matrix_entry(j) = m%matrix_entry(p%row_start(i))
						m%matrix_entry(p%row_start(i)) = d
						p%column_index(j) = p%column_index(p%row_start(i))
						p%column_index(p%row_start(i)) = i
						exit
					end if
				end do
			end do
		end if
	end subroutine sparse_matrix_h_concat
	
	subroutine sparse_matrix_v_concat(m1, m2, m, p)
		type(sparse_matrix), intent(in)::m1
		type(sparse_matrix), intent(in)::m2
		type(sparse_matrix), intent(inout)::m
		type(sparsity_pattern), target, intent(inout)::p
		
		integer::i, j
		real::d
		type(sparsity_pattern), pointer::p1=>NULL()
		type(sparsity_pattern), pointer::p2=>NULL()
		
		p1 => m1%pattern
		p2 => m2%pattern
		if (p1%n_column /= p2%n_column) then
			print *, "in sparse_matrix_v_concat: number of column not match"
			read *
		end if
		call sparsity_pattern_destroy(p)
		call sparse_matrix_destroy(m)

		p%n_row = p1%n_row + p2%n_row
		p%n_column = p1%n_column
		p%n_nonzero = p1%n_nonzero + p2%n_nonzero
		p%n_max_row_length = max(p1%n_max_row_length, p2%n_max_row_length)
		p%n_max_vec_length = p1%n_max_vec_length + p2%n_max_vec_length
		p%is_compressed = .true.
		p%row_start(p%n_row+1) = p%n_nonzero+1
		m%pattern => p

		allocate(p%row_start(1:p%n_row+1))
		allocate(p%column_index(1:p%n_nonzero))
		allocate(m%matrix_entry(1:p%n_nonzero))
		p%row_start(1:p1%n_row+1) = p1%row_start
		p%row_start(p1%n_row+1:p%n_row+1) = p1%n_nonzero + p2%row_start
		p%column_index(1:p1%n_nonzero) = p1%column_index
		p%column_index(p1%n_nonzero+1:p%n_nonzero) = p2%column_index
		m%matrix_entry = (/m1%matrix_entry, m2%matrix_entry/)		

		if (p%n_row == p%n_column) then
			do i=1,p%n_row
				if (p%column_index(p%row_start(i)) == i) then
					cycle
				end if
				do j=p%row_start(i)+1,p%row_start(i+1)-1
					if (p%column_index(j) == i) then
						d = m%matrix_entry(j)
						m%matrix_entry(j) = m%matrix_entry(p%row_start(i))
						m%matrix_entry(p%row_start(i)) = d
						p%column_index(j) = p%column_index(p%row_start(i))
						p%column_index(p%row_start(i)) = i
						exit
					end if
				end do
			end do
		end if
	end subroutine sparse_matrix_v_concat
	
	subroutine sparse_matrix_d_concat(m1, m2, m, p)
		type(sparse_matrix), intent(in)::m1
		type(sparse_matrix), intent(in)::m2
		type(sparse_matrix), intent(inout)::m
		type(sparsity_pattern), target, intent(inout)::p
		
		integer::i, j
		real::d
		type(sparsity_pattern), pointer::p1=>NULL()
		type(sparsity_pattern), pointer::p2=>NULL()
		
		p1 => m1%pattern
		p2 => m2%pattern
		call sparsity_pattern_destroy(p)
		call sparse_matrix_destroy(m)

		p%n_row = p1%n_row + p2%n_row
		p%n_column = p1%n_column + p2%n_column
		p%n_nonzero = p1%n_nonzero + p2%n_nonzero
		p%n_max_row_length = max(p1%n_max_row_length, p2%n_max_row_length)
		p%n_max_vec_length = p1%n_max_vec_length + p2%n_max_vec_length
		p%is_compressed = .true.
		m%pattern => p

		allocate(p%row_start(1:p%n_row+1))
		allocate(p%column_index(1:p%n_nonzero))
		allocate(m%matrix_entry(1:p%n_nonzero))
		p%row_start(1:p1%n_row+1) = p1%row_start
		p%row_start(p1%n_row+1:p%n_row+1) = p1%n_nonzero + p2%row_start
		p%column_index(1:p1%n_nonzero) = p1%column_index
		p%column_index(p1%n_nonzero+1:p%n_nonzero) = p1%n_column + p2%column_index
		m%matrix_entry = (/m1%matrix_entry, m2%matrix_entry/)		

		if (p%n_row == p%n_column) then
			do i=1,p%n_row
				if (p%column_index(p%row_start(i)) == i) then
					cycle
				end if
				do j=p%row_start(i)+1,p%row_start(i+1)-1
					if (p%column_index(j) == i) then
						d = m%matrix_entry(j)
						m%matrix_entry(j) = m%matrix_entry(p%row_start(i))
						m%matrix_entry(p%row_start(i)) = d
						p%column_index(j) = p%column_index(p%row_start(i))
						p%column_index(p%row_start(i)) = i
						exit
					end if
				end do
			end do
		end if
	end subroutine sparse_matrix_d_concat

	subroutine sparsity_pattern_plus(p1, p2, p)
		type(sparsity_pattern), intent(in)::p1
		type(sparsity_pattern), intent(in)::p2
		type(sparsity_pattern), intent(inout)::p

		integer::i, j, m, n
		integer, dimension(1:p1%n_column)::flag
		
		if (p1%n_row /= p2%n_row .or. p1%n_column /= p2%n_column) then
			print *, "in sparsity_pattern_plus: sparsity pattern not match"
			read *
		end if
		if ((.not. p1%is_compressed) .or. (.not. p2%is_compressed)) then
			print *, "in sparsity_pattern_plus: sparsity pattern needed to be compressed"
			read *
		end if
		
		call sparsity_pattern_destroy(p)
		p%n_row = p1%n_row
		p%n_column = p1%n_column
		p%is_compressed = .true.
		p%n_max_row_length = 0
		m = 0
		flag = 0
		do i=1,p%n_row
			n = 0
			do j=p1%row_start(i),p1%row_start(i+1)-1
				if (flag(p1%column_index(j)) /= i) then
					n = n+1
					flag(p1%column_index(j)) = i
				end if
			end do
			do j=p2%row_start(i),p2%row_start(i+1)-1
				if (flag(p2%column_index(j)) /= i) then
					n = n+1
					flag(p2%column_index(j)) = i
				end if
			end do
			p%n_max_row_length = max(p%n_max_row_length, n)
			m = m+n
		end do
		p%n_nonzero = m
		p%n_max_vec_length = m
		allocate(p%row_start(1:p%n_row+1))
		allocate(p%column_index(1:m))
		n = 1
		flag = 0
		p%row_start(1) = 1
		do i=1,p%n_row
			do j=p1%row_start(i),p1%row_start(i+1)-1
				if (flag(p1%column_index(j)) /= i) then
					p%column_index(n) = p1%column_index(j)
					n = n+1
					flag(p1%column_index(j)) = i
				end if
			end do
			do j=p2%row_start(i),p2%row_start(i+1)-1
				if (flag(p2%column_index(j)) /= i) then
					p%column_index(n) = p2%column_index(j)
					n = n+1
					flag(p2%column_index(j)) = i
				end if
			end do
			p%row_start(i+1) = n
		end do

		if (p%n_row == p%n_column) then
			do i=1,p%n_row
				if (p%column_index(p%row_start(i)) == i) then
					cycle
				end if
				do j=p%row_start(i)+1,p%row_start(i+1)-1
					if (p%column_index(j) == i) then
						p%column_index(j) = p%column_index(p%row_start(i))
						p%column_index(p%row_start(i)) = i
						exit
					end if
				end do
			end do
		end if
	end subroutine sparsity_pattern_plus
			
	subroutine sparse_matrix_plus(m1, m2, m, p)
		type(sparse_matrix), intent(in)::m1
		type(sparse_matrix), intent(in)::m2
		type(sparse_matrix), intent(inout)::m
		type(sparsity_pattern), target, optional, intent(inout)::p
		
		integer::i, j, n1, n2
		real::d
		integer, dimension(1:m1%pattern%n_column)::flag
		real, dimension(1:m1%pattern%n_column)::row_entry
		type(sparsity_pattern), pointer::p1=>NULL()
		type(sparsity_pattern), pointer::p2=>NULL()
		
		p1 => m1%pattern
		p2 => m2%pattern
		if (.not. present(p)) then
			call sparse_matrix_initialize(m, p1)
			m%matrix_entry = m1%matrix_entry + m2%matrix_entry
			return
		end if
		if (p1%n_row /= p2%n_row .or. p1%n_column /= p2%n_column) then
			print *, "in sparsity_matrix_plus: sparse matrix dimension not match"
			read *
		end if

		call sparsity_pattern_destroy(p)
		call sparse_matrix_destroy(m)
		p%n_row = p1%n_row
		p%n_column = p1%n_column
		p%is_compressed = .true.
		p%n_max_row_length = 0
		n1 = 0
		flag = 0
		do i=1,p%n_row
			n2 = 0
			do j=p1%row_start(i),p1%row_start(i+1)-1
				if (flag(p1%column_index(j)) /= i) then
					n2 = n2+1
					flag(p1%column_index(j)) = i
				end if
			end do
			do j=p2%row_start(i),p2%row_start(i+1)-1
				if (flag(p2%column_index(j)) /= i) then
					n2 = n2+1
					flag(p2%column_index(j)) = i
				end if
			end do
			p%n_max_row_length = max(p%n_max_row_length, n2)
			n1 = n1+n2
		end do
		p%n_nonzero = n1
		p%n_max_vec_length = n1
		allocate(p%row_start(1:p%n_row+1))
		allocate(p%column_index(1:n1))
		allocate(m%matrix_entry(1:n1))
		n2 = 1
		flag = 0
		p%row_start(1) = 1
		row_entry = 0.0e0
		do i=1,p%n_row
			do j=p1%row_start(i),p1%row_start(i+1)-1
				if (flag(p1%column_index(j)) /= i) then
					p%column_index(n2) = p1%column_index(j)
					n2 = n2+1
					flag(p1%column_index(j)) = i
				end if
				row_entry(p1%column_index(j)) = row_entry(p1%column_index(j)) + m1%matrix_entry(j)
			end do
			do j=p2%row_start(i),p2%row_start(i+1)-1
				if (flag(p2%column_index(j)) /= i) then
					p%column_index(n2) = p2%column_index(j)
					n2 = n2+1
					flag(p2%column_index(j)) = i
				end if
				row_entry(p2%column_index(j)) = row_entry(p2%column_index(j)) + m2%matrix_entry(j)
			end do
			p%row_start(i+1) = n2
			do j=p%row_start(i),p%row_start(i+1)-1
				m%matrix_entry(j) = row_entry(p%column_index(j))
				row_entry(p%column_index(j)) = 0.0e0
			end do
		end do

		if (p%n_row == p%n_column) then
			do i=1,p%n_row
				if (p%column_index(p%row_start(i)) == i) then
					cycle
				end if
				do j=p%row_start(i)+1,p%row_start(i+1)-1
					if (p%column_index(j) == i) then
						d = m%matrix_entry(j)
						m%matrix_entry(j) = m%matrix_entry(p%row_start(i))
						m%matrix_entry(p%row_start(i)) = d
						p%column_index(j) = p%column_index(p%row_start(i))
						p%column_index(p%row_start(i)) = i
						exit
					end if
				end do
			end do
		end if
	end subroutine sparse_matrix_plus

	subroutine sparse_matrix_minus(m1, m2, m, p)
		type(sparse_matrix), intent(in)::m1
		type(sparse_matrix), intent(in)::m2
		type(sparse_matrix), intent(inout)::m
		type(sparsity_pattern), target, optional, intent(inout)::p
		
		integer::i, j, n1, n2
		real::d
		integer, dimension(1:m1%pattern%n_column)::flag
		real, dimension(1:m1%pattern%n_column)::row_entry
		type(sparsity_pattern), pointer::p1=>NULL()
		type(sparsity_pattern), pointer::p2=>NULL()
		
		p1 => m1%pattern
		p2 => m2%pattern
		if (.not. present(p)) then
			call sparse_matrix_initialize(m, p1)
			m%matrix_entry = m1%matrix_entry - m2%matrix_entry
			return
		end if
		if (p1%n_row /= p2%n_row .or. p1%n_column /= p2%n_column) then
			print *, "in sparsity_matrix_plus: sparse matrix dimension not match"
			read *
		end if

		call sparsity_pattern_destroy(p)
		call sparse_matrix_destroy(m)
		p%n_row = p1%n_row
		p%n_column = p1%n_column
		p%is_compressed = .true.
		p%n_max_row_length = 0
		n1 = 0
		flag = 0
		do i=1,p%n_row
			n2 = 0
			do j=p1%row_start(i),p1%row_start(i+1)-1
				if (flag(p1%column_index(j)) /= i) then
					n2 = n2+1
					flag(p1%column_index(j)) = i
				end if
			end do
			do j=p2%row_start(i),p2%row_start(i+1)-1
				if (flag(p2%column_index(j)) /= i) then
					n2 = n2+1
					flag(p2%column_index(j)) = i
				end if
			end do
			p%n_max_row_length = max(p%n_max_row_length, n2)
			n1 = n1+n2
		end do
		p%n_nonzero = n1
		p%n_max_vec_length = n1
		allocate(p%row_start(1:p%n_row+1))
		allocate(p%column_index(1:n1))
		allocate(m%matrix_entry(1:n1))
		n2 = 1
		flag = 0
		p%row_start(1) = 1
		row_entry = 0.0e0
		do i=1,p%n_row
			do j=p1%row_start(i),p1%row_start(i+1)-1
				if (flag(p1%column_index(j)) /= i) then
					p%column_index(n2) = p1%column_index(j)
					n2 = n2+1
					flag(p1%column_index(j)) = i
				end if
				row_entry(p1%column_index(j)) = row_entry(p1%column_index(j)) + m1%matrix_entry(j)
			end do
			do j=p2%row_start(i),p2%row_start(i+1)-1
				if (flag(p2%column_index(j)) /= i) then
					p%column_index(n2) = p2%column_index(j)
					n2 = n2+1
					flag(p2%column_index(j)) = i
				end if
				row_entry(p2%column_index(j)) = row_entry(p2%column_index(j)) - m2%matrix_entry(j)
			end do
			p%row_start(i+1) = n2
			do j=p%row_start(i),p%row_start(i+1)-1
				m%matrix_entry(j) = row_entry(p%column_index(j))
				row_entry(p%column_index(j)) = 0.0e0
			end do
		end do

		if (p%n_row == p%n_column) then
			do i=1,p%n_row
				if (p%column_index(p%row_start(i)) == i) then
					cycle
				end if
				do j=p%row_start(i)+1,p%row_start(i+1)-1
					if (p%column_index(j) == i) then
						d = m%matrix_entry(j)
						m%matrix_entry(j) = m%matrix_entry(p%row_start(i))
						m%matrix_entry(p%row_start(i)) = d
						p%column_index(j) = p%column_index(p%row_start(i))
						p%column_index(p%row_start(i)) = i
						exit
					end if
				end do
			end do
		end if
	end subroutine sparse_matrix_minus
	
	subroutine sparsity_pattern_sub_pattern(p, p1, rows, cols)
		type(sparsity_pattern), intent(in)::p
		type(sparsity_pattern), intent(inout)::p1
		integer, dimension(:), intent(in)::rows
		integer, dimension(:), intent(in)::cols
		
		integer::i, j, k
		type(sparsity_pattern)::p2
		
		p2%n_row = size(rows)
		p2%n_column = p%n_column
		p2%is_compressed = .true.
		allocate(p2%row_start(1:p2%n_row+1))
		p2%n_max_row_length = 0
		p2%row_start(1) = 1
		do i=1,p2%n_row
			j = rows(i)
			k = p%row_start(j+1) - p%row_start(j)
			p2%n_max_row_length = max(p2%n_max_row_length, k)
			p2%row_start(i+1) = p2%row_start(i) + k
		end do
		p2%n_nonzero = p2%row_start(p2%n_row+1)-1
		p2%n_max_vec_length = p2%n_nonzero
		allocate(p2%column_index(1:p2%n_nonzero))
		do i=1,p1%n_row
			j = rows(i)
			p2%column_index(p2%row_start(i):p2%row_start(i+1)-1) = p%column_index(p%row_start(j):p%row_start(j+1)-1)
		end do
		call sparsity_pattern_transpose(p2, p1)
		call sparsity_pattern_destroy(p2)

		p2%n_row = size(cols)
		p2%n_column = p1%n_column
		p2%is_compressed = .true.
		allocate(p2%row_start(1:p2%n_row+1))
		p2%n_max_row_length = 0
		p2%row_start(1) = 1
		do i=1,p2%n_row
			j = cols(i)
			k = p1%row_start(j+1) - p1%row_start(j)
			p2%n_max_row_length = max(p2%n_max_row_length, k)
			p2%row_start(i+1) = p2%row_start(i) + k
		end do
		p2%n_nonzero = p2%row_start(p2%n_row+1)-1
		p2%n_max_vec_length = p2%n_nonzero
		allocate(p2%column_index(1:p2%n_nonzero))
		do i=1,p1%n_row
			j = cols(i)
			p2%column_index(p2%row_start(i):p2%row_start(i+1)-1) = p1%column_index(p1%row_start(j):p1%row_start(j+1)-1)
		end do
		call sparsity_pattern_transpose(p2, p1)
		call sparsity_pattern_destroy(p2)
	end subroutine sparsity_pattern_sub_pattern		

	subroutine sparse_matrix_sub_matrix(m, m1, p1, rows, cols)
		type(sparse_matrix), intent(in)::m
		type(sparse_matrix), intent(inout)::m1
		type(sparsity_pattern), intent(inout)::p1
		integer, dimension(:), intent(in)::rows
		integer, dimension(:), intent(in)::cols
		
		integer::i, j, k
		type(sparsity_pattern), pointer::p=>NULL()
		type(sparsity_pattern), target::p2
		type(sparse_matrix)::m2
		
		p => m%pattern
		
		p2%n_row = size(rows)
		p2%n_column = p%n_column
		p2%is_compressed = .true.
		allocate(p2%row_start(1:p2%n_row+1))
		p2%n_max_row_length = 0
		p2%row_start(1) = 1
		do i=1,p2%n_row
			j = rows(i)
			k = p%row_start(j+1) - p%row_start(j)
			p2%n_max_row_length = max(p2%n_max_row_length, k)
			p2%row_start(i+1) = p2%row_start(i) + k
		end do
		p2%n_nonzero = p2%row_start(p2%n_row+1)-1
		p2%n_max_vec_length = p2%n_nonzero
		m2%pattern => p2
		allocate(p2%column_index(1:p2%n_nonzero))
		allocate(m2%matrix_entry(1:p2%n_nonzero))
		do i=1,p1%n_row
			j = rows(i)
			p2%column_index(p2%row_start(i):p2%row_start(i+1)-1) = p%column_index(p%row_start(j):p%row_start(j+1)-1)
			m2%matrix_entry(p2%row_start(i):p2%row_start(i+1)-1) = m%matrix_entry(p%row_start(j):p%row_start(j+1)-1)
		end do
		call sparse_matrix_transpose(m2, m1, p1)
		call sparse_matrix_destroy(m2)
		call sparsity_pattern_destroy(p2)

		p2%n_row = size(cols)
		p2%n_column = p1%n_column
		p2%is_compressed = .true.
		allocate(p2%row_start(1:p2%n_row+1))
		p2%n_max_row_length = 0
		p2%row_start(1) = 1
		do i=1,p2%n_row
			j = cols(i)
			k = p1%row_start(j+1) - p1%row_start(j)
			p2%n_max_row_length = max(p2%n_max_row_length, k)
			p2%row_start(i+1) = p2%row_start(i) + k
		end do
		p2%n_nonzero = p2%row_start(p2%n_row+1)-1
		p2%n_max_vec_length = p2%n_nonzero
		m2%pattern = p2
		allocate(p2%column_index(1:p2%n_nonzero))
		allocate(m2%matrix_entry(1:p2%n_nonzero))
		do i=1,p1%n_row
			j = cols(i)
			p2%column_index(p2%row_start(i):p2%row_start(i+1)-1) = p1%column_index(p1%row_start(j):p1%row_start(j+1)-1)
			m2%matrix_entry(p2%row_start(i):p2%row_start(i+1)-1) = m1%matrix_entry(p1%row_start(j):p1%row_start(j+1)-1)
		end do
		call sparse_matrix_transpose(m2, m1, p1)
		call sparse_matrix_destroy(m2)
		call sparsity_pattern_destroy(p2)
	end subroutine sparse_matrix_sub_matrix		

	subroutine sparsity_pattern_print_gnuplot(p, filename)
		type(sparsity_pattern), intent(in)::p
		character(*), intent(in)::filename
		
		integer::i, j, fid
		
		fid = 10
		open(unit=fid, file=filename//".gnuplot", action="write")
		do i=1,p%n_row
			do j=p%row_start(i), p%row_start(i+1)-1
				write(unit=fid, fmt=*) p%column_index(j), -i
			end do
		end do
		close(unit=fid)
	end subroutine sparsity_pattern_print_gnuplot
	
	subroutine sparse_matrix_print(m, filename)
		type(sparse_matrix), intent(in)::m
		character(*), intent(in)::filename
		
		integer::i, j, fid
		
		fid = 10
		open(unit=fid, file=filename, action="write")
		do i=1,m%pattern%n_row
			do j=m%pattern%row_start(i), m%pattern%row_start(i+1)-1
				write(unit=fid, fmt=*) i, m%pattern%column_index(j), m%matrix_entry(j)
			end do
		end do
		close(unit=fid)
	end subroutine sparse_matrix_print

END MODULE SPARSE_MATRIX_MODULE

!
! end of file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


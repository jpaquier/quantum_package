 subroutine lapack_diag_complex(eigvalues,eigvectors,H,nmax,n)
  implicit none
  integer :: n,nmax
  complex*16, intent(in)        :: H(nmax,nmax)
  complex*16, intent(out)       :: eigvectors(nmax,n)
  double precision, intent(out) :: eigvalues(n)
  complex*16, allocatable   :: A(:,:)
  complex*16, allocatable   :: work(:)
  integer, allocatable      :: iwork(:)
  integer                   :: lwork, info, i,j,l,k, liwork
  double precision          :: rwork(max(1, 3*n-2))
  character*1 :: jobz 
  jobz='V'
  !if jobz = "N" then it computes only the eigenvalues 
  lwork = 2*n*n + 6*n+ 1
  allocate(A(nmax,n),work(lwork),iwork(liwork))
  A=H
  lwork = 2*n*n + 6*n+ 1
  liwork = 5*n + 3
  call zheev(jobz, 'U', n, A, nmax, eigvalues, work, lwork, rwork, info )               
  if (info < 0) then
    print *, irp_here, ': zheev: the ',-info,'-th argument had an illegal value'
    stop 2
  endif
  eigvectors = A
  deallocate (A,work,iwork)
 end



 subroutine ortho_canonical_complex(overlap,LDA,N,C,LDC,m)
  implicit none
  BEGIN_DOC
  ! Compute C_new=C_old.U.s^-1/2 canonical orthogonalization.
  !
  ! overlap : overlap matrix 
  !
  ! LDA : leftmost dimension of overlap array
  !
  ! N : Overlap matrix is NxN (array is (LDA,N) )
  !
  ! C : Coefficients of the vectors to orthogonalize. On exit,
  !     orthogonal vectors
  !
  ! LDC : leftmost dimension of C
  !
  ! m : Coefficients matrix is MxN, ( array is (LDC,N) )
  !
  END_DOC
  
  integer, intent(in)            :: lda, ldc, n
  integer, intent(out)           :: m
  double precision, intent(in)   :: overlap(lda,n)
  double precision, intent(inout) :: C(ldc,n)
  double precision, allocatable  :: U(:,:)
  double precision, allocatable  :: Vt(:,:)
  double precision, allocatable  :: D(:)
  double precision, allocatable  :: S(:,:)
  !DIR$ ATTRIBUTES ALIGN : 64    :: U, Vt, D
  integer                        :: info, i, j
  
  if (n < 2) then
    return
  endif

  allocate (U(ldc,n), Vt(lda,n), D(n), S(lda,n))

  call svd(overlap,lda,U,ldc,D,Vt,lda,n,n)

  D(:) = dsqrt(D(:))
  m=n
  do i=1,n
    if ( D(i) >= 1.d-6 ) then
      D(i) = 1.d0/D(i)
    else
      m = i-1
      print *,  'Removed Linear dependencies below:', 1.d0/D(m)
      exit
    endif
  enddo
  do i=m+1,n
    D(i) = 0.d0
  enddo

  do i=1,m
    if ( D(i) >= 1.d5 ) then
      print *,  'Warning: Basis set may have linear dependence problems'
    endif
  enddo

  do j=1,n
    do i=1,n
      S(i,j) = U(i,j)*D(j)
    enddo
  enddo

  do j=1,n
    do i=1,n
      U(i,j) = C(i,j)
    enddo
  enddo
  
  call dgemm('N','N',n,n,n,1.d0,U,size(U,1),S,size(S,1),0.d0,C,size(C,1))
  deallocate (U, Vt, D, S)
  
 end
 

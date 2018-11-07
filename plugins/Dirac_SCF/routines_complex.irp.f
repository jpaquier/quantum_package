 subroutine lapack_diag_complex(eigvalues,eigvectors,H,nmax,n)
  BEGIN_DOC
  !Diagonalize a complex matrix H
  !Gives its eigenvalues and eigenvectors
  !Dimension of H is nmax x n 
  END_DOC 
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


 subroutine svd_complex(A,LDA,U,LDU,D,Vt,LDVt,rwork,m,n)
  implicit none
  BEGIN_DOC
  !Compute A = U.D.Vt
  !LDx : leftmost dimension of x
  !Dimension of A is m x n
  END_DOC
  integer, intent(in)             :: LDA, LDU, LDVt, m, n
  complex*16, intent(in)          :: A(LDA,n)
  complex*16, intent(out)         :: U(LDU,m)
  complex*16, intent(out)         :: Vt(LDVt,n)
  double precision, intent(out)   :: D(min(m,n))
  complex*16, allocatable         :: work(:)
  integer                         :: info, lwork, i, j, k
  complex*16, allocatable         :: A_tmp(:,:)
  double precision, intent(out)   :: rwork(5*min(m,n))
  allocate (A_tmp(LDA,n))      
  A_tmp = A                    
  ! Find optimal size for temp arrays
  allocate(work(1))
  lwork = -1
  call zgesvd('A','A', m, n, A_tmp, LDA,                             &
      D, U, LDU, Vt, LDVt, work, lwork, rwork, info)
  lwork = int(work(1))
  deallocate(work)
  allocate(work(lwork))
  call zgesvd('A','A', m, n, A_tmp, LDA,                             &
      D, U, LDU, Vt, LDVt, work, lwork, rwork, info)
  deallocate(work,A_tmp)
  if (info /= 0) then
    print *,  info, ': SVD failed'
    stop
  endif
 end  

 subroutine ortho_canonical_complex(dirac_overlap,LDA,N,C,LDC,m)
  implicit none
  BEGIN_DOC
  !Compute C_new=C_old.U.s^-1/2 canonical orthogonalization.
  !overlap : overlap matrix 
  !LDA : leftmost dimension of overlap array
  !N : Overlap matrix is NxN (array is (LDA,N) )
  !C : Coefficients of the vectors to orthogonalize. On exit, orthogonal vectors
  !LDC : leftmost dimension of C
  !m : Coefficients matrix is MxN, ( array is (LDC,N) )
  END_DOC
  integer, intent(in)            :: lda, ldc, n
  integer, intent(out)           :: m
  complex*16, intent(in)         :: dirac_overlap(lda,n)
  complex*16, intent(inout)      :: C(ldc,n)
  complex*16, allocatable        :: U(:,:)
  complex*16, allocatable        :: Vt(:,:)
  double precision, allocatable  :: D(:)
  double precision, allocatable  :: rwork(:)
  complex*16, allocatable        :: S(:,:)
  integer                        :: info, i, j
  if (n < 2) then
    return
  endif
  allocate (U(ldc,n), Vt(lda,n), D(n), rwork(5*n),S(lda,n))
  call svd_complex(dirac_overlap,lda,U,ldc,D,Vt,lda,rwork,n,n)
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
  call zgemm('N','N',n,n,n,(1.d0,0.d0),U,size(U,1),S,size(S,1),(0.d0,0.d0),C,size(C,1))
  deallocate (U, Vt, D, S)
 end

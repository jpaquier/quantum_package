subroutine lapack_diag_complex(eigvalues,eigvectors,H,nmax,n)
 implicit none
 integer :: n,nmax
 double complex, intent(in)  :: H(nmax,nmax)
 double complex, intent(out) :: eigvectors(nmax,n)
 double complex, intent(out) :: eigvalues(n)

 double complex,allocatable   :: A(:,:)
 double complex,allocatable   :: work(:)
 integer         ,allocatable   :: iwork(:)
 integer                        :: lwork, info, i,j,l,k, liwork
 double precision :: rwork(max(1, 3*n-2))
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

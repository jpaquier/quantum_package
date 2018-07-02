program pouet
 implicit none
 double precision :: test_bart,test,test_bart2
 read_wf = .True.
! call routine
! call test_matrx_svd
call test_rho2_bart(test_bart)
call test_rho2_bart2(test_bart2)
call test_rho2_bourrin(test)
print*,' '
print*,' '
print*,'***********Error*******'
print*,'exact  bourin  = ',test
print*,'test dipole    = ',test_bart2
print*,'test couple    =',test_bart
print*,'Erreur couple  = ',dabs(test_bart - test)
print*,'Erreur dipole  = ',dabs(test_bart2 - test)
end


subroutine test_rho2_bourrin(test)
 implicit none
 double precision, intent(out) :: test
 integer :: j,k,l,istate,n_k_loc,m
 double precision :: r(3),rho2
 double precision, allocatable :: rho2_ap(:)
 double precision :: two_dm_in_r
 double precision :: wall_1, wall_2
 allocate(rho2_ap(N_states))

do istate = 1, N_states
 test = 0.d0
 r(1) = 0.d0
 r(2) = 0.d0
 r(3) = 0.d0
!rho2 = two_dm_in_r(r,r,istate)
 call wall_time(wall_1)
  do j = 1, nucl_num
   do k = 1, n_points_radial_grid  -1
    do l = 1, n_points_integration_angular
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     rho2 = two_dm_in_r(r,r,istate)
!    stop
     test += rho2 * final_weight_functions_at_grid_points(l,k,j)
    enddo
   enddo
  enddo
!print*,'test       = ',test
 enddo
 call wall_time(wall_2)
 print*,'wall time bourrin = ',wall_2 - wall_1
end

subroutine test_rho2_bart(test_bart)
 implicit none
 double precision, intent(out) :: test_bart
 integer :: j,k,l,istate,n_k_loc,m
 double precision :: r(3),rho2
 double precision, allocatable :: rho2_ap(:)
 double precision :: two_dm_in_r
 double precision :: wall_1, wall_2
 allocate(rho2_ap(N_states))

 do istate = 1, N_states
  r(1) = 0.d0
  r(2) = 0.d0
  r(3) = 0.d0
  call on_top_pair_density_approx(r,rho2_ap)
  test_bart= 0.d0 
  call wall_time(wall_1)
   do j = 1, nucl_num
    do k = 1, n_points_radial_grid  -1
     do l = 1, n_points_integration_angular
      r(1) = grid_points_per_atom(1,l,k,j)
      r(2) = grid_points_per_atom(2,l,k,j)
      r(3) = grid_points_per_atom(3,l,k,j)
      call on_top_pair_density_approx(r,rho2_ap)
      test_bart += rho2_ap(istate) * final_weight_functions_at_grid_points(l,k,j) 
     enddo
    enddo
   enddo
! print*,'test Barth = ',test_bart
 enddo
 call wall_time(wall_2)
 print*,'wall time couple  = ',wall_2 - wall_1
end

subroutine test_rho2_bart2(test_bart2)
 implicit none
 double precision, intent(out) :: test_bart2
 integer :: j,k,l,istate,n_k_loc,m
 double precision :: r(3),rho2
 double precision, allocatable :: rho2_dip(:)
 double precision :: two_dm_in_r
 double precision :: wall_1, wall_2
 allocate(rho2_dip(N_states))

 do istate = 1, N_states
  r(1) = 0.d0
  r(2) = 0.d0
  r(3) = 0.d0
  call on_top_pair_density_approx_dip(r,rho2_dip)
  test_bart2= 0.d0
  call wall_time(wall_1)
   do j = 1, nucl_num
    do k = 1, n_points_radial_grid  -1
     do l = 1, n_points_integration_angular
      r(1) = grid_points_per_atom(1,l,k,j)
      r(2) = grid_points_per_atom(2,l,k,j)
      r(3) = grid_points_per_atom(3,l,k,j)
      call on_top_pair_density_approx_dip(r,rho2_dip)
      test_bart2 += rho2_dip(istate) * final_weight_functions_at_grid_points(l,k,j)
     enddo
    enddo
   enddo
! print*,'test Barth dip = ',test_bart2
 enddo
 call wall_time(wall_2)
 print*,'wall time dipole = ',wall_2 - wall_1
end

subroutine test_matrx
 implicit none
 double precision, allocatable :: mat(:,:),eigvec_sym(:,:),eigval_sym(:),mat_save(:,:)
 double precision :: tmp
 integer :: n,i,j
 n = 4
 allocate(mat(n,n),mat_save(n,n),eigvec_sym(n,n),eigval_sym(n))
  mat = 0.d0 
  do i = 1, n 
   mat(i,i) = dble(i)
   do j = 1, n
    if(i==j)cycle
   call RANDOM_NUMBER(tmp)
    !tmp = 0.d0
    mat(i,j) = dble(i+j) + tmp * 0.8d0
   enddo
  enddo
  mat_save = mat
  do i = 1, n
   write(*,'(100(F16.10,X))') mat(i,:)
  enddo

  
  call lapack_diagd(eigval_sym,eigvec_sym,mat_save,n,n)
  print*,'n = ',n
  do i = 1, n
   print*,'eigv = ',eigval_sym(i) 
  enddo

!!!!!!!!!!!!!! part that call the non symmetric eigenvalue problem 
  logical, external :: select_dgees
  character*1 :: JOBVS 
  JOBVS = 'V'
  character*1 :: SORT 
  SORT = 'S'

  integer :: lda 
  lda = n 

  integer :: SDIM ! output : number of eigenvalues 

  double precision, allocatable :: wr(:), wi(:) ! output : real and imaginary part of the eigenvalues 
  allocate(wr(n),wi(n))


  integer :: LDVS 
  LDVS = n
  double precision, allocatable :: vs(:,:) ! orthogonal Schur matrix 
  allocate(vs(LDVS,n))
  
  integer :: LWORK
  LWORK = max(1,3*N) 
 
  double precision, allocatable :: WORK(:)
  allocate(WORK(MAX(1,LWORK)))
  
  logical, allocatable ::  BWORK(:)
  allocate(BWORK(n))
  
  integer :: info

  call dgees(JOBVS,SORT,SELECT_dgees, N , mat, lda, SDIM, WR, WI, VS, LDVS, WORK, LWORK, BWORK, INFO)
  print*,'wouwouwouwouwouw'
  print*,'SDIM = ',SDIM
  print*,' '
  print*,'*************Shur eigval******** '
  do i = 1, SDIM
   print*,'wr(i),wi(i) = ',wr(i),wi(i)
  enddo
 
  print*,' '
  print*,'**********Shur Matrix*****************' 
  do i = 1, n
   write(*,'(100(F16.10,X))') mat(i,:)
  enddo

 
 deallocate(mat)
end


logical function select_dgees(a,b)
 implicit none
 double precision, intent(in) :: a,b
 if(dabs(a).ge.0.d-8 .and. dabs(b).ge. 0.-8)then
  select_dgees = .True. 
 else
  select_dgees = .False.
 endif

end

subroutine test_matrx_svd
 implicit none
 double precision, allocatable :: mat(:,:),eigvec_sym(:,:),eigval_sym(:),mat_save(:,:)
 double precision :: tmp
 double precision, allocatable :: vec_tmp(:)
 integer :: n,i,j
 n = 4
 allocate(vec_tmp(n))
 allocate(mat(n,n),mat_save(n,n),eigvec_sym(n,n),eigval_sym(n))
  mat = 0.d0 
  do i = 1, n 
   call RANDOM_NUMBER(tmp)
   vec_tmp(i) = tmp
   mat(i,i) = dble(i)
   do j = 1, n
    if(i==j)cycle
   call RANDOM_NUMBER(tmp)
    mat(i,j) = dble(i+j) + tmp * 1.8d0
   enddo
  enddo
  mat_save = mat
  do i = 1, n
   write(*,'(100(F16.10,X))') mat(i,:)
  enddo
  double precision :: accu
  accu = 0.d0
  do i = 1, n
   do j = 1, n
    accu += vec_tmp(i) * vec_tmp(j) * mat(j,i)
   enddo
  enddo

  call lapack_diagd(eigval_sym,eigvec_sym,mat_save,n,n)
  print*,'n = ',n
  do i = 1, n
   print*,'eigv = ',eigval_sym(i) 
  enddo
  double precision, allocatable   :: U(:,:)
  double precision, allocatable   :: Vt(:,:)
  double precision, allocatable   :: D(:)
  LDU = n
  LDVt = n
  m = n 
  allocate (U(LDU,m),Vt(LDVt,n),D(min(m,n)))
  call svd(mat,n,U,LDU,D,Vt,LDVt,m,n)


  ! U(i,k)  = <k|i> where <k| is the kth left eigenvector
  ! Vt(k,i) = <i|k> where |k> is the kth right eigenvector
  ! D(k)    =  kth eigenvalue 
  integer :: k,LDU,LDVt,m
  double precision :: accu_2
  double precision :: l_k,r_k
  accu_2 = 0.d0
  do k = 1, n
   print*,'D(k) = ',D(k) 
   ! computes the kth Left and Right eigenvector 
   l_k = 0.d0
   r_k = 0.d0
   do i = 1, n
    l_k += U(i,k) * vec_tmp(i)
    r_k += Vt(k,i) * vec_tmp(i) 
   enddo
   accu_2 += l_k * r_k * D(k)
  enddo
  print*,'accu2 = ',accu_2 
  print*,'accu  = ',accu 

end

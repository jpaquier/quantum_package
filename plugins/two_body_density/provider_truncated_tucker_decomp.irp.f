subroutine twobodydump
  implicit none
  character*(128) :: output
  integer :: i_unit_output,getUnitAndOpen
  output=trim(ezfio_filename)//'.twobody_DUMP'
  i_unit_output = getUnitAndOpen(output,'w')

  integer :: i,j,k,l

  double precision :: integral
  PROVIDE two_bod_alpha_beta_mo_transposed

  do l=1,mo_tot_num
   do k=1,mo_tot_num
    do j=1,mo_tot_num
     do i=1,mo_tot_num
      integral = two_bod_alpha_beta_mo_transposed(i,j,k,l,1)
      if (dabs(integral) > 1d-16) then
        write(i_unit_output,*) integral, i,k,j,l
      else
       integral=0d0
       write(i_unit_output,*) integral, i,k,j,l
      endif
     enddo
    enddo
   enddo
  enddo
end

subroutine tucker_decomposition 
 implicit none
 integer :: i,j,k,l,jkl

 double precision :: integral
 double precision, allocatable :: mat_i(:,:),mat_j(:,:),mat_k(:,:),mat_l(:,:)
 allocate(mat_i(mo_tot_num,mo_tot_num**3),mat_j(mo_tot_num,mo_tot_num**3),mat_k(mo_tot_num,mo_tot_num**3),mat_l(mo_tot_num,mo_tot_num**3))
!!!!!!!unfoldage!!!!!!!

 do i = 1,mo_tot_num
  jkl = 0
  do l=1,mo_tot_num
   do k=1,mo_tot_num
    do j=1,mo_tot_num
     jkl= jkl+1
     mat_i(i,jkl)= two_bod_alpha_beta_mo_transposed(i,j,k,l,1)
    enddo
   enddo
  enddo
 enddo

 do j = 1,mo_tot_num
  jkl = 0
  do l=1,mo_tot_num
   do k=1,mo_tot_num
    do i=1,mo_tot_num
     jkl= jkl+1
     mat_j(j,jkl)= two_bod_alpha_beta_mo_transposed(i,j,k,l,1)
    enddo
   enddo
  enddo
 enddo

 do k = 1,mo_tot_num
  jkl = 0
  do l=1,mo_tot_num
   do j=1,mo_tot_num
    do i=1,mo_tot_num
     jkl= jkl+1
     mat_k(k,jkl)= two_bod_alpha_beta_mo_transposed(i,j,k,l,1)
    enddo
   enddo
  enddo
 enddo

 do l = 1,mo_tot_num
  jkl = 0
  do k=1,mo_tot_num
   do j=1,mo_tot_num
    do i=1,mo_tot_num
     jkl= jkl+1
     mat_l(l,jkl)= two_bod_alpha_beta_mo_transposed(i,j,k,l,1)
    enddo
   enddo
  enddo
 enddo

!!!!!!!test SVD!!!!!!!
 
 double precision, allocatable :: u_i(:,:),vt_i(:,:),D_i(:)
 allocate(u_i(mo_tot_num,mo_tot_num),vt_i(mo_tot_num**3,mo_tot_num**3),D_i(mo_tot_num))
 call svd(mat_i,size(mat_i,1),u_i,size(u_i,1),D_i,vt_i,size(vt_i,1),size(mat_i,1),size(mat_i,2))

 double precision, allocatable :: u_j(:,:),vt_j(:,:),D_j(:)
 allocate(u_j(mo_tot_num,mo_tot_num),vt_j(mo_tot_num**3,mo_tot_num**3),D_j(mo_tot_num))
 call svd(mat_j,size(mat_j,1),u_j,size(u_j,1),D_j,vt_j,size(vt_j,1),size(mat_j,1),size(mat_j,2))

 double precision, allocatable :: u_k(:,:),vt_k(:,:),D_k(:)
 allocate(u_k(mo_tot_num,mo_tot_num),vt_k(mo_tot_num**3,mo_tot_num**3),D_k(mo_tot_num))
 call svd(mat_k,size(mat_k,1),u_k,size(u_k,1),D_k,vt_k,size(vt_k,1),size(mat_k,1),size(mat_k,2))

 double precision, allocatable :: u_l(:,:),vt_l(:,:),D_l(:)
 allocate(u_l(mo_tot_num,mo_tot_num),vt_l(mo_tot_num**3,mo_tot_num**3),D_l(mo_tot_num))
 call svd(mat_l,size(mat_l,1),u_l,size(u_l,1),D_l,vt_l,size(vt_l,1),size(mat_l,1),size(mat_l,2))

 print*, 'D_i(k)=',D_i(1),D_i(2),D_i(3),D_i(4),D_i(5),D_i(6)
 print*, 'D_j(k)=',D_j(1),D_j(2),D_j(3),D_j(4),D_j(5),D_j(6)
 print*, 'D_k(k)=',D_k(1),D_k(2),D_k(3),D_k(4),D_k(5),D_k(6)
 print*, 'D_l(k)=',D_l(1),D_l(2),D_l(3),D_l(4),D_l(5),D_l(6)

!!!!!!!test construction g!!!!!!!
!!!!!!!transposer  
 double precision, allocatable :: u_i_t(:,:)
 allocate(u_i_t(mo_tot_num,mo_tot_num))

 double precision, allocatable :: u_j_t(:,:)
 allocate(u_j_t(mo_tot_num,mo_tot_num))

 double precision, allocatable :: u_k_t(:,:)
 allocate(u_k_t(mo_tot_num,mo_tot_num))

 double precision, allocatable :: u_l_t(:,:)
 allocate(u_l_t(mo_tot_num,mo_tot_num))

 do i = 1,mo_tot_num
  do j = 1,mo_tot_num
   u_i_t(i,j)=u_i(j,i)
   u_j_t(i,j)=u_j(j,i)
   u_k_t(i,j)=u_k(j,i)
   u_l_t(i,j)=u_l(j,i)
  enddo
 enddo


 double precision, allocatable :: g(:,:,:,:)
 allocate(g(mo_tot_num,mo_tot_num,mo_tot_num,mo_tot_num))
integer :: r1,r2,r3,r4

 g = 0d0

 do i = 1,mo_tot_num
  do j = 1,mo_tot_num
   do k = 1,mo_tot_num
    do l =1,mo_tot_num
     do r1 = 1,mo_tot_num
      do r2 = 1,mo_tot_num
       do r3 = 1,mo_tot_num
        do r4 = 1,mo_tot_num
         g(i,j,k,l) += two_bod_alpha_beta_mo_transposed(r1,r2,r3,r4,1) * u_i_t(i,r1) * u_j_t(j,r2) * u_k_t(k,r3) * u_l_t(l,r4)
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo

!!!!!!!!!!tucker truc!!!!!

 double precision, allocatable :: n2_bis(:,:,:,:)
 allocate(n2_bis(mo_tot_num,mo_tot_num,mo_tot_num,mo_tot_num))

 n2_bis= 0d0

 do i = 1,mo_tot_num
  do j = 1,mo_tot_num
   do k = 1,mo_tot_num
    do l =1,mo_tot_num
     do r1 =1,mo_tot_num
      do r2 = 1,mo_tot_num
       do r3 = 1,mo_tot_num
        do r4 = 1,mo_tot_num
         n2_bis(i,j,k,l) += g(r1,r2,r3,r4) * u_i(i,r1) * u_j(j,r2) * u_k(k,r3) * u_l(l,r4)
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo


 double precision :: accu
 accu= 0d0
 do i = 1,mo_tot_num
  do j = 1,mo_tot_num
   do k = 1,mo_tot_num
    do l =1,mo_tot_num
      accu +=  dabs(n2_bis(i,j,k,l) - two_bod_alpha_beta_mo_transposed(i,j,k,l,1))
    enddo
   enddo
  enddo
 enddo


 print*, '**************'
 print*, 'accu=', accu
 print*, '**************'
end


subroutine truncated_tucker_decomposition 
 implicit none
 integer :: i,j,k,l,jkl

 double precision :: integral
 double precision, allocatable :: mat_i(:,:),mat_j(:,:),mat_k(:,:),mat_l(:,:)
 allocate(mat_i(mo_tot_num,mo_tot_num**3),mat_j(mo_tot_num,mo_tot_num**3),mat_k(mo_tot_num,mo_tot_num**3),mat_l(mo_tot_num,mo_tot_num**3))
!!!!!!!unfoldage!!!!!!!

 do i = 1,mo_tot_num
  jkl = 0
  do l=1,mo_tot_num
   do k=1,mo_tot_num
    do j=1,mo_tot_num
     jkl= jkl+1
     mat_i(i,jkl)= two_bod_alpha_beta_mo_transposed(i,j,k,l,1)
    enddo
   enddo
  enddo
 enddo

 do j = 1,mo_tot_num
  jkl = 0
  do l=1,mo_tot_num
   do k=1,mo_tot_num
    do i=1,mo_tot_num
     jkl= jkl+1
     mat_j(j,jkl)= two_bod_alpha_beta_mo_transposed(i,j,k,l,1)
    enddo
   enddo
  enddo
 enddo

!do k = 1,mo_tot_num
! jkl = 0
! do l=1,mo_tot_num
!  do j=1,mo_tot_num
!   do i=1,mo_tot_num
!    jkl= jkl+1
!    mat_k(k,jkl)= two_bod_alpha_beta_mo_transposed(i,j,k,l,1)
!   enddo
!  enddo
! enddo
!enddo

!do l = 1,mo_tot_num
! jkl = 0
! do k=1,mo_tot_num
!  do j=1,mo_tot_num
!   do i=1,mo_tot_num
!    jkl= jkl+1
!    mat_l(l,jkl)= two_bod_alpha_beta_mo_transposed(i,j,k,l,1)
!   enddo
!  enddo
! enddo
!enddo

!!!!!!!test SVD!!!!!!!
 
 double precision, allocatable :: u_i(:,:),vt_i(:,:),D_i(:)
 allocate(u_i(mo_tot_num,mo_tot_num),vt_i(mo_tot_num**3,mo_tot_num**3),D_i(mo_tot_num))
 call svd(mat_i,size(mat_i,1),u_i,size(u_i,1),D_i,vt_i,size(vt_i,1),size(mat_i,1),size(mat_i,2))

 double precision, allocatable :: u_j(:,:),vt_j(:,:),D_j(:)
 allocate(u_j(mo_tot_num,mo_tot_num),vt_j(mo_tot_num**3,mo_tot_num**3),D_j(mo_tot_num))
 call svd(mat_j,size(mat_j,1),u_j,size(u_j,1),D_j,vt_j,size(vt_j,1),size(mat_j,1),size(mat_j,2))

!double precision, allocatable :: u_k(:,:),vt_k(:,:),D_k(:)
!allocate(u_k(mo_tot_num,mo_tot_num),vt_k(mo_tot_num**3,mo_tot_num**3),D_k(mo_tot_num))
!call svd(mat_k,size(mat_k,1),u_k,size(u_k,1),D_k,vt_k,size(vt_k,1),size(mat_k,1),size(mat_k,2))

!double precision, allocatable :: u_l(:,:),vt_l(:,:),D_l(:)
!allocate(u_l(mo_tot_num,mo_tot_num),vt_l(mo_tot_num**3,mo_tot_num**3),D_l(mo_tot_num))
!call svd(mat_l,size(mat_l,1),u_l,size(u_l,1),D_l,vt_l,size(vt_l,1),size(mat_l,1),size(mat_l,2))

!!!!!!!Selection valeur propre!!!!!!!
 double precision :: threshinou
 integer :: n_eigen_i,n_eigen_j
 threshinou= 1d-5
 print*,'************************'
 print*,'threshinou    =',threshinou
 print*,'************************'

 n_eigen_i = 1 
 do while ( (dabs(D_i(n_eigen_i)) .gt. threshinou) .AND. (n_eigen_i .lt. mo_tot_num)  )
  n_eigen_i += 1
  print*,n_eigen_i,D_i(n_eigen_i)
 enddo
 print*,'************************'
 print*,'n_eigen_i    =',n_eigen_i
 print*,'************************'


 n_eigen_j = 1 
 do while ((dabs(D_j(n_eigen_j)) .gt. threshinou) .AND. (n_eigen_j .lt. mo_tot_num))
  n_eigen_j += 1
  print*,n_eigen_j,D_j(n_eigen_j)
 enddo
 print*,'************************'
 print*,'n_eigen_j    =',n_eigen_j
 print*,'************************'
!!!!!!!Creation _truncated_U_i _truncated_U_j!!!!!!!

 double precision, allocatable :: truncated_u_i(:,:)
 allocate(truncated_u_i(mo_tot_num,n_eigen_i))

 double precision, allocatable :: truncated_u_j(:,:)
 allocate(truncated_u_j(mo_tot_num,n_eigen_j))

 do i=1,mo_tot_num
  do j=1,n_eigen_i
   truncated_u_i(i,j) = u_i(i,j)
  enddo
 enddo

 do i=1,mo_tot_num
  do j=1,n_eigen_j
   truncated_u_j(i,j) = u_j(i,j)
  enddo
 enddo

!!!!!!!test construction g!!!!!!!
!!!!!!!transposer  
 double precision, allocatable :: truncated_u_i_t(:,:)
 allocate(truncated_u_i_t(n_eigen_i,mo_tot_num))

 double precision, allocatable :: truncated_u_j_t(:,:)
 allocate(truncated_u_j_t(n_eigen_j,mo_tot_num))

!double precision, allocatable :: u_k_t(:,:)
!allocate(u_k_t(mo_tot_num,mo_tot_num))

!double precision, allocatable :: u_l_t(:,:)
!allocate(u_l_t(mo_tot_num,mo_tot_num))

 do i = 1,n_eigen_i
  do j = 1,mo_tot_num
   truncated_u_i_t(i,j)=u_i(j,i)
  enddo
 enddo

 do i = 1,n_eigen_j
  do j = 1,mo_tot_num
   truncated_u_j_t(i,j)=u_j(j,i)
  enddo
 enddo


 double precision, allocatable :: g(:,:,:,:)
 allocate(g(n_eigen_i,n_eigen_j,n_eigen_i,n_eigen_j))
 integer :: r1,r2,r3,r4

 g = 0d0

 do i = 1,n_eigen_i
  do j = 1,n_eigen_j
   do k = 1,n_eigen_i
    do l =1,n_eigen_j
     do r1 = 1,mo_tot_num
      do r2 = 1,mo_tot_num
       do r3 = 1,mo_tot_num
        do r4 = 1,mo_tot_num
         g(i,j,k,l) += two_bod_alpha_beta_mo_transposed(r1,r2,r3,r4,1) * truncated_u_i_t(i,r1) * truncated_u_j_t(j,r2) * truncated_u_i_t(k,r3) * truncated_u_j_t(l,r4)
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo

!!!!!!!!!!tucker truc!!!!!

!double precision, allocatable :: n2_bis(:,:,:,:)
!allocate(n2_bis(mo_tot_num,mo_tot_num,mo_tot_num,mo_tot_num))

!n2_bis= 0d0

!do i = 1,mo_tot_num
! do j = 1,mo_tot_num
!  do k = 1,mo_tot_num
!   do l =1,mo_tot_num
!    do r1 =1,n_eigen_i
!     do r2 = 1,n_eigen_j
!      do r3 = 1,n_eigen_i
!       do r4 = 1,n_eigen_j
!        n2_bis(i,j,k,l) += g(r1,r2,r3,r4) * truncated_u_i(i,r1) * truncated_u_j(j,r2) * truncated_u_i(k,r3) * truncated_u_j(l,r4)
!       enddo
!      enddo
!     enddo
!    enddo
!   enddo
!  enddo
! enddo
!enddo


!double precision :: accu
!accu= 0d0
!do i = 1,mo_tot_num
! do j = 1,mo_tot_num
!  do k = 1,mo_tot_num
!   do l =1,mo_tot_num
!     accu +=  dabs(n2_bis(i,j,k,l) - two_bod_alpha_beta_mo_transposed(i,j,k,l,1))
!   enddo
!  enddo
! enddo
!enddo


!print*, '**************'
!print*, 'accu=', accu
!print*, '**************'

!test sur l'energie!!!!!
 double precision :: get_mo_bielec_integral_ijkl_r3
 double precision, allocatable :: integrals_ij(:,:)
 allocate(integrals_ij(mo_tot_num,mo_tot_num))


 double precision :: n2_bis_ijkl,E_cor_tot_normal,E_cor_tot_tucker

 n2_bis_ijkl= 0d0
 E_cor_tot_normal = 0d0
 E_cor_tot_tucker = 0d0


 do i = 1,mo_tot_num
  do j = 1,mo_tot_num
   call get_mo_bielec_integrals_ijkl_r3_ij(i,j,mo_tot_num,integrals_ij,mo_integrals_ijkl_r3_map)
   do k = 1,mo_tot_num
    do l =1,mo_tot_num
     n2_bis_ijkl = 0d0

     do r1 =1,n_eigen_i
      do r2 = 1,n_eigen_j
       do r3 = 1,n_eigen_i
        do r4 = 1,n_eigen_j
         n2_bis_ijkl += g(r1,r2,r3,r4) * truncated_u_i(i,r1) * truncated_u_j(j,r2) * truncated_u_i(k,r3) * truncated_u_j(l,r4)
        enddo
       enddo
      enddo
     enddo

     E_cor_tot_normal += two_bod_alpha_beta_mo_transposed(k,l,j,i,1) * integrals_ij(l,k)
     E_cor_tot_tucker += n2_bis_ijkl * integrals_ij(l,k)
    enddo
   enddo
  enddo
 enddo


 double precision :: accu
 accu= E_cor_tot_normal-E_cor_tot_tucker 


 print*, '**************'
 print*, 'accu=', accu
 print*, '**************'

end

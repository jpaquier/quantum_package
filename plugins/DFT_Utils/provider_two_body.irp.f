 BEGIN_PROVIDER [double precision,trace_abs_max,(N_states)]
&BEGIN_PROVIDER [double precision,trace_abs_max_couple,(mo_tot_num,mo_tot_num,N_states)]
 implicit none
 integer :: i,j,k,l,istate
 double precision, allocatable :: mat(:,:),D(:),U(:,:),Vt(:,:) 
 double precision :: trace_1,trace_mat_abs
 integer :: t,LDU,LDVt,m
 allocate(mat(mo_tot_num,mo_tot_num),U(mo_tot_num,mo_tot_num),Vt(mo_tot_num,mo_tot_num),D(mo_tot_num))
 LDU = mo_tot_num 
 LDVt = mo_tot_num
 trace_abs_max=0
 trace_abs_max_couple=0.d0
 do istate = 1, N_states
  do i = 1, mo_tot_num ! loop over the first electron 
   do j = 1, mo_tot_num ! loop over the second electron 
    ! get matrix 
    do k = 1, mo_tot_num
     do l = 1, mo_tot_num
      mat(l,k) = two_bod_alpha_beta_mo_transposed(l,k,j,i,istate)
     enddo
    enddo
    call svd(mat,mo_tot_num,U,LDU,D,Vt,LDVt,mo_tot_num,mo_tot_num)
    ! D(k)    =  kth eigenvalue 
    do t = 1, mo_tot_num
      trace_abs_max_couple(j,i,istate) += dabs(D(t))  
     ! print*,'D(t)    = ',D(t)
    enddo
    trace_abs_max(istate)+= trace_abs_max_couple(j,i,istate)
   enddo
  enddo
 enddo
 deallocate(mat,U,Vt,D)
 END_PROVIDER


 BEGIN_PROVIDER [integer,n_couple,(N_states)]
&BEGIN_PROVIDER [integer,n_couple_max]
&BEGIN_PROVIDER [integer,list_couple,(mo_tot_num**2,2,N_states)]
&BEGIN_PROVIDER [integer,n_k_dens,(mo_tot_num**2,N_states)]
&BEGIN_PROVIDER [integer,list_k_dens,(mo_tot_num,mo_tot_num**2,N_states)]
 implicit none
 integer :: i,j,k,l,istate,n_total
 double precision, allocatable :: mat(:,:),D(:),U(:,:),Vt(:,:)
 double precision :: trace_1,trace_mat_abs
 integer :: t,LDU,LDVt,m
 allocate(mat(mo_tot_num,mo_tot_num),U(mo_tot_num,mo_tot_num),Vt(mo_tot_num,mo_tot_num),D(mo_tot_num))
 LDU = mo_tot_num
 LDVt = mo_tot_num
 n_couple_max=0 
 list_couple=0
 print*,'*****'
 print*,'thr_eig_2dm    = ',thr_eig_2dm
 print*,'thr_couple_2dm = ',thr_couple_2dm
 print*,'N orbital ',mo_tot_num
 do istate = 1, N_states 
  n_couple(istate)=0
  do i = 1, mo_tot_num ! loop over the first electron 
   do j = 1, mo_tot_num ! loop over the second electron 
    trace_mat_abs=0.d0
    ! get matrix 
    do k = 1, mo_tot_num
     do l = 1, mo_tot_num
      mat(l,k) = two_bod_alpha_beta_mo_transposed(l,k,j,i,istate)
     enddo
    enddo
    call svd(mat,mo_tot_num,U,LDU,D,Vt,LDVt,mo_tot_num,mo_tot_num)
    do t = 1, mo_tot_num
     trace_mat_abs += dabs(D(t))
    enddo 
    ! If trace > thr_couple_2dm: diagonalize the matrix 
    if(dabs((trace_mat_abs)/trace_abs_max(istate)) .ge. thr_couple_2dm)then
     n_couple(istate) += 1
     list_couple(n_couple,1,istate)=i
     list_couple(n_couple,2,istate)=j
     n_k_dens(n_couple,istate)=0   
     ! Selection of the eigenvalue > thr_eig_2dm
     do k = 1, mo_tot_num
     if(dabs(D(k)/trace_abs_max_couple(j,i,istate)) .ge. thr_eig_2dm)then
       n_k_dens(n_couple,istate)+=1
       list_k_dens(n_k_dens(n_couple,istate),n_couple,istate) = k
     endif
     enddo
    endif
   enddo
  enddo
  if(n_couple(istate) .gt. n_couple_max)then
   n_couple_max= n_couple(istate)
  endif 
  print*,'***********COUPLE**********'
  print*,'INFO ON THE TWO BODY DM PRE DIAGONALIZATION '
  print*,'max(N_couple) = ',n_couple_max
  print*,'N couple etat = ',n_couple(istate)
  print*,'N couple MAX  = ',mo_tot_num**2
  n_total = 0
  do i = 1, n_couple(istate)
   n_total += n_k_dens(i,istate)
  enddo
  print*,'N total       = ',n_total
  print*,'N operation   = ',n_total * mo_tot_num
  print*,'N operation max=',mo_tot_num **4
  print*,'*************'
 enddo
 deallocate(mat,U,Vt,D)
 END_PROVIDER


 BEGIN_PROVIDER [double precision,lambda_k,(mo_tot_num,n_couple_max,N_states)]
&BEGIN_PROVIDER [double precision,psi_k_couple_l,(mo_tot_num,mo_tot_num,n_couple_max,N_states)]
&BEGIN_PROVIDER [double precision,psi_k_couple_r,(mo_tot_num,mo_tot_num,n_couple_max,N_states)]
 implicit none
 integer :: i,j,k,l,m,t,istate,LDU,LDVt
 double precision, allocatable :: mat(:,:),D(:),U(:,:),Vt(:,:)
 allocate(mat(mo_tot_num,mo_tot_num),U(mo_tot_num,mo_tot_num),Vt(mo_tot_num,mo_tot_num),D(mo_tot_num))
 LDU = mo_tot_num
 LDVt = mo_tot_num
 do istate = 1, N_states
  do m = 1, n_couple(istate) ! loop over the first electron 
   ! get matrix 
   i=list_couple(m,1,istate)
   j=list_couple(m,2,istate)
   do k = 1, mo_tot_num
    do l = 1, mo_tot_num
     mat(l,k) = two_bod_alpha_beta_mo_transposed(l,k,j,i,istate)
    enddo
   enddo
   call svd(mat,mo_tot_num,U,LDU,D,Vt,LDVt,mo_tot_num,mo_tot_num)
   ! U(i,k)  = <k|i> where <k| is the kth left eigenvector
   ! Vt(k,i) = <i|k> where |k> is the kth right eigenvector
   ! D(k)    =  kth eigenvalue 
   ! SAve eigenvalues & right and left eigenvetors
   do k = 1, n_k_dens(m,istate)
    lambda_k(k,m,istate)= D(list_k_dens(k,m,istate))
    do t=1,mo_tot_num
     psi_k_couple_l(t,k,m,istate)=U(t,list_k_dens(k,m,istate))
     psi_k_couple_r(t,k,m,istate)=Vt(list_k_dens(k,m,istate),t)
    enddo
   enddo
  enddo
 enddo
 deallocate(mat,U,Vt,D)
 END_PROVIDER


subroutine on_top_pair_density_approx(r,rho2_ap)
 implicit none
 double precision, intent(in)  :: r(3)
 double precision, intent(out) :: rho2_ap(N_states)
 double precision :: psi_temp_l,psi_temp_r,tmp
 double precision :: mos_array(mo_tot_num)
 integer :: i,j,k,l,t,m,istate
 call give_all_mos_at_r(r,mos_array)
 psi_temp_l = 0.d0
 psi_temp_r = 0.d0
 do istate = 1, N_states
  rho2_ap(istate) = 0.d0
  do m = 1, n_couple(istate) ! loop over the first electron 
   i=list_couple(m,1,istate)
   j=list_couple(m,2,istate)
   tmp = 0.d0
   do k = 1, n_k_dens(m,istate) 
    psi_temp_l=0.d0
    psi_temp_r=0.d0
    do t=1,mo_tot_num
     psi_temp_r+=psi_k_couple_r(t,k,m,istate)*mos_array(t)
     psi_temp_l+=psi_k_couple_l(t,k,m,istate)*mos_array(t)
    enddo
    tmp += lambda_k(k,m,istate)*psi_temp_l*psi_temp_r
   enddo
   rho2_ap(istate) += mos_array(i)*mos_array(j)*tmp
  enddo
 enddo
end

 BEGIN_PROVIDER [double precision,trace_abs_max_dip,(N_states)]
&BEGIN_PROVIDER [double precision,trace_abs_max_couple_dip,(mo_tot_num,mo_tot_num,N_states)]
 implicit none
 integer :: i,j,k,l,istate
 double precision, allocatable :: mat(:,:),eigvec(:,:),eigval(:)
 double precision :: trace_1,trace_mat_abs
 allocate(mat(mo_tot_num,mo_tot_num),eigvec(mo_tot_num,mo_tot_num),eigval(mo_tot_num))
 trace_abs_max_dip=0.d0
 trace_abs_max_couple_dip=0.d0
 do istate = 1, N_states
  do i = 1, mo_tot_num ! loop over the first electron 
   do j = 1, mo_tot_num ! loop over the second electron 
    ! get matrix 
    do k = 1, mo_tot_num
     do l = 1, mo_tot_num
      mat(l,k) = two_bod_alpha_beta_mo(l,k,j,i,istate)
     enddo
    enddo
    call lapack_diagd(eigval,eigvec,mat,mo_tot_num,mo_tot_num)
    do k = 1, mo_tot_num
     trace_abs_max_couple_dip(j,i,istate) += dabs(eigval(k))
    enddo
    trace_abs_max_dip(istate)+= trace_abs_max_couple_dip(j,i,istate)
   enddo
  enddo
 enddo
 deallocate(mat,eigvec,eigval)
 END_PROVIDER

 BEGIN_PROVIDER [integer,n_couple_dip,(N_states)]
&BEGIN_PROVIDER [integer,n_couple_max_dip]
&BEGIN_PROVIDER [integer,list_couple_dip,(mo_tot_num**2,2,N_states)]
&BEGIN_PROVIDER [integer,n_k_dens_dip,(mo_tot_num**2,N_states)]
&BEGIN_PROVIDER [integer,list_k_dens_dip,(mo_tot_num,mo_tot_num**2,N_states)]
 implicit none
 integer :: i,j,k,l,istate,n_total
 double precision, allocatable :: mat(:,:),eigvec(:,:),eigval(:)
 double precision :: trace_1,trace_mat_abs
 allocate(mat(mo_tot_num,mo_tot_num),eigvec(mo_tot_num,mo_tot_num),eigval(mo_tot_num))
 n_couple_max_dip=0
 list_couple_dip=0
! print*,'*****'
! print*,'thr_eig_2dm    = ',thr_eig_2dm
! print*,'thr_couple_2dm = ',thr_couple_2dm
! print*,'N orbital ',mo_tot_num
 do istate = 1, N_states
  n_couple_dip(istate)=0
  do i = 1, mo_tot_num ! loop over the first electron 
   do j = 1, mo_tot_num ! loop over the second electron 
    trace_mat_abs=0.d0
    ! get matrix 
    do k = 1, mo_tot_num
     do l = 1, mo_tot_num
      mat(l,k) = two_bod_alpha_beta_mo(l,k,j,i,istate)
     enddo
    enddo
    call lapack_diagd(eigval,eigvec,mat,mo_tot_num,mo_tot_num)
    do k = 1, mo_tot_num
     trace_mat_abs += dabs(eigval(k)) 
    enddo
    ! If trace > thr_couple_2dm: diagonalize the matrix 
    if(dabs((trace_mat_abs)/trace_abs_max_dip(istate)) .gt. thr_couple_2dm)then
     n_couple_dip(istate) += 1
     list_couple_dip(n_couple_dip,1,istate)=i
     list_couple_dip(n_couple_dip,2,istate)=j
     n_k_dens_dip(n_couple_dip,istate)=0
     ! Selection of the eigenvalue > thr_eig_2dm
     do k = 1, mo_tot_num
     !print*,'k          =  ',k
     !print*,'eig/trace  =  ', dabs(eigval(k)/trace_abs_max_couple(j,i,istate))
     !print*,'Thr        =  ',thr_eig_2dm 
     if(dabs(eigval(k)/trace_abs_max_couple_dip(j,i,istate)) .gt. thr_eig_2dm)then
       n_k_dens_dip(n_couple_dip,istate)+=1
       list_k_dens_dip(n_k_dens_dip(n_couple_dip,istate),n_couple_dip,istate) = k
     endif
     enddo
    endif
   enddo
  enddo
  if(n_couple_dip(istate) .gt. n_couple_max_dip)then
   n_couple_max_dip= n_couple_dip(istate)
  endif
  print*,'************DIPOLE************'
  print*,'INFO ON THE TWO BODY DM PRE DIAGONALIZATION '
  print*,'max(N_couple) = ',n_couple_max_dip
  print*,'N couple etat = ',n_couple_dip(istate)
  print*,'N couple MAX  = ',mo_tot_num**2
  n_total = 0
  do i = 1, n_couple_dip(istate)
   n_total += n_k_dens_dip(i,istate)
  enddo
  print*,'N total       = ',n_total
  print*,'N operation   = ',n_total * mo_tot_num
  print*,'N operation max=',mo_tot_num **4
  print*,'*************'
 enddo
 deallocate(mat,eigvec,eigval)
 END_PROVIDER


 BEGIN_PROVIDER [double precision,lambda_k_dip,(mo_tot_num,n_couple_max_dip,N_states)]
&BEGIN_PROVIDER [double precision,psi_k_couple_dip,(mo_tot_num,mo_tot_num,n_couple_max_dip,N_states)]
 implicit none
 integer :: i,j,k,l,m,t,istate
 double precision, allocatable :: mat(:,:),eigvec(:,:),eigval(:)
 allocate(mat(mo_tot_num,mo_tot_num),eigvec(mo_tot_num,mo_tot_num),eigval(mo_tot_num))
 do istate = 1, N_states
  do m = 1, n_couple_dip(istate) ! loop over the first electron 
   ! get matrix 
   i=list_couple_dip(m,1,istate)
   j=list_couple_dip(m,2,istate)
   do k = 1, mo_tot_num
    do l = 1, mo_tot_num
     mat(l,k) = two_bod_alpha_beta_mo(l,k,j,i,istate)
    enddo
   enddo
   call lapack_diagd(eigval,eigvec,mat,mo_tot_num,mo_tot_num)
   ! SAve eigenvalues & eigenvetors
   do k = 1, n_k_dens_dip(m,istate)
    lambda_k_dip(k,m,istate)= eigval(list_k_dens_dip(k,m,istate))
    do t=1,mo_tot_num
     psi_k_couple_dip(t,k,m,istate)=eigvec(t,list_k_dens_dip(k,m,istate))
    enddo
   enddo
  enddo
 enddo
 deallocate(mat,eigvec,eigval)
 END_PROVIDER

subroutine on_top_pair_density_approx_dip(r,rho2_dip)
 implicit none
 double precision, intent(in)  :: r(3)
 double precision, intent(out) :: rho2_dip(N_states)
 double precision :: psi_temp,tmp
 double precision :: mos_array(mo_tot_num)
 integer :: i,j,k,l,t,m,istate
 call give_all_mos_at_r(r,mos_array)
 psi_temp = 0.d0
 do istate = 1, N_states
  rho2_dip(istate) = 0.d0
  do m = 1, n_couple_dip(istate) ! loop over the first electron 
   i=list_couple_dip(m,1,istate)
   j=list_couple_dip(m,2,istate)
   tmp = 0.d0
   do k = 1, n_k_dens_dip(m,istate)
    psi_temp=0.d0
    do t=1,mo_tot_num
     psi_temp+=psi_k_couple_dip(t,k,m,istate)*mos_array(t)
    enddo
    tmp += lambda_k_dip(k,m,istate)*psi_temp*psi_temp
   enddo
   rho2_dip(istate) += mos_array(i)*mos_array(j)*tmp
  enddo
 enddo
end

 BEGIN_PROVIDER [double precision, thr_eig_2dm]
&BEGIN_PROVIDER [double precision, thr_couple_2dm]
 implicit none
 thr_couple_2dm = thr_ontop_approx * 1e-1 
 thr_eig_2dm = thr_ontop_approx
 END_PROVIDER



 BEGIN_PROVIDER [integer, couple_to_array, (mo_tot_num,mo_tot_num)]
&BEGIN_PROVIDER [integer, couple_to_array_reverse, (mo_tot_num*mo_tot_num,2)]
 implicit none
 integer :: i,j,compt
 compt = 0
 do i = 1, mo_tot_num ! loop over the first electron 
  do j = 1, mo_tot_num ! loop over the second electron 
   ! get matrix 
   compt += 1
   couple_to_array(j,i) = compt
   couple_to_array_reverse(compt,1) = i 
   couple_to_array_reverse(compt,2) = j 
  enddo
 enddo

 END_PROVIDER 


 BEGIN_PROVIDER [double precision,E_cor_tot,(N_states)]
&BEGIN_PROVIDER [double precision,E_cor_couple_sorted,(mo_tot_num*mo_tot_num,N_states)]
&BEGIN_PROVIDER [integer,E_cor_couple_sorted_order,(mo_tot_num*mo_tot_num,N_states)]
 implicit none
 integer :: i,j,k,l,istate,t,icouple
 integer, allocatable :: order_loc(:)
 double precision :: get_mo_bielec_integral_ijkl_r3
 double precision, allocatable :: E_cor_coupl(:,:)
 allocate(E_cor_coupl(mo_tot_num*mo_tot_num,2),order_loc(mo_tot_num*mo_tot_num))
 E_cor_tot = 0.d0
 do istate = 1, N_states
   do i = 1, mo_tot_num ! loop over the first electron 
    do j = 1, mo_tot_num ! loop over the second electron 
     ! get matrix 
     icouple = couple_to_array(j,i)
     order_loc(icouple)= icouple
     E_cor_coupl(icouple,:) = 0.d0
     do k = 1, mo_tot_num
      do l = 1, mo_tot_num
       E_cor_coupl(icouple,1) -= dabs(two_bod_alpha_beta_mo_transposed(l,k,j,i,istate) * get_mo_bielec_integral_ijkl_r3(l,k,j,i,mo_integrals_ijkl_r3_map))
       E_cor_coupl(icouple,2) +=     (two_bod_alpha_beta_mo_transposed(l,k,j,i,istate) * get_mo_bielec_integral_ijkl_r3(l,k,j,i,mo_integrals_ijkl_r3_map))
      enddo
     enddo
!    print*,'couple '
!    print*,couple_to_array_reverse(icouple,1),couple_to_array_reverse(icouple,2)
!    print*,E_cor_coupl(icouple,1)
     E_cor_tot(istate) += E_cor_coupl(icouple,2)
    enddo
   enddo
   call dsort(E_cor_coupl(1,1),order_loc,mo_tot_num*mo_tot_num)
   do t = 1, mo_tot_num*mo_tot_num
    E_cor_couple_sorted(t,istate)= E_cor_coupl(order_loc(t),2)
    E_cor_couple_sorted_order(t,istate)= order_loc(t)
   enddo
  enddo
 deallocate(E_cor_coupl,order_loc)
 END_PROVIDER



 BEGIN_PROVIDER [integer,n_couple_ec,(N_states)]
&BEGIN_PROVIDER [integer,n_couple_max_ec]
&BEGIN_PROVIDER [double precision ,E_cor_couple, (N_states)]
 implicit none
 integer :: m,istate
 double precision :: E_cor_loc
 n_couple_max_ec = 0 
 do istate = 1, N_states
  n_couple_ec(istate) = 0
  m = 1
  n_couple_ec(istate) += 1
  E_cor_loc = E_cor_couple_sorted(m,istate)
  do while (dabs(E_cor_loc - E_cor_tot(istate))/dabs(E_cor_tot(istate)) .gt. thr_couple_2dm  )
   m += 1
   n_couple_ec(istate) += 1
   E_cor_loc += E_cor_couple_sorted(m,istate)
  enddo
  E_cor_couple(istate) = E_cor_loc
 enddo
 n_couple_max_ec = maxval(n_couple_ec)
 print*,'n_ couple max    = ',n_couple_ec(1)
END_PROVIDER 


 BEGIN_PROVIDER [integer,identity_eig,(2,mo_tot_num*n_couple_max_ec,N_states)]
&BEGIN_PROVIDER [integer,identity_eig_reverse,(mo_tot_num,n_couple_max_ec,N_states)]
 implicit none
 integer :: m,t,istate, i_eigen
 i_eigen = 0
 do istate = 1, N_states
  do m = 1, n_couple_ec(istate)
   do t = 1, mo_tot_num
    i_eigen += 1
    identity_eig(1,i_eigen,istate) = m
    identity_eig(2,i_eigen,istate) = t
    identity_eig_reverse(t,m,istate) = i_eigen
   enddo
  enddo
 enddo

END_PROVIDER





 BEGIN_PROVIDER [integer,n_k_tot,(N_states)]
&BEGIN_PROVIDER [integer,n_k_tot_max]
&BEGIN_PROVIDER [double precision,ec_eigen_sorted,(mo_tot_num*n_couple_max_ec,N_states)]
&BEGIN_PROVIDER [integer,ec_eigen_sorted_order,(mo_tot_num*n_couple_max_ec,N_states)]
 implicit none
 integer :: i_eigen,s,l,k,i,j,LDU,LDVt,t,q,r,m,istate
 double precision :: E_cor_loc,k_ec_tmp,get_mo_bielec_integral_ijkl_r3
 double precision, allocatable :: mat(:,:),D(:),U(:,:),Vt(:,:),ec_eigen(:,:)
 integer, allocatable :: order_loc(:)
 allocate(order_loc(mo_tot_num*n_couple_max_ec),ec_eigen(mo_tot_num*n_couple_max_ec,2),mat(mo_tot_num,mo_tot_num),U(mo_tot_num,mo_tot_num),Vt(mo_tot_num,mo_tot_num),D(mo_tot_num))
 LDU = mo_tot_num
 LDVt = mo_tot_num
 do istate = 1, N_states
  i_eigen=0
  E_cor_loc = 0.d0
  do m = 1, n_couple_ec(istate)
   i = couple_to_array_reverse(E_cor_couple_sorted_order(m,istate),1)
   j = couple_to_array_reverse(E_cor_couple_sorted_order(m,istate),2)
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
   do t = 1, mo_tot_num  
    i_eigen += 1
    order_loc(i_eigen) = i_eigen
    k_ec_tmp = 0.d0
   do q =1, mo_tot_num
    do r =1, mo_tot_num
      k_ec_tmp += U(q,t)*Vt(t,r)*get_mo_bielec_integral_ijkl_r3(r,q,j,i,mo_integrals_ijkl_r3_map)
    enddo
   enddo
  ec_eigen(i_eigen,1) = -abs(D(t) * k_ec_tmp)
  ec_eigen(i_eigen,2) = D(t) * k_ec_tmp 
  E_cor_loc += ec_eigen(i_eigen,2)
   enddo
  enddo
  call dsort(ec_eigen(1,1),order_loc,mo_tot_num*n_couple_max_ec)
  E_cor_loc = 0.d0
  n_k_tot(istate) = 1
  ec_eigen_sorted(n_k_tot(istate),istate)= ec_eigen(order_loc(n_k_tot(istate)),2)
  E_cor_loc += ec_eigen_sorted(n_k_tot(istate),istate)
  ec_eigen_sorted_order(n_k_tot(istate),istate)= order_loc(n_k_tot(istate)) 
  do while (dabs(E_cor_loc - E_cor_couple(istate))/dabs(E_cor_couple(istate)) .gt. thr_eig_2dm )
   n_k_tot(istate) += 1
   ec_eigen_sorted(n_k_tot(istate),istate)= ec_eigen(order_loc(n_k_tot(istate)),2)
   ec_eigen_sorted_order(n_k_tot(istate),istate)= order_loc(n_k_tot(istate))
   E_cor_loc += ec_eigen_sorted(n_k_tot(istate),istate) 
  enddo
  double precision :: per_k,per_c
  print*,'n couples    =',n_couple_ec(istate)
  print*,'n couples max=',mo_tot_num**2
  print*,'n k          =',n_k_tot(istate)
  print*,'n k max      =',mo_tot_num**3
  per_c = (dble(n_couple_ec(istate))/(dble(mo_tot_num)**2))*100
  per_k = (dble(n_k_tot(istate))/(dble(mo_tot_num)**3))*100
  print*,'*******************************'
  print*,'Percentage of selected couples = ',per_c
  print*,'Percentage of selected of k    = ',per_k
  print*,'*******************************'
 enddo
 n_k_tot_max = maxval(n_k_tot)
 END_PROVIDER
 
 BEGIN_PROVIDER [integer,n_k_by_couple,(mo_tot_num,mo_tot_num,N_states)]
 implicit none
 integer :: istate,i_eigen,m,i,j 
 n_k_by_couple=0
 do istate = 1, N_states
  do i_eigen= 1,n_k_tot(istate) 
   m = identity_eig(1,ec_eigen_sorted_order(i_eigen,istate),istate)
   i = couple_to_array_reverse(E_cor_couple_sorted_order(m,istate),1)
   j = couple_to_array_reverse(E_cor_couple_sorted_order(m,istate),2)
  
   n_k_by_couple(j,i,istate) += 1 
  enddo
 enddo
 END_PROVIDER


 BEGIN_PROVIDER [double precision,lambda_k_ec,(2,n_k_tot_max,N_states)]
&BEGIN_PROVIDER [double precision,psi_k_couple_l_ec,(mo_tot_num,n_k_tot_max,N_states)]
&BEGIN_PROVIDER [double precision,psi_k_couple_r_ec,(mo_tot_num,n_k_tot_max,N_states)]
 implicit none
 integer :: i_eigen,i_eigen_temp,s,l,k,i,j,LDU,LDVt,t,q,r,m,istate
 integer, allocatable :: order_loc(:)
 double precision, allocatable :: mat(:,:),D(:),U(:,:),Vt(:,:),ec_eigen(:,:)
 double precision :: E_cor_loc,k_ec_tmp,get_mo_bielec_integral_ijkl_r3
 allocate(order_loc(mo_tot_num),ec_eigen(mo_tot_num,2),mat(mo_tot_num,mo_tot_num),U(mo_tot_num,mo_tot_num),Vt(mo_tot_num,mo_tot_num),D(mo_tot_num))
 LDU = mo_tot_num
 LDVt = mo_tot_num
 do istate = 1, N_states
  i_eigen=0
  do m = 1, n_couple_ec(istate)
   i = couple_to_array_reverse(E_cor_couple_sorted_order(m,istate),1)
   j = couple_to_array_reverse(E_cor_couple_sorted_order(m,istate),2)
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
   i_eigen_temp = 0
   order_loc = 0
   do t = 1, mo_tot_num
   i_eigen_temp += 1
   order_loc(i_eigen_temp) = i_eigen_temp
   k_ec_tmp = 0.d0
   do q =1, mo_tot_num
    do r =1, mo_tot_num
      k_ec_tmp += U(q,t)*Vt(t,r)*get_mo_bielec_integral_ijkl_r3(r,q,j,i,mo_integrals_ijkl_r3_map)
    enddo
   enddo
   ec_eigen(i_eigen_temp,1) = -abs(D(t) * k_ec_tmp)
   ec_eigen(i_eigen_temp,2) = D(t) * k_ec_tmp
   enddo
   call dsort(ec_eigen(1,1),order_loc,mo_tot_num)
   !blablablabla
   do k = 1, n_k_by_couple(j,i,istate)
    i_eigen += 1
    lambda_k_ec(1,i_eigen,istate)= D(order_loc(k))
    lambda_k_ec(2,i_eigen,istate)= m
    do t=1,mo_tot_num
     psi_k_couple_l_ec(t,i_eigen,istate)=U(t,order_loc(k))
     psi_k_couple_r_ec(t,i_eigen,istate)=Vt(order_loc(k),t)
    enddo
   enddo
  enddo
 enddo
END_PROVIDER


subroutine on_top_pair_density_thresh_ec(rho2_ana)
 implicit none
 double precision, intent(out) :: rho2_ana(N_states)
 double precision :: tmp,tmp2,get_mo_bielec_integral_ijkl_r3
 integer :: i,j,k,t,m,s,istate
 double precision :: wall_1, wall_2
 do istate = 1, N_states
   rho2_ana(istate) = 0.d0
   tmp = 0.d0
   call wall_time(wall_1)
   do k = 1, n_k_tot(istate)
    m=lambda_k_ec(2,k,istate)
    i = couple_to_array_reverse(E_cor_couple_sorted_order(m,istate),1)
    j = couple_to_array_reverse(E_cor_couple_sorted_order(m,istate),2)
    tmp2 = 0.d0
    do t=1,mo_tot_num
     do s=1,mo_tot_num
     tmp2 +=psi_k_couple_l_ec(t,k,istate)*psi_k_couple_r_ec(s,k,istate)*get_mo_bielec_integral_ijkl_r3(s,t,j,i,mo_integrals_ijkl_r3_map)
     enddo
    enddo
    tmp += lambda_k_ec(1,k,istate)*tmp2
   enddo
   rho2_ana(istate) += tmp
 enddo
 call wall_time(wall_2)
 print*,'wall time selected k analitycal integrals = ',wall_2 - wall_1
end

double precision function two_dm_in_r_k_selected(r1,r2,istate)
 implicit none
 integer, intent(in) :: istate
 double precision, intent(in) :: r1(3),r2(3)
 double precision, allocatable :: mos_array_r1(:), mos_array_r2(:)
 integer :: i,j,k,l,m,t
 double precision :: tmp,psi_temp_l,psi_temp_r
 allocate(mos_array_r2(mo_tot_num), mos_array_r1(mo_tot_num))
 call give_all_mos_at_r(r1,mos_array_r1)
 call give_all_mos_at_r(r2,mos_array_r2)
 two_dm_in_r_k_selected = 0.d0
  tmp = 0.d0
  do k = 1, n_k_tot(istate)
   m=lambda_k_ec(2,k,istate)
   i = couple_to_array_reverse(E_cor_couple_sorted_order(m,istate),1)
   j = couple_to_array_reverse(E_cor_couple_sorted_order(m,istate),2)
   psi_temp_l=0.d0
   psi_temp_r=0.d0
   do t=1,mo_tot_num
    psi_temp_r+=psi_k_couple_r_ec(t,k,istate)*mos_array_r1(t)
    psi_temp_l+=psi_k_couple_l_ec(t,k,istate)*mos_array_r2(t)
   enddo
   two_dm_in_r_k_selected += lambda_k_ec(1,k,istate)*psi_temp_l*psi_temp_r*mos_array_r1(i)*mos_array_r1(j)
  enddo
end


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

 BEGIN_PROVIDER [double precision, thr_eig_2dm]
&BEGIN_PROVIDER [double precision, thr_couple_2dm]
 implicit none
 BEGIN_DOC
 ! Determine from the general 'thr_ontop_approx' threshold the two thresholds used to approximated the on top. The first one (thr_couple_2dm) impact the number of selected electron couple and the second one (thr_eig_2dm) impact the total number of svd eigenvalue taken into account to compute the on top
 END_DOC
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
 BEGIN_DOC
 ! Compute the correlation energy by electron couple and rank the results by descending order in E_cor_couple_sorted
 END_DOC
 integer :: i,j,k,l,istate,t,icouple
 integer, allocatable :: order_loc(:)
 double precision :: get_mo_bielec_integral_ijkl_r3
 double precision, allocatable :: E_cor_coupl(:,:)
 allocate(E_cor_coupl(mo_tot_num*mo_tot_num,2),order_loc(mo_tot_num*mo_tot_num))

 double precision, allocatable :: integrals_ij(:,:)
 allocate(integrals_ij(mo_tot_num,mo_tot_num))
 double precision :: cpu0,cpu1
 E_cor_tot = 0.d0
 cpu0 = dabs(two_bod_alpha_beta_mo_transposed(1,1,1,1,1) * get_mo_bielec_integral_ijkl_r3(1,1,1,1,mo_integrals_ijkl_r3_map))
 call cpu_time(cpu0)
 do istate = 1, N_states
   do i = 1, mo_tot_num ! loop over the first electron 
    do j = 1, mo_tot_num ! loop over the second electron 
     ! get matrix 
     icouple = couple_to_array(j,i)
     order_loc(icouple)= icouple
     E_cor_coupl(icouple,:) = 0.d0
     call get_mo_bielec_integrals_ijkl_r3_ij(i,j,mo_tot_num,integrals_ij,mo_integrals_ijkl_r3_map)
     do k = 1, mo_tot_num
      do l = 1, mo_tot_num
       E_cor_coupl(icouple,1) -= dabs(two_bod_alpha_beta_mo_transposed(l,k,j,i,istate) * integrals_ij(l,k))
       E_cor_coupl(icouple,2) +=     (two_bod_alpha_beta_mo_transposed(l,k,j,i,istate) * integrals_ij(l,k))
      enddo
     enddo
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
 deallocate(integrals_ij)
 call cpu_time(cpu1)
 print*,'Time to provide E_cor_tot   = ',cpu1-cpu0
 END_PROVIDER



 BEGIN_PROVIDER [integer,n_couple_ec,(N_states)]
&BEGIN_PROVIDER [integer,n_couple_max_ec]
&BEGIN_PROVIDER [double precision ,E_cor_couple, (N_states)]
 implicit none
 BEGIN_DOC
 ! Compute the number of elecron couple by state needed to reach the 'thr_couple_2dm' threshold 
 END_DOC
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
 !print*,'n_ couple max    = ',n_couple_ec(1)
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
&BEGIN_PROVIDER [double precision,lambda_k_ec,(2,mo_tot_num*n_couple_max_ec,N_states)]
&BEGIN_PROVIDER [double precision,psi_k_couple_l_ec,(mo_tot_num,mo_tot_num*n_couple_max_ec,N_states)]
&BEGIN_PROVIDER [double precision,psi_k_couple_r_ec,(mo_tot_num,mo_tot_num*n_couple_max_ec,N_states)]
 implicit none
 BEGIN_DOC
 ! SVD computation of the matrices extracted from the n selected couple. And ranking by decreasing order of all the solution in ec_eigen_sorted
 END_DOC 
 integer :: i_eigen,s,l,k,i,j,LDU,LDVt,t,q,r,m,istate,w
 double precision :: E_cor_loc,k_ec_tmp
 double precision, allocatable :: mat(:,:),D(:),U(:,:),Vt(:,:),ec_eigen(:,:)
 integer, allocatable :: order_loc(:)
 allocate(order_loc(mo_tot_num*n_couple_max_ec),ec_eigen(mo_tot_num*n_couple_max_ec,2),mat(mo_tot_num,mo_tot_num),U(mo_tot_num,mo_tot_num),Vt(mo_tot_num,mo_tot_num),D(mo_tot_num))
 double precision, allocatable :: integrals_ij(:,:)
 double precision :: svd1,svd2,svdtot
 double precision :: average1,average2,averagetot
 double precision :: sort1,sort2
 double precision :: cpu1,cpu2
 integer :: lwork_opt
 allocate(integrals_ij(mo_tot_num,mo_tot_num))
 LDU = mo_tot_num
 LDVt = mo_tot_num
 svdtot = 0.d0
 
 call cpu_time(cpu1)
 do istate = 1, N_states
  i_eigen=0
  E_cor_loc = 0.d0
  m = 1
   i = couple_to_array_reverse(E_cor_couple_sorted_order(m,istate),1)
   j = couple_to_array_reverse(E_cor_couple_sorted_order(m,istate),2)
   do k = 1, mo_tot_num
    do l = 1, mo_tot_num
     mat(l,k) = two_bod_alpha_beta_mo_transposed(l,k,j,i,istate)
    enddo
   enddo  
   call find_optimal_lwork_svd(mat,mo_tot_num,U,LDU,D,Vt,LDVt,mo_tot_num,mo_tot_num,lwork_opt)
   print*,"lwork_opt     =",lwork_opt
  do m = 1, n_couple_ec(istate)
   i = couple_to_array_reverse(E_cor_couple_sorted_order(m,istate),1)
   j = couple_to_array_reverse(E_cor_couple_sorted_order(m,istate),2)
   do k = 1, mo_tot_num
    do l = 1, mo_tot_num
     mat(l,k) = two_bod_alpha_beta_mo_transposed(l,k,j,i,istate)
    enddo
   enddo  
   call cpu_time(svd1)
   call svd_lwork_in(mat,mo_tot_num,U,LDU,D,Vt,LDVt,mo_tot_num,mo_tot_num,lwork_opt)
   call cpu_time(svd2)
   svdtot += dabs(svd2-svd1)
   ! U(i,k)  = <k|i> where <k| is the kth left eigenvector
   ! Vt(k,i) = <i|k> where |k> is the kth right eigenvector
   ! D(k)    =  kth eigenvalue 
   ! SAve eigenvalues & right and left eigenvetors
   call cpu_time(average1)
   call get_mo_bielec_integrals_ijkl_r3_ij(i,j,mo_tot_num,integrals_ij,mo_integrals_ijkl_r3_map) 
   
   do t = 1, mo_tot_num  
    i_eigen += 1
    order_loc(i_eigen) = i_eigen
    k_ec_tmp = 0.d0
    do q =1, mo_tot_num
     do r =1, mo_tot_num
        k_ec_tmp += U(q,t)*Vt(t,r)*integrals_ij(r,q)
     enddo
    enddo
    ec_eigen(i_eigen,1) = -abs(D(t) * k_ec_tmp)
    ec_eigen(i_eigen,2) = D(t) * k_ec_tmp
    lambda_k_ec(1,i_eigen,istate)= D(t)
    lambda_k_ec(2,i_eigen,istate)= m
    do w=1,mo_tot_num
     psi_k_couple_l_ec(w,i_eigen,istate)=U(w,t)
     psi_k_couple_r_ec(w,i_eigen,istate)=Vt(t,w)
    enddo
   enddo
   call cpu_time(average2)
   averagetot += dabs(average2-average1)
  enddo
  call cpu_time(sort1)
  call dsort(ec_eigen(1,1),order_loc,mo_tot_num*n_couple_max_ec)
  call cpu_time(sort2)
  E_cor_loc = 0.d0
  n_k_tot(istate) = 1
  ec_eigen_sorted(n_k_tot(istate),istate)= ec_eigen(order_loc(n_k_tot(istate)),2)
  E_cor_loc += ec_eigen_sorted(n_k_tot(istate),istate)
  ec_eigen_sorted_order(n_k_tot(istate),istate)= order_loc(n_k_tot(istate)) 
  do while (dabs(E_cor_loc - E_cor_couple(istate))/dabs(E_cor_couple(istate)) .ge. max(thr_eig_2dm,1.d-14) )
   n_k_tot(istate) += 1
   ec_eigen_sorted(n_k_tot(istate),istate)= ec_eigen(order_loc(n_k_tot(istate)),2)
   ec_eigen_sorted_order(n_k_tot(istate),istate)= order_loc(n_k_tot(istate))
   E_cor_loc += ec_eigen_sorted(n_k_tot(istate),istate) 
  enddo

!!!!!Print selection!!!!
  double precision :: per_k,per_c
  print*,'n couples selected =',n_couple_ec(istate)
  print*,'n couples max      =',mo_tot_num**2
  print*,'n k selected       =',n_k_tot(istate)
  print*,'n k max            =',mo_tot_num**3
  per_c = (dble(n_couple_ec(istate))/(dble(mo_tot_num)**2))*100
  per_k = (dble(n_k_tot(istate))/(dble(mo_tot_num)**3))*100
  print*,'*******************************'
  print*,'Percentage of selected couples = ',per_c
  print*,'Percentage of selected of k    = ',per_k
  print*,'*******************************'
  print*,'time for all SVD                       = ',svdtot
  print*,'time for all averages                  = ',averagetot
  print*,'time for sorting                       = ',sort2-sort1
 call cpu_time(cpu2)
  print*,'total time to provide all eigenvectors = ',cpu2-cpu1
  svdtot = svdtot/dble(n_couple_ec(istate))
  averagetot = averagetot/dble(n_couple_ec(istate))
  print*,'Average time per SVD                   = ',svdtot
  print*,'Average time per average               = ',averagetot
 enddo
 n_k_tot_max = maxval(n_k_tot)
deallocate(order_loc,mat,D,U,Vt,ec_eigen,integrals_ij)
 END_PROVIDER
 


subroutine on_top_pair_density_thresh_ec(rho2_ana)
 implicit none
 BEGIN_DOC
 !Compute the on top with analytical bielectronic integrals
 END_DOC
 double precision, intent(out) :: rho2_ana(N_states)
 double precision :: tmp,tmp2
 integer :: i,j,k,t,m,s,istate
 double precision :: wall_1, wall_2
 double precision, allocatable :: integrals_ij(:,:)
 allocate(integrals_ij(mo_tot_num,mo_tot_num))
 do istate = 1, N_states
   rho2_ana(istate) = 0.d0
   tmp = 0.d0
   call wall_time(wall_1)
   do k = 1, n_k_tot(istate)
    m=lambda_k_ec(2,ec_eigen_sorted_order(k,istate),istate)
    i = couple_to_array_reverse(E_cor_couple_sorted_order(m,istate),1)
    j = couple_to_array_reverse(E_cor_couple_sorted_order(m,istate),2)
    tmp2 = 0.d0
    call get_mo_bielec_integrals_ijkl_r3_ij(i,j,mo_tot_num,integrals_ij,mo_integrals_ijkl_r3_map)
    do t=1,mo_tot_num
     do s=1,mo_tot_num
     tmp2 +=psi_k_couple_l_ec(t,ec_eigen_sorted_order(k,istate),istate)*psi_k_couple_r_ec(s,ec_eigen_sorted_order(k,istate),istate)*integrals_ij(s,t)
     enddo
    enddo
    tmp += lambda_k_ec(1,ec_eigen_sorted_order(k,istate),istate)*tmp2
   enddo
   rho2_ana(istate) += tmp
 enddo
 call wall_time(wall_2)
 print*,'wall time selected k analitycal integrals = ',wall_2 - wall_1
deallocate(integrals_ij)
end

double precision function two_dm_in_r_k_selected(r1,r2,istate)
 implicit none
 BEGIN_DOC
 !Compute the on top locally  
 END_DOC
 integer, intent(in) :: istate
 double precision, intent(in) :: r1(3),r2(3)
 double precision, allocatable :: mos_array_r1(:), mos_array_r2(:)
 integer :: i,j,k,l,m,t
 double precision :: tmp,psi_temp_l,psi_temp_r
 double precision :: u_dot_v
 allocate(mos_array_r2(mo_tot_num), mos_array_r1(mo_tot_num))
 call give_all_mos_at_r(r1,mos_array_r1)
 call give_all_mos_at_r(r2,mos_array_r2)
 two_dm_in_r_k_selected = 0.d0
  tmp = 0.d0
  do k = 1, n_k_tot(istate)
   m = lambda_k_ec(2,ec_eigen_sorted_order(k,istate),istate)
   i = couple_to_array_reverse(E_cor_couple_sorted_order(m,istate),1)
   j = couple_to_array_reverse(E_cor_couple_sorted_order(m,istate),2)
   psi_temp_l=0.d0
   psi_temp_r=0.d0
   psi_temp_r =  u_dot_v(psi_k_couple_r_ec(1,ec_eigen_sorted_order(k,istate),istate),mos_array_r1,mo_tot_num)
   psi_temp_l =  u_dot_v(psi_k_couple_l_ec(1,ec_eigen_sorted_order(k,istate),istate),mos_array_r2,mo_tot_num)
   two_dm_in_r_k_selected += lambda_k_ec(1,ec_eigen_sorted_order(k,istate),istate)*psi_temp_l*psi_temp_r*mos_array_r1(i)*mos_array_r1(j)
  enddo
deallocate(mos_array_r2,mos_array_r1)
end


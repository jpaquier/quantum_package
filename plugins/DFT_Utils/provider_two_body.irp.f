
 BEGIN_PROVIDER [integer,thr_couple]
&BEGIN_PROVIDER [integer,thr_eig]
 thr_eig = 1.d-18
 thr_couple = 1.d-18

 END_PROVIDER


 BEGIN_PROVIDER [integer,n_couple,(N_states)]
&BEGIN_PROVIDER [integer,n_couple_max]
&BEGIN_PROVIDER [integer,list_couple,(mo_tot_num**2,2,N_states)]
&BEGIN_PROVIDER [integer,n_k_dens,(mo_tot_num**2,N_states)]
 implicit none
 integer :: i,j,k,l,istate
 double precision, allocatable :: mat(:,:),eigvec(:,:),eigval(:)
 double precision :: trace_1,trace_mat_abs
 allocate(mat(mo_tot_num,mo_tot_num),eigvec(mo_tot_num,mo_tot_num),eigval(mo_tot_num))
 n_couple_max=0 
 list_couple=0
 do istate = 1, N_states 
 n_couple(istate)=0
 do i = 1, mo_tot_num ! loop over the first electron 
  do j = 1, mo_tot_num ! loop over the second electron 
   ! get matrix 
   do k = 1, mo_tot_num
    do l = 1, mo_tot_num
     mat(l,k) = two_bod_alpha_beta_mo_transposed(l,k,j,i,istate)
    enddo
    trace_mat_abs += dabs(mat(k,k))
   enddo
   ! If trace > thr_couple: diagonalize the matrix 
   if(dabs(trace_mat_abs) .gt. thr_couple)then
     n_couple(istate) += 1
     list_couple(n_couple,1,istate)=i
     list_couple(n_couple,2,istate)=j
     call lapack_diagd(eigval,eigvec,mat,mo_tot_num,mo_tot_num)
     n_k_dens(n_couple,istate)=0   
     ! Selection of the eigenvalue > thr_eig
       do k = 1, mo_tot_num
          if(dabs(eigval(k)) .gt. thr_eig)then
          n_k_dens(n_couple,istate)+=1
          endif
       enddo
   endif
  enddo
 enddo
 if(n_couple(istate) .gt. n_couple_max)then
 n_couple_max= n_couple(istate)
 endif
 enddo
 deallocate(mat,eigvec,eigval)
 END_PROVIDER


 BEGIN_PROVIDER [double precision,lambda_k,(mo_tot_num,mo_tot_num**2,N_states)]
&BEGIN_PROVIDER [double precision,psi_k_couple,(mo_tot_num,mo_tot_num,n_couple_max,N_states)]
 implicit none
 integer :: i,j,k,l,m,n_k_loc,t,istate
 double precision, allocatable :: mat(:,:),eigvec(:,:),eigval(:)
 allocate(mat(mo_tot_num,mo_tot_num),eigvec(mo_tot_num,mo_tot_num),eigval(mo_tot_num))
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
   call lapack_diagd(eigval,eigvec,mat,mo_tot_num,mo_tot_num)
   n_k_loc=0
   ! SAve eigenvalues & eigenvetors
   do k = 1, mo_tot_num
    if(n_k_loc .lt. n_k_dens(m,istate) .AND. dabs(eigval(k)) .gt. thr_eig)then
     n_k_loc+=1
     lambda_k(n_k_loc,m,istate)= eigval(k)
     do t=1,mo_tot_num
      psi_k_couple(t,n_k_loc,m,istate)=eigvec(t,k)
     enddo
    endif
   enddo
 enddo
 enddo
 deallocate(mat,eigvec,eigval)
 END_PROVIDER


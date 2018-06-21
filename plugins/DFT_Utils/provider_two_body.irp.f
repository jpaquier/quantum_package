
 BEGIN_PROVIDER [integer,thr_couple]
&BEGIN_PROVIDER [integer,thr_eig]
 thr_eig = 1.d-22
 thr_couple = 1.d-22

 END_PROVIDER


 BEGIN_PROVIDER [integer,n_couple,(N_states)]
&BEGIN_PROVIDER [integer,n_couple_max]
&BEGIN_PROVIDER [integer,list_couple,(mo_tot_num**2,2,N_states)]
&BEGIN_PROVIDER [integer,n_k_dens,(mo_tot_num**2,N_states)]
&BEGIN_PROVIDER [integer,list_k_dens,(mo_tot_num,mo_tot_num**2,N_states)]
 implicit none
 integer :: i,j,k,l,istate
 double precision, allocatable :: mat(:,:),eigvec(:,:),eigval(:)
 double precision :: trace_1,trace_mat_abs
 allocate(mat(mo_tot_num,mo_tot_num),eigvec(mo_tot_num,mo_tot_num),eigval(mo_tot_num))
 n_couple_max=0 
 list_couple=0
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
          list_k_dens(n_k_dens(n_couple,istate),n_couple,istate) = k
         endif
       enddo
    !print*,'N vecteur propre pour couple ',n_k_dens(n_couple,istate)
    endif
   enddo
  enddo
  if(n_couple(istate) .gt. n_couple_max)then
   n_couple_max= n_couple(istate)
  endif
!! print*,'N couple max ',n_couple_max
!! print*,'N couple etat0  ',n_couple(istate)
 enddo
 deallocate(mat,eigvec,eigval)
 END_PROVIDER


 BEGIN_PROVIDER [double precision,lambda_k,(mo_tot_num,n_couple_max,N_states)]
&BEGIN_PROVIDER [double precision,psi_k_couple,(mo_tot_num,mo_tot_num,n_couple_max,N_states)]
 implicit none
 integer :: i,j,k,l,m,t,istate
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
   ! SAve eigenvalues & eigenvetors
   do k = 1, n_k_dens(m,istate)
    lambda_k(k,m,istate)= eigval(list_k_dens(k,m,istate))
    do t=1,mo_tot_num
     psi_k_couple(t,k,m,istate)=eigvec(t,list_k_dens(k,m,istate))
    enddo
   enddo
  enddo
 enddo
 deallocate(mat,eigvec,eigval)
 END_PROVIDER


subroutine on_top_pair_density_approx(r,rho2_ap)
 implicit none
 double precision, intent(in)  :: r(3)
 double precision, intent(out) :: rho2_ap(N_states)
 double precision :: psi_temp
 double precision :: mos_array(mo_tot_num)
 integer :: i,j,k,l,t,m,istate
 call give_all_mos_at_r(r,mos_array)
 psi_temp = 0.d0
 do istate = 1, N_states
  rho2_ap(istate) = 0.d0
  do m = 1, n_couple(istate) ! loop over the first electron 
   i=list_couple(m,1,istate)
   j=list_couple(m,2,istate)
    do k = 1, n_k_dens(m,istate) 
     psi_temp=0.d0
     do t=1,mo_tot_num
      psi_temp+=psi_k_couple(t,k,m,istate)*mos_array(t)
     enddo
     rho2_ap(istate) += mos_array(i)*mos_array(j)*lambda_k(k,m,istate)*psi_temp*psi_temp
   enddo
  enddo
 enddo
end


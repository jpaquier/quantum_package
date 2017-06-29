BEGIN_PROVIDER [double precision, dressing_ref_hamiltonian, (n_det_ref,n_det_ref,N_states)]
 implicit none
 integer :: i,j,k,l
 integer :: ii,jj,istate
 double precision :: hij,sec_order,H_ref(N_det_ref),hik,hkl
 integer          :: idx(0:N_det_ref)
 double precision :: accu_negative,accu_positive,phase
 integer :: degree_exc_ionic,degree_exc_neutral,exc(0:2,2,2)
 dressing_ref_hamiltonian = 0.d0
 accu_negative = 0.d0
 accu_positive = 0.d0
 integer :: h1,p1,h2,p2,s1,s2
 do istate = 1, N_states
   do i = 1, N_det_non_ref
    call filter_connected_i_H_psi0(psi_ref,psi_non_ref(1,1,i),N_int,N_det_ref,idx)
    H_ref = 0.d0
    do ii=1,idx(0)
      k = idx(ii)
      !DEC$ FORCEINLINE
      call i_H_j(psi_ref(1,1,k),psi_non_ref(1,1,i),N_int,hij)
      H_ref(k)  = hij
    enddo
    do ii= 1, idx(0)
     k = idx(ii)
     hik = H_ref(k) * lambda_special(istate,i)
     do jj = 1, idx(0)
      l = idx(jj)
      dressing_ref_hamiltonian(k,l,istate) += hik * H_ref(l)
     enddo
    enddo
   enddo
 enddo
END_PROVIDER 

BEGIN_PROVIDER [double precision, hamiltonian_total_dressed, (n_det_ref,n_det_ref,N_states)]
 implicit none
 integer :: i,j,k
 do k = 1, N_states
  do i = 1, N_det_ref
   do j = 1, N_det_ref
    hamiltonian_total_dressed(j,i,k) = dressing_ref_hamiltonian(j,i,k) + ref_hamiltonian_matrix(j,i)
   enddo
  enddo
 enddo

END_PROVIDER 



 BEGIN_PROVIDER [ double precision, lambda_special, (N_states, N_det_non_ref) ]
&BEGIN_PROVIDER [ integer, lambda_special_pt2, (0:psi_det_size) ]
&BEGIN_PROVIDER [ integer, lambda_special_kept, (0:psi_det_size) ]
  implicit none
  BEGIN_DOC
  ! cm/<Psi_0|H|D_m> or perturbative 1/Delta_E(m)
  END_DOC
  integer :: i,k
  double precision               :: ihpsi_current(N_states)
  integer                        :: i_pert_count
  double precision               :: hii, lambda_pert
  integer                        :: N_lambda_special_pt2, N_lambda_special_pt3
  
  i_pert_count = 0
  lambda_special = 0.d0
  N_lambda_special_pt2 = 0
  N_lambda_special_pt3 = 0
  lambda_special_pt2(0) = 0
  lambda_special_kept(0) = 0

  do i=1,N_det_non_ref
    call i_h_psi(psi_non_ref(1,1,i), psi_ref, psi_ref_coef, N_int, N_det_ref,&
        size(psi_ref_coef,1), N_states,ihpsi_current)
    call i_H_j(psi_non_ref(1,1,i),psi_non_ref(1,1,i),N_int,hii)
    do k=1,N_states
      if (ihpsi_current(k) == 0.d0) then
        ihpsi_current(k) = 1.d-32
      endif
!      lambda_special(k,i) = psi_non_ref_coef(i,k)/ihpsi_current(k) 
      lambda_special(k,i) = min(-1.d-32,psi_non_ref_coef(i,k)/ihpsi_current(k) )
      lambda_pert = 1.d0 / (psi_ref_energy_diagonalized(k)-hii)
      if (lambda_pert / lambda_special(k,i)  < 0.5d0) then
        ! Ignore lamdba
        i_pert_count += 1
        lambda_special(k,i) = 0.d0
        if (lambda_special_pt2(N_lambda_special_pt2) /= i) then
          N_lambda_special_pt2 += 1
          lambda_special_pt2(N_lambda_special_pt2) = i
        endif
      else
        ! Keep lamdba
        if (lambda_special_kept(N_lambda_special_pt3) /= i) then
          N_lambda_special_pt3 += 1
          lambda_special_kept(N_lambda_special_pt3) = i
        endif
      endif
    enddo
  enddo
  lambda_special_pt2(0) = N_lambda_special_pt2
  lambda_special_kept(0) = N_lambda_special_pt3
  print*,'N_det_non_ref = ',N_det_non_ref
  print*,'psi_coef_ref_ratio = ',psi_ref_coef(2,1)/psi_ref_coef(1,1)
  print*,'lambda max = ',maxval(dabs(lambda_special))
  print*,'Number of ignored determinants = ',i_pert_count  

END_PROVIDER

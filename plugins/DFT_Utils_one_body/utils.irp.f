
subroutine GGA_type_functionals(r,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, & 
                                ex,vx_rho_a,vx_rho_b,vx_grad_rho_a_2,vx_grad_rho_b_2,vx_grad_rho_a_b, &  
                                ec,vc_rho_a,vc_rho_b,vc_grad_rho_a_2,vc_grad_rho_b_2,vc_grad_rho_a_b )
 implicit none
 double precision, intent(in)  :: r(3),rho_a(N_states),rho_b(N_states),grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
 double precision, intent(out) :: ex(N_states),vx_rho_a(N_states),vx_rho_b(N_states),vx_grad_rho_a_2(N_states),vx_grad_rho_b_2(N_states),vx_grad_rho_a_b(N_states)
 double precision, intent(out) :: ec(N_states),vc_rho_a(N_states),vc_rho_b(N_states),vc_grad_rho_a_2(N_states),vc_grad_rho_b_2(N_states),vc_grad_rho_a_b(N_states)
 integer          :: istate
 double precision :: r2(3),dr2(3), local_potential,r12,dx2,mu
 do istate = 1, N_states
  if(exchange_functional.EQ."short_range_PBE")then
   call ex_pbe_sr(mu_erf,rho_a(istate),rho_b(istate),grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),ex(istate),vx_rho_a(istate),vx_rho_b(istate),vx_grad_rho_a_2(istate),vx_grad_rho_b_2(istate),vx_grad_rho_a_b(istate))
  else if(exchange_functional.EQ."None")then
   ex = 0.d0
   vx_rho_a = 0.d0
   vx_rho_b = 0.d0
   vx_grad_rho_a_2 = 0.d0
   vx_grad_rho_a_b = 0.d0
   vx_grad_rho_b_2 = 0.d0
  else 
   print*, 'Exchange functional required does not exist ...'
   print*,'exchange_functional',exchange_functional
   stop
  endif

  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
  if(correlation_functional.EQ."short_range_PBE")then
  ! convertion from (alpha,beta) formalism to (closed, open) formalism
   call rho_ab_to_rho_oc(rho_a(istate),rho_b(istate),rhoo,rhoc)
   call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),sigmaoo,sigmacc,sigmaco)

   call ec_pbe_sr(mu_erf,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,ec(istate),vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo)

   call v_rho_oc_to_v_rho_ab(vrhoo,vrhoc,vc_rho_a(istate),vc_rho_b(istate))
   call v_grad_rho_oc_to_v_grad_rho_ab(vsigmaoo,vsigmacc,vsigmaco,vc_grad_rho_a_2(istate),vc_grad_rho_b_2(istate),vc_grad_rho_a_b(istate))
  else if(correlation_functional.EQ."None")then
   ec = 0.d0
   vc_rho_a = 0.d0
   vc_rho_b = 0.d0
   vc_grad_rho_a_2 = 0.d0
   vc_grad_rho_a_b = 0.d0
   vc_grad_rho_b_2 = 0.d0
  else 
   print*, 'Correlation functional required does not exist ...'
   print*, 'correlation_functional',correlation_functional
   stop
  endif
 enddo
end


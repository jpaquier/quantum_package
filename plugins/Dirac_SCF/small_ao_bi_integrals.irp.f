 double precision function small_ao_bielec_integral(i,j,k,l)
  implicit none
  BEGIN_DOC
  !  integral of the AO basis <ik|jl> or (ij|kl)
  !     i(r1) j(r1) 1/r12 k(r2) l(r2)
  END_DOC
  integer,intent(in)             :: i,j,k,l
  integer                        :: p,q,r,s
  double precision               :: I_center(3),J_center(3),K_center(3),L_center(3)
  integer                        :: num_i,num_j,num_k,num_l,dim1,I_power(3),J_power(3),K_power(3),L_power(3)
  double precision               :: integral
  include 'Utils/constants.include.F'
  double precision               :: P_new(0:max_dim,3),P_center(3),fact_p,pp
  double precision               :: Q_new(0:max_dim,3),Q_center(3),fact_q,qq
  integer                        :: iorder_p(3), iorder_q(3)
  double precision               :: ao_bielec_integral_schwartz_accel
! if (ao_prim_num(i) * ao_prim_num(j) * ao_prim_num(k) * ao_prim_num(l) > 1024 ) then
!  ao_bielec_integral = ao_bielec_integral_schwartz_accel(i,j,k,l)
!  return
! endif
  dim1 = n_pt_max_integrals
  num_i = small_ao_nucl(i)
  num_j = small_ao_nucl(j)
  num_k = small_ao_nucl(k)
  num_l = small_ao_nucl(l)
  small_ao_bielec_integral = 0.d0
  if (num_i /= num_j .or. num_k /= num_l .or. num_j /= num_k)then
   do p = 1, 3
    I_power(p) = small_ao_power(i,p)
    J_power(p) = small_ao_power(j,p)
    K_power(p) = small_ao_power(k,p)
    L_power(p) = small_ao_power(l,p)
    I_center(p) = nucl_coord(num_i,p)
    J_center(p) = nucl_coord(num_j,p)
    K_center(p) = nucl_coord(num_k,p)
    L_center(p) = nucl_coord(num_l,p)
   enddo
  double precision               :: coef1, coef2, coef3, coef4
  double precision               :: p_inv,q_inv
  double precision               :: general_primitive_integral
   coef1 = small_ao_coef_normalized(i)
   coef2 = coef1*small_ao_coef_normalized(j)
   call give_explicit_poly_and_gaussian(P_new,P_center,pp,fact_p,iorder_p,&
        small_ao_expo(i),small_ao_expo(j),                 &
        I_power,J_power,I_center,J_center,dim1)
   p_inv = 1.d0/pp
   coef3 = coef2*small_ao_coef_normalized(k)
   coef4 = coef3*small_ao_coef_normalized(l)
   call give_explicit_poly_and_gaussian(Q_new,Q_center,qq,fact_q,iorder_q,&
        small_ao_expo(k),small_ao_expo(l),             &
        K_power,L_power,K_center,L_center,dim1)
   q_inv = 1.d0/qq
   integral = general_primitive_integral(dim1,              &
              P_new,P_center,fact_p,pp,p_inv,iorder_p,             &
              Q_new,Q_center,fact_q,qq,q_inv,iorder_q)
   small_ao_bielec_integral +=+  coef4 * integral
  else
   do p = 1, 3
    I_power(p) = small_ao_power(i,p)
    J_power(p) = small_ao_power(j,p)
    K_power(p) = small_ao_power(k,p)
    L_power(p) = small_ao_power(l,p)
   enddo
  double  precision              :: ERI
   coef1 = small_ao_coef_normalized(i)
   coef2 = coef1*small_ao_coef_normalized(j)
   coef3 = coef2*small_ao_coef_normalized(k)
   coef4 = coef3*small_ao_coef_normalized(l)
   integral = ERI(                                          &
              small_ao_expo(i),small_ao_expo(j),small_ao_expo(k),small_ao_expo(l),&
              I_power(1),J_power(1),K_power(1),L_power(1),         &
              I_power(2),J_power(2),K_power(2),L_power(2),         &
              I_power(3),J_power(3),K_power(3),L_power(3))
   small_ao_bielec_integral += coef4 * integral
  endif
 end


program test_dm
 implicit none
 read_wf = .True.
 touch read_wf 
 call print_laplacian
end

subroutine print_laplacian
 implicit none
 double precision :: r(3)
 integer :: i,nx,istate
 double precision :: dx,xmax
 istate = 1
 nx = 10
 xmax = 3.d0
 dx = xmax/dble(nx)
 r = 0.d0
 write(33,*)'r(1) , two_dm_laplacian / two_dm - two_dm_laplacian_HF / two_dm_HF , two_dm , two_dm_HF , two_dm_laplacian , two_dm_laplacian_HF , dm_a+dm_b'
 do i = 1, nx
  double precision :: two_dm,two_dm_laplacian,total_dm,two_dm_HF,two_dm_laplacian_HF,total_dm_HF,dm_a,dm_b 
  call dm_dft_alpha_beta_at_r(r,dm_a,dm_b)
  call spherical_averaged_two_dm_at_second_order(r,0.d0,istate,two_dm,two_dm_laplacian,total_dm)
  call spherical_averaged_two_dm_HF_at_second_order(r,0.d0,istate,two_dm_HF,two_dm_laplacian_HF,total_dm_HF)
 !corr_hole_2 = (two_dm_laplacian - two_dm_laplacian_HF * two_dm/total_dm_HF)/two_dm_HF
  two_dm = max(two_dm,1.d-15)
  two_dm_HF = max(two_dm_HF,1.d-15)
  write(33,'(100(F16.8,X))')r(1),two_dm_laplacian / two_dm - two_dm_laplacian_HF / two_dm_HF,two_dm,two_dm_HF,two_dm_laplacian,two_dm_laplacian_HF,dm_a+dm_b
!! approximated polynom 
! mu0 =  dpi * (two_dm_laplacian / two_dm - two_dm_laplacian_HF / two_dm_HF)
! mu0 = max(mu0,1.d-15)
! alpha = 1.d0 + 2.d0/(dsqrt(dacos(-1.d0)) * mu0)
! mu = mu0*alpha
!! new version with exact polynom 
! alpha_bis = (two_dm_laplacian / two_dm - two_dm_laplacian_HF / two_dm_HF)
! alpha_bis = max(alpha_bis,1.d-15)
! beta = 2.d0/(3.d0*dacos(-1.d0))
! delta = 2.d0/dsqrt(dacos(-1.d0))
! mu = 1.d0/(2.d0*beta)*(alpha_bis + dsqrt(alpha_bis*alpha_bis + 4.d0 * alpha_bis * beta * delta))
  r(1) += dx
 enddo


end


!subroutine test_energy_dm
!integer :: i,j,k,l
!double precision :: energy_2,get_mo_bielec_integral
!energy_2 = 0.d0
!do istate = 1, N_states
! do i = 1, mo_tot_num
!  do j = 1, mo_tot_num 
!   do k = 1, mo_tot_num
!    do l = 1, mo_tot_num
!     energy_2 += get_mo_bielec_integral(k,i,l,j,mo_integrals_map) * two_bod_alpha_beta_mo(k,l,j,i,istate)
!    enddo
!   enddo
!  enddo
! enddo
!enddo
!integer :: sze
!sze = mo_tot_num
!double precision, allocatable :: out_array(:,:)
!double precision :: integral
!allocate(out_array(sze,sze))
!do istate = 1, N_states
! do i = 1, mo_tot_num
!  do j = 1, mo_tot_num 
!   call get_mo_bielec_integrals_ij(i,j,sze,out_array,mo_map_integrals) 
!   do k = 1, mo_tot_num
!    do l = 1, mo_tot_num
!     integral = get_mo_bielec_integral(i,j,k,l,mo_integrals_map)
!     if(dabs(integral - out_array(k,l)).gt.1.d-10)then
!      print*,'i,j,k,l',i,j,k,l 
!      print*,dabs(integral - out_array(k,l)), integral, out_array(k,l)
!     endif
!    enddo
!   enddo
!  enddo
! enddo
!enddo
!
!print*,'psi_energy_bielec = ',psi_energy_bielec 
!print*,'test              = ',energy_2


!end

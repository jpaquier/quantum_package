program pouet
 read_wf = .True.
 touch read_wf
 !call test_delta      
 call test_delta_derivative 
end



subroutine test_delta  
 implicit none
 integer :: i,j,k,l,m,n
 print*,'*****************$******************'
 print*,'ECMD LDA standart =',Energy_c_md_LDA 
 print*,'ECMD LDA barth    =',Energy_c_md_LDA_barth
end


subroutine test_delta_derivative
 implicit none
 integer :: i,j,k,l,m,n
 double precision :: r(3),rs,rsplus,rsmoins,drs,accu
 double precision :: d_delta_barth_bis,delta_barth,d_delta_barth,delta_plus,delta_moins
 double precision :: aos_array(ao_num)
 double precision :: dn,rhot,rhos,rhoa,rhob,xi,mu,pi,denominator_delta,wignerseitz_radius,d_wignerseitz_radius

!!!!!!!!!!!delta2!!!!!!!!
 double precision :: accu_2,d_delta_barth_bis_2,delta_2,delta_plus_delta_2,delta_moins_delta_2,d_1st_deltaterm
!!!!!!!!!!!delta3!!!!!!!!
 double precision :: accu_3,d_delta_barth_bis_3,delta_3,delta_plus_delta_3,delta_moins_delta_3,d_2nd_deltaterm 
!!!!!!!!!!!delta4!!!!!!!!
 double precision :: accu_4,d_delta_barth_bis_4,delta_4,delta_plus_delta_4,delta_moins_delta_4,d_3rd_deltaterm
!!!!!!!!!!!delta5!!!!!!!!
 double precision :: accu_5,d_delta_barth_bis_5,delta_5,delta_plus_delta_5,delta_moins_delta_5,d_4th_deltaterm
!!!!!!!!!!!delta6!!!!!!!!
 double precision :: accu_6,d_delta_barth_bis_6,delta_6,delta_plus_delta_6,delta_moins_delta_6,d_5th_deltaterm


 pi=dacos(-1.d0)
 mu=mu_erf

 do n = 1, 16
!dn = 10d0**(-n)
 accu = 0d0
 accu_2= 0d0
 accu_3 = 0d0
 accu_5 = 0d0
 accu_4 = 0d0
 accu_6 = 0d0
  do j = 1,nucl_num
   do k = 1, n_points_radial_grid -1
    do l = 1 , n_points_integration_angular
     r(1) = grid_points_per_atom(1,l,k,j)
     r(2) = grid_points_per_atom(2,l,k,j)
     r(3) = grid_points_per_atom(3,l,k,j)
     call dm_dft_alpha_beta_and_all_aos_at_r(r,rhoa,rhob,aos_array)
     rhot = rhoa + rhob
     rhos = rhoa - rhob
     
     dn = rhot*10d0**(-n)

     xi = (rhoa-rhob)/(rhoa+rhob)
     rs = wignerseitz_radius(rhot) 
     drs= d_wignerseitz_radius(rhot)

     rsplus = wignerseitz_radius(rhot+dn)
     rsmoins = wignerseitz_radius(rhot-dn) 

!!!!!!!!!!!!!!!!!testsss delta2!!!!!!!!!
     
    delta_plus_delta_2 = delta_2(rsplus)*mu**2/denominator_delta(rsplus,xi,mu)
    delta_moins_delta_2 = delta_2(rsmoins)*mu**2/denominator_delta(rsmoins,xi,mu)
     
    d_delta_barth_bis_2= (delta_plus_delta_2-delta_moins_delta_2)/(2.d0 * dn)
 
    accu_2 += dabs(d_delta_barth_bis_2 - d_1st_deltaterm(rs,mu,xi)*drs)*final_weight_functions_at_grid_points(l,k,j)

!!!!!!!!!!!!!!!!!testsss delta3!!!!!!!!!

    delta_plus_delta_3 = delta_3(rsplus,xi)*mu**3/denominator_delta(rsplus,xi,mu)
    delta_moins_delta_3= delta_3(rsmoins,xi)*mu**3/denominator_delta(rsmoins,xi,mu)


    d_delta_barth_bis_3= (delta_plus_delta_3-delta_moins_delta_3)/(2.d0 * dn)

    accu_3 += dabs(d_delta_barth_bis_3-d_2nd_deltaterm(rs,mu,xi)*drs)*final_weight_functions_at_grid_points(l,k,j)


!!!!!!!!!!!!!!!!!TEST delta_4!!!!!!!!!
 
    delta_plus_delta_4 = delta_4(rsplus,xi)*mu**4/denominator_delta(rsplus,xi,mu)
    delta_moins_delta_4= delta_4(rsmoins,xi)*mu**4/denominator_delta(rsmoins,xi,mu)
         
 
    d_delta_barth_bis_4= (delta_plus_delta_4-delta_moins_delta_4)/(2.d0 * dn)
 
    accu_4 += dabs(d_delta_barth_bis_4-d_3rd_deltaterm(rs,mu,xi)*drs)*final_weight_functions_at_grid_points(l,k,j)



!!!!!!!!!!!!!!!!!TEST delta_5!!!!!!!!!

    delta_plus_delta_5 = delta_5(rsplus,xi)*mu**5/denominator_delta(rsplus,xi,mu)
    delta_moins_delta_5= delta_5(rsmoins,xi)*mu**5/denominator_delta(rsmoins,xi,mu)



    d_delta_barth_bis_5= (delta_plus_delta_5-delta_moins_delta_5)/(2.d0 * dn)

    accu_5 += dabs(d_delta_barth_bis_5-d_4th_deltaterm(rs,mu,xi)*drs)*final_weight_functions_at_grid_points(l,k,j)

!!!!!!!!!!!!!!!!!TEST delta_6!!!!!!!!!

    delta_plus_delta_6 = delta_6(rsplus,xi)*mu**6/denominator_delta(rsplus,xi,mu)
    delta_moins_delta_6= delta_6(rsmoins,xi)*mu**6/denominator_delta(rsmoins,xi,mu)
 
    d_delta_barth_bis_6= (delta_plus_delta_6-delta_moins_delta_6)/(2.d0 * dn)
 
    accu_6 += dabs(d_delta_barth_bis_6-d_5th_deltaterm(rs,mu,xi)*drs)*final_weight_functions_at_grid_points(l,k,j)
 
!!!!!!!!!!!!!!!!!TOTAL TEST!!!!!!!!!

    delta_plus = delta_barth(rsplus,xi,mu)
    delta_moins = delta_barth(rsmoins,xi,mu)

    d_delta_barth_bis = (delta_plus - delta_moins) /(2.d0 * dn)

    accu += dabs(d_delta_barth_bis - d_delta_barth(rs,xi,mu)*drs)*rhot*final_weight_functions_at_grid_points(l,k,j)

    enddo
   enddo
  enddo
  !print*,'dn=n*',10d0**(-n),accu_2,accu_3,accu_4,accu_5,accu_6
  print*,'dn=n*',10d0**(-n),accu
 enddo

end

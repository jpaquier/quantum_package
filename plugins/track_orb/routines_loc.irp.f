subroutine loc_cele_routine
     implicit none
     integer id1,i_atom,shift,shift_h
     parameter (id1=300)
     character*1 jobz,uplo
     character*64 file1,file2
     character*72 string(id1,8),cdum
     double precision               :: cmo(id1,id1,1),cmoref(id1,id1,1),newcmo(id1,id1,1)
     double precision               :: s(id1,id1,1),dum,ddum(id1,id1),ovl(id1,id1)
     double precision               :: w(id1),work(3*id1),t(id1,id1),wi(id1,id1)
     integer n,i,j,k,l,nmo(8),isym,nsym,idum,nrot(8),irot(id1,8)
     integer ipiv(id1),info,lwork
     logical *1 z54
     integer                        :: index_rot(1000,1)
     double precision               :: accu_norm

     z54=.false.
     accu_norm = 0.d0
     do i =1,mo_tot_num
       accu_norm += dabs(mo_overlap(i,i))
     enddo
     nsym = 1
     nmo(1) = mo_tot_num
     cmo = 0.d0
     do isym=1,nsym
       do i=1,nmo(isym)
         do j = 1, ao_num
           cmo(j,i,isym) = mo_coef(j,i)
         enddo
       enddo
     enddo
     do isym=1,nsym
       do j=1,mo_tot_num
         do i=1,ao_num
           newcmo(i,j,isym)=cmo(i,j,isym)
         enddo
       enddo
     enddo
     nrot(1) = n_act_orb   ! number of orbitals to be localized
     
     
     
     
     cmoref = 0.d0
     irot = 0
     
     do i = 1, n_act_orb
      irot(i,1) = list_act(i)
      do j = 1, ao_num
       cmoref(j,i,1)   = mo_coef_begin_iteration(j,irot(i,1))   !
      enddo
     enddo
     
     
     
     do i = 1, ao_num
       do j =1, ao_num
         s(i,j,1) =  ao_overlap(i,j)
       enddo
     enddo
     !Now big loop over symmetry
     
     
     
     do isym=1,nsym
       
       if (nrot(isym).eq.0) cycle
       do j=1,nrot(isym)
         do i=1,ao_num
           ddum(i,j)=0.d0
           do k=1,ao_num
             ddum(i,j)=ddum(i,j)+s(i,k,isym)*cmo(k,irot(j,isym),isym)
           enddo
         enddo
       enddo
       
       
       
       do i=1,nrot(isym)
         do j=1,nrot(isym)
           ovl(i,j)=0.d0
           do k=1,ao_num
             ovl(i,j)=ovl(i,j)+cmoref(k,i,isym)*ddum(k,j)
           enddo
         enddo
       enddo
       call maxovl_no_print(nrot(isym),nrot(isym),ovl,t,wi)
       do i=1,nrot(isym)
         do j=1,ao_num
           !         write (6,*) 'isym,',isym,nrot(isym),nmo(isym)
           newcmo(j,irot(i,isym),isym)=0.d0
           do k=1,nrot(isym)
             newcmo(j,irot(i,isym),isym)=newcmo(j,irot(i,isym),isym) + cmo(j,irot(k,isym),isym)*t(k,i)
           enddo
         enddo
       enddo
       
     enddo !big loop over symmetry
     
     10 format (4E19.12)
     
     
     !  Now we copyt the newcmo into the mo_coef
     
     mo_coef = 0.d0
     do isym=1,nsym
       do i=1,nmo(isym)
         do j = 1, ao_num
           mo_coef(j,i) = newcmo(j,i,isym)
         enddo
       enddo
     enddo
     !      pause
     
     
     ! we say that it hase been touched, and valid and that everything that
     ! depends on mo_coef must not be reprovided
     accu_norm = 0.d0
     do i =1,mo_tot_num
       accu_norm += dabs(mo_overlap(i,i))
     enddo
     
end

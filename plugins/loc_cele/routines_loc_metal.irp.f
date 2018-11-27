subroutine loc_cele_routine_metal
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
     nrot(1) = n_orb_loc ! number of orbitals to be localized
     
     
     integer                        :: index_rot(1000,1)
     
     
     cmoref = 0.d0
     irot = 0
     
  do i = 1, n_orb_loc
   irot(i,1) = index_loc_orb(i)
  enddo
  do i = 1, n_orb_loc
   do j = 1, ao_num
    cmoref(j,i,1)   = guess_loc_mos(j,i)
   enddo
  enddo

  do i = 1, nrot(1)
    print*,'irot(i,1) = ',irot(i,1)
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
       
       
       
      integer :: m
       do i=1,nrot(isym) ! mos
         do j=1,nrot(isym)! ref vector 
           ovl(i,j)=0.d0
            do m = 1, ao_num
             do k=1,ao_num
               ovl(j,i)=ovl(j,i)+cmoref(k,j,isym)*cmo(m,irot(i,isym),isym) * s(k,m,isym)
             enddo
            enddo
         enddo
       enddo
       call maxovl(nrot(isym),nrot(isym),ovl,t,wi)
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
     double precision               :: accu_norm
!    touch mo_coef
!    print*,'after  = '
     accu_norm = 0.d0
     do i =1,mo_tot_num
       accu_norm += dabs(mo_overlap(i,i))
     enddo
     print*, 'accu_norm = ',accu_norm
!    print*,'accu_norm = ',accu_norm
     ! We call the routine that saves mo_coef in the ezfio format
!    call save_mos
     
end

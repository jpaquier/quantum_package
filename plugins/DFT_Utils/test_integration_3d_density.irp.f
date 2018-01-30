program pouet
 print*,'coucou'
 read_wf = .True.
 touch read_wf
!print*,'m_knowles = ',m_knowles
 call routine3
end

subroutine routine3
 implicit none
 integer :: i,j,k,l
 double precision :: accu
 accu = 0.d0
 print*, 'energy_x',energy_x
 print*, 'energy_c',energy_c
end

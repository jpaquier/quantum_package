program pouet
 implicit none
 read_wf = .true. 
 touch read_wf
 call truncated_tucker_decomposition 
 !call tucker_decomposition
end  

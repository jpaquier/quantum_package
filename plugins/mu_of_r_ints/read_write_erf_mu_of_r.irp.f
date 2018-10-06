BEGIN_PROVIDER [ logical, read_ao_integrals_mu_of_r ]
&BEGIN_PROVIDER [ logical, read_ao_integrals_mu_of_r_left ]
&BEGIN_PROVIDER [ logical, read_mo_integrals_mu_of_r ]
&BEGIN_PROVIDER [ logical, write_ao_integrals_mu_of_r ]
&BEGIN_PROVIDER [ logical, write_ao_integrals_mu_of_r_left ]
&BEGIN_PROVIDER [ logical, write_mo_integrals_mu_of_r ]

 BEGIN_DOC
! One level of abstraction for disk_access_ao_integrals_mu_of_r and disk_access_mo_integrals_mu_of_r
 END_DOC
implicit none

    if (disk_access_ao_integrals_mu_of_r.EQ.'Read') then
        read_ao_integrals_mu_of_r =  .True.
        write_ao_integrals_mu_of_r = .False.

    else if  (disk_access_ao_integrals_mu_of_r.EQ.'Write') then
        read_ao_integrals_mu_of_r = .False.
        write_ao_integrals_mu_of_r =  .True.
    
    else if (disk_access_ao_integrals_mu_of_r.EQ.'None') then
        read_ao_integrals_mu_of_r = .False.
        write_ao_integrals_mu_of_r = .False.

    else
        print *, 'bielec_integrals_mu_of_r/disk_access_ao_integrals_mu_of_r has a wrong type'
        stop 1

    endif

    if (disk_access_ao_integrals_mu_of_r_left.EQ.'Read') then
        read_ao_integrals_mu_of_r_left =  .True.
        write_ao_integrals_mu_of_r_left = .False.

    else if  (disk_access_ao_integrals_mu_of_r_left.EQ.'Write') then
        read_ao_integrals_mu_of_r_left = .False.
        write_ao_integrals_mu_of_r_left =  .True.
    
    else if (disk_access_ao_integrals_mu_of_r_left.EQ.'None') then
        read_ao_integrals_mu_of_r_left = .False.
        write_ao_integrals_mu_of_r_left = .False.

    else
        print *, 'bielec_integrals_mu_of_r/disk_access_ao_integrals_mu_of_r_left has a wrong type'
        stop 1

    endif

    if (disk_access_mo_integrals_mu_of_r.EQ.'Read') then
        read_mo_integrals_mu_of_r =  .True.
        write_mo_integrals_mu_of_r = .False.

    else if  (disk_access_mo_integrals_mu_of_r.EQ.'Write') then
        read_mo_integrals_mu_of_r = .False.
        write_mo_integrals_mu_of_r =  .True.

    else if (disk_access_mo_integrals_mu_of_r.EQ.'None') then
        read_mo_integrals_mu_of_r = .False.
        write_mo_integrals_mu_of_r = .False.

    else
        print *, 'bielec_integrals_mu_of_r/disk_access_mo_integrals_mu_of_r has a wrong type'
        stop 1

    endif

END_PROVIDER

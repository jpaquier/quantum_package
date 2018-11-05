#!/bin/bash
# specify the QP folder 
QP=$QP_ROOT
# sourcing the quantum_package.rc file
. ${QP}/quantum_package.rc



ezfio=$1
if [ "$ezfio" = "--help" ]; then
 echo "                             "
 echo "                             "
 echo "                             "
 echo "********            *********"
 echo "This script runs a RSH-CIPSI calculation starting from a ezfio folder"
 echo "Here are the ASSUMTIONS"
 echo "******** ASSUMPTIONS *********"
 echo "1) The mo_class should be specified in the EZFIO folder"
 echo "2) The wave function stored in the EZFIO folder, whatever it is, will be used as the starting point to run the RSH calculation"
 echo "********            *********"
 echo "                             "
 echo "Here are the INPUTS ARGUMENT for the script"
 echo "******** INPUTS OF THE SCRIPT *********"
 echo "\$1 : ezfio = EZFIO folder on which you are going to run the script"
 echo "     no orbitals would be changed, it will only modify the wave function stored in the EZFIO folder"
 echo "\$2 : functional = EXCHANGE / CORRELATION functionals to be used in RS-DFT calculation"
 echo "     can be LDA or PBE"
 echo "\$3 : mu = SPLITTING of the INTERACTION to be used in RS-DFT calculation "
 echo "      0 <=> pure DFT calculation, Infinity <=> pure WFT calculation "
 echo "\$4 : pt2max  : maximum value of the PT2 for the CIPSI calculation (note that it is with the effective hamiltonian so it can be self-consistent) "
 echo "      0 <=> FCI calculation, any bigger value is selected CI"
 echo "\$5 : thresh : value of the convergence of the energy for the self-consistent CIPSI calculation for a given set of Slater determinants"  
 echo "      a value of 0.0000001 is in general more than sufficient for regular calculations"
 echo "********            *********"
 echo "                             "
 echo "                             "
 exit
fi
if [ -f $ezfio ]; then
   echo "File $ezfio exists."
else
   echo "Input EZFIO folder does not exists !"
   echo "Folder $ezfio does not exist."
   echo "stopping ..."
   echo "......"
   echo "......"
   echo "......"
fi
exit
# define the exchange / correlation functionals to be used in RS-DFT calculation
functional=$2
# splitting of the interaction to be used in RS-DFT calculation 
mu=$3
# maximum value of the PT2 for the CIPSI calculation (note that it is with the effective hamiltonian so it can be self-consistent)
pt2max=$4
# value of the convergence of the energy for the self-consistent CIPSI calculation at a given number of determinant
thresh=$5


################################################## CREATION OF THE EZFIO FOLDER ##########################################################

################################################## RUNNING THE RS-KS-DFT CALCULATION #####################################################
# set the exchange / correlation functionals 
echo "short_range_${functional}" > ${ezfio}/dft_keywords/exchange_functional
echo "short_range_${functional}" > ${ezfio}/dft_keywords/correlation_functional
# set the mu value of the splitting of the bi-electronic interaction
echo  $mu                        > ${ezfio}/dft_keywords/mu_erf 

################################################## RUNNING THE SELF-CONSISTENT CIPSI-RS-DFT CALCULATION  #################################
# set the maximum value of the PT2 for CIPSI calculation 
echo $pt2max > ${ezfio}/perturbation/pt2_max

# ####### INITIALIZATION OF THE RS-DFT CALCULATION : CIPSI WITH AN EFFECTIVE HAMILTONIAN BUILT WITH THE RS-KS DENSITY ################## #
# specify that you use the wave function stored in the EZFIO (i.e. RS_KS) to build the density used in the construction of the effective short-range potential 
echo "WFT"  > ${ezfio}/dft_keywords/density_for_dft
qp_edit -c ${ezfio}

# write the effective Hamiltonian containing long-range interaction and short-range effective potential to be diagonalized in a self-consistent way
qp_run write_integrals_restart_dft_no_ecmd ${ezfio} | tee rsdft-0

# save the RS-KS one-body density for the damping on the density 
qp_run save_one_body_dm ${ezfio} 
# specify that you will do some damping on the density later on 
echo "damping_rs_dft"  > ${ezfio}/dft_keywords/density_for_dft
# specify the damping factor on the density : 0 == no update of the density, 1 == full update of the density 
echo "0.75"            > ${ezfio}/dft_keywords/damping_for_rs_dft

for i in {1..3}
do
#  run the CIPSI calculation with the effective Hamiltonian already stored in the EZFIO folder 
   qp_run fci_zmq ${ezfio} | tee H2O-fci-$i
   # run 
   EV=0

   for j in {1..100}
   do
      # write the new effective Hamiltonian with the damped density (and the current density to be damped with the next density)
      qp_run write_integrals_restart_dft_no_ecmd ${ezfio} | tee rsdft-${i}-${j}
      # value of the variational RS-DFT energy 
      EV_new=`grep "TOTAL ENERGY        =" rsdft-${i}-${j} | cut -d "=" -f 2`
      # rediagonalize the new effective Hamiltonian to obtain a new wave function and a new density 
      qp_run diag_restart_save_lowest_state ${ezfio} | tee diag-${i}-${j}
      # checking the convergence
      DE=`echo "${EV} - ${EV_new}" | bc`
      DE=`echo "print abs(${DE})" | python `
      CONV=`echo "${DE} > ${thresh}" | bc`
      echo $j $EV $EV_new $DE $thresh $CONV >> data_conv_${i}
      if [ "$CONV" -eq "0" ]; then
        break
      fi
      DE=0
      DE1=0
      EV=$EV_new
    done
    qp_run write_integrals_restart_dft_no_ecmd ${ezfio} | tee rsdft-${i}-final
done

#!/bin/bash
# specify the QP folder 
QP=$QP_ROOT
# sourcing the quantum_package.rc file
. ${QP}/quantum_package.rc

# define the xyz file to be used
xyzfile=H2O
# define the name of the ezfio folder 
ezfio=${xyzfile}.ezfio
# define the basis set for the oxygen and hydrogen atoms 
basis_O=cc-pcvtz
basis_H=cc-pvtz
# define the exchange / correlation functionals to be used in RS-DFT calculation
functional=PBE
# splitting of the interaction to be used in RS-DFT calculation 
mu=0.5
# maximum value of the PT2 for the CIPSI calculation (note that it is with the effective hamiltonian so it can be self-consistent)
pt2max="0.0001"
# value of the convergence of the energy for the self-consistent CIPSI calculation at a given number of determinant
thresh=0.00000001 # 10^-10


################################################## CREATION OF THE EZFIO FOLDER ##########################################################
# create the ezfio folder 
qp_create_ezfio_from_xyz -b "O:${basis_O} | H:${basis_H}" ${xyzfile}.xyz -o $ezfio
# check that everything is correct in input file
qp_edit -c ${ezfio}

################################################## RUNNING THE RS-KS-DFT CALCULATION #####################################################
# set the exchange / correlation functionals 
echo "short_range_${functional}" > ${ezfio}/dft_keywords/exchange_functional
echo "short_range_${functional}" > ${ezfio}/dft_keywords/correlation_functional
# set the mu value of the splitting of the bi-electronic interaction
echo  $mu                        > ${ezfio}/dft_keywords/mu_erf 
qp_edit -c ${ezfio}

# run a range separated KS DFT calculation 
qp_run RS_KS_SCF ${ezfio} | tee H2O_RS_KS_SCF.out

################################################## RUNNING THE SELF-CONSISTENT CIPSI-RS-DFT CALCULATION  #################################
# set the maximum value of the PT2 for CIPSI calculation 
echo $pt2max > ${ezfio}/perturbation/pt2_max

# ####### INITIALIZATION OF THE RS-DFT CALCULATION : CIPSI WITH AN EFFECTIVE HAMILTONIAN BUILT WITH THE RS-KS DENSITY ################## #
# specify that you use the wave function stored in the EZFIO (i.e. RS_KS) to build the density used in the construction of the effective short-range potential 
echo "WFT"  > ${ezfio}/dft_keywords/density_for_dft
# write the effective Hamiltonian containing long-range interaction and short-range effective potential to be diagonalized in a self-consistent way
qp_run write_integrals_restart_dft_no_ecmd ${ezfio} | tee ${ezfio}_rsdft-0

# save the RS-KS one-body density for the damping on the density 
qp_run save_one_body_dm ${ezfio} 
# specify that you will do some damping on the density later on 
echo "damping_rs_dft"  > ${ezfio}/dft_keywords/density_for_dft
# specify the damping factor on the density : 0 == no update of the density, 1 == full update of the density 
echo "0.75"            > ${ezfio}/dft_keywords/damping_for_rs_dft

for i in {1..3}
do
#  run the CIPSI calculation with the effective Hamiltonian already stored in the EZFIO folder 
   qp_run fci_zmq ${ezfio} | tee ${ezfio}/fci-$i
   # run 
   EV=0

   echo "#" iter evar old     evar new    deltae      threshold  >> ${ezfio}_data_conv_${i}
   for j in {1..100}
   do
      # write the new effective Hamiltonian with the damped density (and the current density to be damped with the next density)
      qp_run write_integrals_restart_dft_no_ecmd ${ezfio} | tee ${ezfio}/rsdft-${i}-${j}
      # value of the variational RS-DFT energy 
      EV_new=`grep "TOTAL ENERGY        =" ${ezfio}/rsdft-${i}-${j} | cut -d "=" -f 2`
      # rediagonalize the new effective Hamiltonian to obtain a new wave function and a new density 
      qp_run diag_restart_save_lowest_state ${ezfio} | tee ${ezfio}/diag-${i}-${j}
      # checking the convergence
      DE=`echo "${EV_new} - ${EV}" | bc`
      DEabs=`echo "print abs(${DE})" | python `
      CONV=`echo "print ${DEabs} < ${thresh}" | python`
      echo $j $EV $EV_new $DE $DEabs $thresh >> ${ezfio}_data_conv_${i}
      if [ "$CONV" = "True" ]; then
        break
      fi
      EV=$EV_new
    done
    qp_run write_integrals_restart_dft_no_ecmd ${ezfio} | tee rsdft-${i}-final
done

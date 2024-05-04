#!/bin/csh -f
#TEST-PROGRAM sander
#TEST-DESCRIP TO_BE_DEtermined
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

set sander = "${AMBERHOME}/bin/sander"
if( $?TESTsander ) then
    set sander = $TESTsander
endif

set can_do_mpi=0
if( `echo $sander | grep "MPI"` != '') then
 echo 'Using MPI version of sander'
 set can_do_mpi=1
endif

set mpi=0
if( ! $?DO_PARALLEL ) then
  set DO_PARALLEL=' '
  echo "Running test using file-based exchange (DO_PARALLEL not set)"
  echo "Parallelism treated within TeraChem"
else
  if( $can_do_mpi == 1 ) then
    echo "Running test using MPI data exchange (DO_PARALLEL is set)"
    echo "Parallelism treated within TeraChem"
    set mpi=1
    set numprocs = `$DO_PARALLEL echo | wc -l`
    if ( $numprocs > 1 ) then
	echo "This test requires >> mpirun -np 1 <<"
	echo "You specified >> mpirun -np $numprocs <<"
	exit 0
    endif
  else
    echo "Can not do an MPI run using serial version of sander."
    echo "Unset DO_PARALLEL to perform serial test"
    exit(1)
  endif
endif

cat > mdin <<EOF
NR QM/MM/ CAM-B3LYP production run correlated (dt=0.5fs, no shake)
 &cntrl
  imin   = 0,           !no minimization
  irest  = 1,           !restart
  ntx    = 5,           !coordinates and velocities are read in
  cut    = 999.0,       !non-bonded interactions cutoff
  dt     = 0.0005,      !0.5fs time step
  ntb    = 0,           !no periodicity and PME off!
  ntt    = 1,           !turn off thermostat
  tautp  = 1.0,         !time constant for heat coupling 1.0ps
  tempi  = 300.0,       !initial temperature
  temp0  = 300.0,       !temp = 300K
  ntpr   = 50,          !print details to log every step
  ntwx   = 200,         !write coordinates to mdcrd every 500 steps (1000 Snaps)
  ntwr   = 200,         !write restart file at last step
  nstlim = 100000,    	!run for 50 ps  
  nscm   = 2000,	!Remove COM every #nscm steps
  ifqnt  = 1,		!Turn on QM/MM dynamics
 /
 &ewald 
 dsum_tol = 0.000001	!width of the direct sum part of the Ewald sum
 /
 &qmmm
  qmmask    = ':1',	!Chomophore Residue-id
  qm_theory = 'dftb3',	!Calling External QM Softwater
  qmmm_int  = 1,	!Electrostatic Embedding
  qmcharge = 0,		!Charge on the System
  spin = 1,		!Multiplicity
 /
EOF


set output = qmmm.mdout
set restrt = qmmm2.restrt
rm -f $output


touch dummy
$DO_PARALLEL $sander -O -p ../solute_solvent_drop.prmtop -c ../equil.ncrst -o $output -r $restrt -x gs_nbdnh2_dmso.mdcrd -inf qmmm.info< dummy || goto error

# We do this due to rounding errors when reading from a file in the file-based 
# data exchange vs receiving the data through MPI
${AMBERHOME}/test/dacdif -a 0.0003  $output.save $output
${AMBERHOME}/test/dacdif -a 0.0003  $restrt.save $restrt

#/bin/rm -f mdin mdinfo mdcrd dummy inpfile.xyz ptchrg.xyz startfile tc_job* old.tc*
#/bin/rm -r scr
exit(0)

error:
echo "  ${0}:  Program error"
exit(1)



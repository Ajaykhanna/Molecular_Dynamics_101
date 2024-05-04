#!/bin/csh -f
#TEST-PROGRAM sander
#TEST-DESCRIP TO_BE_DEtermined
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

set sander = "${AMBERHOME}/bin/sander"
if( $?TESTsander ) then
    set sander = $TESTsander
endif

/home/akhanna2/data/test_amber_modifications/amber20/test/check_GAUSSIAN.x
if( $status > 0) then
  exit(0)
endif

if( $?DO_PARALLEL ) then
   echo "Running $TESTsander in parallel on "`$DO_PARALLEL echo | wc -l`" cores."
else set DO_PARALLEL = ''
endif

# check on how many CPU core Gaussian shall be running
if ( $?GAU_NCPUS) then
    if ( $GAU_NCPUS > 1 ) then
        echo "GAUSSIAN will run on $GAU_NCPUS CPU cores"
    endif
else
    set GAU_NCPUS = 1
    echo "GAUSSIAN will run on 1 CPU core"
    echo "(environment variable GAU_NCPUS not set)"
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
  qm_theory = 'EXTERN',	!Calling External QM Softwater
  qmmm_int  = 1,	!Electrostatic Embedding
  qmcharge = 0,		!Charge on the System
  spin = 1,		!Multiplicity
 /
 &gau
  method   = 'CAM-B3LYP',
  basis    = '6-31G(d)',
  num_threads = 32,
  mem         = 100GB, 
 /
EOF

set output = qmmm.mdout
set restrt = qmmm2.restrt
rm -f $output

touch dummy
$DO_PARALLEL $sander -O -p ../solute_solvent_drop.prmtop -c ../equil.ncrst -o $output -r $restrt -x p_nitroaniline_gs.mdcrd -inf qmmm.info< dummy || goto error

${AMBERHOME}/test/dacdif -t 1 $output.save $output
${AMBERHOME}/test/dacdif -a 0.00002 $restrt.save $restrt

/bin/rm -f dummy mdin mdinfo mdcrd restrt inpfile*.xyz gau_job* old.gau_job* fort*
/bin/rm -rf 000
exit(0)

error:
echo "  ${0}:  Program error"
exit(1)



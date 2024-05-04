Example: 7BenzylAmino-NBD in DMSO
Initial Geometry Optimization and Level of Theory for both solvent and solute in Gaussian 16
# Optimization
opt=(calcall,tight) freq=(savenm,hpmodes) cam-b3lyp/6-31g(d)
scrf=(iefpcm, solvent=dmso) nosymm geom=connectivity int=ultrafine

# RESP Charge
HF/6-31G* nosymm iop(6/33=2) pop(chelpg,regular)
scrf=(iefpcm, solvent=dmso) nosymm geom=AllCheck guess=Read

Forcefield parameterization:
Antechamber --> Charges: RESP --> GAFF2
Radius of Droplet 40A

AIMD Level of Theory:
&qmmm
  qmmask    = ':1',     !Chomophore Residue-id
  qm_theory = 'EXTERN', !Calling External QM Softwater
  qmmm_int  = 1,        !Electrostatic Embedding
  qmcharge = 0,         !Charge on the System
  spin = 1,             !Multiplicity
 /
 &tc
  method   = 'camb3lyp',!QM Level of Theory
  basis    = '6-31g*',  !Basis-Set
  ngpus    = 1,         !#GPUs
/


Disclaimer: The files in the repo are only for tutorial purposes. Please check them carefully before performing any tests or experiments with them.
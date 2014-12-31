#This model performs multi-algorithm simulation using Spatiocyte next reaction
#method with lattice-based particle reaction-diffusion process.
#A is a HD species with molecules assumed to be homogeneously distributed in the
#compartment. B is a individually diffused species by the particle simulator.
#This model executes the reaction A + B -> C using Spatiocyte next reaction #method.

Stepper SpatiocyteStepper(SS) { VoxelRadius 4.4e-9; }
System System(/) {
  StepperID SS;
  Variable Variable(GEOMETRY) { Value 0; }
  Variable Variable(LENGTHX) { Value 5e-7; }
  Variable Variable(LENGTHY) { Value 5e-7; }
  Variable Variable(LENGTHZ) { Value 5e-7; }
  Variable Variable(VACANT) { Value 0; }
  Variable Variable(A) {
      Value 1000;
      Name "HD"; }
  Variable Variable(B) { Value 1500; }
  Variable Variable(C) { Value 0; }
  Process MoleculePopulateProcess(populate) {
    VariableReferenceList [_ Variable:/:B]; }
  Process DiffusionProcess(diffuseB) {
    VariableReferenceList [_ Variable:/:B];
    D 1e-12; } # m^2/s
  Process DiffusionProcess(diffuseC) {
    VariableReferenceList [_ Variable:/:C];
    D 1e-12; } # m^2/s
  Process VisualizationLogProcess(visualize) {
    VariableReferenceList [_ Variable:/:B]
                          [_ Variable:/:C];
    LogInterval 5e-5; } # s
  Process IteratingLogProcess(logiter) {
    VariableReferenceList [_ Variable:/:A]
                          [_ Variable:/:B]
                          [_ Variable:/:C];
    LogEnd 0.03;
    LogInterval 1e-4; }
  Process SpatiocyteNextReactionProcess(reaction) {
    VariableReferenceList [_ Variable:/:A -1]
                          [_ Variable:/:B -1]
                          [_ Variable:/:C 1];
    k 1e-20; }
}

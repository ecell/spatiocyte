Stepper SpatiocyteStepper(SS) { 
  VoxelRadius 4.4e-9;
  SearchVacant 0; }
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
  Variable Variable(B) { Value 2000; }
  Variable Variable(C) { Value 0; }
  Process MoleculePopulateProcess(populate) {
    VariableReferenceList [_ Variable:/:B]; }
  Process DiffusionProcess(diffuseB) {
    VariableReferenceList [_ Variable:/:B];
    D 16e-12; } # m^2/s
  Process DiffusionProcess(diffuseC) {
    VariableReferenceList [_ Variable:/:C];
    D 16e-12; } # m^2/s
  Process VisualizationLogProcess(visualize) {
    VariableReferenceList [_ Variable:/:B]
                          [_ Variable:/:C];
    LogInterval 1e-3; } # s
  Process IteratingLogProcess(logiter) {
    VariableReferenceList [_ Variable:/:A]
                          [_ Variable:/:B]
                          [_ Variable:/:C];
    LogEnd 10;
    LogStart 0;
    LogInterval 1e-3; }
  Process SpatiocyteNextReactionProcess(reaction) {
    VariableReferenceList [_ Variable:/:A -1]
                          [_ Variable:/:B -1]
                          [_ Variable:/:C 1];
    k 1e-20; }
}

Filaments = 13
RotateAngle = 0
MTRadius = 12.5e-9
VoxelRadius = 0.4e-8
KinesinRadius = 0.4e-8
totalMTLength = 0.6e-6

theSimulator.createStepper('SpatiocyteStepper', 'SS').VoxelRadius = VoxelRadius
theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = 0.7e-6
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = 0.15e-6
theSimulator.createEntity('Variable', 'Variable:/:LENGTHZ').Value = 0.15e-6

theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:Kinesin').Value = 90
theSimulator.createEntity('Variable', 'Variable:/:MTKinesin' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:MTKinesinATP' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:Tubulin' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:TubulinM' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:TubulinP' ).Value = 0

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:populate')
populator.VariableReferenceList = [['_', 'Variable:/:MTKinesin']]

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:populateK')
populator.VariableReferenceList = [['_', 'Variable:/:Kinesin']]

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:detachPlus')
react.VariableReferenceList = [['_', 'Variable:/:MTKinesin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:TubulinP','-1']]
react.VariableReferenceList = [['_', 'Variable:/:Kinesin','1']]
react.VariableReferenceList = [['_', 'Variable:/:TubulinP','1']]
react.p = 1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:explicitAttach')
react.VariableReferenceList = [['_', 'Variable:/:Kinesin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:Tubulin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:MTKinesin','1']]
react.k = 2.5863133e-18

diffuse = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseKinesin')
diffuse.VariableReferenceList = [['_', 'Variable:/:Kinesin']]
diffuse.D = 4e-12

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:detach')
react.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP','-1']]
react.VariableReferenceList = [['_', 'Variable:/:Tubulin','1']]
react.VariableReferenceList = [['_', 'Variable:/:Kinesin','1']]
react.SearchVacant = 1
react.k = 2.5

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:hydrolysis')
react.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP','-1']]
react.VariableReferenceList = [['_', 'Variable:/:MTKinesin','1']]
react.SearchVacant = 1
react.k = 100000

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:phosphorylate')
react.VariableReferenceList = [['_', 'Variable:/:MTKinesin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP','1']]
react.SearchVacant = 1
react.k = 14500

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:ratchet')
react.VariableReferenceList = [['_', 'Variable:/:MTKinesin','-1']]
react.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP','1']]
react.VariableReferenceList = [['_', 'Variable:/:Tubulin','1']]
react.BindingSite = 1
react.k = 55000

diffuse = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePlus')
diffuse.VariableReferenceList = [['_', 'Variable:/:MTKinesin']]
diffuse.VariableReferenceList = [['_', 'Variable:/:Tubulin', '1']]
diffuse.D = 0.04e-10

coord = theSimulator.createEntity('CoordinateLogProcess', 'Process:/:coord')
coord.VariableReferenceList = [['_', 'Variable:/:Tubulin', '-1' ]]
coord.VariableReferenceList = [['_', 'Variable:/:TubulinM']]
coord.VariableReferenceList = [['_', 'Variable:/:TubulinP']]
coord.VariableReferenceList = [['_', 'Variable:/:Kinesin']]
coord.VariableReferenceList = [['_', 'Variable:/:MTKinesin']]
coord.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP']]
coord.LogInterval = 8e-6

visualLogger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:visualLogger')
visualLogger.VariableReferenceList = [['_', 'Variable:/:Tubulin']]
visualLogger.VariableReferenceList = [['_', 'Variable:/:TubulinM']]
visualLogger.VariableReferenceList = [['_', 'Variable:/:TubulinP']]
visualLogger.VariableReferenceList = [['_', 'Variable:/:Kinesin']]
visualLogger.VariableReferenceList = [['_', 'Variable:/:MTKinesin']]
visualLogger.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP']]
visualLogger.LogInterval = 8e-6

Microtubule = theSimulator.createEntity('MicrotubuleProcess', 'Process:/:Microtubule')
Microtubule.OriginX = 0
Microtubule.OriginY = 0
Microtubule.OriginZ = 0
Microtubule.RotateX = 0
Microtubule.RotateY = 0
Microtubule.RotateZ = RotateAngle
Microtubule.Radius = MTRadius
Microtubule.SubunitRadius = KinesinRadius
Microtubule.Length = totalMTLength
Microtubule.Filaments = Filaments
Microtubule.Periodic = 0
Microtubule.VariableReferenceList = [['_', 'Variable:/:MTKinesin' ]]
Microtubule.VariableReferenceList = [['_', 'Variable:/:MTKinesinATP' ]]
Microtubule.VariableReferenceList = [['_', 'Variable:/:Tubulin' , '-1']]
Microtubule.VariableReferenceList = [['_', 'Variable:/:TubulinM' , '-2']]
Microtubule.VariableReferenceList = [['_', 'Variable:/:TubulinP' , '-3']]

run(0.01)


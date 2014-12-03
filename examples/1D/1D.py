
theSimulator.createStepper('SpatiocyteStepper', 'SS').VoxelRadius = 0.4e-8
theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = 0.7e-6
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = 0.15e-6
theSimulator.createEntity('Variable', 'Variable:/:LENGTHZ').Value = 0.15e-6

theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:A').Value = 10
theSimulator.createEntity('Variable', 'Variable:/:A_Filament' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:Subunit' ).Value = 0

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:populate')
populator.VariableReferenceList = [['_', 'Variable:/:A_Filament']]

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:populateK')
populator.VariableReferenceList = [['_', 'Variable:/:A']]

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:explicitAttach')
react.VariableReferenceList = [['_', 'Variable:/:A','-1']]
react.VariableReferenceList = [['_', 'Variable:/:Subunit','-1']]
react.VariableReferenceList = [['_', 'Variable:/:A_Filament','1']]
react.k = 2.5863133e-18

diffuse = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseA')
diffuse.VariableReferenceList = [['_', 'Variable:/:A']]
diffuse.D = 4e-12

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:detach')
react.VariableReferenceList = [['_', 'Variable:/:A_Filament','-1']]
react.VariableReferenceList = [['_', 'Variable:/:Subunit','1']]
react.VariableReferenceList = [['_', 'Variable:/:A','1']]
react.SearchVacant = 1
react.k = 10000

diffuse = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePlus')
diffuse.VariableReferenceList = [['_', 'Variable:/:A_Filament']]
diffuse.VariableReferenceList = [['_', 'Variable:/:Subunit', '1']]
diffuse.D = 0.04e-10

coord = theSimulator.createEntity('CoordinateLogProcess', 'Process:/:coord')
coord.VariableReferenceList = [['_', 'Variable:/:Subunit']]
coord.VariableReferenceList = [['_', 'Variable:/:A']]
coord.VariableReferenceList = [['_', 'Variable:/:A_Filament']]
coord.LogInterval = 8e-6

visualLogger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:visualLogger')
visualLogger.VariableReferenceList = [['_', 'Variable:/:Subunit']]
visualLogger.VariableReferenceList = [['_', 'Variable:/:A']]
visualLogger.VariableReferenceList = [['_', 'Variable:/:A_Filament']]
visualLogger.LogInterval = 8e-6

Filament = theSimulator.createEntity('FilamentProcess', 'Process:/:Filament')
Filament.OriginX = 0
Filament.OriginY = 0
Filament.OriginZ = 0
Filament.RotateX = 0
Filament.RotateY = 0
Filament.RotateZ = 0.2
Filament.SubunitRadius = 0.4e-8
Filament.Length = 0.6e-6
Filament.Periodic = 1
Filament.VariableReferenceList = [['_', 'Variable:/:A_Filament' ]]
Filament.VariableReferenceList = [['_', 'Variable:/:Subunit' , '-1']]

run(0.01)


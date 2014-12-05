
theSimulator.createStepper('SpatiocyteStepper', 'SS').VoxelRadius = 0.4e-8
theSimulator.rootSystem.StepperID = 'SS'
theSimulator.createEntity('Variable', 'Variable:/:LENGTHX').Value = 0.7e-6
theSimulator.createEntity('Variable', 'Variable:/:LENGTHY').Value = 0.15e-6
theSimulator.createEntity('Variable', 'Variable:/:LENGTHZ').Value = 0.15e-6

theSimulator.createEntity('Variable', 'Variable:/:VACANT')
theSimulator.createEntity('Variable', 'Variable:/:A').Value = 20
theSimulator.createEntity('Variable', 'Variable:/:A_Filament' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:A_FilamentATP' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:Subunit' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:SubunitM' ).Value = 0
theSimulator.createEntity('Variable', 'Variable:/:SubunitP' ).Value = 0

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:populate')
populator.VariableReferenceList = [['_', 'Variable:/:A_Filament']]

populator = theSimulator.createEntity('MoleculePopulateProcess', 'Process:/:populateK')
populator.VariableReferenceList = [['_', 'Variable:/:A']]

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:detachPlus')
react.VariableReferenceList = [['_', 'Variable:/:A_Filament','-1']]
react.VariableReferenceList = [['_', 'Variable:/:SubunitP','-1']]
react.VariableReferenceList = [['_', 'Variable:/:A','1']]
react.VariableReferenceList = [['_', 'Variable:/:SubunitP','1']]
react.p = 1

react = theSimulator.createEntity('DiffusionInfluencedReactionProcess', 'Process:/:explicitAttach')
react.VariableReferenceList = [['_', 'Variable:/:A','-1']]
react.VariableReferenceList = [['_', 'Variable:/:Subunit','-1']]
react.VariableReferenceList = [['_', 'Variable:/:A_Filament','1']]
react.p = 1

diffuse = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffuseA')
diffuse.VariableReferenceList = [['_', 'Variable:/:A']]
diffuse.D = 4e-12

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:detach')
react.VariableReferenceList = [['_', 'Variable:/:A_FilamentATP','-1']]
react.VariableReferenceList = [['_', 'Variable:/:Subunit','1']]
react.VariableReferenceList = [['_', 'Variable:/:A','1']]
react.SearchVacant = 1
react.k = 150

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:hydrolysis')
react.VariableReferenceList = [['_', 'Variable:/:A_FilamentATP','-1']]
react.VariableReferenceList = [['_', 'Variable:/:A_Filament','1']]
react.SearchVacant = 1
react.k = 100000

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:phosphorylate')
react.VariableReferenceList = [['_', 'Variable:/:A_Filament','-1']]
react.VariableReferenceList = [['_', 'Variable:/:A_FilamentATP','1']]
react.SearchVacant = 1
react.k = 14500

react = theSimulator.createEntity('SpatiocyteTauLeapProcess', 'Process:/:ratchet')
react.VariableReferenceList = [['_', 'Variable:/:A_Filament','-1']]
react.VariableReferenceList = [['_', 'Variable:/:A_FilamentATP','1']]
react.VariableReferenceList = [['_', 'Variable:/:Subunit','1']]
react.BindingSite = 1
react.k = 55000

diffuse = theSimulator.createEntity('DiffusionProcess', 'Process:/:diffusePlus')
diffuse.VariableReferenceList = [['_', 'Variable:/:A_Filament']]
diffuse.VariableReferenceList = [['_', 'Variable:/:Subunit', '1']]
diffuse.D = 0.04e-10

coord = theSimulator.createEntity('CoordinateLogProcess', 'Process:/:coord')
coord.VariableReferenceList = [['_', 'Variable:/:Subunit']]
coord.VariableReferenceList = [['_', 'Variable:/:SubunitM']]
coord.VariableReferenceList = [['_', 'Variable:/:SubunitP']]
coord.VariableReferenceList = [['_', 'Variable:/:A']]
coord.VariableReferenceList = [['_', 'Variable:/:A_Filament']]
coord.VariableReferenceList = [['_', 'Variable:/:A_FilamentATP']]
coord.LogInterval = 8e-6

visualLogger = theSimulator.createEntity('VisualizationLogProcess', 'Process:/:visualLogger')
visualLogger.VariableReferenceList = [['_', 'Variable:/:Subunit']]
visualLogger.VariableReferenceList = [['_', 'Variable:/:SubunitM']]
visualLogger.VariableReferenceList = [['_', 'Variable:/:SubunitP']]
visualLogger.VariableReferenceList = [['_', 'Variable:/:A']]
visualLogger.VariableReferenceList = [['_', 'Variable:/:A_Filament']]
visualLogger.VariableReferenceList = [['_', 'Variable:/:A_FilamentATP']]
visualLogger.LogInterval = 8e-6

visualLogger = theSimulator.createEntity('MicroscopyTrackingProcess', 'Process:/:microLogger')
visualLogger.VariableReferenceList = [['_', 'Variable:/:A', '1']]
visualLogger.VariableReferenceList = [['_', 'Variable:/:A', '-1']]
visualLogger.LogInterval = 1e-4

Filament = theSimulator.createEntity('FilamentProcess', 'Process:/:Filament')
Filament.OriginX = 0
Filament.OriginY = 0
Filament.OriginZ = 0
Filament.RotateX = 0
Filament.RotateY = 0
Filament.RotateZ = 0
Filament.SubunitRadius = 0.4e-8
Filament.Length = 0.6e-6
Filament.Periodic = 0
Filament.VariableReferenceList = [['_', 'Variable:/:A_Filament' ]]
Filament.VariableReferenceList = [['_', 'Variable:/:A_FilamentATP' ]]
Filament.VariableReferenceList = [['_', 'Variable:/:Subunit' , '-1']]
Filament.VariableReferenceList = [['_', 'Variable:/:SubunitM' , '-2']]
Filament.VariableReferenceList = [['_', 'Variable:/:SubunitP' , '-3']]

run(0.05)


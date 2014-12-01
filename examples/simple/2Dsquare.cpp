#include <libecs/libecs.hpp>
#include <libecs/Model.hpp>
#include <libecs/Entity.hpp>
#include <libecs/Variable.hpp>
#include <libecs/Process.hpp>
#include <libecs/SpatiocyteCommon.hpp>


libecs::Variable& createVariable(libecs::Model& model,
                                 const libecs::String& fullid,
                                 const double& value)
{
  libecs::Entity& entity(*model.createEntity("Variable",
                                             libecs::FullID(fullid)));
  entity.loadProperty("Value", libecs::Polymorph(value));
  return dynamic_cast<libecs::Variable&>(entity);
}

libecs::Process& createProcess(libecs::Model& model,
                               const libecs::String& class_name,
                               const libecs::String& fullid)
{
  return dynamic_cast<libecs::Process&>(
                  *model.createEntity(class_name, libecs::FullID(fullid)));
}

libecs::System& createSystem(libecs::Model& model,
                             const libecs::String& fullid,
                             const libecs::String& stepperID)
{

  libecs::Entity& entity(*model.createEntity("System",
                                             libecs::FullID(fullid)));
  entity.loadProperty("StepperID", libecs::Polymorph(stepperID));
  return dynamic_cast<libecs::System&>(entity);
}

int main()
{
  libecs::initialize();
  libecs::Model& 
    model(*(new libecs::Model(*libecs::createDefaultModuleMaker())));
  model.setup();
  libecs::Stepper& stepper(*model.createStepper("SpatiocyteStepper", "SS"));
  stepper.setProperty("VoxelRadius", libecs::Polymorph(4.4e-9)); 
  stepper.setProperty("SearchVacant", libecs::Polymorph(libecs::Integer(0))); 
  model.getRootSystem()->setProperty("StepperID", libecs::Polymorph("SS"));
  createVariable(model, "Variable:/:GEOMETRY", CUBOID);
  createVariable(model, "Variable:/:LENGTHX", 1e-8);
  createVariable(model, "Variable:/:LENGTHY", 1e-6);
  createVariable(model, "Variable:/:LENGTHZ", 1e-6);
  createVariable(model, "Variable:/:XYPLANE", REMOVE_BOTH);
  createVariable(model, "Variable:/:XZPLANE", REMOVE_BOTH);
  createVariable(model, "Variable:/:YZPLANE", REMOVE_LOWER);
  createVariable(model, "Variable:/:VACANT", 0);

  createSystem(model, "System:/:Surface", "SS");
  createVariable(model, "Variable:/Surface:DIMENSION", 2);
  createVariable(model, "Variable:/Surface:VACANT", 0);
  createVariable(model, "Variable:/Surface:A", 500);
  createVariable(model, "Variable:/Surface:As", 0); 
  
  libecs::Process& vis(createProcess(model, "VisualizationLogProcess",
                                     "Process:/:logger"));
  vis.registerVariableReference("_", libecs::String("Variable:/Surface:A"), 0);
  vis.registerVariableReference("_", libecs::String("Variable:/Surface:As"), 0);
  vis.loadProperty("LogInterval", libecs::Polymorph(0.01));

  libecs::Process& pop(createProcess(model, "MoleculePopulateProcess",
                                     "Process:/:pop"));
  pop.registerVariableReference("_", libecs::String("Variable:/Surface:A"), 0);
  pop.registerVariableReference("_", libecs::String("Variable:/Surface:As"), 0);

  libecs::Process& r1(createProcess(model, "DiffusionInfluencedReactionProcess",
                                    "Process:/:reaction1"));
  r1.registerVariableReference("_", libecs::String("Variable:/Surface:A"), -1);
  r1.registerVariableReference("_", libecs::String("Variable:/Surface:A"), -1);
  r1.registerVariableReference("_", libecs::String("Variable:/Surface:As"), 1);
  r1.registerVariableReference("_", libecs::String("Variable:/Surface:As"), 1);
  r1.loadProperty("p", libecs::Polymorph(0.0001));

  libecs::Process& r2(createProcess(model, "DiffusionInfluencedReactionProcess",
                                    "Process:/:reaction2"));
  r2.registerVariableReference("_", libecs::String("Variable:/Surface:A"), -1);
  r2.registerVariableReference("_", libecs::String("Variable:/Surface:As"), -1);
  r2.registerVariableReference("_", libecs::String("Variable:/Surface:As"), 1);
  r2.registerVariableReference("_", libecs::String("Variable:/Surface:As"), 1);
  r2.loadProperty("p", libecs::Polymorph(1.0));

  libecs::Process& r3(createProcess(model, "SpatiocyteNextReactionProcess",
                                    "Process:/:reaction3"));
  r3.registerVariableReference("_", libecs::String("Variable:/Surface:As"), -1);
  r3.registerVariableReference("_", libecs::String("Variable:/Surface:A"), 1);
  r3.loadProperty("Deoligomerize", libecs::Polymorph(libecs::Integer(6)));
  r3.loadProperty("Rates", libecs::Polymorph(
                  boost::make_tuple(243.0, 81.0, 27.0, 9.0, 3.0, 1.0)));

  libecs::Process& dif(createProcess(model, "DiffusionProcess",
                                     "Process:/:diffuseA"));
  dif.registerVariableReference("_", libecs::String("Variable:/Surface:A"), 0);
  dif.loadProperty("D", libecs::Polymorph(1e-13));

  model.initialize();
  double t(model.getCurrentTime());
  while(t < 100)
    { 
      model.step();
      t = model.getCurrentTime();
    }

  delete &model;
  libecs::finalize(); 
}

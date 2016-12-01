#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "Configuration.h"
#include "Simulation.h"
#include "TurbulentSimulation.h"
#include "parallelManagers/PetscParallelConfiguration.h"
#include "MeshsizeFactory.h"
#include <iomanip>

int main (int argc, char *argv[]) {

    // Parallelization related. Initialize and identify
    // ---------------------------------------------------
    int rank;   // This processor's identifier
    int nproc;  // Number of processors in the group
    PetscInitialize(&argc, &argv, "petsc_commandline_arg", PETSC_NULL);
    MPI_Comm_size(PETSC_COMM_WORLD, &nproc);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    std::cout << "Rank: " << rank << ", Nproc: " << nproc << std::endl;
    //----------------------------------------------------


    // read configuration and store information in parameters object
    Configuration configuration(argv[1]);
    Parameters parameters;
    configuration.loadParameters(parameters);
    PetscParallelConfiguration parallelConfiguration(parameters);
    MeshsizeFactory::getInstance().initMeshsize(parameters);
    FlowField *flowField = NULL;
    Simulation *simulation = NULL;

    #ifdef DEBUG
    std::cout << "Processor " << parameters.parallel.rank << " with index ";
    std::cout << parameters.parallel.indices[0] << ",";
    std::cout << parameters.parallel.indices[1] << ",";
    std::cout << parameters.parallel.indices[2];
    std::cout <<    " is computing the size of its subdomain and obtains ";
    std::cout << parameters.parallel.localSize[0] << ", ";
    std::cout << parameters.parallel.localSize[1] << " and ";
    std::cout << parameters.parallel.localSize[2];
    std::cout << ". Left neighbour: " << parameters.parallel.leftNb;
    std::cout << ", right neighbour: " << parameters.parallel.rightNb;
    std::cout << std::endl;
    std::cout << "Min. meshsizes: " << parameters.meshsize->getDxMin() << ", " << parameters.meshsize->getDyMin() << ", " << parameters.meshsize->getDzMin() << std::endl;
    #endif

    // initialise simulation
    if (parameters.simulation.type=="turbulence"){
      // TODO WS2: initialise turbulent flow field and turbulent simulation object
      if(rank==0){ std::cout << "Start RANS simulation in " << parameters.geometry.dim << "D" << std::endl; }
      flowField = new TurbFlowField(parameters);
      if(flowField == NULL){ handleError(1, "flowField==NULL!"); }
      simulation = new TurbulentSimulation(parameters, *(dynamic_cast<TurbFlowField *>(flowField)));
    //   handleError(1,"Turbulence currently not supported yet!");
    } else if (parameters.simulation.type=="dns"){
      if(rank==0){ std::cout << "Start DNS simulation in " << parameters.geometry.dim << "D" << std::endl; }
      flowField = new FlowField(parameters);
      if(flowField == NULL){ handleError(1, "flowField==NULL!"); }
      simulation = new Simulation(parameters,*flowField);
    } else {
      handleError(1, "Unknown simulation type! Currently supported: dns, turbulence");
    }
    // call initialization of simulation (initialize flow field)
    if(simulation == NULL){ handleError(1, "simulation==NULL!"); }
    simulation->initializeFlowField();
    //flowField->getFlags().show();

    FLOAT time = 0.0;
    FLOAT timeStdOut=parameters.stdOut.interval;
    int timeSteps = 0;

    // Checks when the VTK file is to be outputted
    FLOAT timeVTKOut = parameters.vtk.interval;
    // VTK: plot initial state
    simulation->plotVTK(timeSteps);
    std::cout << "Plotting VTK file at time: " << time << std::endl << std::endl;

    // time loop
    while (time < parameters.simulation.finalTime){

      simulation->solveTimestep();

      time += parameters.timestep.dt;

      // std-out: terminal info
      if ( (rank==0) && (timeStdOut <= time) ){
          std::cout << "Current time: " << time << "\ttimestep: " <<
                        parameters.timestep.dt << std::endl;
          timeStdOut += parameters.stdOut.interval;
      }

      // VTK: trigger VTK output
      // NOTE: Only sequential code.
      // Hence the check for rank == 0. In case parallel code is run, only rank 0 will write.
      if ( (rank == 0) && (timeVTKOut <= time) ){
          simulation->plotVTK(timeSteps);
          std::cout << "Plotting VTK file at time: " << time << std::endl << std::endl;
          timeVTKOut += parameters.vtk.interval;
      }

      timeSteps++;
    }

    // VTK: plot final output
    // if ( rank == 0 ) {
    //     simulation->plotVTK(timeSteps);
    //     std::cout << "Plotting VTK file at time: " << time << std::endl << std::endl;
    // }

    delete simulation; simulation=NULL;
    delete flowField;  flowField= NULL;

    PetscFinalize();
}

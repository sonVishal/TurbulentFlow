#ifndef _TURBULENT_SIMULATION_H_
#define _TURBULENT_SIMULATION_H_

#define BINARY_TURB_VTK

#include "Simulation.h"
#include "TurbFlowField.h"
#include "stencils/MixingLengthStencil.h"
#ifdef BINARY_TURB_VTK
	#include "stencils/VTKBinaryTurbStencil.h"
#else
	#include "stencils/TurbVTKStencil.h"
#endif

#include "stencils/TurbFGHStencil.h"
#include "stencils/TurbLPmodelStencil.h"
#include "stencils/MaxTurbViscosityStencil.h"

class TurbulentSimulation : public Simulation {
private:
    // TODO: Add other members such as TurbFGH when done
    TurbFlowField &_turbFlowField;

	#ifdef BINARY_TURB_VTK
    	VTKBinaryTurbStencil _turbVTKStencil;
	#else
		TurbVTKStencil _turbVTKStencil;
	#endif

    FieldIterator<TurbFlowField> _turbVTKIterator;

    TurbFGHStencil _turbFGHStencil;
    FieldIterator<TurbFlowField> _turbFGHIterator;

    TurbLPmodel _turbLPmodelStencil;
    FieldIterator<TurbFlowField> _turbLPmodelIterator;

    MaxTurbViscosityStencil _maxTurbViscosityStencil;
    FieldIterator<TurbFlowField> _maxTurbViscosityFieldIterator;
    //GlobalBoundaryIterator<TurbFlowField> _minTurbViscosityBoundaryIterator;

    GlobalBoundaryIterator<TurbFlowField> _wallTurbViscosityIterator;
public:
    TurbulentSimulation(Parameters &parameters, TurbFlowField &turbFlowField) :
        Simulation(parameters, turbFlowField),
        _turbFlowField(turbFlowField),
		#ifdef BINARY_TURB_VTK
				_turbVTKStencil(parameters),
		#else
				_turbVTKStencil(parameters),
		#endif
        _turbVTKIterator(turbFlowField,parameters,_turbVTKStencil),
        _turbFGHStencil(parameters),
        _turbFGHIterator(turbFlowField,parameters,_turbFGHStencil),
        _turbLPmodelStencil(parameters),
        _turbLPmodelIterator(turbFlowField,parameters,_turbLPmodelStencil),
        _maxTurbViscosityStencil(parameters),
        _maxTurbViscosityFieldIterator(turbFlowField,parameters,_maxTurbViscosityStencil,2,-1),
        //_minTurbViscosityBoundaryIterator(turbFlowField,parameters,_minTurbViscosityStencil),
        _wallTurbViscosityIterator(_globalBoundaryFactory.getGlobalBoundaryTurbViscosityIterator(_turbFlowField)){}
    ~TurbulentSimulation() {}
    void initializeFlowField() {
        Simulation::initializeFlowField();

        // Initialize the mixing length
        MixingLengthStencil mixingLengthStencil(_parameters);
        FieldIterator<TurbFlowField> mixingLengthIterator(_turbFlowField,_parameters,mixingLengthStencil,1,0);
        mixingLengthIterator.iterate();
    }
    void solveTimestep() {
        // compute viscosity
        _turbLPmodelIterator.iterate();
        // communicate turbulent viscosity values
        _petscParallelManager.communicateViscosity();
		// update wall values
        _wallTurbViscosityIterator.iterate();

        // Do the viscosity first so that the min viscosity is not 0 while setting up dt
        // This is equivalent as changing order in a cyclic manner does not change the algorithm
        // determine and set max. timestep which is allowed in this simulation
        setTimeStep();
        // compute fgh
        _turbFGHIterator.iterate();
        // set global boundary values
        _wallFGHIterator.iterate();
        // compute the right hand side
        _rhsIterator.iterate();
        // solve for pressure
        _solver.solve();
        // communicate pressure values
        _petscParallelManager.communicatePressure();
        // compute velocity
        _velocityIterator.iterate();
    	// set obstacle boundaries
    	_obstacleIterator.iterate();
        // communicate velocity values
    	_petscParallelManager.communicateVelocity();
        // Iterate for velocities on the boundary
        _wallVelocityIterator.iterate();
    }
    void plotVTK(int timeStep) {

	#ifdef BINARY_TURB_VTK
    	_turbVTKIterator.iterate();
    	_turbVTKStencil.write(_turbFlowField, timeStep);
	#else
		if (_turbVTKStencil.openFile(timeStep)) {
			_turbVTKIterator.iterate();
			_turbVTKStencil.write(_turbFlowField, timeStep);
		} else {
			std::cout << "ERROR: Plotting VTK file at time: " << timeStep << " FAILED!" << std::endl;
			handleError(1,"Could not open the file for writing!");
		}
	#endif
    }
private:
    void setTimeStep() {
        FLOAT localMin, globalMin;
        assertion(_parameters.geometry.dim == 2 || _parameters.geometry.dim == 3);
        FLOAT factor = 1.0/(_parameters.meshsize->getDxMin() * _parameters.meshsize->getDxMin()) +
                       1.0/(_parameters.meshsize->getDyMin() * _parameters.meshsize->getDyMin());

        // determine maximum velocity
        _maxUStencil.reset();
        _maxUFieldIterator.iterate();
        _maxUBoundaryIterator.iterate();

        // determine the maximum turbulent viscosity (we store nu + nu_t as turbulent viscosity)
        _maxTurbViscosityStencil.reset();
        _maxTurbViscosityFieldIterator.iterate();
        //_minTurbViscosityBoundaryIterator.iterate();

        assertion(_maxTurbViscosityStencil.getMaxValue() > 0.0);

        if (_parameters.geometry.dim == 3) {
          factor += 1.0/(_parameters.meshsize->getDzMin() * _parameters.meshsize->getDzMin());
          _parameters.timestep.dt = 1.0 / _maxUStencil.getMaxValues()[2];
        } else {
          _parameters.timestep.dt = 1.0 / _maxUStencil.getMaxValues()[0];
        }

        localMin = std::min(_parameters.timestep.dt,
                                          std::min(std::min(1.0/(2*factor*_maxTurbViscosityStencil.getMaxValue()),
                                          1.0 / _maxUStencil.getMaxValues()[0]),
                                          1.0 / _maxUStencil.getMaxValues()[1]));

        // Here, we select the type of operation before compiling. This allows to use the correct
        // data type for MPI. Not a concern for small simulations, but useful if using heterogeneous
        // machines.

        globalMin = MY_FLOAT_MAX;
        MPI_Allreduce(&localMin, &globalMin, 1, MY_MPI_FLOAT, MPI_MIN, PETSC_COMM_WORLD);

        _parameters.timestep.dt = globalMin;
        _parameters.timestep.dt *= _parameters.timestep.tau;
    }
};

#endif

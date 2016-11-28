#ifndef _TURBULENT_SIMULATION_H_
#define _TURBULENT_SIMULATION_H_
#include "Simulation.h"
#include "TurbFlowField.h"
#include "stencils/DistToWallStencil.h"
#include "stencils/TurbVTKStencil.h"
#include "stencils/TurbFGHStencil.h"
#include "stencils/TurbLPmodelStencil.h"
class TurbulentSimulation : public Simulation {
private:
    // TODO: Add other members such as TurbFGH when done
    TurbFlowField &_turbFlowField;
    TurbVTKStencil _turbVTKStencil;
    FieldIterator<TurbFlowField> _turbVTKIterator;

    TurbFGHStencil _turbFGHStencil;
    FieldIterator<TurbFlowField> _turbFGHIterator;

    TurbLPmodel _turbLPmodelStencil;
    FieldIterator<TurbFlowField> _turbLPmodelIterator;
public:
    TurbulentSimulation(Parameters &parameters, TurbFlowField &turbFlowField) :
        Simulation(parameters, turbFlowField),
        _turbFlowField(turbFlowField),
        _turbVTKStencil(parameters),
        _turbVTKIterator(turbFlowField,parameters,_turbVTKStencil),
        _turbFGHStencil(parameters),
        _turbFGHIterator(turbFlowField,parameters,_turbFGHStencil),
        _turbLPmodelStencil(parameters),
        _turbLPmodelIterator(turbFlowField,parameters,_turbLPmodelStencil) {}
    ~TurbulentSimulation() {}
    void initializeFlowField() {
        Simulation::initializeFlowField();
        // Initialize the distance to wall
        DistToWallStencil distToWallStencil(_parameters);
        FieldIterator<TurbFlowField> distToWallIterator(_turbFlowField,_parameters,distToWallStencil);
        distToWallIterator.iterate();
    }
    void solveTimestep() {
        // TODO
        handleError(1,"TODO");

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
        // TODO WS2: communicate pressure values
        // compute velocity
        _velocityIterator.iterate();
        // compute viscosity
        _turbLPmodelIterator.iterate();
    	// set obstacle boundaries
    	_obstacleIterator.iterate();
        // TODO WS2: communicate velocity values
        // Iterate for velocities on the boundary
        _wallVelocityIterator.iterate();
        // TODO: Wall viscosity iterator??
    }
    void plotVTK(int timeStep) {
        if (_turbVTKStencil.openFile(timeStep)) {
            _turbVTKIterator.iterate();
            _turbVTKStencil.write(_turbFlowField, timeStep);
        } else {
            std::cout << "ERROR: Plotting VTK file at time: " << timeStep << " FAILED!" << std::endl;
            handleError(1,"Could not open the file for writing!");
        }
    }
private:
    void setTimeStep() {
        // TODO
        handleError(1,"TODO");
    }
};

#endif

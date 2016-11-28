#ifndef _TURBULENT_SIMULATION_H_
#define _TURBULENT_SIMULATION_H_
#include "Simulation.h"
#include "TurbFlowField.h"
#include "stencils/DistToWallStencil.h"
#include "stencils/TurbVTKStencil.h"
class TurbulentSimulation : public Simulation {
private:
    // TODO: Add other members such as TurbFGH when done
    TurbFlowField &_turbFlowField;
    TurbVTKStencil _turbVTKStencil;
    FieldIterator<TurbFlowField> _turbVTKIterator;
public:
    TurbulentSimulation(Parameters &parameters, TurbFlowField &turbFlowField) :
        Simulation(parameters, turbFlowField),
        _turbFlowField(turbFlowField),
        _turbVTKStencil(parameters),
        _turbVTKIterator(turbFlowField,parameters,_turbVTKStencil) {}
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

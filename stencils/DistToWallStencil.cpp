#include "DistToWallStencil.h"

DistToWallStencil::DistToWallStencil ( const Parameters & parameters ) :
    FieldStencil<TurbFlowField> ( parameters ) {}


void DistToWallStencil::apply ( TurbFlowField & flowField,  int i, int j ){
    const int obstacle = flowField.getFlags().getValue(i, j);

    if ((obstacle & OBSTACLE_SELF) == 0){ // If this is a fluid cell
        if (_parameters.simulation.scenario == "cavity") {
            /* code */
        } else if (_parameters.simulation.scenario == "channel") {
            /* code */
        } else {
            handleError(1, "This block should not be reached.\nSomething went wrong with the scenario.\n");
        }
    }
}


void DistToWallStencil::apply ( TurbFlowField & flowField, int i, int j, int k ){
    const int obstacle = flowField.getFlags().getValue(i, j, k);

    if ((obstacle & OBSTACLE_SELF) == 0){ // If this is a fluid cell
        if (_parameters.simulation.scenario == "cavity") {
            /* code */
        } else if (_parameters.simulation.scenario == "channel") {
            /* code */
        } else {
            handleError(1, "This block should not be reached.\nSomething went wrong with the scenario.\n");
        }
    }
}

void DistToWallStencil::getCellCenter( int i, int j, FLOAT* cellCenter ) {
    FLOAT cellX   = _parameters.meshsize->getPosX(i,j);
    FLOAT cellY   = _parameters.meshsize->getPosY(i,j);
    FLOAT cellDx  = _parameters.meshsize->getDx(i,j);
    FLOAT cellDy  = _parameters.meshsize->getDy(i,j);
    cellCenter[0] = cellX + cellDx/2.0;
    cellCenter[1] = cellY + cellDy/2.0;
}

void DistToWallStencil::getCellCenter( int i, int j, int k, FLOAT* cellCenter ) {
    FLOAT cellX   = _parameters.meshsize->getPosX(i,j,k);
    FLOAT cellY   = _parameters.meshsize->getPosY(i,j,k);
    FLOAT cellZ   = _parameters.meshsize->getPosZ(i,j,k);
    FLOAT cellDx  = _parameters.meshsize->getDx(i,j,k);
    FLOAT cellDy  = _parameters.meshsize->getDy(i,j,k);
    FLOAT cellDz  = _parameters.meshsize->getDz(i,j,k);
    cellCenter[0] = cellX + cellDx/2.0;
    cellCenter[1] = cellY + cellDy/2.0;
    cellCenter[2] = cellZ + cellDz/2.0;
}

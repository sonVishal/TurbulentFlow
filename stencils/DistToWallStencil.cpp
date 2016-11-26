#include "DistToWallStencil.h"

DistToWallStencil::DistToWallStencil ( const Parameters & parameters ) :
    FieldStencil<TurbFlowField> ( parameters ) {}


void DistToWallStencil::apply ( TurbFlowField & flowField,  int i, int j ){
    const int obstacle = flowField.getFlags().getValue(i, j);
    FLOAT cellCenterPos[2];
    getCellCenter(i, j, cellCenterPos);
    if ((obstacle & OBSTACLE_SELF) == 0){ // If this is a fluid cell
        if (_parameters.simulation.scenario == "cavity") { // Cavity scenario means top is always not a wall
            // Min is obtained by checking left, right and bottom.
            flowField.getDistanceToWall().getScalar(i, j) = std::min(cellCenterPos[1], std::min(_parameters.geometry.lengthX-cellCenterPos[0], cellCenterPos[0]));
        } else if (_parameters.simulation.scenario == "channel") { // Channel means left and right are always not a wall
            // Min is obtained by checking top and bottom
            if (cellCenterPos[0] < _parameters.bfStep.xRatio*_parameters.geometry.lengthX) { // Step present and cell above step
                flowField.getDistanceToWall().getScalar(i, j) = std::min(_parameters.geometry.lengthY-cellCenterPos[1], _parameters.bfStep.yRatio*_parameters.geometry.lengthY);
            } else { // Step present and cell beyond step or step not present since xRatio = 0.0
                flowField.getDistanceToWall().getScalar(i, j) = std::min(cellCenterPos[1], _parameters.geometry.lengthY-cellCenterPos[1]);
            } // Not required to check in x direction since we do not have walls in x direction
        } else {
            handleError(1, "This block should not be reached.\nSomething went wrong with the scenario.\n");
        }
    }
}


void DistToWallStencil::apply ( TurbFlowField & flowField, int i, int j, int k ){
    const int obstacle = flowField.getFlags().getValue(i, j, k);
    FLOAT cellCenterPos[3];
    getCellCenter(i, j, k, cellCenterPos);
    if ((obstacle & OBSTACLE_SELF) == 0){ // If this is a fluid cell
        if (_parameters.simulation.scenario == "cavity") {// Cavity scenario means top is always not a wall
            // Min is obtained by checking left, right, front, back and bottom.
            FLOAT minX = std::min(cellCenterPos[0], _parameters.geometry.lengthX-cellCenterPos[0]);
            FLOAT minZ = std::min(cellCenterPos[2], _parameters.geometry.lengthZ-cellCenterPos[2]);
            flowField.getDistanceToWall().getScalar(i, j, k) = std::min(cellCenterPos[1], std::min(minX, minZ));
        } else if (_parameters.simulation.scenario == "channel") { // Channel means left and right are always not a wall
            FLOAT min1;
            // Min is obtained by checking top, bottom, front and back
            // Check for the step first
            if (cellCenterPos[0] < _parameters.bfStep.xRatio*_parameters.geometry.lengthX) { // Step present and cell above step
                min1 = std::min(_parameters.geometry.lengthY-cellCenterPos[1], _parameters.bfStep.yRatio*_parameters.geometry.lengthY);
            } else { // Step present and cell beyond step or step not present since xRatio = 0.0
                min1 = std::min(cellCenterPos[1], _parameters.geometry.lengthY-cellCenterPos[1]);
            } // Not required to check in x direction since we do not have walls in x direction
            // Once the step is checked we check for front and back
            flowField.getDistanceToWall().getScalar(i, j, k) = std::min(min1,std::min(cellCenterPos[2],_parameters.geometry.lengthZ-cellCenterPos[2]));
        } else {
            handleError(1, "This block should not be reached.\nSomething went wrong with the scenario.\n");
        }
    }
}

// Get the coordinates of the current cell center
void DistToWallStencil::getCellCenter( int i, int j, FLOAT* cellCenter ) {
    FLOAT cellX   = _parameters.meshsize->getPosX(i,j);
    FLOAT cellY   = _parameters.meshsize->getPosY(i,j);
    FLOAT cellDx  = _parameters.meshsize->getDx(i,j);
    FLOAT cellDy  = _parameters.meshsize->getDy(i,j);
    cellCenter[0] = cellX + cellDx/2.0;
    cellCenter[1] = cellY + cellDy/2.0;
}

// Get the coordinates of the current cell center
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

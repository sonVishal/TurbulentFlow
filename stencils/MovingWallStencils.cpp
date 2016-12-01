#include "MovingWallStencils.h"

MovingWallVelocityStencil::MovingWallVelocityStencil ( const Parameters & parameters ) :
    BoundaryStencil<FlowField> ( parameters ) {}

// 2D stencils

void MovingWallVelocityStencil::applyLeftWall ( FlowField & flowField, int i, int j ){

    flowField.getVelocity().getVector(i, j)[0] = _parameters.walls.vectorLeft[0];
    flowField.getVelocity().getVector(i, j)[1] = 2 * _parameters.walls.vectorLeft[1]
                                               - flowField.getVelocity().getVector(i+1, j)[1];
}


void MovingWallVelocityStencil::applyRightWall ( FlowField & flowField, int i, int j ){

    flowField.getVelocity().getVector(i-1, j)[0] = _parameters.walls.vectorRight[0];
    flowField.getVelocity().getVector(i, j)[1] =
        2 * _parameters.walls.vectorRight[1] -
        flowField.getVelocity().getVector(i-1, j)[1];
}


void MovingWallVelocityStencil::applyBottomWall ( FlowField & flowField, int i, int j ){

    flowField.getVelocity().getVector(i, j)[0] = 2 * _parameters.walls.vectorBottom[0]
                                               - flowField.getVelocity().getVector(i, j+1)[0];
    flowField.getVelocity().getVector(i, j)[1] = _parameters.walls.vectorBottom[1];
}


void MovingWallVelocityStencil::applyTopWall ( FlowField & flowField, int i, int j ){

    flowField.getVelocity().getVector(i, j)[0] =
        2 * _parameters.walls.vectorTop[0] -
        flowField.getVelocity().getVector(i, j-1)[0];
    flowField.getVelocity().getVector(i, j-1)[1] = _parameters.walls.vectorTop[1];
}


// 3D stencils

void MovingWallVelocityStencil::applyLeftWall ( FlowField & flowField, int i, int j, int k ){

    flowField.getVelocity().getVector(i, j, k)[0] = _parameters.walls.vectorLeft[0];
    flowField.getVelocity().getVector(i, j, k)[1] = 2 * _parameters.walls.vectorLeft[1]
                                               - flowField.getVelocity().getVector(i+1, j, k)[1];
    flowField.getVelocity().getVector(i, j, k)[2] = 2 * _parameters.walls.vectorLeft[2]
                                               - flowField.getVelocity().getVector(i+1, j, k)[2];
}


void MovingWallVelocityStencil::applyRightWall ( FlowField & flowField, int i, int j , int k ){

    flowField.getVelocity().getVector(i-1, j, k)[0] = _parameters.walls.vectorRight[0];
    flowField.getVelocity().getVector(i, j, k)[1] =
        2 * _parameters.walls.vectorRight[1] -
        flowField.getVelocity().getVector(i-1, j, k)[1];
    flowField.getVelocity().getVector(i, j, k)[2] =
        2 * _parameters.walls.vectorRight[2] -
        flowField.getVelocity().getVector(i-1, j, k)[2];
}


void MovingWallVelocityStencil::applyBottomWall ( FlowField & flowField, int i, int j, int k ){

    flowField.getVelocity().getVector(i, j, k)[0] = 2 * _parameters.walls.vectorBottom[0]
                                               - flowField.getVelocity().getVector(i, j+1, k)[0];
    flowField.getVelocity().getVector(i, j, k)[1] = _parameters.walls.vectorBottom[1];
    flowField.getVelocity().getVector(i, j, k)[2] = 2 * _parameters.walls.vectorBottom[2]
                                               - flowField.getVelocity().getVector(i, j+1, k)[2];
}


void MovingWallVelocityStencil::applyTopWall ( FlowField & flowField, int i, int j, int k ){

    flowField.getVelocity().getVector(i, j, k)[0] =
        2 * _parameters.walls.vectorTop[0] -
        flowField.getVelocity().getVector(i, j-1, k)[0];
    flowField.getVelocity().getVector(i, j-1, k)[1] = _parameters.walls.vectorTop[1];
    flowField.getVelocity().getVector(i, j, k)[2] =
        2 * _parameters.walls.vectorTop[2] -
        flowField.getVelocity().getVector(i, j-1, k)[2];
}


void MovingWallVelocityStencil::applyFrontWall ( FlowField & flowField, int i, int j, int k ){

    flowField.getVelocity().getVector(i, j, k)[0] = 2 * _parameters.walls.vectorFront[0]
        - flowField.getVelocity().getVector(i, j, k+1)[0];
    flowField.getVelocity().getVector(i, j, k)[1] = 2 * _parameters.walls.vectorFront[1]
        - flowField.getVelocity().getVector(i, j, k+1)[1];
    flowField.getVelocity().getVector(i, j, k)[2] = _parameters.walls.vectorFront[2];
}


void MovingWallVelocityStencil::applyBackWall ( FlowField & flowField, int i, int j, int k ){

    flowField.getVelocity().getVector(i, j, k)[0] = 2 * _parameters.walls.vectorBack[0]
        - flowField.getVelocity().getVector(i, j, k-1)[0];
    flowField.getVelocity().getVector(i, j, k)[1] = 2 * _parameters.walls.vectorBack[1]
        - flowField.getVelocity().getVector(i, j, k-1)[1];
    flowField.getVelocity().getVector(i, j, k-1)[2] = _parameters.walls.vectorBack[2];
}


MovingWallFGHStencil::MovingWallFGHStencil ( const Parameters & parameters ) :
    BoundaryStencil<FlowField> ( parameters ) {}

// 2D stencils

void MovingWallFGHStencil::applyLeftWall ( FlowField & flowField, int i, int j ){
    flowField.getFGH().getVector(i, j)[0] = _parameters.walls.vectorLeft[0];
}


void MovingWallFGHStencil::applyRightWall ( FlowField & flowField, int i, int j ){
    flowField.getFGH().getVector(i-1, j)[0] = _parameters.walls.vectorRight[0];
}


void MovingWallFGHStencil::applyBottomWall ( FlowField & flowField, int i, int j ){
    flowField.getFGH().getVector(i, j)[1] = _parameters.walls.vectorBottom[1];
}


void MovingWallFGHStencil::applyTopWall ( FlowField & flowField, int i, int j ){
    flowField.getFGH().getVector(i, j-1)[1] = _parameters.walls.vectorTop[1];
}


// 3D stencils

void MovingWallFGHStencil::applyLeftWall ( FlowField & flowField, int i, int j, int k ){
    flowField.getFGH().getVector(i, j, k)[0] = _parameters.walls.vectorLeft[0];
}


void MovingWallFGHStencil::applyRightWall ( FlowField & flowField, int i, int j , int k ){
    flowField.getFGH().getVector(i-1, j, k)[0] = _parameters.walls.vectorRight[0];
}


void MovingWallFGHStencil::applyBottomWall ( FlowField & flowField, int i, int j, int k ){
    flowField.getFGH().getVector(i, j, k)[1] = _parameters.walls.vectorBottom[1];
}


void MovingWallFGHStencil::applyTopWall ( FlowField & flowField, int i, int j, int k ){
    flowField.getFGH().getVector(i, j-1, k)[1] = _parameters.walls.vectorTop[1];
}


void MovingWallFGHStencil::applyFrontWall ( FlowField & flowField, int i, int j, int k ){
    flowField.getFGH().getVector(i, j, k)[2] = _parameters.walls.vectorFront[2];
}


void MovingWallFGHStencil::applyBackWall ( FlowField & flowField, int i, int j, int k ){
    flowField.getFGH().getVector(i, j, k-1)[2] = _parameters.walls.vectorBack[2];
}

MovingWallTurbViscosityStencil::MovingWallTurbViscosityStencil ( const Parameters & parameters ) :
    BoundaryStencil<TurbFlowField> ( parameters ) {}

// 2D stencils

void MovingWallTurbViscosityStencil::applyLeftWall ( TurbFlowField & flowField, int i, int j ){

    flowField.getTurbViscosity().getScalar(i, j) = -flowField.getTurbViscosity().getScalar(i+1, j);
}


void MovingWallTurbViscosityStencil::applyRightWall ( TurbFlowField & flowField, int i, int j ){

    flowField.getTurbViscosity().getScalar(i, j) = -flowField.getTurbViscosity().getScalar(i-1, j);
}


void MovingWallTurbViscosityStencil::applyBottomWall ( TurbFlowField & flowField, int i, int j ){

    flowField.getTurbViscosity().getScalar(i, j) = flowField.getTurbViscosity().getScalar(i, j+1);
}


void MovingWallTurbViscosityStencil::applyTopWall ( TurbFlowField & flowField, int i, int j ){

    flowField.getTurbViscosity().getScalar(i, j) = flowField.getTurbViscosity().getScalar(i, j-1);
}


// 3D stencils

void MovingWallTurbViscosityStencil::applyLeftWall ( TurbFlowField & flowField, int i, int j, int k ){

    flowField.getTurbViscosity().getScalar(i, j, k) = -flowField.getTurbViscosity().getScalar(i+1, j, k);
}


void MovingWallTurbViscosityStencil::applyRightWall ( TurbFlowField & flowField, int i, int j , int k ){

    flowField.getTurbViscosity().getScalar(i, j, k) = -flowField.getTurbViscosity().getScalar(i-1, j, k);
}


void MovingWallTurbViscosityStencil::applyBottomWall ( TurbFlowField & flowField, int i, int j, int k ){

    flowField.getTurbViscosity().getScalar(i, j, k) = -flowField.getTurbViscosity().getScalar(i, j+1, k);
}


void MovingWallTurbViscosityStencil::applyTopWall ( TurbFlowField & flowField, int i, int j, int k ){

    flowField.getTurbViscosity().getScalar(i, j, k) = -flowField.getTurbViscosity().getScalar(i, j-1, k);
}


void MovingWallTurbViscosityStencil::applyFrontWall ( TurbFlowField & flowField, int i, int j, int k ){

    flowField.getTurbViscosity().getScalar(i, j, k) = -flowField.getTurbViscosity().getScalar(i, j, k+1);
}


void MovingWallTurbViscosityStencil::applyBackWall ( TurbFlowField & flowField, int i, int j, int k ){

    flowField.getTurbViscosity().getScalar(i, j, k) = -flowField.getTurbViscosity().getScalar(i, j, k-1);
}

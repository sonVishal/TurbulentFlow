#include "TurbFlowField.h"

TurbFlowField::TurbFlowField ( int Nx, int Ny ) :
    FlowField( Nx, Ny ),
    _turbViscosity ( ScalarField ( Nx + 3, Ny + 3 ) ),
    _distToWall ( ScalarField ( Nx + 3, Ny + 3 ) )
{ }


TurbFlowField::TurbFlowField ( int Nx, int Ny, int Nz ) :
    FlowField( Nx, Ny, Nz ),
    _turbViscosity ( ScalarField ( Nx + 3, Ny + 3, Nz + 3 ) ),
    _distToWall ( ScalarField ( Nx + 3, Ny + 3, Nz + 3 ) )
{ }


TurbFlowField::TurbFlowField (const Parameters & parameters) :
    FlowField( parameters ),
    // The size member variable is initialized by the parent class FlowField
    _turbViscosity(parameters.geometry.dim==2?ScalarField(_size_x + 3, _size_y + 3):
            ScalarField(_size_x + 3, _size_y + 3, _size_z + 3)),
    _distToWall(parameters.geometry.dim==2?ScalarField(_size_x + 3, _size_y + 3):
            ScalarField(_size_x + 3, _size_y + 3, _size_z + 3))
{ }

TurbFlowField::~TurbFlowField () {

}

ScalarField & TurbFlowField::getTurbViscosity () {
    return _turbViscosity;
}

ScalarField & TurbFlowField::getDistToWall() {
    return _distToWall;
}

void TurbFlowField::getPressureVelocityAndTurbVisc(FLOAT &pressure, FLOAT &turbViscosity, FLOAT* const velocity,  int i, int j){
    FLOAT * v_here = getVelocity().getVector(i, j);
    FLOAT * v_left = getVelocity().getVector(i-1, j);
    FLOAT * v_down = getVelocity().getVector(i, j-1);

    velocity[0] = ( v_here[0] + v_left[0] ) / 2;
    velocity[1] = ( v_here[1] + v_down[1] ) / 2;

    pressure = getPressure().getScalar(i,j);
    turbViscosity = getTurbViscosity().getScalar(i,j);
}

void TurbFlowField::getPressureVelocityAndTurbVisc(FLOAT &pressure, FLOAT &turbViscosity, FLOAT* const velocity, int i, int j, int k){
    FLOAT * v_here = getVelocity().getVector(i, j, k);
    FLOAT * v_left = getVelocity().getVector(i-1, j, k);
    FLOAT * v_down = getVelocity().getVector(i, j-1, k);
    FLOAT * v_back = getVelocity().getVector(i, j, k-1);

    velocity[0] = ( v_here[0] + v_left[0] ) / 2;
    velocity[1] = ( v_here[1] + v_down[1] ) / 2;
    velocity[2] = ( v_here[2] + v_back[2] ) / 2;

    pressure = getPressure().getScalar(i,j,k);
    turbViscosity = getTurbViscosity().getScalar(i,j,k);
}

void TurbFlowField::getDistanceToWall(FLOAT &distToWall, int i, int j) {
    distToWall = getDistToWall().getScalar(i,j);
}

void TurbFlowField::getDistanceToWall(FLOAT &distToWall, int i, int j, int k) {
    distToWall = getDistToWall().getScalar(i,j,k);
}

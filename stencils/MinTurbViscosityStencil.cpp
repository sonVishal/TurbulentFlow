#include "MinTurbViscosityStencil.h"
#include <algorithm>
#include <math.h>
#include <limits>


MinTurbViscosityStencil::MinTurbViscosityStencil (const Parameters & parameters) :
    FieldStencil<TurbFlowField> (parameters), BoundaryStencil<TurbFlowField> (parameters) {
    reset();
}

void MinTurbViscosityStencil::apply (TurbFlowField & flowField, int i, int j){
    cellMinValue(flowField, i, j);
}

void MinTurbViscosityStencil::apply (TurbFlowField & flowField, int i, int j, int k){
    cellMinValue(flowField, i, j, k);
}

void MinTurbViscosityStencil::applyLeftWall   ( TurbFlowField & flowField, int i, int j ){
    cellMinValue(flowField, i, j);
}

void MinTurbViscosityStencil::applyRightWall  ( TurbFlowField & flowField, int i, int j ){
    cellMinValue(flowField, i, j);
}

void MinTurbViscosityStencil::applyBottomWall ( TurbFlowField & flowField, int i, int j ){
    cellMinValue(flowField, i, j);
}

void MinTurbViscosityStencil::applyTopWall    ( TurbFlowField & flowField, int i, int j ){
    cellMinValue(flowField, i, j);
}

void MinTurbViscosityStencil::applyLeftWall   ( TurbFlowField & flowField, int i, int j, int k ){
    cellMinValue(flowField, i, j, k);
}

void MinTurbViscosityStencil::applyRightWall  ( TurbFlowField & flowField, int i, int j, int k ){
    cellMinValue(flowField, i, j, k);
}

void MinTurbViscosityStencil::applyBottomWall ( TurbFlowField & flowField, int i, int j, int k ){
    cellMinValue(flowField, i, j, k);
}

void MinTurbViscosityStencil::applyTopWall    ( TurbFlowField & flowField, int i, int j, int k ){
    cellMinValue(flowField, i, j, k);
}

void MinTurbViscosityStencil::applyFrontWall  ( TurbFlowField & flowField, int i, int j, int k ){
    cellMinValue(flowField, i, j, k);
}

void MinTurbViscosityStencil::applyBackWall   ( TurbFlowField & flowField, int i, int j, int k ){
    cellMinValue(flowField, i, j, k);
}


void MinTurbViscosityStencil::cellMinValue(TurbFlowField & flowField, int i, int j){
    FLOAT turbViscosity = flowField.getTurbViscosity().getScalar(i, j);
    if (fabs(turbViscosity) < _minValue) {
        _minValue = turbViscosity;
    }
}

void MinTurbViscosityStencil::cellMinValue(TurbFlowField & flowField, int i, int j, int k){
    FLOAT turbViscosity = flowField.getTurbViscosity().getScalar(i, j, k);
    if (fabs(turbViscosity) < _minValue) {
        _minValue = turbViscosity;
    }
}

void MinTurbViscosityStencil::reset () {
    _minValue = std::numeric_limits<FLOAT>::max();
}

const FLOAT MinTurbViscosityStencil::getMinValue() const{
    return _minValue;
}

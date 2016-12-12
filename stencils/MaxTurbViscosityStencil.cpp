#include "MaxTurbViscosityStencil.h"
#include <algorithm>
#include <math.h>
#include <limits>


MaxTurbViscosityStencil::MaxTurbViscosityStencil (const Parameters & parameters) :
    FieldStencil<TurbFlowField> (parameters), BoundaryStencil<TurbFlowField> (parameters) {
    reset();
}

void MaxTurbViscosityStencil::apply (TurbFlowField & flowField, int i, int j){
    cellMaxValue(flowField, i, j);
}

void MaxTurbViscosityStencil::apply (TurbFlowField & flowField, int i, int j, int k){
    cellMaxValue(flowField, i, j, k);
}

void MaxTurbViscosityStencil::applyLeftWall   ( TurbFlowField & flowField, int i, int j ){
    cellMaxValue(flowField, i, j);
}

void MaxTurbViscosityStencil::applyRightWall  ( TurbFlowField & flowField, int i, int j ){
    cellMaxValue(flowField, i, j);
}

void MaxTurbViscosityStencil::applyBottomWall ( TurbFlowField & flowField, int i, int j ){
    cellMaxValue(flowField, i, j);
}

void MaxTurbViscosityStencil::applyTopWall    ( TurbFlowField & flowField, int i, int j ){
    cellMaxValue(flowField, i, j);
}

void MaxTurbViscosityStencil::applyLeftWall   ( TurbFlowField & flowField, int i, int j, int k ){
    cellMaxValue(flowField, i, j, k);
}

void MaxTurbViscosityStencil::applyRightWall  ( TurbFlowField & flowField, int i, int j, int k ){
    cellMaxValue(flowField, i, j, k);
}

void MaxTurbViscosityStencil::applyBottomWall ( TurbFlowField & flowField, int i, int j, int k ){
    cellMaxValue(flowField, i, j, k);
}

void MaxTurbViscosityStencil::applyTopWall    ( TurbFlowField & flowField, int i, int j, int k ){
    cellMaxValue(flowField, i, j, k);
}

void MaxTurbViscosityStencil::applyFrontWall  ( TurbFlowField & flowField, int i, int j, int k ){
    cellMaxValue(flowField, i, j, k);
}

void MaxTurbViscosityStencil::applyBackWall   ( TurbFlowField & flowField, int i, int j, int k ){
    cellMaxValue(flowField, i, j, k);
}


void MaxTurbViscosityStencil::cellMaxValue(TurbFlowField & flowField, int i, int j){
    FLOAT turbViscosity = flowField.getTurbViscosity().getScalar(i, j);
    _maxValue = std::max(_maxValue,turbViscosity);
}

void MaxTurbViscosityStencil::cellMaxValue(TurbFlowField & flowField, int i, int j, int k){
    FLOAT turbViscosity = flowField.getTurbViscosity().getScalar(i, j, k);
    _maxValue = std::max(_maxValue,turbViscosity);
}

void MaxTurbViscosityStencil::reset () {
    _maxValue = 0.0;
}

FLOAT MaxTurbViscosityStencil::getMaxValue() const{
    return _maxValue;
}

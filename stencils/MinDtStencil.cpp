#include "MinDtStencil.h"
#include <algorithm>
#include <math.h>
#include <limits>


MinDtStencil::MinDtStencil (const Parameters & parameters) :
    FieldStencil<TurbFlowField> (parameters) {
    reset();
}

void MinDtStencil::apply (TurbFlowField & flowField, int i, int j){
    cellMinValue(flowField, i, j);
}

void MinDtStencil::apply (TurbFlowField & flowField, int i, int j, int k){
    cellMinValue(flowField, i, j, k);
}

void MinDtStencil::cellMinValue(TurbFlowField & flowField, int i, int j){
    FLOAT turbViscosity = flowField.getTurbViscosity().getScalar(i,j);
    FLOAT dx = _parameters.meshsize->getDx(i,j);
    FLOAT dy = _parameters.meshsize->getDy(i,j);
    FLOAT funcVal = 1/(2*turbViscosity*(1/dx/dx + 1/dy/dy));
    _minValue = std::min(_minValue,funcVal);
}

void MinDtStencil::cellMinValue(TurbFlowField & flowField, int i, int j, int k){
    FLOAT turbViscosity = flowField.getTurbViscosity().getScalar(i,j,k);
    FLOAT dx = _parameters.meshsize->getDx(i,j,k);
    FLOAT dy = _parameters.meshsize->getDy(i,j,k);
    FLOAT dz = _parameters.meshsize->getDz(i,j,k);
    FLOAT funcVal = 1/(2*turbViscosity*(1/dx/dx + 1/dy/dy + 1/dz/dz));
    _minValue = std::min(_minValue,funcVal);
}

void MinDtStencil::reset () {
    _minValue = MY_FLOAT_MAX;
}

FLOAT MinDtStencil::getMinValue() const{
    return _minValue;
}

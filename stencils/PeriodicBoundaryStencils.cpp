#include "PeriodicBoundaryStencils.h"

PeriodicBoundaryVelocityStencil::PeriodicBoundaryVelocityStencil(const Parameters & parameters):
    BoundaryStencil<FlowField>(parameters) {}

// 2D problem

void PeriodicBoundaryVelocityStencil::applyLeftWall(FlowField & flowField, int i, int j){
    flowField.getVelocity().getVector(0, j)[0] =
        flowField.getVelocity().getVector(flowField.getNx(), j)[0];
    flowField.getVelocity().getVector(1, j)[1] =
        flowField.getVelocity().getVector(flowField.getNx()+1, j)[1];
}

void PeriodicBoundaryVelocityStencil::applyRightWall(FlowField & flowField, int i, int j){
    flowField.getVelocity().getVector(flowField.getNx()+2, j)[0] =
        flowField.getVelocity().getVector(2, j)[0];
    flowField.getVelocity().getVector(flowField.getNx()+2, j)[1] =
        flowField.getVelocity().getVector(2, j)[1];
}

void PeriodicBoundaryVelocityStencil::applyBottomWall(FlowField & flowField, int i, int j){
    flowField.getVelocity().getVector(i, 0)[1] =
        flowField.getVelocity().getVector(i, flowField.getNy())[1];
    flowField.getVelocity().getVector(i, 1)[0] =
        flowField.getVelocity().getVector(i, flowField.getNy()+1)[0];
}

void PeriodicBoundaryVelocityStencil::applyTopWall(FlowField & flowField, int i, int j){
    flowField.getVelocity().getVector(i, flowField.getNy()+2)[1] =
        flowField.getVelocity().getVector(i, 2)[1];
    flowField.getVelocity().getVector(i, flowField.getNy()+2)[0] =
        flowField.getVelocity().getVector(i, 2)[0];
}

// 3D Problem

void PeriodicBoundaryVelocityStencil::applyLeftWall(FlowField & flowField, int i, int j, int k){
    flowField.getVelocity().getVector(0, j, k)[0] =
        flowField.getVelocity().getVector(flowField.getNx(), j, k)[0];
    flowField.getVelocity().getVector(1, j, k)[1] =
        flowField.getVelocity().getVector(flowField.getNx()+1, j, k)[1];
    flowField.getVelocity().getVector(1, j, k)[2] =
        flowField.getVelocity().getVector(flowField.getNx()+1, j, k)[2];
}

void PeriodicBoundaryVelocityStencil::applyRightWall(FlowField & flowField, int i, int j, int k){
    flowField.getVelocity().getVector(flowField.getNx()+2, j, k)[0] =
        flowField.getVelocity().getVector(2, j, k)[0];
    flowField.getVelocity().getVector(flowField.getNx()+2, j, k)[1] =
        flowField.getVelocity().getVector(2, j, k)[1];
    flowField.getVelocity().getVector(flowField.getNx()+2, j, k)[2] =
        flowField.getVelocity().getVector(2, j, k)[2];
}

void PeriodicBoundaryVelocityStencil::applyBottomWall(FlowField & flowField, int i, int j, int k){
    flowField.getVelocity().getVector(i, 1, k)[0] =
        flowField.getVelocity().getVector(i, flowField.getNy()+1, k)[0];
    flowField.getVelocity().getVector(i, 0, k)[1] =
        flowField.getVelocity().getVector(i, flowField.getNy(), k)[1];
    flowField.getVelocity().getVector(i, 1, k)[2] =
        flowField.getVelocity().getVector(i, flowField.getNy()+1, k)[2];
}

void PeriodicBoundaryVelocityStencil::applyTopWall(FlowField & flowField, int i, int j, int k){
    flowField.getVelocity().getVector(i, flowField.getNy()+2, k)[0] =
        flowField.getVelocity().getVector(i, 2, k)[0];
    flowField.getVelocity().getVector(i, flowField.getNy()+2, k)[1] =
        flowField.getVelocity().getVector(i, 2, k)[1];
    flowField.getVelocity().getVector(i, flowField.getNy()+2, k)[2] =
        flowField.getVelocity().getVector(i, 2, k)[2];
}
void PeriodicBoundaryVelocityStencil::applyFrontWall(FlowField & flowField, int i, int j, int k){
    flowField.getVelocity().getVector(i, j, 1)[0] =
        flowField.getVelocity().getVector(i, j, flowField.getNz()+1)[0];
    flowField.getVelocity().getVector(i, j, 1)[1] =
        flowField.getVelocity().getVector(i, j, flowField.getNz()+1)[1];
    flowField.getVelocity().getVector(i, j, 0)[2] =
        flowField.getVelocity().getVector(i, j, flowField.getNz())[2];
}

void PeriodicBoundaryVelocityStencil::applyBackWall(FlowField & flowField, int i, int j, int k){
    flowField.getVelocity().getVector(i, j, flowField.getNz()+2)[0] =
        flowField.getVelocity().getVector(i, j, 2)[0];
    flowField.getVelocity().getVector(i, j, flowField.getNz()+2)[1] =
        flowField.getVelocity().getVector(i, j, 2)[1];
    flowField.getVelocity().getVector(i, j, flowField.getNz()+2)[2] =
        flowField.getVelocity().getVector(i, j, 2)[2];
}

PeriodicBoundaryTurbViscosityStencil::PeriodicBoundaryTurbViscosityStencil(const Parameters & parameters):
    BoundaryStencil<TurbFlowField>(parameters) {}

// 2D problem

void PeriodicBoundaryTurbViscosityStencil::applyLeftWall(TurbFlowField & flowField, int i, int j){
    flowField.getTurbViscosity().getScalar(0, j) =
        flowField.getTurbViscosity().getScalar(flowField.getNx(), j);
}

void PeriodicBoundaryTurbViscosityStencil::applyRightWall(TurbFlowField & flowField, int i, int j){
    flowField.getTurbViscosity().getScalar(flowField.getNx()+2, j) =
        flowField.getTurbViscosity().getScalar(2, j);
}

void PeriodicBoundaryTurbViscosityStencil::applyBottomWall(TurbFlowField & flowField, int i, int j){
    flowField.getTurbViscosity().getScalar(i, 0) =
        flowField.getTurbViscosity().getScalar(i, flowField.getNy());
}

void PeriodicBoundaryTurbViscosityStencil::applyTopWall(TurbFlowField & flowField, int i, int j){
    flowField.getTurbViscosity().getScalar(i, flowField.getNy()+2) =
        flowField.getTurbViscosity().getScalar(i, 2);
}

// 3D Problem

void PeriodicBoundaryTurbViscosityStencil::applyLeftWall(TurbFlowField & flowField, int i, int j, int k){
    flowField.getTurbViscosity().getScalar(0, j, k) =
        flowField.getTurbViscosity().getScalar(flowField.getNx(), j, k);
}

void PeriodicBoundaryTurbViscosityStencil::applyRightWall(TurbFlowField & flowField, int i, int j, int k){
    flowField.getTurbViscosity().getScalar(flowField.getNx()+2, j, k) =
        flowField.getTurbViscosity().getScalar(2, j, k);
}

void PeriodicBoundaryTurbViscosityStencil::applyBottomWall(TurbFlowField & flowField, int i, int j, int k){
    flowField.getTurbViscosity().getScalar(i, 0, k) =
        flowField.getTurbViscosity().getScalar(i, flowField.getNy(), k);
}

void PeriodicBoundaryTurbViscosityStencil::applyTopWall(TurbFlowField & flowField, int i, int j, int k){
    flowField.getTurbViscosity().getScalar(i, flowField.getNy()+2, k) =
        flowField.getTurbViscosity().getScalar(i, 2, k);
}
void PeriodicBoundaryTurbViscosityStencil::applyFrontWall(TurbFlowField & flowField, int i, int j, int k){
    flowField.getTurbViscosity().getScalar(i, j, 0) =
        flowField.getTurbViscosity().getScalar(i, j, flowField.getNz());
}

void PeriodicBoundaryTurbViscosityStencil::applyBackWall(TurbFlowField & flowField, int i, int j, int k){
    flowField.getTurbViscosity().getScalar(i, j, flowField.getNz()+2) =
        flowField.getTurbViscosity().getScalar(i, j, 2);
}


PeriodicBoundaryFGHStencil::PeriodicBoundaryFGHStencil(const Parameters & parameters):
    BoundaryStencil<FlowField>(parameters) {}

// 2D problem
void PeriodicBoundaryFGHStencil::applyLeftWall(FlowField & flowField, int i, int j){}
void PeriodicBoundaryFGHStencil::applyRightWall(FlowField & flowField, int i, int j){}
void PeriodicBoundaryFGHStencil::applyBottomWall(FlowField & flowField, int i, int j){}
void PeriodicBoundaryFGHStencil::applyTopWall(FlowField & flowField, int i, int j){}

// 3D Problem
void PeriodicBoundaryFGHStencil::applyLeftWall(FlowField & flowField, int i, int j, int k){}
void PeriodicBoundaryFGHStencil::applyRightWall(FlowField & flowField, int i, int j, int k){}
void PeriodicBoundaryFGHStencil::applyBottomWall(FlowField & flowField, int i, int j, int k){}
void PeriodicBoundaryFGHStencil::applyTopWall(FlowField & flowField, int i, int j, int k){}
void PeriodicBoundaryFGHStencil::applyFrontWall(FlowField & flowField, int i, int j, int k){}
void PeriodicBoundaryFGHStencil::applyBackWall(FlowField & flowField, int i, int j, int k){}

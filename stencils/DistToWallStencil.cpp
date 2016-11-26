#include "DistToWallStencil.h"

DistToWallStencil::DistToWallStencil ( const Parameters & parameters ) : FieldStencil<TurbFlowField> ( parameters ) {}


void DistToWallStencil::apply ( TurbFlowField & flowField,  int i, int j ){
    // TODO
    handleError(1, "TODO Implement DistToWallStencil\n");
}


void DistToWallStencil::apply ( TurbFlowField & flowField, int i, int j, int k ){
    // TODO
    handleError(1, "TODO Implement DistToWallStencil\n");
}

/// TODO boundary values???
#include "TurbLPmodelStencil.h"

TurbLPmodel::TurbLPmodel ( const Parameters & parameters ) : FieldStencil<TurbFlowField> ( parameters ) {}

void TurbLPmodel::apply( TurbFlowField & flowField, int i, int j ){

    loadLocalVelocity2D(  flowField, _localVelocity, i, j);
    loadLocalMeshsize2D(_parameters, _localMeshsize, i, j);

    ///not included into  StencilFunctions because wall distance is stored in (turbulent)flowfield and that's not
    /// included in StencilFuncitons.h (only local..)
    FLOAT dudx_ =  dudx(_localVelocity,_localMeshsize);
    FLOAT dvdy_ =  dvdy(_localVelocity,_localMeshsize);
    FLOAT dudy_ =  dudy(_localVelocity,_localMeshsize);
    FLOAT dvdx_ =  dvdx(_localVelocity,_localMeshsize);

    FLOAT l_mix = flowField.getDistanceToWall().getScalar(i, j);

    flowField.setTurbulentViscosity(l_mix * l_mix *
            std::sqrt( 2*( dudx_ * dudx_  +  dvdy_ * dvdy_  +  2 * (dudy_ + dvdx_) * (dudy_ + dvdx_) ) ),i,j);

}

void TurbLPmodel::apply ( TurbFlowField & flowField, int i, int j, int k ){

    loadLocalVelocity3D(  flowField, _localVelocity, i, j, k);
    loadLocalMeshsize3D(_parameters, _localMeshsize, i, j, k);

    FLOAT dudx_ =  dudx(_localVelocity,_localMeshsize);
    FLOAT dudy_ =  dudy(_localVelocity,_localMeshsize);
    FLOAT dudz_ =  dudz(_localVelocity,_localMeshsize);
    FLOAT dvdx_ =  dvdx(_localVelocity,_localMeshsize);
    FLOAT dvdy_ =  dvdy(_localVelocity,_localMeshsize);
    FLOAT dvdz_ =  dvdz(_localVelocity,_localMeshsize);
    FLOAT dwdx_ =  dwdx(_localVelocity,_localMeshsize);
    FLOAT dwdy_ =  dwdy(_localVelocity,_localMeshsize);
    FLOAT dwdz_ =  dwdz(_localVelocity,_localMeshsize);

    FLOAT l_mix = flowField.getDistanceToWall().getScalar(i, j);

    flowField.setTurbulentViscosity(l_mix * l_mix *
        std::sqrt( 2*( dudx_ * dudx_  +  dvdy_ * dvdy_  +  dwdz_ * dwdz_ +
            2*( (dudy_ + dvdx_) * (dudy_ + dvdx_) + (dudz_ + dwdx_) * (dudz_ + dwdx_) + (dvdz_ + dwdy_) * (dvdz_ + dwdy_) ) ) ),i,j,k);
}

// void TurbLPmodel::getMixingLength( TurbFlowField & flowField, FLOAT& l_mix, int i, int j ) {
//     ///is boundary layer thickness going to be calculated for each time step? (probably not)
//     /// is boundary layer going to be calculated lecaly? how for BFS and LDC? (probably not) How to decide which value is chosen?
//     l_mix = std::min(flowField.getDistToWall().getScalar(i,j)*_parameters.turbulence.kappa, 0.09*_parameters.turbulence.bdLayerThickness);
//     /// make l_mix property that stores all l_mix (scalar field) in case BL stays constant => determine only once
//     /// to get wall distance: once stencil that contains the BoundaryIterator => Iterator iterates over all cells and
//     /// for each cell all boundary cells are iterated (compare x,y,z, vales => min distance).
// }
//
// void TurbLPmodel::getMixingLength( TurbFlowField & flowField, FLOAT& l_mix, int i, int j, int k ) {
//     ///is boundary layer thickness going to be calculated for each time step? (probably not)
//     /// is boundary layer going to be calculated lecaly? how for BFS and LDC? (probably not) How to decide which value is chosen?
//     l_mix = std::min(flowField.getDistToWall().getScalar(i,j,k)*_parameters.turbulence.kappa, 0.09*_parameters.turbulence.bdLayerThickness);
//     /// make l_mix property that stores all l_mix (scalar field) in case BL stays constant => determine only once
//     /// to get wall distance: once stencil that contains the BoundaryIterator => Iterator iterates over all cells and
//     /// for each cell all boundary cells are iterated (compare x,y,z, vales => min distance).
// }


/// not part of LPStencil...just as model
// void FGHStencil::apply ( FlowField & flowField, int i, int j, int k ){
//     // The same as in 2D, with slight modifications
//
//     const int obstacle = flowField.getFlags().getValue(i, j, k);
//
//     FLOAT * const values = flowField.getFGH().getVector(i,j,k);
//
//     if ((obstacle & OBSTACLE_SELF) == 0){   // If the cell is fluid
//
//         loadLocalVelocity3D(  flowField, _localVelocity, i, j, k);
//         loadLocalMeshsize3D(_parameters, _localMeshsize, i, j, k);
//
//         if ((obstacle & OBSTACLE_RIGHT) == 0) { // If the right cell is fluid
//             values [0] = computeF3D(_localVelocity, _localMeshsize, _parameters, _parameters.timestep.dt);
//         }
//         if ((obstacle & OBSTACLE_TOP) == 0) {
//             values [1] = computeG3D(_localVelocity, _localMeshsize, _parameters, _parameters.timestep.dt);
//         }
//         if ((obstacle & OBSTACLE_BACK) == 0) {
//             values [2] = computeH3D(_localVelocity, _localMeshsize, _parameters, _parameters.timestep.dt);
//         }
//     }
// }

/// TODO boundary values???
#include "LPmodelStencil.h"
#include "StencilFunctions.h"
//#include "Definitions.h"

LPmodel::LPmodel ( const Parameters & parameters ) : FieldStencil<FlowField> ( parameters ) {} ///TODO flowfield -> Turbulent flowfield

void LPmodel::apply( FlowField & flowField, int i, int j ){  ///TODO flowfield -> Turbulent flowfield

    loadLocalVelocity2D(  flowField, _localVelocity, i, j);
    loadLocalMeshsize2D(_parameters, _localMeshsize, i, j);

    ///not included into  StencilFunctions because wall distance is stored in (turbulent)flowfield and that's not
    /// included in StencilFuncitons.h (only local..)
    dudx_ =  dudx(_localVelocity,_localMeshsize);
    dvdy_ =  dvdy(_localVelocity,_localMeshsize);
    dudy_ =  dudy(_localVelocity,_localMeshsize);
    dvdx_ =  dvdx(_localVelocity,_localMeshsize);

    flowfieldTODO.getTurbulentViscosity().getScalar(i,j) = l_mix * l_mix *
            std::sqrt( 2*( dudx_ * dudx_  +  dvdy_ * dvdy_  +  2 * (dudy_ + dvdx_) * (dudy_ + dvdx_) ) );

}

void LPmodel::apply ( FlowField & flowField, int i, int j, int k ){ ///TODO flowfield -> Turbulent flowfield

    loadLocalVelocity3D(  flowField, _localVelocity, i, j, k);
    loadLocalMeshsize3D(_parameters, _localMeshsize, i, j, k);

    dudx_ =  dudx(_localVelocity,_localMeshsize);
    dudy_ =  dudy(_localVelocity,_localMeshsize);
    dudz_ =  dudz(_localVelocity,_localMeshsize);
    dvdx_ =  dvdx(_localVelocity,_localMeshsize);
    dvdy_ =  dvdy(_localVelocity,_localMeshsize);
    dvdz_ =  dvdz(_localVelocity,_localMeshsize);
    dwdx_ =  dwdx(_localVelocity,_localMeshsize);
    dwdy_ =  dwdy(_localVelocity,_localMeshsize);
    dwdz_ =  dwdz(_localVelocity,_localMeshsize);

    l_mix = flowfieldTODO.getWallDistance().getScalar(i,j,k);

    flowfieldTODO.getTurbulentViscosity().getScalar(i,j,k) = l_mix * l_mix *
        std::sqrt( 2*( dudx_ * dudx_  +  dvdy_ * dvdy_  +  dwdz * dwdz +
            2*( (dudy_ + dvdx_) * (dudy_ + dvdx_) + (dudz + dwdx) * (dudz + dwdx) + (dvdz + dwdy) * (dvdz + dwdy) ));
}

void LPmodel::getMixingLength(){
    ///is boundary layer thickness going to be calculated for each time step? (probably not)
    /// is boundary layer going to be calculated lecaly? how for BFS and LDC? (probably not) How to decide which value is chosen?
    l_mix = std::min(flowfieldTODO.getWallDistance().getScalar(i,j)*kappa, 0.09*BOUNDARYLAYERTHICKNESS);
    /// make l_mix property that stores all l_mix (scalar field) in case BL stays constant => determine only once
    /// to get wall distance: once stencil that contains the BoundaryIterator => Iterator iterates over all cells and
    /// for each cell all boundary cells are iterated (compare x,y,z, vales => min distance).
}

/// not part of LPStencil...just as model
void FGHStencil::apply ( FlowField & flowField, int i, int j, int k ){
    // The same as in 2D, with slight modifications

    const int obstacle = flowField.getFlags().getValue(i, j, k);

    FLOAT * const values = flowField.getFGH().getVector(i,j,k);

    if ((obstacle & OBSTACLE_SELF) == 0){   // If the cell is fluid

        loadLocalVelocity3D(  flowField, _localVelocity, i, j, k);
        loadLocalMeshsize3D(_parameters, _localMeshsize, i, j, k);

        if ((obstacle & OBSTACLE_RIGHT) == 0) { // If the right cell is fluid
            values [0] = computeF3D(_localVelocity, _localMeshsize, _parameters, _parameters.timestep.dt);
        }
        if ((obstacle & OBSTACLE_TOP) == 0) {
            values [1] = computeG3D(_localVelocity, _localMeshsize, _parameters, _parameters.timestep.dt);
        }
        if ((obstacle & OBSTACLE_BACK) == 0) {
            values [2] = computeH3D(_localVelocity, _localMeshsize, _parameters, _parameters.timestep.dt);
        }
    }
}

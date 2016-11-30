/// TODO boundary values???
#include "TurbLPmodelStencil.h"

TurbLPmodel::TurbLPmodel ( const Parameters & parameters ) : FieldStencil<TurbFlowField> ( parameters ) {}

void TurbLPmodel::apply( TurbFlowField & flowField, int i, int j ){

    // TODO: get mixing length. For now use distToWall*kappa
    FLOAT l_mix = flowField.getDistanceToWall().getScalar(i, j)*_parameters.turbulence.kappa;

    FLOAT tensorProd = 0.0;

    getShearStressTensorProduct(flowField, tensorProd, i, j);

    flowField.setTurbulentViscosity(l_mix * l_mix * tensorProd + _parameters.flow.viscosity, i, j);
}

void TurbLPmodel::apply ( TurbFlowField & flowField, int i, int j, int k ){

    // TODO: get mixing length. For now use distToWall*kappa
    FLOAT l_mix = flowField.getDistanceToWall().getScalar(i, j, k)*_parameters.turbulence.kappa;

    FLOAT tensorProd = 0.0;

    getShearStressTensorProduct(flowField, tensorProd, i, j, k);

    flowField.setTurbulentViscosity(l_mix * l_mix * tensorProd + _parameters.flow.viscosity, i, j, k);
}

void TurbLPmodel::getShearStressTensorProduct( TurbFlowField flowField,
    FLOAT &prod, int i, int j )
{
    loadLocalVelocity2D(  flowField, _localVelocity, i, j);
    loadLocalMeshsize2D( _parameters, _localMeshsize, i, j);

    FLOAT dudx_ =  dudx( _localVelocity, _localMeshsize );
    FLOAT dvdy_ =  dvdy( _localVelocity, _localMeshsize );
    FLOAT dudy_ =  dudy( _localVelocity, _localMeshsize );
    FLOAT dvdx_ =  dvdx( _localVelocity, _localMeshsize );

    prod = std::sqrt( 2*( dudx_ * dudx_  +  dvdy_ * dvdy_  +  2 * (dudy_ + dvdx_) * (dudy_ + dvdx_) ) );

}

void TurbLPmodel::getShearStressTensorProduct( TurbFlowField flowField,
    FLOAT &prod, int i, int j, int k )
{
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

    prod = std::sqrt( 2*( dudx_ * dudx_  +  dvdy_ * dvdy_  +  dwdz_ * dwdz_ +
        2*( (dudy_ + dvdx_) * (dudy_ + dvdx_) + (dudz_ + dwdx_) * (dudz_ + dwdx_) + (dvdz_ + dwdy_) * (dvdz_ + dwdy_) ) ) );
}

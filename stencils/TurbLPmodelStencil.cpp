#include "TurbLPmodelStencil.h"

TurbLPmodel::TurbLPmodel ( const Parameters & parameters ) : FieldStencil<TurbFlowField> ( parameters ) {}

void TurbLPmodel::apply( TurbFlowField & flowField, int i, int j ){
    // Don't check for fluid cell as we will get base visocity for obstacle cells which makes min turb viscosity non-zero : look at dt
    // TODO: get mixing length. For now use distToWall*kappa
    FLOAT l_mix = flowField.getMixingLength().getScalar(i, j);

    FLOAT tensorProd = 0.0;

    getShearStressTensorProduct(flowField, tensorProd, i, j);

    flowField.getTurbViscosity().getScalar(i, j) = l_mix * l_mix * tensorProd + 1/_parameters.flow.Re;
}

void TurbLPmodel::apply ( TurbFlowField & flowField, int i, int j, int k ){
    // Don't check for fluid cell as we will get base visocity for obstacle cells which makes min turb viscosity non-zero : look at dt
    // TODO: get mixing length. For now use distToWall*kappa
    FLOAT l_mix = flowField.getMixingLength().getScalar(i, j, k);

    FLOAT tensorProd = 0.0;

    getShearStressTensorProduct(flowField, tensorProd, i, j, k);

    flowField.getTurbViscosity().getScalar(i, j, k) = l_mix * l_mix * tensorProd + 1/_parameters.flow.Re;
}

void TurbLPmodel::getShearStressTensorProduct( TurbFlowField &flowField,
    FLOAT &prod, int i, int j )
{
    loadLocalVelocity2D(  flowField, _localVelocity, i, j);
    loadLocalMeshsize2D( _parameters, _localMeshsize, i, j);

    FLOAT dudx_ =  dudx( _localVelocity, _localMeshsize );
    FLOAT dvdy_ =  dvdy( _localVelocity, _localMeshsize );
    FLOAT dudy_ =  dudy( _localVelocity, _localMeshsize );
    FLOAT dvdx_ =  dvdx( _localVelocity, _localMeshsize );


    prod = std::sqrt( 2*dudx_ * dudx_  +  2*dvdy_ * dvdy_   +  (dudy_ + dvdx_) * (dudy_ + dvdx_) );

}

void TurbLPmodel::getShearStressTensorProduct( TurbFlowField &flowField,
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

    prod = std::sqrt( 2*(dudx_ * dudx_  +  dvdy_ * dvdy_  +  dwdz_ * dwdz_) +
        (dudy_ + dvdx_) * (dudy_ + dvdx_) + (dudz_ + dwdx_) * (dudz_ + dwdx_) + (dvdz_ + dwdy_) * (dvdz_ + dwdy_)   );
}

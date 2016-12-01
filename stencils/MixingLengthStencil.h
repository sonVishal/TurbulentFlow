#ifndef _STENCIL_MIX_LEN_H_
#define _STENCIL_MIX_LEN_H_

#include "../TurbFlowField.h"
#include "../Stencil.h"
#include "../Parameters.h"

class MixingLengthStencil : public FieldStencil<TurbFlowField> {
private:
    typedef FLOAT(MixingLengthStencil::*mixLenFuncPtr2D)(FLOAT, FLOAT, int,int);
    typedef FLOAT(MixingLengthStencil::*mixLenFuncPtr3D)(FLOAT, FLOAT, int,int,int);

    mixLenFuncPtr2D _mixLenFuncPtr2D[4];
    mixLenFuncPtr3D _mixLenFuncPtr3D[4];
    FLOAT _stepCornerPos[2];
    void getCellCenter( int i, int j, FLOAT* cellCenter);
    void getCellCenter( int i, int j, int k, FLOAT* cellCenter);

    FLOAT mixLenPrandtl2D(FLOAT distToWall, FLOAT meanVelocity, int i, int j);
    FLOAT mixLenPrandtl3D(FLOAT distToWall, FLOAT meanVelocity, int i, int j, int k);

    FLOAT mixLenIgnoreBdLayer2D(FLOAT distToWall, FLOAT meanVelocity, int i, int j);
    FLOAT mixLenIgnoreBdLayer3D(FLOAT distToWall, FLOAT meanVelocity, int i, int j, int k);

    FLOAT mixLenLaminarPlate2D(FLOAT distToWall, FLOAT meanVelocity, int i, int j);
    FLOAT mixLenLaminarPlate3D(FLOAT distToWall, FLOAT meanVelocity, int i, int j, int k);

    FLOAT mixLenTurbulentPlate2D(FLOAT distToWall, FLOAT meanVelocity, int i, int j);
    FLOAT mixLenTurbulentPlate3D(FLOAT distToWall, FLOAT meanVelocity, int i, int j, int k);
public:

    MixingLengthStencil ( const Parameters & parameters );

    /** Apply the stencil in 2D
     *
     * Performs the operation of the stencil in a single position given by the indexes.
     * @param flowField State of the flow
     * @param i Index in the x direction
     * @param j Index in the y direction
     */
    void apply ( TurbFlowField & flowField, int i, int j );

    /** Apply the stencil in 3D
     *
     * @param flowField State of the flow
     * @param i Index in the x direction
     * @param j Index in the y direction
     * @param k Index in the z direction
     */
    void apply ( TurbFlowField & flowField, int i, int j, int k );
};
#endif

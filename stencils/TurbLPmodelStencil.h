#ifndef _LPMODEL_STENCIL_H_
#define _LPMODEL_STENCIL_H_

#include "../Stencil.h"
#include "../Parameters.h"
#include "../TurbFlowField.h"
#include "StencilFunctions.h"
// #include "Definitions.h"

/** Stencil to compute the velocity once the pressure has been found.
 */
class TurbLPmodel : public FieldStencil<TurbFlowField> {

    private:
        // A local velocity variable that will be used to approximate derivatives. Size matches 3D
        // case, but can be used for 2D as well.
        FLOAT _localVelocity [ 27 * 3 ];
        // local meshsize
        FLOAT _localMeshsize [ 27 * 3 ];

        void getShearStressTensorProduct( TurbFlowField flowField,
            FLOAT &prod, int i, int j );

        void getShearStressTensorProduct( TurbFlowField flowField,
            FLOAT &prod, int i, int j, int k );

    public:

        /** Constructor
         * @param parameters Parameters of the problem
         */
        TurbLPmodel(const Parameters & parameters);

        /** Apply the stencil in 2D
         * @param flowField Flow field information
         * @param i Position in the X direction
         * @param j Position in the Y direction
         */
        void apply ( TurbFlowField & flowField, int i, int j );

        /** Apply the stencil in 3D
         * @param flowField Flow field information
         * @param i Position in the X direction
         * @param j Position in the Y direction
         * @param k Position in the Z direction
         */
        void apply ( TurbFlowField & flowField, int i, int j, int k );

        // /** get the mixing length
        //  *
        //  */
        // void getMixingLength( TurbFlowField & flowField, FLOAT& l_mix, int i, int j );
        // void getMixingLength( TurbFlowField & flowField, FLOAT& l_mix, int i, int j, int k );
};



#endif

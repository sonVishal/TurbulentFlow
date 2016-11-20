#ifndef _LPMODEL_STENCIL_H_
#define _LPMODEL_STENCIL_H_

#include "../Stencil.h"
#include "../Parameters.h"
#include "../FlowField.h"

/** Stencil to compute the velocity once the pressure has been found.
 */
class LPmodel : public FieldStencil<FlowField> {    ///TODO flowfield -> Turbulent flowfield

    private:
        // A local velocity variable that will be used to approximate derivatives. Size matches 3D
        // case, but can be used for 2D as well.
        FLOAT _localVelocity [ 27 * 3 ];
        // local meshsize
        FLOAT _localMeshsize [ 27 * 3 ];

    public:

        /** Constructor
         * @param parameters Parameters of the problem
         */
        LPmodel(const Parameters & parameters);

        /** Apply the stencil in 2D
         * @param flowField Flow field information
         * @param i Position in the X direction
         * @param j Position in the Y direction
         */
        void apply ( FlowField & flowField, int i, int j ); ///TODO flowfield -> Turbulent flowfield

        /** Apply the stencil in 3D
         * @param flowField Flow field information
         * @param i Position in the X direction
         * @param j Position in the Y direction
         * @param k Position in the Z direction
         */
        void apply ( FlowField & flowField, int i, int j, int k ); ///TODO flowfield -> Turbulent flowfield

        /** get the mixing length
         *
         */
        void getMixingLength();
};



#endif

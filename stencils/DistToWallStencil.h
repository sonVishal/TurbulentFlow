#ifndef _STENCIL_DIST_TO_WALL_
#define _STENCIL_DIST_TO_WALL_

#include "../TurbFlowField.h"
#include "../Stencil.h"
#include "../Parameters.h"

class DistToWallStencil : public FieldStencil<TurbFlowField> {
private:

public:

    DistToWallStencil ( const Parameters & parameters );

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

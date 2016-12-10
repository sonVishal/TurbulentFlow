#ifndef _MIN_DT_STENCIL_H_
#define _MIN_DT_STENCIL_H_

#include "../Stencil.h"
#include "../Parameters.h"
#include "../TurbFlowField.h"


/** this class computes the maximum value of max(velocity)/meshsize for all grid cells.
 *  Originally, one would compute the max. velocity only and adapt it with the meshsize afterwards.
 *  This, however, becomes inconsistent when dealing with non-uniform, e.g. stretched, meshes, since
 *  the meshsize may be different for every grid cell. We therefore determine the max(velocity)/meshsize
 *  and synchronise this value over whole computational domain.
 *  @author Philipp Neumann
 */
class MinDtStencil : public FieldStencil<TurbFlowField> {

    private:

        FLOAT _minValue;

        /** Sets the maximum value to the value of the cell if it surpasses the current one.
         *
         * 2D version of the function
         * @param flowField Flow field
         * @param i Position in the X direction.
         * @param j Position in the Y direction.
         */
        void cellMinValue(TurbFlowField & flowField, int i, int j);

        /** Sets the maximum value to the value of the cell if it surpasses the current one.
         *
         * 3D version of the function
         * @param flowField Flow field
         * @param i Position in the X direction.
         * @param j Position in the Y direction.
         * @param k Position in the Z direction.
         */
        void cellMinValue(TurbFlowField & flowField, int i, int j, int k);

    public:

        /** Constructor
         *
         * @param parameters Parameters of the problem
         */
        MinDtStencil (const Parameters & parameters);

        //@ brief Body iterations
        //@{
        void apply (TurbFlowField & flowField, int i, int j);
        void apply (TurbFlowField & flowField, int i, int j, int k);
        //@}

        /** Resets the maximum values to zero before computing the timestep
         */
        void reset ();

        /** Returns the array with the maximum modules of the components of the velocity,
         *  divided by the respective local meshsize
         */
        const FLOAT getMinValue() const;
};

#endif

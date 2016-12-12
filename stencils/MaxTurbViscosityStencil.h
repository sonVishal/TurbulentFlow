#ifndef _MAX_NU_T_STENCIL_H_
#define _MAX_NU_T_STENCIL_H_

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
class MaxTurbViscosityStencil : public FieldStencil<TurbFlowField>, public BoundaryStencil<TurbFlowField> {

    private:

        FLOAT _maxValue;

        /** Sets the maximum value to the value of the cell if it surpasses the current one.
         *
         * 2D version of the function
         * @param flowField Flow field
         * @param i Position in the X direction.
         * @param j Position in the Y direction.
         */
        void cellMaxValue(TurbFlowField & flowField, int i, int j);

        /** Sets the maximum value to the value of the cell if it surpasses the current one.
         *
         * 3D version of the function
         * @param flowField Flow field
         * @param i Position in the X direction.
         * @param j Position in the Y direction.
         * @param k Position in the Z direction.
         */
        void cellMaxValue(TurbFlowField & flowField, int i, int j, int k);

    public:

        /** Constructor
         *
         * @param parameters Parameters of the problem
         */
        MaxTurbViscosityStencil (const Parameters & parameters);

        //@ brief Body iterations
        //@{
        void apply (TurbFlowField & flowField, int i, int j);
        void apply (TurbFlowField & flowField, int i, int j, int k);
        //@}

        //@ brief Boundary iterations for the 2D problem
        //@param flowField Flow field with the state of the fluid
        //@param i Position in the X direction
        //@param j Position in the Y direction
        //@{
        void applyLeftWall   ( TurbFlowField & flowField, int i, int j );
        void applyRightWall  ( TurbFlowField & flowField, int i, int j );
        void applyBottomWall ( TurbFlowField & flowField, int i, int j );
        void applyTopWall    ( TurbFlowField & flowField, int i, int j );
        //@}

        //@ brief Boundary iterations for the 3D problem
        //@param flowField Flow field with the state of the fluid
        //@param i Position in the X direction
        //@param j Position in the Y direction
        //@param k Position in the Z direction
        //@{
        void applyLeftWall   ( TurbFlowField & flowField, int i, int j, int k );
        void applyRightWall  ( TurbFlowField & flowField, int i, int j, int k );
        void applyBottomWall ( TurbFlowField & flowField, int i, int j, int k );
        void applyTopWall    ( TurbFlowField & flowField, int i, int j, int k );
        void applyFrontWall  ( TurbFlowField & flowField, int i, int j, int k );
        void applyBackWall   ( TurbFlowField & flowField, int i, int j, int k );
        //@}

        /** Resets the maximum values to zero before computing the timestep
         */
        void reset ();

        /** Returns the array with the maximum modules of the components of the velocity,
         *  divided by the respective local meshsize
         */
        FLOAT getMaxValue() const;
};

#endif

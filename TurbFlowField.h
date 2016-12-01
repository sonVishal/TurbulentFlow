#ifndef _TURB_FLOW_FIELD_H_
#define _TURB_FLOW_FIELD_H_

#include "FlowField.h"
#include "DataStructures.h"
#include "Parameters.h"

/** Turbulent Flow field
 *
 * Class intended to contain the state of the domain for the turbulent simulation.
 */
class TurbFlowField : public FlowField {

    private:

        ScalarField _turbViscosity; //! Scalar field representing the turbulent viscosity (= nu + nu_t)

        ScalarField _mixingLength; //! Scalar field representing the mixing length
        // Mixing length is computed here as it depends on static parameters
        // NOTE: Works only for Channel flow since we know the input velocity profile
        // Otherwise, for a general case we need to communicate the mean velocity at each x location
        //  and use that to compute the Re_x term in the boundary layer thickness. Sorry!

    public:

        /** Constructor for the 2D flow field
         *
         * Constructor for the flow field. Allocates all the fields and sets
         * the sizes. Currently, this contructor is only used for testing purposes.
         *
         * @param Nx Size of the fuild domain (non-ghost cells), in the X direction
         * @param Ny Size of the fuild domain (non-ghost cells), in the Y direction
         */
        TurbFlowField ( int Nx, int Ny );

        /** Constructor for the 3D flow field
         *
         * Constructor for the flow field. Allocates all the fields and sets
         * the sizes. Currently, this contructor is only used for testing purposes.
         *
         * @param Nx Size of the fuild domain (non-ghost cells), in the X direction
         * @param Ny Size of the fuild domain (non-ghost cells), in the Y direction
         * @param Nz Size of the fuild domain (non-ghost cells), in the Z direction
         */
        TurbFlowField (int Nx, int Ny, int Nz);

        /** Constructs a field from parameters object
         *
         * Constructs a field from a parameters object, so that it dimensionality can be defined in
         * the configuration file.
         *
         * @param parameters Parameters object with geometric information
         */
        TurbFlowField (const Parameters & parameters);

        /** Destructor
        */
        ~TurbFlowField ();

        /** Get the turbulent viscosity
         * @return Scalar field with the turbulent viscosity
         */
        ScalarField & getTurbViscosity ();

        /** Get the Mixing length
         * @return Scalar field with the mixing length
         */
        ScalarField & getMixingLength ();

        void getPressureVelocityAndTurbVisc(FLOAT &pressure, FLOAT &turbViscosity, FLOAT* const velocity, int i, int j);
        void getPressureVelocityAndTurbVisc(FLOAT &pressure, FLOAT &turbViscosity, FLOAT* const velocity, int i, int j, int k);
};

#endif

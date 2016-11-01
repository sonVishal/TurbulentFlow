#ifndef _VTK_STENCIL_H_
#define _VTK_STENCIL_H_

#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../FlowField.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

/** Stencil for writting VTK files
 *
 * When iterated with, creates a VTK file.
 */
class VTKStencil : public FieldStencil<FlowField> {
    private:

        /** File stream variable used for writing the VTK file
         */
        std::ofstream _outputFileHandle;

        /** String stream variable used for storing pressure values at current time
         */
        std::stringstream _pressure;

        /** String stream variable used for storing velocity values at current time
         */
        std::stringstream _velocity;

        /** Size of the mesh in X,Y, and Z direction
         */
        int _sizeX, _sizeY, _sizeZ;

        /** Writes the header for the VTK file as well as the coordinates
         */
        void writeHeaderAndCoords();

    public:

        /** Constructor
         *
         * @param prefix String with the prefix of the name of the VTK files
         */
        VTKStencil ( const Parameters& parameters );

        /** 2D operation for one position
         *
         * @param flowField State of the flow field
         * @param i Position in the x direction
         * @param j Position in the y direction
         */
        void apply ( FlowField & flowField, int i, int j );

        /** 3D operation for one position
         *
         * @param flowField State of the flow field
         * @param i Position in the x direction
         * @param j Position in the y direction
         * @param k Position in the z direction
         */
        void apply ( FlowField & flowField, int i, int j, int k );

        /** Writes the information to the file
         * @param flowField Flow field to be written
         */
        void write ( FlowField & flowField, int timeStep );

};

#endif

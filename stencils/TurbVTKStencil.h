#ifndef _TURB_VTK_STENCIL_H_
#define _TURB_VTK_STENCIL_H_

#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../TurbFlowField.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

/** Stencil for writting VTK files
 *
 * When iterated with, creates a VTK file.
 */
class TurbVTKStencil : public FieldStencil<TurbFlowField> {
    private:

        /** File stream variable used for writing the VTK file
         */
        std::ofstream _outputFile;

        /** String stream variable used for storing pressure values at current time
         * As per my knowledge the stream buffers are internally allocated on the heap
         */
        std::ostringstream _pressureStream;

        /** String stream variable used for storing turbulent viscosity at current time
         * As per my knowledge the stream buffers are internally allocated on the heap
         */
        std::ostringstream _turbViscosityStream;

        /** String stream variable used for storing distance to wall at current time
         * As per my knowledge the stream buffers are internally allocated on the heap
         */
        std::ostringstream _turbDistToWall;

        /** String stream variable used for storing velocity values at current time
         * As per my knowledge the stream buffers are internally allocated on the heap
         */
        std::ostringstream _velocityStream;

        /** Local size of the mesh in X,Y, and Z direction
         */
        int _localSize[3];

        /** Writes the header for the VTK file as well as the coordinates
         */
        void writeHeaderAndCoords();

    public:

        /** Constructor
         *
         * @param prefix String with the prefix of the name of the VTK files
         */
        TurbVTKStencil ( const Parameters& parameters );

        /** 2D operation for one position
         *
         * @param flowField State of the flow field
         * @param i Position in the x direction
         * @param j Position in the y direction
         */
        void apply ( TurbFlowField & flowField, int i, int j );

        /** 3D operation for one position
         *
         * @param flowField State of the flow field
         * @param i Position in the x direction
         * @param j Position in the y direction
         * @param k Position in the z direction
         */
        void apply ( TurbFlowField & flowField, int i, int j, int k );

        /** Writes the information to the file
         * @param flowField Flow field to be written
         */
        void write ( TurbFlowField & flowField, int timeStep );

        /** Opens the file for writing and returns true if it was opened otherwise false
         * @param timeStep used to create the file name for current time step
         */
        bool openFile ( int timeStep );

};

#endif

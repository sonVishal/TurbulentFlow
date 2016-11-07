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
        std::ofstream _outputFile;

        /** String stream variable used for storing pressure values at current time
         * Put on the heap since it will become large
         */
        std::stringstream* _pressureStream;

        /** String stream variable used for storing velocity values at current time
         * Put on the heap since it will become large
         */
        std::stringstream* _velocityStream;

        /** Flag indicating if the string streams were allocated memory
         */
        bool _streamMemoryAllocFlag;

        /** Local size of the mesh in X,Y, and Z direction
         */
        const int* _localSize;

        /** Writes the header for the VTK file as well as the coordinates
         */
        void writeHeaderAndCoords();

    public:

        /** Constructor
         *
         * @param prefix String with the prefix of the name of the VTK files
         */
        VTKStencil ( const Parameters& parameters );

        /** Destructor
         */
        ~VTKStencil();

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

        /** Opens the file for writing and returns true if it was opened otherwise false
         * @param timeStep used to create the file name for current time step
         */
        bool openFile ( int timeStep );

};

#endif

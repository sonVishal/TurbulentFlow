#include "VTKStencil.h"

/** Constructor
 *
 * @param prefix String with the prefix of the name of the VTK files
 */
VTKStencil::VTKStencil ( const Parameters & parameters ) : FieldStencil<FlowField> ( parameters ) { }

/** 2D operation for one position
 *
 * @param flowField State of the flow field
 * @param i Position in the x direction
 * @param j Position in the y direction
 */
void VTKStencil::apply ( FlowField & flowField, int i, int j ) {

}

/** 3D operation for one position
 *
 * @param flowField State of the flow field
 * @param i Position in the x direction
 * @param j Position in the y direction
 * @param k Position in the z direction
 */
void VTKStencil::apply ( FlowField & flowField, int i, int j, int k ) {

}

void VTKStencil::writeHeader( std::ofstream& fileHandle ) {
    fileHandle << "# vtk DataFile Version 2.0\n" << "I need something to put here\n";
    fileHandle << "ASCII\n" << std::endl;
}

/** Writes the information to the file
 * @param flowField Flow field to be written
 */
void VTKStencil::write ( FlowField & flowField, int timeStep ) {
    std::stringstream fileName;
    fileName << _parameters.vtk.prefix << '_' << timeStep << ".vtk";
    std::ofstream vtkFileHandle (fileName.str().c_str());
    if(vtkFileHandle.is_open()) {
        writeHeader(vtkFileHandle);
        vtkFileHandle.close();
    }
}

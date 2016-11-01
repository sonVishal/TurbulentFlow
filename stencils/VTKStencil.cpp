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

void VTKStencil::writeHeaderAndCoords() {
    _outputFileHandle << "# vtk DataFile Version 2.0\n" << "I need something to put here\n";
    _outputFileHandle << "ASCII\n" << std::endl;

    const int* localSize = _parameters.parallel.localSize;

    int sizeXp1 = localSize[0]+1;
    int sizeYp1 = localSize[1]+1;
    int sizeZp1 = localSize[2]+1;

    _outputFileHandle << "DATASET STRUCTURED_GRID\n" << "DIMENSIONS " << sizeXp1 << ' ' << sizeYp1 << ' ' << sizeZp1 << std::endl;
    _outputFileHandle << "POINTS " << sizeXp1*sizeYp1*sizeZp1 << " float" << std::endl;

    FLOAT x = 0.0, y = 0.0, z = 0.0;

    const int* firstCorner = _parameters.parallel.firstCorner;

    sizeXp1 += firstCorner[0];
    sizeYp1 += firstCorner[1];
    sizeZp1 += firstCorner[2];

    if (sizeZp1 == 1) {
        for (int j = firstCorner[1]; j < sizeYp1; j++) {
            y = _parameters.meshsize->getPosY(firstCorner[0],j);
            for (int i = firstCorner[0]; i < sizeXp1; i++) {
                x = _parameters.meshsize->getPosX(i,j);
                _outputFileHandle << x << ' ' << y << ' ' << z << std::endl;
            }
        }
    } else {
        for (int k = firstCorner[2]; k < sizeZp1; k++) {
            z = _parameters.meshsize->getPosZ(firstCorner[0],firstCorner[1],k);
            for (int j = firstCorner[1]; j < sizeYp1; j++) {
                y = _parameters.meshsize->getPosY(firstCorner[0],j,k);
                for (int i = firstCorner[0]; i < sizeXp1; i++) {
                    x = _parameters.meshsize->getPosX(i,j,k);
                    _outputFileHandle << x << ' ' << y << ' ' << z << std::endl;
                }
            }
        }
    }
}

/** Writes the information to the file
 * @param flowField Flow field to be written
 */
void VTKStencil::write ( FlowField & flowField, int timeStep ) {
    std::stringstream fileName;
    fileName << _parameters.vtk.prefix << '_' << timeStep << ".vtk";
    _outputFileHandle.open(fileName.str().c_str());
    if(_outputFileHandle.is_open()) {
        writeHeaderAndCoords();
        _outputFileHandle.close();
    }
}

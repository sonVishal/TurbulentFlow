#include "VTKStencil.h"

/** Constructor
 *
 * @param prefix String with the prefix of the name of the VTK files
 */
VTKStencil::VTKStencil ( const Parameters & parameters ) : FieldStencil<FlowField> ( parameters ) {
    _numLocalCells = parameters.parallel.localSize[0]*parameters.parallel.localSize[1];
    _numLocalCells *= (parameters.parallel.localSize[2] == 0)?1:parameters.parallel.localSize[2];
}

/** 2D operation for one position
 *
 * @param flowField State of the flow field
 * @param i Position in the x direction
 * @param j Position in the y direction
 */
void VTKStencil::apply ( FlowField & flowField, int i, int j ) {
    if (i == 1 || j == 1) {
        return;
    }
    if ((flowField.getFlags().getValue(i,j) & OBSTACLE_SELF) == 0) { // Fluid cell
        FLOAT pressure;
        FLOAT velocity[2];
        flowField.getPressureAndVelocity(pressure, velocity, i, j);
        _pressure << pressure << std::endl;
        _velocity << velocity[0] << ' ' << velocity[1] << ' ' << 0 << std::endl;
    } else {
        _pressure << 0 << std::endl;
        _velocity << 0 << ' ' << 0 << ' ' << 0 << std::endl;
    }
}

/** 3D operation for one position
 *
 * @param flowField State of the flow field
 * @param i Position in the x direction
 * @param j Position in the y direction
 * @param k Position in the z direction
 */
void VTKStencil::apply ( FlowField & flowField, int i, int j, int k ) {
    if (i == 1 || j == 1 || k == 1) {
        return;
    }
    if ((flowField.getFlags().getValue(i,j,k) & OBSTACLE_SELF) == 0) { // Fluid cell
        FLOAT pressure;
        FLOAT velocity[3];
        flowField.getPressureAndVelocity(pressure, velocity, i, j, k);
        _pressure << pressure << std::endl;
        _velocity << velocity[0] << ' ' << velocity[1] << ' ' << velocity[2] << std::endl;
    } else {
        _pressure << 0 << std::endl;
        _velocity << 0 << ' ' << 0 << ' ' << 0 << std::endl;
    }
}

void VTKStencil::writeHeaderAndCoords() {
    _outputFileHandle << "# vtk DataFile Version 2.0\n" << "I need something to put here\n";
    _outputFileHandle << "ASCII\n" << std::endl;

    static const int* localSize = _parameters.parallel.localSize;
    static const int* firstCorner = _parameters.parallel.firstCorner;

    static const int sizeXp1 = localSize[0]-firstCorner[0]+1;
    static const int sizeYp1 = localSize[1]-firstCorner[1]+1;
    static const int sizeZp1 = localSize[2]-firstCorner[2]+1;

    _outputFileHandle << "DATASET STRUCTURED_GRID\n" << "DIMENSIONS " << sizeXp1 << ' ' << sizeYp1 << ' ' << sizeZp1 << std::endl;
    _outputFileHandle << "POINTS " << sizeXp1*sizeYp1*sizeZp1 << " float" << std::endl;

    FLOAT x = 0.0, y = 0.0, z = 0.0;

    if (_parameters.geometry.dim == 2) {
        for (int j = 2; j < sizeYp1+2; j++) {
            y = _parameters.meshsize->getPosY(2,j);
            for (int i = 2; i < sizeXp1+2; i++) {
                x = _parameters.meshsize->getPosX(i,j);
                _outputFileHandle << x << ' ' << y << ' ' << z << std::endl;
            }
        }
    } else {
        for (int k = 2; k < sizeZp1+2; k++) {
            z = _parameters.meshsize->getPosZ(2,2,k);
            for (int j = 2; j < sizeYp1+2; j++) {
                y = _parameters.meshsize->getPosY(2,j,k);
                for (int i = 2; i < sizeXp1+2; i++) {
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
    _outputFileHandle.open(fileName.str().c_str(),std::ios::out);
    if(_outputFileHandle.is_open()) {
        writeHeaderAndCoords();
        _outputFileHandle << "\nCELL_DATA " << _numLocalCells << std::endl;
        _outputFileHandle << "SCALARS pressure float 1\n";
        _outputFileHandle << "LOOKUP_TABLE default\n";
        _outputFileHandle << _pressure.str();
        _outputFileHandle << "\nVECTORS velocity float\n";
        _outputFileHandle << _velocity.str();
        _outputFileHandle.close();
    }
    _pressure.str(std::string());
    _velocity.str(std::string());
}

#include "VTKStencil.h"

/** Constructor
 *
 * @param prefix String with the prefix of the name of the VTK files
 */
VTKStencil::VTKStencil ( const Parameters & parameters ) : FieldStencil<FlowField> ( parameters ) {
    // Get the local size of the mesh
    _sizeX = parameters.parallel.localSize[0];
    _sizeY = parameters.parallel.localSize[1];
    _sizeZ = parameters.parallel.localSize[2];

    _pressure << std::setprecision(6);
    _velocity << std::setprecision(6);
}

/** 2D operation for one position
 *
 * @param flowField State of the flow field
 * @param i Position in the x direction
 * @param j Position in the y direction
 */
void VTKStencil::apply ( FlowField & flowField, int i, int j ) {
    // Since the field iterator includes the ghost layers we skip them
    if (i == 1 || j == 1) {
        return;
    }

    // Check if fluid cell and then print out to the stringstreams
    if ((flowField.getFlags().getValue(i,j) & OBSTACLE_SELF) == 0) {
        FLOAT pressure;
        FLOAT velocity[2];
        flowField.getPressureAndVelocity(pressure, velocity, i, j);
        _pressure << pressure << std::endl;
        _velocity << velocity[0] << ' ' << velocity[1] << ' ' << 0 << std::endl;
    } else {
    // Otherwise output 0 values
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
    // Since the field iterator includes the ghost layers we skip them
    if (i == 1 || j == 1 || k == 1) {
        return;
    }

    // Check if fluid cell and then print out to the stringstreams
    if ((flowField.getFlags().getValue(i,j,k) & OBSTACLE_SELF) == 0) {
        FLOAT pressure;
        FLOAT velocity[3];
        flowField.getPressureAndVelocity(pressure, velocity, i, j, k);
        _pressure << pressure << std::endl;
        _velocity << velocity[0] << ' ' << velocity[1] << ' ' << velocity[2] << std::endl;
    } else {
    // Otherwise output 0 values
        _pressure << 0 << std::endl;
        _velocity << 0 << ' ' << 0 << ' ' << 0 << std::endl;
    }
}

/** Writes the header for the VTK file as well as the coordinates
 */
void VTKStencil::writeHeaderAndCoords() {
    // This is the header
    _outputFileHandle << "# vtk DataFile Version 2.0\n" << "I need something to put here\n";
    _outputFileHandle << "ASCII\n" << std::endl;

    // Write the header for the coordinates
    // Number of coordinates is 1 more than the number of cells in that direction
    _outputFileHandle << "DATASET STRUCTURED_GRID\n" << "DIMENSIONS ";
    _outputFileHandle << _sizeX+1 << ' ' << _sizeY+1 << ' ' << _sizeZ+1 << std::endl;
    _outputFileHandle << "POINTS " << (_sizeX+1)*(_sizeY+1)*(_sizeZ+1) << " float" << std::endl;

    // Temporary variables to avoid calling function many times
    FLOAT y = 0.0, z = 0.0;

    // Call for loop depending on dimension
    // Start from 2 as we skip the ghost cells (both lower and upper)
    if (_parameters.geometry.dim == 2) {
        for (int j = 2; j < _sizeY+3; j++) {
            y = _parameters.meshsize->getPosY(2,j);
            for (int i = 2; i < _sizeX+3; i++) {
                _outputFileHandle << _parameters.meshsize->getPosX(i,j);
                _outputFileHandle << ' ' << y << ' ' << z << std::endl;
            }
        }
    } else {
        for (int k = 2; k < _sizeZ+3; k++) {
            z = _parameters.meshsize->getPosZ(2,2,k);
            for (int j = 2; j < _sizeY+3; j++) {
                y = _parameters.meshsize->getPosY(2,j,k);
                for (int i = 2; i < _sizeX+3; i++) {
                    _outputFileHandle << _parameters.meshsize->getPosX(i,j,k);
                    _outputFileHandle << ' ' << y << ' ' << z << std::endl;
                }
            }
        }
    }
}

/** Writes the information to the file
 * @param flowField Flow field to be written
 */
void VTKStencil::write ( FlowField & flowField, int timeStep ) {
    // fileName stores the file name as per XML file
    // it adds current timestep + .vtk extension to it
    std::stringstream fileName;
    fileName << _parameters.vtk.prefix << '_' << timeStep << ".vtk";

    // Open the file for writing out
    _outputFileHandle.open(fileName.str().c_str(),std::ios::out);
    // If file is opened write the stuff
    if(_outputFileHandle.is_open()) {
        // Write the header and coordinates
        writeHeaderAndCoords();
        // Write the header for pressure field
        _outputFileHandle << "\nCELL_DATA " << _sizeX*_sizeY*(_sizeZ==0?1:_sizeZ) << std::endl;
        _outputFileHandle << "SCALARS pressure float 1\n";
        _outputFileHandle << "LOOKUP_TABLE default\n";
        // Write the pressure
        _outputFileHandle << _pressure.str();
        // Write the header for velocity field
        _outputFileHandle << "\nVECTORS velocity float\n";
        // Write the velocity
        _outputFileHandle << _velocity.str();
        // Close the file
        _outputFileHandle.close();
    }
    // Clear the stringstream buffers
    _pressure.str(std::string());
    _velocity.str(std::string());
}

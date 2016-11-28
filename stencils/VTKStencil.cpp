#include "VTKStencil.h"

/** Constructor
 *
 * @param prefix String with the prefix of the name of the VTK files
 */
VTKStencil::VTKStencil ( const Parameters & parameters ) : FieldStencil<FlowField> ( parameters ) {
    // Get the local size and first corner of this subdomain
    _localSize[0] = parameters.parallel.localSize[0];
    _localSize[1] = parameters.parallel.localSize[1];
    _localSize[2] = parameters.parallel.localSize[2];
    if (parameters.geometry.dim == 2) {
        _localSize[2] = 0;
    }

    _pressureStream << std::setprecision(6);
    _velocityStream << std::setprecision(6);
}

/** 2D operation for one position
 *
 * @param flowField State of the flow field
 * @param i Position in the x direction
 * @param j Position in the y direction
 */
void VTKStencil::apply ( FlowField & flowField, int i, int j ) {
    // Since the field iterator includes the left and bottom ghost layers we skip them
    if (i == 1 || j == 1) {
        return;
    }

    // Check if fluid cell and then print out to the stringstreams
    if ((flowField.getFlags().getValue(i,j) & OBSTACLE_SELF) == 0) {
        FLOAT pressure;
        FLOAT velocity[2];
        flowField.getPressureAndVelocity(pressure, velocity, i, j);
        _pressureStream << pressure << std::endl;
        _velocityStream << velocity[0] << ' ' << velocity[1] << ' ' << 0 << std::endl;
    } else {
    // Otherwise output 0 values
        _pressureStream << 0 << std::endl;
        _velocityStream << 0 << ' ' << 0 << ' ' << 0 << std::endl;
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
    // Since the field iterator include the left, bottom and front ghost layers we skip them
    if (i == 1 || j == 1 || k == 1) {
        return;
    }

    // Check if fluid cell and then print out to the stringstreams
    if ((flowField.getFlags().getValue(i,j,k) & OBSTACLE_SELF) == 0) {
        FLOAT pressure;
        FLOAT velocity[3];
        flowField.getPressureAndVelocity(pressure, velocity, i, j, k);
        _pressureStream << pressure << std::endl;
        _velocityStream << velocity[0] << ' ' << velocity[1] << ' ' << velocity[2] << std::endl;
    } else {
    // Otherwise output 0 values
        _pressureStream << 0 << std::endl;
        _velocityStream << 0 << ' ' << 0 << ' ' << 0 << std::endl;
    }
}

/** Writes the header for the VTK file as well as the coordinates
 */
void VTKStencil::writeHeaderAndCoords() {
    // This is the header
    _outputFile << "# vtk DataFile Version 2.0\n" << "I need something to put here\n";
    _outputFile << "ASCII\n" << std::endl;

    // Write the header for the coordinates
    // Number of coordinates is 1 more than the number of cells in that direction
    _outputFile << "DATASET STRUCTURED_GRID\n" << "DIMENSIONS ";
    _outputFile << _localSize[0]+1 << ' ' << _localSize[1]+1 << ' ' << _localSize[2]+1 << std::endl;
    _outputFile << "POINTS " << (_localSize[0]+1)*(_localSize[1]+1)*(_localSize[2]+1) << " float" << std::endl;

    // Temporary variables to avoid calling function many times
    FLOAT y = 0.0, z = 0.0;

    // Call for loop depending on dimension
    // Start from 2 as we skip the ghost cells (both lower and upper)
    if (_parameters.geometry.dim == 2) {
        for (int j = 2; j < _localSize[1]+3; j++) {
            y = _parameters.meshsize->getPosY(2,j);
            for (int i = 2; i < _localSize[0]+3; i++) {
                _outputFile << _parameters.meshsize->getPosX(i,j);
                _outputFile << ' ' << y << ' ' << z << std::endl;
            }
        }
    } else {
        for (int k = 2; k < _localSize[2]+3; k++) {
            z = _parameters.meshsize->getPosZ(2,2,k);
            for (int j = 2; j < _localSize[1]+3; j++) {
                y = _parameters.meshsize->getPosY(2,j,k);
                for (int i = 2; i < _localSize[0]+3; i++) {
                    _outputFile << _parameters.meshsize->getPosX(i,j,k);
                    _outputFile << ' ' << y << ' ' << z << std::endl;
                }
            }
        }
    }
}

/** Opens the file for writing and returns true if it was opened otherwise false
 * returns false if file was not opened or memory was not allocated for
 * pressure and velocity streams.
 * @param timeStep used to create the file name for current time step
 */
bool VTKStencil::openFile( int timeStep ) {
    // fileName stores the file name as per XML file
    // it adds current timestep + .vtk extension to it
    std::stringstream fileName;
    fileName << _parameters.vtk.prefix << '_' << timeStep << ".vtk";

    // Open the file for writing out
    _outputFile.open(fileName.str().c_str(),std::ios::out);

    return _outputFile.is_open();
}

/** Writes the information to the file
 * @param flowField Flow field to be written
 */
void VTKStencil::write ( FlowField & flowField, int timeStep ) {
    // Write the header and coordinates
    writeHeaderAndCoords();
    // Write the header for pressure field
    _outputFile << "\nCELL_DATA " << _localSize[0]*_localSize[1]*(_localSize[2]==0?1:_localSize[2]) << std::endl;
    _outputFile << "SCALARS pressure float 1\n";
    _outputFile << "LOOKUP_TABLE default\n";
    // Write the pressure
    _outputFile << _pressureStream.str();
    // Write the header for velocity field
    _outputFile << "\nVECTORS velocity float\n";
    // Write the velocity
    _outputFile << _velocityStream.str();
    // Close the file
    _outputFile.close();

    // Clear the stringstream buffers
    _pressureStream.str(std::string());
    _velocityStream.str(std::string());
}

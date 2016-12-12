#ifndef _VTK_BINARY_TURB_STENCIL_H_
#define _VTK_BINARY_TURB_STENCIL_H_

#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../TurbFlowField.h"
#include <string>
#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>


typedef unsigned char byte;


class VTKBinaryTurbStencil : public FieldStencil<TurbFlowField> {


    public:

        /** Constructor
         *
         * @param prefix String with the prefix of the name of the VTK files
         */
        VTKBinaryTurbStencil ( const Parameters & parameters );

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

        void populateBytes(float tmp, byte* byte1, byte* byte2, byte* byte3, byte* byte4);
        ~VTKBinaryTurbStencil(); // destructor

	private:
        std::string _file_path;
		//Stores header data
    	std::string _header;
		//Stores coordinate data
		std::string _coord_data;
		//Stores cell data scalars
		std::string _cd_scalar_p;
		//Stores cell data vectors
		std::string _cd_velo;
		//Stores cell data scalars
		std::string _cd_scalar_vis;
		//Stores cell data vectors
		std::string _cd_scalar_wall;

		//buffer for coordinates
		unsigned int coord_b_len;
		byte *coord_b;

		//buffer for pressure
		unsigned int pres_b_len;
		unsigned int cur_pres_pos;
		byte *pres_b;

		//buffer for velocity
		unsigned int velo_b_len;
		unsigned int cur_velo_pos;
		byte *velo_b;

		//buffer for viscosity
		unsigned int vis_b_len;
		unsigned int cur_vis_pos;
		byte *vis_b;

		//buffer for distance to wall
		unsigned int wall_b_len;
		unsigned int cur_wall_pos;
		byte *wall_b;


};

#endif

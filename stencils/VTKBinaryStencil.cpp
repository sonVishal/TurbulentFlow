#include "VTKBinaryStencil.h"


#define VTKBINARYDEBUG

VTKBinaryStencil::VTKBinaryStencil ( const Parameters & parameters ): FieldStencil<FlowField> (parameters){
	const int dir_err = system("mkdir -p binary_vtk");
	if (-1 == dir_err) {
			perror("Cannot create vtk directory!");
			exit(1);
	}

	std::stringstream ss("");
	ss<<"binary_vtk/"<<parameters.vtk.prefix<<"_Re"<<parameters.flow.Re<<"_"<<parameters.geometry.sizeX<<"x"
			<<parameters.geometry.sizeY<<"x"<<parameters.geometry.sizeZ;

	this->_file_path = ss.str();
	ss.str("");


	int x=parameters.parallel.localSize[0];
	int y=parameters.parallel.localSize[1];
	int z=parameters.parallel.localSize[2];


	//Initializing each variable
	ss<<"# vtk DataFile Version 2.0\nMy fancy data\nBINARY\n"<<
			"DATASET STRUCTURED_GRID\nDIMENSIONS "<<(x+1)<<" "<<(y+1)<<" "<<(z+1)<<"\n"<<
			"POINTS "<<((x+1)*(y+1)*(z+1))<<" float\n";
	this->_header= ss.str();
	ss.str("");


	/*
	 * Coordinates
	 */

	//Amount of points in the simulation *3(coordinate x,y,z) *4(because a float can be stores in4 bytes)
	this->coord_b_len = (x+1)*(y+1)*(z+1)*4*3;
	//holds the coordination values
	this->coord_b = new byte[coord_b_len];


	unsigned int buf_counter=0;
	//Fill in values into buffer
	for(int k=1;k<z+2;k++){
		for(int j=1;j<y+2;j++){
			for(int i=1;i<x+2;i++){
				FLOAT x = parameters.parallel.firstCorner[0]+parameters.meshsize->getPosX(i,j,k)+parameters.meshsize->getDx(i,j,k);
				FLOAT y = parameters.parallel.firstCorner[1]+parameters.meshsize->getPosY(i,j,k)+parameters.meshsize->getDy(i,j,k);
				FLOAT z = parameters.parallel.firstCorner[2]+parameters.meshsize->getPosZ(i,j,k)+parameters.meshsize->getDz(i,j,k);

				populateBytes((float)x, &coord_b[buf_counter],&coord_b[buf_counter+1],&coord_b[buf_counter+2],&coord_b[buf_counter+3]);
				buf_counter+=4;
				populateBytes((float)y, &coord_b[buf_counter],&coord_b[buf_counter+1],&coord_b[buf_counter+2],&coord_b[buf_counter+3]);
				buf_counter+=4;
				populateBytes((float)z, &coord_b[buf_counter],&coord_b[buf_counter+1],&coord_b[buf_counter+2],&coord_b[buf_counter+3]);
				buf_counter+=4;
			}
		}
	}
	/*
	 * Pressure
	 */

	//amount of cells times 4(Byte) because of float
	this->pres_b_len = x*y*z*4;
	this->cur_pres_pos=0;
	this->pres_b = new byte[pres_b_len];


	ss << "CELL_DATA "<<(x*y*z)<<"\nSCALARS pressure float 1\nLOOKUP_TABLE default\n";
	this->_cd_scalar_p = ss.str();
	ss.str("");

	/*
	 * Velocity
	 */

	//amount of cells times 3(x,y,z) times 4(Bytes) because of float
	this->velo_b_len =x*y*z*3*4;
	this->cur_velo_pos=0;
	this->velo_b = new byte[this->velo_b_len];


	ss << "VECTORS velocity float\n";
	this->_cd_velo = ss.str();
}

void VTKBinaryStencil::apply ( FlowField & flowField, int i, int j ){

}

void VTKBinaryStencil::apply ( FlowField & flowField, int i, int j, int k ){

	if(flowField.getFlags().getValue(i,j,k)==0){
		FLOAT pres;
		FLOAT* velo;

		velo = (FLOAT*)calloc(3,sizeof(velo));
		flowField.getPressureAndVelocity(pres, velo,i,j,k);

		//pressure
		populateBytes((float)pres, &pres_b[cur_pres_pos],&pres_b[cur_pres_pos+1],&pres_b[cur_pres_pos+2],&pres_b[cur_pres_pos+3]);
		cur_pres_pos+=4;

		//velo x-direciton
		populateBytes((float)velo[0], &velo_b[cur_velo_pos],&velo_b[cur_velo_pos+1],&velo_b[cur_velo_pos+2],&velo_b[cur_velo_pos+3]);
		cur_velo_pos+=4;

		//velo y-direciton
		populateBytes((float)velo[1], &velo_b[cur_velo_pos],&velo_b[cur_velo_pos+1],&velo_b[cur_velo_pos+2],&velo_b[cur_velo_pos+3]);
		cur_velo_pos+=4;

		//velo z-direciton
		populateBytes((float)velo[2], &velo_b[cur_velo_pos],&velo_b[cur_velo_pos+1],&velo_b[cur_velo_pos+2],&velo_b[cur_velo_pos+3]);
		cur_velo_pos+=4;


	}else{
		populateBytes(0.0, &pres_b[cur_pres_pos],&pres_b[cur_pres_pos+1],&pres_b[cur_pres_pos+2],&pres_b[cur_pres_pos+3]);
		cur_pres_pos+=4;

		//velo x-direciton
		populateBytes(0.0, &velo_b[cur_velo_pos],&velo_b[cur_velo_pos+1],&velo_b[cur_velo_pos+2],&velo_b[cur_velo_pos+3]);
		cur_velo_pos+=4;

		//velo y-direciton
		populateBytes(0.0, &velo_b[cur_velo_pos],&velo_b[cur_velo_pos+1],&velo_b[cur_velo_pos+2],&velo_b[cur_velo_pos+3]);
		cur_velo_pos+=4;

		//velo z-direciton
		populateBytes(0.0, &velo_b[cur_velo_pos],&velo_b[cur_velo_pos+1],&velo_b[cur_velo_pos+2],&velo_b[cur_velo_pos+3]);
		cur_velo_pos+=4;
	}



}

void VTKBinaryStencil::write ( FlowField & flowField, int timeStep ){
	std::stringstream ss("");
	ss<<this->_file_path<<"_"<<timeStep<<".vtk";
	#ifdef VTKBINARYDEBUG
		std::cout<<"VTKStancil writer "<<ss.str().c_str()<<std::endl;
	#endif

	std::ofstream output (ss.str().c_str(), std::ios::out | std::ios::binary);
	if (output.is_open()) {

		//write header data
		output.write(this->_header.c_str(), _header.size());

		//write coordinate data
		output.write((const char *)coord_b, coord_b_len);

		//write pressure data
		output.write(this->_cd_scalar_p.c_str(), _cd_scalar_p.size());
		output.write((const char *)pres_b, pres_b_len);

		//write velocity data
		output.write(this->_cd_velo.c_str(), _cd_velo.size());
		output.write((const char *)velo_b, velo_b_len);


	}
	else std::cerr<<"Cannot open file "<<ss.str().c_str()<<" !"<<std::endl;

	output.close();

	//reset buffer
	this->cur_pres_pos=0;
	this->cur_velo_pos=0;
}

/**
 * Populate the bytes with the float
 */
inline void VTKBinaryStencil::populateBytes(float tmp, byte* byte1, byte* byte2, byte* byte3, byte* byte4){
	union{
		float tmp_float;
		byte binary[4];
	}convert;

	convert.tmp_float = tmp;

	*byte1=convert.binary[3];
	*byte2=convert.binary[2];
	*byte3=convert.binary[1];
	*byte4=convert.binary[0];
}

VTKBinaryStencil::~VTKBinaryStencil(){
	delete this->coord_b;
	delete this->pres_b;
	delete this->velo_b;
}


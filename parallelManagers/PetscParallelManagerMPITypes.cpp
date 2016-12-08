#include <parallelManagers/PetscParallelManagerMPITypes.h>

//#define DEBUG_PARMNG
//todo bug somewhere in pressure or velocity by spliting in y direction

PetscParallelManagerMPITypes::PetscParallelManagerMPITypes(FlowField & flowField, const Parameters & parameters)
: _parameters(parameters),
  _flowField(flowField){

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	_localSize = _parameters.parallel.localSize;


	int x = _localSize[0];
	int y = _localSize[1];
	int z = _localSize[2];

	int dim = _parameters.geometry.dim;



	if(dim==2){
		//Start locations for the communication
		_send_right_k = _recv_right_k = _send_left_k = _recv_left_k =
		_send_top_k = _recv_top_k = _send_bot_k = _recv_bot_k =
		_send_front_i = _send_front_j = _send_front_k=
		_recv_front_i = _recv_front_j = _recv_front_k =
		_send_back_i = _send_back_j = _send_back_k =
		_recv_back_i=  _recv_back_j = _recv_back_k = 0;

		_send_right_i = _localSize[0];
		_send_right_j = 2;
		_recv_left_i = 0;
		_recv_left_j = 2;

		_send_left_i = 2;
		_send_left_j = 2;
		_recv_right_i = _localSize[0]+2;
		_recv_right_j = 2;

		_send_top_i = 2;
		_send_top_j = _localSize[1];
		_recv_bot_i = 2;
		_recv_bot_j = 0;

		_send_bot_i = 2;
		_send_bot_j = 2;
		_recv_top_i = 2;
		_recv_top_j = _localSize[1]+2;

		/*
		 * Pressure Datatypes
		 */

		//Initialize datatype to communicate to the left and right neighbor
		//Vector to left
		MPI_Type_vector(
				y, //count how many blocks we want to send
				1, //blocklength how long
				x+3, //stride offset to the data
				MY_MPI_FLOAT,	//MPI Datatype
				&_pres_left
		);
		//vector to right
		MPI_Type_vector(
				y, //count how many blocks we want to send
				2, //blocklength how long
				x+3, //stride offset to the data
				MY_MPI_FLOAT,	//MPI Datatype
				&_pres_right
		);


		//Initialize datatype to communicate to the top and bottom
		//Vector to bot
		MPI_Type_vector(
				1, //count how many blocks we want to send
				x, //blocklength how long
				0, //stride offset to the data
				MY_MPI_FLOAT,	//MPI Datatype
				&_pres_bottom
		);
		//Vector to top
		MPI_Type_vector(
				2, //count how many blocks we want to send
				2*x, //blocklength how long
				4, //stride offset to the data
				MY_MPI_FLOAT,	//MPI Datatype
				&_pres_top
		);

		MPI_Type_vector(
				0, //count how many blocks we want to send
				0, //blocklength how long
				0, //stride offset to the data
				MY_MPI_FLOAT,	//MPI Datatype
				&_pres_front
		);
		MPI_Type_vector(
				0, //count how many blocks we want to send
				0, //blocklength how long
				0, //stride offset to the data
				MY_MPI_FLOAT,	//MPI Datatype
				&_pres_back
		);


		/*
		 * Velocity Datatypes
		 */

		//Initialize datatype to communicate to the left and right neighbor
		//Vector to left
		MPI_Type_vector(
				y, //count how many blocks we want to send
				2, //blocklength how long
				(x+3)*2, //stride offset to the data
				MY_MPI_FLOAT,	//MPI Datatype
				&_velo_left
		);
		//vector to right
		MPI_Type_vector(
				y, //count how many blocks we want to send
				4, //blocklength how long
				(x+3)*2, //stride offset to the data
				MY_MPI_FLOAT,	//MPI Datatype
				&_velo_right
		);
		//Initialize datatype to communicate to the top and bottom
		//vector to top
		MPI_Type_vector(
				1, //count how many blocks we want to send
				x*2, //blocklength how long
				0, //stride offset to the data
				MY_MPI_FLOAT,	//MPI Datatype
				&_velo_bottom
		);
		//vector to bot
		MPI_Type_vector(
				2, //count how many blocks we want to send
				x*2, //blocklength how long
				6, //stride offset to the data 4*2
				MY_MPI_FLOAT,	//MPI Datatype
				&_velo_top
		);

		MPI_Type_vector(
				0, //count how many blocks we want to send
				0, //blocklength how long
				0, //stride offset to the data
				MY_MPI_FLOAT,	//MPI Datatype
				&_velo_front
		);
		MPI_Type_vector(
				0, //count how many blocks we want to send
				0, //blocklength how long
				0, //stride offset to the data
				MY_MPI_FLOAT,	//MPI Datatype
				&_velo_back
		);

	}else if(dim==3){
		_send_right_i = x;
		_send_right_j = 0;
		_send_right_k = 0;
		_recv_left_i  = 0;
		_recv_left_j  = 0;
		_recv_left_k  = 0;

		_send_left_i  = 2;
		_send_left_j  = 0;
		_send_left_k  = 0;
		_recv_right_i = x+2;
		_recv_right_j = 0;
		_recv_right_k = 0;

		_send_top_i   = 0;
		_send_top_j   = y;
		_send_top_k   = 0;
		_recv_bot_i   = 0;
		_recv_bot_j   = 0;
		_recv_bot_k   = 0;

		_send_bot_i   = 0;
		_send_bot_j   = 2;
		_send_bot_k   = 0;
		_recv_top_i   = 0;
		_recv_top_j   = y+2;
		_recv_top_k   = 0;

		_send_front_i = 0;
		_send_front_j = 0;
		_send_front_k = 2;
		_recv_back_i = 0;
		_recv_back_j = 0;
		_recv_back_k = z+2;

		_send_back_i = 0;
		_send_back_j = 0;
		_send_back_k = z;
		_recv_front_i = 0;
		_recv_front_j = 0;
		_recv_front_k = 0;



		/*
		 * Pressure Datatypes 3D
		 */
		//Vector to right
		MPI_Type_vector(
				(y+3)*(z+3), //count how many blocks we want to send
				2, //blocklength how long
				x+3, //stride offset to the data
				MY_MPI_FLOAT,	//MPI Datatype
				&_pres_right
		);


		//Vector to left
		MPI_Type_vector(
				(y+3)*(z+3), //count how many blocks we want to send
				1, //blocklength how long
				x+3, //stride offset to the data
				MY_MPI_FLOAT,	//MPI Datatype
				&_pres_left
		);

		//Vector to bot
		MPI_Type_vector(
				z+3, //count how many blocks we want to send
				x+3, //blocklength how long
				(y+3)*(x+3), //stride offset to the data
				MY_MPI_FLOAT,	//MPI Datatype
				&_pres_bottom
		);
		//Vector to top
		MPI_Type_vector(
				z+3, //count how many blocks we want to send
				2*(x+3), //blocklength how long
				(y+3)*(x+3), //stride offset to the data
				MY_MPI_FLOAT,	//MPI Datatype
				&_pres_top
		);

		//Vector to front
		MPI_Type_vector(
				1, //count how many blocks we want to send
				(x+3)*(y+3), //blocklength how long
				0, //stride offset to the data
				MY_MPI_FLOAT,	//MPI Datatype
				&_pres_front
		);
		//Vector to back
		MPI_Type_vector(
				1, //count how many blocks we want to send
				(x+3)*(y+3)*2, //blocklength how long
				0, //stride offset to the data
				MY_MPI_FLOAT,	//MPI Datatype
				&_pres_back
		);


		/*
		 * Velocity Datatypes 3D
		 */

		//Initialize datatype to communicate to the left and right neighbor
		//Vector to left
		MPI_Type_vector(
				(y+3)*(z+3), //count how many blocks we want to send
				3, //blocklength how long
				(x+3)*3, //stride offset to the data
				MY_MPI_FLOAT,	//MPI Datatype
				&_velo_left
		);
		//vector to right
		MPI_Type_vector(
				(y+3)*(z+3), //count how many blocks we want to send
				6, //blocklength how long
				(x+3)*3, //stride offset to the data
				MY_MPI_FLOAT,	//MPI Datatype
				&_velo_right
		);


		//Initialize datatype to communicate to the top and bottom
		//vector to top
		MPI_Type_vector(
				z+3, //count how many blocks we want to send
				(x+3)*3, //blocklength how long
				(y+3)*(x+3)*3, //stride offset to the data
				MY_MPI_FLOAT,	//MPI Datatype
				&_velo_bottom
		);
		//vector to bot
		MPI_Type_vector(
				z+3, //count how many blocks we want to send
				(x+3)*6, //blocklength how long
				(y+3)*(x+3)*3, //stride offset to the data 4*2
				MY_MPI_FLOAT,	//MPI Datatype
				&_velo_top
		);

		//Vector to front
		MPI_Type_vector(
				1, //count how many blocks we want to send
				(x+3)*(y+3)*3, //blocklength how long
				0, //stride offset to the data
				MY_MPI_FLOAT,	//MPI Datatype
				&_velo_front
		);
		//Vector to back
		MPI_Type_vector(
				1, //count how many blocks we want to send
				(x+3)*(y+3)*6, //blocklength how long
				0, //stride offset to the data
				MY_MPI_FLOAT,	//MPI Datatype
				&_velo_back
		);
	}

	MPI_Type_commit(&_pres_left);
	MPI_Type_commit(&_pres_right);
	MPI_Type_commit(&_pres_bottom);
	MPI_Type_commit(&_pres_top);
	MPI_Type_commit(&_pres_front);
	MPI_Type_commit(&_pres_back);

	MPI_Type_commit(&_velo_left);
	MPI_Type_commit(&_velo_right);
	MPI_Type_commit(&_velo_bottom);
	MPI_Type_commit(&_velo_top);
	MPI_Type_commit(&_velo_front);
	MPI_Type_commit(&_velo_back);
}



PetscParallelManagerMPITypes::~PetscParallelManagerMPITypes(){}

void PetscParallelManagerMPITypes::communicatePressure(){


	MPI_Status stat;

	//Communication to the right neighbor
	MPI_Sendrecv(
			&_flowField.getPressure().getScalar(_send_right_i, _send_right_j, _send_right_k),		//sendbuffer in this case memory location
			1,			//sizeof datatype in this case 1
			_pres_right,
			_parameters.parallel.rightNb,
			COMM_RIGHT_PRES,
			&_flowField.getPressure().getScalar(_recv_left_i, _recv_left_j, _recv_left_k),		//resvbuffer in this case memory location
			1,			//size
			_pres_right,
			_parameters.parallel.leftNb,
			COMM_RIGHT_PRES, //recvtag
			MPI_COMM_WORLD,
			&stat			//MPIStatus
	);


	MPI_Sendrecv(
			&_flowField.getPressure().getScalar(_send_left_i, _send_left_j, _send_left_k),		//sendbuffer in this case memory location
			1,			//sizeof datatype in this case 1
			_pres_left,
			_parameters.parallel.leftNb,
			COMM_LEFT_PRES,
			&_flowField.getPressure().getScalar(_recv_right_i, _recv_right_j, _recv_right_k),		//resvbuffer in this case memory location
			1,			//size
			_pres_left,
			_parameters.parallel.rightNb,
			COMM_LEFT_PRES, //recvtag
			MPI_COMM_WORLD,
			&stat			//MPIStatus
	);

	MPI_Sendrecv(
			&_flowField.getPressure().getScalar(_send_top_i, _send_top_j, _send_top_k),		//sendbuffer in this case memory location
			1,			//sizeof datatype in this case 1
			_pres_top,
			_parameters.parallel.topNb,
			COMM_TOP_PRES,
			&_flowField.getPressure().getScalar(_recv_bot_i, _recv_bot_j, _recv_bot_k),		//resvbuffer in this case memory location
			1,			//size
			_pres_top,
			_parameters.parallel.bottomNb,
			COMM_TOP_PRES, //recvtag
			MPI_COMM_WORLD,
			&stat			//MPIStatus
	);


	MPI_Sendrecv(
			&_flowField.getPressure().getScalar(_send_bot_i, _send_bot_j, _send_bot_k),		//sendbuffer in this case memory location
			1,			//sizeof datatype in this case 1
			_pres_bottom,
			_parameters.parallel.bottomNb,
			COMM_BOTTOM_PRES,
			&_flowField.getPressure().getScalar(_recv_top_i, _recv_top_j, _recv_top_k),		//resvbuffer in this case memory location
			1,			//size
			_pres_bottom,
			_parameters.parallel.topNb,
			COMM_BOTTOM_PRES, //recvtag
			MPI_COMM_WORLD,
			&stat			//MPIStatus
	);


	//front back
	MPI_Sendrecv(
			&_flowField.getPressure().getScalar(_send_back_i, _send_back_j, _send_back_k),		//sendbuffer in this case memory location
			1,			//sizeof datatype in this case 1
			_pres_back,
			_parameters.parallel.backNb,
			COMM_BACK_PRES,
			&_flowField.getPressure().getScalar(_recv_front_i, _recv_front_j, _recv_front_k),		//resvbuffer in this case memory location
			1,			//size
			_pres_back,
			_parameters.parallel.frontNb,
			COMM_BACK_PRES, //recvtag
			MPI_COMM_WORLD,
			&stat			//MPIStatus
	);

	//front back
	MPI_Sendrecv(
			&_flowField.getPressure().getScalar(_send_front_i, _send_front_j, _send_front_k),		//sendbuffer in this case memory location
			1,			//sizeof datatype in this case 1
			_pres_front,
			_parameters.parallel.frontNb,
			COMM_FRONT_PRES,
			&_flowField.getPressure().getScalar(_recv_back_i, _recv_back_j, _recv_back_k),		//resvbuffer in this case memory location
			1,			//size
			_pres_front,
			_parameters.parallel.backNb,
			COMM_FRONT_PRES, //recvtag
			MPI_COMM_WORLD,
			&stat			//MPIStatus
	);
}


void PetscParallelManagerMPITypes::communicateVelocity(){

	MPI_Status stat;

	MPI_Sendrecv(
			_flowField.getVelocity().getVector(_send_right_i, _send_right_j, _send_right_k),		//sendbuffer in this case memory location
			1,			//sizeof datatype in this case 1
			_velo_right,
			_parameters.parallel.rightNb,
			COMM_RIGHT_VELO,
			_flowField.getVelocity().getVector(_recv_left_i, _recv_left_j, _recv_left_k),		//resvbuffer in this case memory location
			1,			//size
			_velo_right,
			_parameters.parallel.leftNb,
			COMM_RIGHT_VELO, //recvtag
			MPI_COMM_WORLD,
			&stat			//MPIStatus
	);

	MPI_Sendrecv(
			_flowField.getVelocity().getVector(_send_left_i, _send_left_j, _send_left_k),		//sendbuffer in this case memory location
			1,			//sizeof datatype in this case 1
			_velo_left,
			_parameters.parallel.leftNb,
			COMM_LEFT_VELO,
			_flowField.getVelocity().getVector(_recv_right_i, _recv_right_j, _recv_right_k),		//resvbuffer in this case memory location
			1,			//size
			_velo_left,
			_parameters.parallel.rightNb,
			COMM_LEFT_VELO, //recvtag
			MPI_COMM_WORLD,
			&stat			//MPIStatus
	);

	MPI_Sendrecv(
			_flowField.getVelocity().getVector(_send_top_i, _send_top_j, _send_top_k),		//sendbuffer in this case memory location
			1,			//sizeof datatype in this case 1
			_velo_top,
			_parameters.parallel.topNb,
			COMM_TOP_VELO,
			_flowField.getVelocity().getVector(_recv_bot_i, _recv_bot_j, _recv_bot_k),		//resvbuffer in this case memory location
			1,			//size
			_velo_top,
			_parameters.parallel.bottomNb,
			COMM_TOP_VELO, //recvtag
			MPI_COMM_WORLD,
			&stat			//MPIStatus
	);

	MPI_Sendrecv(
			_flowField.getVelocity().getVector(_send_bot_i, _send_bot_j, _send_bot_k),		//sendbuffer in this case memory location
			1,			//sizeof datatype in this case 1
			_velo_bottom,
			_parameters.parallel.bottomNb,
			COMM_BOTTOM_VELO,
			_flowField.getVelocity().getVector(_recv_top_i, _recv_top_j, _recv_top_k),		//resvbuffer in this case memory location
			1,			//size
			_velo_bottom,
			_parameters.parallel.topNb,
			COMM_BOTTOM_VELO, //recvtag
			MPI_COMM_WORLD,
			&stat			//MPIStatus
	);

	//front back
	MPI_Sendrecv(
			_flowField.getVelocity().getVector(_send_back_i, _send_back_j, _send_back_k),		//sendbuffer in this case memory location
			1,			//sizeof datatype in this case 1
			_velo_back,
			_parameters.parallel.backNb,
			COMM_BACK_VELO,
			_flowField.getVelocity().getVector(_recv_front_i, _recv_front_j, _recv_front_k),		//resvbuffer in this case memory location
			1,			//size
			_velo_back,
			_parameters.parallel.frontNb,
			COMM_BACK_VELO, //recvtag
			MPI_COMM_WORLD,
			&stat			//MPIStatus
	);

	//front back
	MPI_Sendrecv(
			_flowField.getVelocity().getVector(_send_front_i, _send_front_j, _send_front_k),		//sendbuffer in this case memory location
			1,			//sizeof datatype in this case 1
			_velo_front,
			_parameters.parallel.frontNb,
			COMM_FRONT_VELO,
			_flowField.getVelocity().getVector(_recv_back_i, _recv_back_j, _recv_back_k),		//resvbuffer in this case memory location
			1,			//size
			_velo_front,
			_parameters.parallel.backNb,
			COMM_FRONT_VELO, //recvtag
			MPI_COMM_WORLD,
			&stat			//MPIStatus
	);
}


void PetscParallelManagerMPITypes::communicateViscosity(){
	//NOTE! Datatype for pressure is the same for the viscosity thereful its reused!
	MPI_Status stat;

	TurbFlowField & turbField = dynamic_cast<TurbFlowField &>(_flowField);

	//Communication to the right neighbor
	MPI_Sendrecv(
			&turbField.getTurbViscosity().getScalar(_send_right_i, _send_right_j, _send_right_k),		//sendbuffer in this case memory location
			1,			//sizeof datatype in this case 1
			_pres_right,
			_parameters.parallel.rightNb,
			COMM_RIGHT_VIS,
			&turbField.getPressure().getScalar(_recv_left_i, _recv_left_j, _recv_left_k),		//resvbuffer in this case memory location
			1,			//size
			_pres_right,
			_parameters.parallel.leftNb,
			COMM_RIGHT_VIS, //recvtag
			MPI_COMM_WORLD,
			&stat			//MPIStatus
	);


	MPI_Sendrecv(
			&turbField.getPressure().getScalar(_send_left_i, _send_left_j, _send_left_k),		//sendbuffer in this case memory location
			1,			//sizeof datatype in this case 1
			_pres_left,
			_parameters.parallel.leftNb,
			COMM_LEFT_VIS,
			&turbField.getPressure().getScalar(_recv_right_i, _recv_right_j, _recv_right_k),		//resvbuffer in this case memory location
			1,			//size
			_pres_left,
			_parameters.parallel.rightNb,
			COMM_LEFT_VIS, //recvtag
			MPI_COMM_WORLD,
			&stat			//MPIStatus
	);

	MPI_Sendrecv(
			&turbField.getPressure().getScalar(_send_top_i, _send_top_j, _send_top_k),		//sendbuffer in this case memory location
			1,			//sizeof datatype in this case 1
			_pres_top,
			_parameters.parallel.topNb,
			COMM_TOP_VIS,
			&turbField.getPressure().getScalar(_recv_bot_i, _recv_bot_j, _recv_bot_k),		//resvbuffer in this case memory location
			1,			//size
			_pres_top,
			_parameters.parallel.bottomNb,
			COMM_TOP_VIS, //recvtag
			MPI_COMM_WORLD,
			&stat			//MPIStatus
	);


	MPI_Sendrecv(
			&turbField.getPressure().getScalar(_send_bot_i, _send_bot_j, _send_bot_k),		//sendbuffer in this case memory location
			1,			//sizeof datatype in this case 1
			_pres_bottom,
			_parameters.parallel.bottomNb,
			COMM_BOTTOM_VIS,
			&turbField.getPressure().getScalar(_recv_top_i, _recv_top_j, _recv_top_k),		//resvbuffer in this case memory location
			1,			//size
			_pres_bottom,
			_parameters.parallel.topNb,
			COMM_BOTTOM_VIS, //recvtag
			MPI_COMM_WORLD,
			&stat			//MPIStatus
	);


	//front back
	MPI_Sendrecv(
			&turbField.getPressure().getScalar(_send_back_i, _send_back_j, _send_back_k),		//sendbuffer in this case memory location
			1,			//sizeof datatype in this case 1
			_pres_back,
			_parameters.parallel.backNb,
			COMM_BACK_VIS,
			&turbField.getPressure().getScalar(_recv_front_i, _recv_front_j, _recv_front_k),		//resvbuffer in this case memory location
			1,			//size
			_pres_back,
			_parameters.parallel.frontNb,
			COMM_BACK_VIS, //recvtag
			MPI_COMM_WORLD,
			&stat			//MPIStatus
	);

	//front back
	MPI_Sendrecv(
			&turbField.getPressure().getScalar(_send_front_i, _send_front_j, _send_front_k),		//sendbuffer in this case memory location
			1,			//sizeof datatype in this case 1
			_pres_front,
			_parameters.parallel.frontNb,
			COMM_FRONT_VIS,
			&turbField.getPressure().getScalar(_recv_back_i, _recv_back_j, _recv_back_k),		//resvbuffer in this case memory location
			1,			//size
			_pres_front,
			_parameters.parallel.backNb,
			COMM_FRONT_VIS, //recvtag
			MPI_COMM_WORLD,
			&stat			//MPIStatus
	);
}




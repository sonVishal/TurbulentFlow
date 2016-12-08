#ifndef PARALLELMANAGERS_PETSCPARALLELMANAGERMPITYPES_H_
#define PARALLELMANAGERS_PETSCPARALLELMANAGERMPITypes_H_

#include "Iterators.h"
#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../FlowField.h"
#include "../TurbFlowField.h"

#define COMM_LEFT_PRES 0
#define COMM_RIGHT_PRES 1
#define COMM_TOP_PRES 2
#define COMM_BOTTOM_PRES 3
#define COMM_FRONT_PRES 4
#define COMM_BACK_PRES 5

#define COMM_LEFT_VELO 6
#define COMM_RIGHT_VELO 7
#define COMM_TOP_VELO 8
#define COMM_BOTTOM_VELO 9
#define COMM_FRONT_VELO 10
#define COMM_BACK_VELO 11

#define COMM_LEFT_VIS 12
#define COMM_RIGHT_VIS 13
#define COMM_TOP_VIS 14
#define COMM_BOTTOM_VIS 15
#define COMM_FRONT_VIS 16
#define COMM_BACK_VIS 17

//#define DEBUG_PARMNG


class PetscParallelManagerMPITypes {
public:
	PetscParallelManagerMPITypes(FlowField & _flowField, const Parameters & _parameters);
	virtual ~PetscParallelManagerMPITypes();
	void communicatePressure();
	void communicateVelocity();
	void communicateViscosity();


private:
	const Parameters & _parameters;
	FlowField & _flowField;

	const int* _localSize;

	int _send_right_i, _send_right_j, _send_right_k;
	int _recv_right_i, _recv_right_j, _recv_right_k;

	int _send_left_i, _send_left_j, _send_left_k;
	int _recv_left_i, _recv_left_j, _recv_left_k;

	int _send_top_i, _send_top_j, _send_top_k;
	int _recv_top_i, _recv_top_j, _recv_top_k;

	int _send_bot_i, _send_bot_j, _send_bot_k;
	int _recv_bot_i, _recv_bot_j, _recv_bot_k;

	int _send_front_i, _send_front_j, _send_front_k;
	int _recv_front_i, _recv_front_j, _recv_front_k;

	int _send_back_i, _send_back_j, _send_back_k;
	int _recv_back_i, _recv_back_j, _recv_back_k;



	//Datatype to communicate to the left and right rank
	MPI_Datatype _pres_right;
	MPI_Datatype _pres_left;

	//Datatype to communicate to the top and the bottom rank
	MPI_Datatype _pres_top;
	MPI_Datatype _pres_bottom;

	//Datatype to communicate to the front and back rank
	MPI_Datatype _pres_front;
	MPI_Datatype _pres_back;

	//Datatype to communicate to the left and right rank
	MPI_Datatype _velo_right;
	MPI_Datatype _velo_left;

	//Datatype to communicate to the top and the bottom rank
	MPI_Datatype _velo_top;
	MPI_Datatype _velo_bottom;

	//Datatype to communicate to the front and back rank
	MPI_Datatype _velo_front;
	MPI_Datatype _velo_back;

	int rank;
	int size;

};

#endif

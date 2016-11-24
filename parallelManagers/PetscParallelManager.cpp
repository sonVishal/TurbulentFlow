#include <parallelManagers/PetscParallelManager.h>

//#define DEBUG_PARMNG

PetscParallelManager::PetscParallelManager(FlowField & flowField, const Parameters & parameters)
: _parameters(parameters),
  _flowField(flowField),
  _pres_buf_fill(parameters),
  _pres_buf_read(parameters),
  _presFillIterator(flowField,parameters,_pres_buf_fill,2,-1),
  _presReadIterator(flowField,parameters,_pres_buf_read,2,-1),
  _velo_buf_fill(parameters),
  _velo_buf_read(parameters),
  _veloFillIterator(flowField,parameters,_velo_buf_fill,2,-1),
  _veloReadIterator(flowField,parameters,_velo_buf_read,2,-1){

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
}



PetscParallelManager::~PetscParallelManager(){}

void PetscParallelManager::communicatePressure(){
	#ifdef DEBUG_PARMNG
		std::cout<<"Communicate Pressure"<<std::endl;
	#endif

	//Iterate over all communicator cells
	MPI_Barrier(MPI_COMM_WORLD);
		for(int i=0;i<size;i++){
			if(i==rank){
				//std::cout<<"rank "<<rank<<std::endl;
				_presFillIterator.iterate();
			}

			MPI_Barrier(MPI_COMM_WORLD);
		}

	_pres_buf_fill.resetCounter();

	MPI_Status stat;

	//Communication to the right neighbor
	MPI_Sendrecv(
			_pres_buf_fill.getRight()->arr, //sendbuffer
			_pres_buf_fill.getRight()->size,//sizeof buffer;
			MY_MPI_FLOAT,
			_parameters.parallel.rightNb,
			COMM_RIGHT_PRES,
			_pres_buf_read.getLeft()->arr,
			_pres_buf_read.getLeft()->size,
			MY_MPI_FLOAT,
			_parameters.parallel.leftNb,
			COMM_RIGHT_PRES, //recvtag
			MPI_COMM_WORLD,
			&stat			//MPIStatus
	);
	#ifdef DEBUG_PARMNG
		MPI_Barrier(MPI_COMM_WORLD);
		for(int i=0;i<size;i++){
			if(i==rank){
				std::cout<<"Rank "<<rank<<"send to "<<_parameters.parallel.rightNb<<"\n"<<printComBuf(_pres_buf_fill.getRight())<<std::endl;
				std::cout<<"Rank "<<rank<<"received "<<_parameters.parallel.leftNb<<"\n"<<printComBuf(_pres_buf_read.getLeft())<<std::endl;
			}

			MPI_Barrier(MPI_COMM_WORLD);
		}

	#endif


	//Communication to the left neighbor
	MPI_Sendrecv(
			_pres_buf_fill.getLeft()->arr, //sendbuffer
			_pres_buf_fill.getLeft()->size,//sizeof buffer;
			MY_MPI_FLOAT,
			_parameters.parallel.leftNb,
			COMM_LEFT_PRES,
			_pres_buf_read.getRight()->arr,
			_pres_buf_read.getRight()->size,
			MY_MPI_FLOAT,
			_parameters.parallel.rightNb,
			COMM_LEFT_PRES, //recvtag
			MPI_COMM_WORLD,
			&stat			//MPIStatus
	);

	#ifdef DEBUG_PARMNG
		MPI_Barrier(MPI_COMM_WORLD);
		for(int i=0;i<size;i++){
			if(i==rank){
				std::cout<<"Rank "<<rank<<"send to "<<_parameters.parallel.leftNb<<"\n"<<printComBuf(_pres_buf_fill.getLeft())<<std::endl;
				std::cout<<"Rank "<<rank<<"received "<<_parameters.parallel.rightNb<<"\n"<<printComBuf(_pres_buf_read.getRight())<<std::endl;
			}

			MPI_Barrier(MPI_COMM_WORLD);
		}

	#endif

	//communication to the top neighbor
	MPI_Sendrecv(
			_pres_buf_fill.getTop()->arr, //sendbuffer
			_pres_buf_fill.getTop()->size,//sizeof buffer;
			MY_MPI_FLOAT,
			_parameters.parallel.topNb,
			COMM_TOP_PRES,
			_pres_buf_read.getBottom()->arr,
			_pres_buf_read.getBottom()->size,
			MY_MPI_FLOAT,
			_parameters.parallel.bottomNb,
			COMM_TOP_PRES, //recvtag
			MPI_COMM_WORLD,
			&stat			//MPIStatus
	);

	#ifdef DEBUG_PARMNG
		std::cout<<"Rank "<<rank<<"send to "<<_parameters.parallel.topNb<<"\n"<<printComBuf(_pres_buf_fill.getTop())<<std::endl;
		std::cout<<"Rank "<<rank<<"received "<<_parameters.parallel.bottomNb<<"\n"<<printComBuf(_pres_buf_read.getBottom())<<std::endl;
	#endif

	//communication to the bottom neighbor
	MPI_Sendrecv(
			_pres_buf_fill.getBottom()->arr, //sendbuffer
			_pres_buf_fill.getBottom()->size,//sizeof buffer;
			MY_MPI_FLOAT,
			_parameters.parallel.bottomNb,
			COMM_BOTTOM_PRES,
			_pres_buf_read.getTop()->arr,
			_pres_buf_read.getTop()->size,
			MY_MPI_FLOAT,
			_parameters.parallel.topNb,
			COMM_BOTTOM_PRES, //recvtag
			MPI_COMM_WORLD,
			&stat			//MPIStatus
	);

	#ifdef DEBUG_PARMNG
		std::cout<<"Rank "<<rank<<"send to "<<_parameters.parallel.bottomNb<<"\n"<<printComBuf(_pres_buf_fill.getBottom())<<std::endl;
		std::cout<<"Rank "<<rank<<"received "<<_parameters.parallel.topNb<<"\n"<<printComBuf(_pres_buf_read.getTop())<<std::endl;
	#endif




	MPI_Barrier(MPI_COMM_WORLD);
	for(int i=0;i<size;i++){
		if(i==rank){
			//std::cout<<"rank "<<rank<<std::endl;
			_presReadIterator.iterate();
		}

		MPI_Barrier(MPI_COMM_WORLD);
	}
	_pres_buf_read.resetCounter();
}


void PetscParallelManager::communicateVelocity(){
	#ifdef DEBUG_PARMNG
		std::cout<<"Communicate Velocity"<<std::endl;
	#endif


	MPI_Barrier(MPI_COMM_WORLD);
	for(int i=0;i<size;i++){
		if(i==rank){
			//std::cout<<"rank "<<rank<<std::endl;
			_veloFillIterator.iterate();
		}

		MPI_Barrier(MPI_COMM_WORLD);
	}

	_velo_buf_fill.resetCounter();

	MPI_Status stat;

	MPI_Sendrecv(
			_velo_buf_fill.getRight()->arr, //sendbuffer
			_velo_buf_fill.getRight()->size,//sizeof buffer;
			MY_MPI_FLOAT,
			_parameters.parallel.rightNb,
			COMM_RIGHT_VELO,
			_velo_buf_read.getLeft()->arr,
			_velo_buf_read.getLeft()->size,
			MY_MPI_FLOAT,
			_parameters.parallel.leftNb,
			COMM_RIGHT_VELO, //recvtag
			MPI_COMM_WORLD,
			&stat			//MPIStatus
	);

	MPI_Sendrecv(
			_velo_buf_fill.getLeft()->arr, //sendbuffer
			_velo_buf_fill.getLeft()->size,//sizeof buffer;
			MY_MPI_FLOAT,
			_parameters.parallel.leftNb,
			COMM_LEFT_VELO,
			_velo_buf_read.getRight()->arr,
			_velo_buf_read.getRight()->size,
			MY_MPI_FLOAT,
			_parameters.parallel.rightNb,
			COMM_LEFT_VELO, //recvtag
			MPI_COMM_WORLD,
			&stat			//MPIStatus
	);

	MPI_Sendrecv(
			_velo_buf_fill.getTop()->arr, //sendbuffer
			_velo_buf_fill.getTop()->size,//sizeof buffer;
			MY_MPI_FLOAT,
			_parameters.parallel.topNb,
			COMM_TOP_VELO,
			_velo_buf_read.getBottom()->arr,
			_velo_buf_read.getBottom()->size,
			MY_MPI_FLOAT,
			_parameters.parallel.bottomNb,
			COMM_TOP_VELO, //recvtag
			MPI_COMM_WORLD,
			&stat			//MPIStatus
	);

	MPI_Sendrecv(
			_velo_buf_fill.getBottom()->arr, //sendbuffer
			_velo_buf_fill.getBottom()->size,//sizeof buffer;
			MY_MPI_FLOAT,
			_parameters.parallel.bottomNb,
			COMM_BOTTOM_VELO,
			_velo_buf_read.getTop()->arr,
			_velo_buf_read.getTop()->size,
			MY_MPI_FLOAT,
			_parameters.parallel.topNb,
			COMM_BOTTOM_VELO, //recvtag
			MPI_COMM_WORLD,
			&stat			//MPIStatus
	);



	MPI_Barrier(MPI_COMM_WORLD);
	for(int i=0;i<size;i++){
		if(i==rank){
			//std::cout<<"rank "<<rank<<std::endl;
			_veloReadIterator.iterate();
		}

		MPI_Barrier(MPI_COMM_WORLD);
	}

	_velo_buf_read.resetCounter();
}




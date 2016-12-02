#ifndef PARALLELMANAGERS_PETSCPARALLELMANAGER_H_
#define PARALLELMANAGERS_PETSCPARALLELMANAGER_H_

#include "stencils/PressureBufferFillStencil.h"
#include "stencils/PressureBufferReadStencil.h"
#include "stencils/VelocityBufferFillStencil.h"
#include "stencils/VelocityBufferReadStencil.h"
#include "Iterators.h"
#include "stencils/CommBuffer.h"

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

//#define DEBUG_PARMNG


class PetscParallelManager {
public:
	PetscParallelManager(FlowField & _flowField, const Parameters & _parameters);
	virtual ~PetscParallelManager();
	void communicatePressure();
	void communicateVelocity();


private:
	const Parameters & _parameters;
	FlowField & _flowField;
	PressureBufferFillStencil _pres_buf_fill;
	PressureBufferReadStencil _pres_buf_read;
	ParallelBoundaryIterator<FlowField> _presFillIterator;
	ParallelBoundaryIterator<FlowField> _presReadIterator;

	VelocityBufferFillStencil _velo_buf_fill;
	VelocityBufferReadStencil _velo_buf_read;
	ParallelBoundaryIterator<FlowField> _veloFillIterator;
	ParallelBoundaryIterator<FlowField> _veloReadIterator;




	int rank;
	int size;

	std::string printComBuf(ComBuf* buffer){
		std::stringstream ss;
		ss<<"Buffer size "<<buffer->size<<std::endl;
		ss<<"Buffer counter"<<buffer->counter<<std::endl;

		for(unsigned int i=0;i<buffer->size;i++){
			ss<<buffer->arr[i]<<" ";
			if((i+1)%10==0) ss<<std::endl;
		}
		ss<<std::endl;

		return ss.str();
	}

};

#endif /* PARALLELMANAGERS_PETSCPARALLELMANAGER_H_ */

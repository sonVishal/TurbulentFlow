#include <stencils/PressureBufferReadStencil.h>

//#define DEBUG_PRESREAD

PressureBufferReadStencil::PressureBufferReadStencil(const Parameters & parameters) : BoundaryStencil<FlowField> ( parameters ){
	unsigned int value_size=(unsigned int)sizeof(FLOAT);
	int x = parameters.parallel.localSize[0]+3;
	int y = parameters.parallel.localSize[1]+3;
	int z = parameters.parallel.localSize[2]+3;

	//3d
	if(parameters.geometry.dim==2){
		//todo reserve exactly the buffers which are needed
		left_recv.size=2*y;
		right_recv.size=y;
		top_recv.size=x;
		bottom_recv.size=2*x;
		back_recv.size=0;
		front_recv.size=0;


		back_recv.arr=NULL;
		front_recv.arr=NULL;

	}else{
		right_recv.size=y*z;
		left_recv.size=2*y*z;

		top_recv.size=x*z;
		bottom_recv.size=2*x*z;

		front_recv.size=x*y;
		back_recv.size=2*x*y;

		front_recv.arr =(FLOAT*)calloc(front_recv.size,value_size);
		back_recv.arr =(FLOAT*)calloc(back_recv.size,value_size);

	}
	//Resets all counters
	resetCounter();

	left_recv.arr=(FLOAT*)calloc(left_recv.size,value_size);
	right_recv.arr=(FLOAT*)calloc(right_recv.size,value_size);
	top_recv.arr=(FLOAT*)calloc(top_recv.size,value_size);
	bottom_recv.arr=(FLOAT*)calloc(bottom_recv.size,value_size);
}
PressureBufferReadStencil::~PressureBufferReadStencil(){
	free(left_recv.arr);
	free(right_recv.arr);
	free(top_recv.arr);
	free(bottom_recv.arr);
	free(front_recv.arr);
	free(back_recv.arr);
}


void PressureBufferReadStencil::applyLeftWall(FlowField & flowField, int i, int j){
	#ifdef DEBUG_PRESREAD
		std::cout<<"Left Wall"<<std::endl;
		std::cout<<i<<" "<<j<<std::endl;
		std::cout<<i<<" "<<j<<std::endl;
		if(left_recv.arr==NULL) std::cerr<<"Left Wall was not initialized!"<<std::endl;
	#endif
//	flowField.getPressure().getScalar(i-2,j)=left_recv.arr[left_recv.counter++];
//	flowField.getPressure().getScalar(i-1,j)=left_recv.arr[left_recv.counter++];
		int x=i-2;
		int y=j-2;
	flowField.getPressure().getScalar(i-2,j)=left_recv.arr[x+2+y];
	x++;
	flowField.getPressure().getScalar(i-1,j)=left_recv.arr[x+2*y];
}
void PressureBufferReadStencil::applyRightWall(FlowField & flowField, int i, int j){
	#ifdef DEBUG_PRESREAD
		std::cout<<"Right Wall"<<std::endl;
		std::cout<<i<<" "<<j<<std::endl;
		if(right_recv.arr==NULL) std::cerr<<"Right Wall was not initialized!"<<std::endl;
	#endif
//	flowField.getPressure().getScalar(i+1,j)=right_recv.arr[right_recv.counter++];
	flowField.getPressure().getScalar(i+1,j)=right_recv.arr[j-2];
}
void PressureBufferReadStencil::applyBottomWall(FlowField & flowField, int i, int j){
	#ifdef DEBUG_PRESREAD
		std::cout<<"Bottom Wall"<<std::endl;
		std::cout<<i<<" "<<j<<std::endl;
		std::cout<<i<<" "<<j<<std::endl;
		if(bottom_recv.arr==NULL) std::cerr<<"Bottom Wall was not initialized!"<<std::endl;
	#endif
	//flowField.getPressure().getScalar(i,j-2)=bottom_recv.arr[bottom_recv.counter++];
	//flowField.getPressure().getScalar(i,j-1)=bottom_recv.arr[bottom_recv.counter++];
	int x=i-2;
	int y=_parameters.parallel.localSize[0];
	flowField.getPressure().getScalar(i,j-2)=bottom_recv.arr[x];
	flowField.getPressure().getScalar(i,j-1)=bottom_recv.arr[x+y];
}
void PressureBufferReadStencil::applyTopWall(FlowField & flowField, int i, int j){
	#ifdef DEBUG_PRESREAD
		std::cout<<"Top Wall"<<std::endl;
		std::cout<<i<<" "<<j<<std::endl;
		if(top_recv.arr==NULL) std::cerr<<"Top Wall was not initialized!"<<std::endl;
	#endif
	//flowField.getPressure().getScalar(i,j+1)=top_recv.arr[top_recv.counter++];
	flowField.getPressure().getScalar(i,j+1)=top_recv.arr[i-2];
}


void PressureBufferReadStencil::applyLeftWall   (FlowField & flowField, int i, int j, int k){
	#ifdef DEBUG_PRESREAD
		std::cout<<"Left Wall"<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k<<std::endl;
		if(left_recv.arr==NULL) std::cerr<<"Left Wall was not initialized!"<<std::endl;
	#endif
	flowField.getPressure().getScalar(i-2,j,k)=left_recv.arr[left_recv.counter++];
	flowField.getPressure().getScalar(i-2,j,k)=left_recv.arr[left_recv.counter++];
}
void PressureBufferReadStencil::applyRightWall  (FlowField & flowField, int i, int j, int k){
	#ifdef DEBUG_PRESREAD
		std::cout<<"Right Wall"<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k<<std::endl;
		if(right_recv.arr==NULL) std::cerr<<"Right Wall was not initialized!"<<std::endl;
	#endif
	flowField.getPressure().getScalar(i+1,j,k)=right_recv.arr[right_recv.counter++];
}
void PressureBufferReadStencil::applyBottomWall (FlowField & flowField, int i, int j, int k){
	#ifdef DEBUG_PRESREAD
		std::cout<<"Bottom Wall"<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k<<std::endl;
		if(bottom_recv.arr==NULL) std::cerr<<"Bottom Wall was not initialized!"<<std::endl;
	#endif
	flowField.getPressure().getScalar(i,j-2,k)=bottom_recv.arr[bottom_recv.counter++];
	flowField.getPressure().getScalar(i,j-1,k)=bottom_recv.arr[bottom_recv.counter++];
}
void PressureBufferReadStencil::applyTopWall    (FlowField & flowField, int i, int j, int k){
	#ifdef DEBUG_PRESREAD
		std::cout<<"Top Wall"<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k<<std::endl;
		if(top_recv.arr==NULL) std::cerr<<"Top Wall was not initialized!"<<std::endl;
	#endif
	flowField.getPressure().getScalar(i,j+1,k)=top_recv.arr[top_recv.counter++];
}
void PressureBufferReadStencil::applyFrontWall  (FlowField & flowField, int i, int j, int k){
	#ifdef DEBUG_PRESREAD
		std::cout<<"Front Wall"<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k<<std::endl;
		if(front_recv.arr==NULL) std::cerr<<"Front Wall was not initialized!"<<std::endl;
	#endif
	flowField.getPressure().getScalar(i,j,k-2)=bottom_recv.arr[bottom_recv.counter++];
	flowField.getPressure().getScalar(i,j,k-1)=bottom_recv.arr[bottom_recv.counter++];
}
void PressureBufferReadStencil::applyBackWall   (FlowField & flowField, int i, int j, int k){
	#ifdef DEBUG_PRESREAD
		std::cout<<"Back Wall"<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k<<std::endl;
		if(back_recv.arr==NULL) std::cerr<<"Back Wall was not initialized!"<<std::endl;
	#endif
	flowField.getPressure().getScalar(i,j,k) = top_recv.arr[top_recv.counter++];
}

/*
 * This functions needs to be triggered when all data is stored in the buffers
 */
void PressureBufferReadStencil::resetCounter(){
	left_recv.counter=0;
	right_recv.counter=0;
	top_recv.counter=0;
	bottom_recv.counter=0;
	front_recv.counter=0;
	back_recv.counter=0;
}

ComBuf* PressureBufferReadStencil::getLeft(){return &left_recv;}
ComBuf* PressureBufferReadStencil::getRight(){return &right_recv;}
ComBuf* PressureBufferReadStencil::getTop(){return &top_recv;}
ComBuf* PressureBufferReadStencil::getBottom(){return &bottom_recv;}
ComBuf* PressureBufferReadStencil::getFront(){return &front_recv;}
ComBuf* PressureBufferReadStencil::getBack(){return &back_recv;}

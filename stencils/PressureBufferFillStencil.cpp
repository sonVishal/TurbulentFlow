#include <stencils/PressureBufferFillStencil.h>

//#define DEBUG_PRESFILL

PressureBufferFillStencil::PressureBufferFillStencil ( const Parameters & parameters ) : BoundaryStencil<FlowField> ( parameters ){
	unsigned int value_size=(unsigned int)sizeof(FLOAT);
	int x = parameters.parallel.localSize[0]+3;
	int y = parameters.parallel.localSize[1]+3;
	int z = parameters.parallel.localSize[2]+3;

	if(parameters.geometry.dim==2){
		//todo reserve exactly the buffers which are needed
		//Set the size of the arrays for 2d
		left_send.size = y;
		right_send.size = 2*y;
		top_send.size = 2*x;
		bottom_send.size = x;
		front_send.size = 0;
		back_send.size = 0;

		front_send.arr=NULL;
		back_send.arr=NULL;

	}else{
		//Set the size of the arrays for 2d
		left_send.size = y*z;
		right_send.size = 2*y*z;
		top_send.size = 2*x*z;
		bottom_send.size = x*z;
		front_send.size = 2*x*y;
		back_send.size = x*y;


		front_send.arr=(FLOAT*)calloc(front_send.size, value_size);
		back_send.arr=(FLOAT*)calloc(back_send.size, value_size);
	}

	//Reset all counters
	resetCounter();

	left_send.arr=(FLOAT*)calloc(left_send.size, value_size);
	right_send.arr=(FLOAT*)calloc(right_send.size, value_size);
	top_send.arr=(FLOAT*)calloc(top_send.size, value_size);
	bottom_send.arr=(FLOAT*)calloc(bottom_send.size, value_size);
}

PressureBufferFillStencil::~PressureBufferFillStencil(){

	free(left_send.arr); // sends one columns
	free(right_send.arr); // sends two columns
	free(top_send.arr); //sends one row
	free(bottom_send.arr); //sends two rows

	//3d
	free(front_send.arr);
	free(back_send.arr);

}


void PressureBufferFillStencil::applyLeftWall   (FlowField & flowField, int i, int j){
	#ifdef DEBUG_PRESFILL
		std::cout<<"Left Wall"<<std::endl;
		std::cout<<i<<" "<<j<<std::endl;
	#endif
	//left_send.arr[left_send.counter++]=flowField.getPressure().getScalar(i,j);
	left_send.arr[j-2]=flowField.getPressure().getScalar(i,j);
}
void PressureBufferFillStencil::applyRightWall  (FlowField & flowField, int i, int j){
	#ifdef DEBUG_PRESFILL
		std::cout<<"Right Wall"<<std::endl;
		std::cout<<i-1<<" "<<j<<std::endl;
		std::cout<<i<<" "<<j<<std::endl;
#endif
//		right_send.arr[right_send.counter++]=flowField.getPressure().getScalar(i-1,j);
//		right_send.arr[right_send.counter++]=flowField.getPressure().getScalar(i,j);

	int x=i-_parameters.parallel.localSize[0];
	int y=j-2;
	right_send.arr[2*y+x]=flowField.getPressure().getScalar(i-1,j);
	x++;
	right_send.arr[2*y+x]=flowField.getPressure().getScalar(i,j);
}
void PressureBufferFillStencil::applyBottomWall (FlowField & flowField, int i, int j){
	#ifdef DEBUG_PRESFILL
		std::cout<<"Bottom Wall"<<std::endl;
		std::cout<<i<<" "<<j<<std::endl;
	#endif
//		bottom_send.arr[bottom_send.counter++]=flowField.getPressure().getScalar(i,j);

		bottom_send.arr[i-2]=flowField.getPressure().getScalar(i,j);
}
void PressureBufferFillStencil::applyTopWall    (FlowField & flowField, int i, int j){
	#ifdef DEBUG_PRESFILL
		std::cout<<"Top Wall"<<std::endl;
		std::cout<<i<<" "<<j-1<<std::endl;
		std::cout<<i<<" "<<j<<std::endl;
	#endif
//	top_send.arr[top_send.counter++]=flowField.getPressure().getScalar(i,j-1);
//	top_send.arr[top_send.counter++]=flowField.getPressure().getScalar(i,j);
	int x=i-2;
	int y=_parameters.parallel.localSize[0];
	top_send.arr[x]=flowField.getPressure().getScalar(i,j-1);

	top_send.arr[x+y]=flowField.getPressure().getScalar(i,j);

}


void PressureBufferFillStencil::applyLeftWall   (FlowField & flowField, int i, int j, int k){
	#ifdef DEBUG_PRESFILL
		std::cout<<"Left Wall"<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k<<std::endl;
	#endif
	left_send.arr[left_send.counter++]=flowField.getPressure().getScalar(i,j,k);
}
void PressureBufferFillStencil::applyRightWall  (FlowField & flowField, int i, int j, int k){
	#ifdef DEBUG_PRESFILL
		std::cout<<"Right Wall"<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k<<std::endl;
	#endif
	right_send.arr[right_send.counter++]=flowField.getPressure().getScalar(i-1,j,k);
	right_send.arr[right_send.counter++]=flowField.getPressure().getScalar(i,j,k);
}
void PressureBufferFillStencil::applyBottomWall (FlowField & flowField, int i, int j, int k){
	#ifdef DEBUG_PRESFILL
		std::cout<<"Bottom Wall"<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k<<std::endl;
	#endif
	bottom_send.arr[bottom_send.counter++]=flowField.getPressure().getScalar(i,j,k);
}
void PressureBufferFillStencil::applyTopWall    (FlowField & flowField, int i, int j, int k){
	#ifdef DEBUG_PRESFILL
		std::cout<<"Top Wall"<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k<<std::endl;
	#endif
	top_send.arr[top_send.counter++]=flowField.getPressure().getScalar(i,j-1,k);
	top_send.arr[top_send.counter++]=flowField.getPressure().getScalar(i,j,k);
}
void PressureBufferFillStencil::applyFrontWall  (FlowField & flowField, int i, int j, int k){
	#ifdef DEBUG_PRESFILL
		std::cout<<"Front Wall"<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k<<std::endl;
	#endif
	front_send.arr[front_send.counter++]=flowField.getPressure().getScalar(i,j,k);
}
void PressureBufferFillStencil::applyBackWall   (FlowField & flowField, int i, int j, int k){
	#ifdef DEBUG_PRESFILL
		std::cout<<"Back Wall"<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k<<std::endl;
	#endif
	back_send.arr[back_send.counter++]=flowField.getPressure().getScalar(i,j,k-1);
	back_send.arr[back_send.counter++]=flowField.getPressure().getScalar(i,j,k);
}

/*
 * This functions needs to be triggered when all data is stored in the buffers
 */
void PressureBufferFillStencil::resetCounter(){
	left_send.counter=0;
	right_send.counter=0;
	top_send.counter=0;
	bottom_send.counter=0;
	front_send.counter=0;
	back_send.counter=0;
}


ComBuf* PressureBufferFillStencil::getLeft(){return &left_send;}
ComBuf* PressureBufferFillStencil::getRight(){return &right_send;}
ComBuf* PressureBufferFillStencil::getTop(){return &top_send;}
ComBuf* PressureBufferFillStencil::getBottom(){return &bottom_send;}
ComBuf* PressureBufferFillStencil::getFront(){return &front_send;}
ComBuf* PressureBufferFillStencil::getBack(){return &back_send;}

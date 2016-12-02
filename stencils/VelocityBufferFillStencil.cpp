#include <stencils/VelocityBufferFillStencil.h>

//#define DEBUG_VELOFILL

VelocityBufferFillStencil::VelocityBufferFillStencil ( const Parameters & parameters ) : BoundaryStencil<FlowField> ( parameters ){
	unsigned int value_size=(unsigned int)sizeof(FLOAT);
	int x = parameters.parallel.localSize[0]+3;
	int y = parameters.parallel.localSize[1]+3;
	int z = parameters.parallel.localSize[2]+3;

	if(parameters.geometry.dim==2){
		//todo reserve exactly the buffers which are needed
		//Set the size of the arrays for 2d
		left_send.size = 2*y;
		right_send.size = 2*2*y;
		top_send.size = 2*2*x;
		bottom_send.size = 2*x;
		front_send.size = 0;
		back_send.size = 0;

		front_send.arr=NULL;
		back_send.arr=NULL;

	}else{
		//Set the size of the arrays for 2d
		left_send.size = 3*y*z;
		right_send.size = 3*2*y*z;
		top_send.size = 3*2*x*z;
		bottom_send.size = 3*x*z;
		front_send.size = 3*2*x*y;
		back_send.size = 3*x*y;


		front_send.arr=(FLOAT*)calloc(front_send.size,value_size);
		back_send.arr=(FLOAT*)calloc(back_send.size,value_size);
	}

	//Reset all counters
	resetCounter();


	left_send.arr=(FLOAT*)calloc(left_send.size,value_size);
	right_send.arr=(FLOAT*)calloc(right_send.size,value_size);
	top_send.arr=(FLOAT*)calloc(top_send.size,value_size);
	bottom_send.arr=(FLOAT*)calloc(bottom_send.size,value_size);
}

VelocityBufferFillStencil::~VelocityBufferFillStencil(){

	free(left_send.arr); // sends one columns
	free(right_send.arr); // sends two columns
	free(top_send.arr); //sends one row
	free(bottom_send.arr); //sends two rows

	//3d
	free(front_send.arr);
	free(back_send.arr);

}


void VelocityBufferFillStencil::applyLeftWall   (FlowField & flowField, int i, int j){
	#ifdef DEBUG_VELOFILL
		std::cout<<"Left Wall"<<std::endl;
		std::cout<<i<<" "<<j<<std::endl;
	#endif

	left_send.arr[left_send.counter++]=flowField.getVelocity().getVector(i,j)[0];
	left_send.arr[left_send.counter++]=flowField.getVelocity().getVector(i,j)[1];
}
void VelocityBufferFillStencil::applyRightWall  (FlowField & flowField, int i, int j){
	#ifdef DEBUG_VELOFILL
		std::cout<<"Right Wall"<<std::endl;
		std::cout<<i-1<<" "<<j<<std::endl;
		std::cout<<i<<" "<<j<<std::endl;
	#endif
	right_send.arr[right_send.counter++]=flowField.getVelocity().getVector(i-1,j)[0];
	right_send.arr[right_send.counter++]=flowField.getVelocity().getVector(i-1,j)[1];
	right_send.arr[right_send.counter++]=flowField.getVelocity().getVector(i,j)[0];
	right_send.arr[right_send.counter++]=flowField.getVelocity().getVector(i,j)[1];


}
void VelocityBufferFillStencil::applyBottomWall (FlowField & flowField, int i, int j){
	#ifdef DEBUG_VELOFILL
		std::cout<<"Bottom Wall"<<std::endl;
		std::cout<<i<<" "<<j<<std::endl;
	#endif

	bottom_send.arr[bottom_send.counter++]=flowField.getVelocity().getVector(i,j)[0];
	bottom_send.arr[bottom_send.counter++]=flowField.getVelocity().getVector(i,j)[1];


}
void VelocityBufferFillStencil::applyTopWall    (FlowField & flowField, int i, int j){
	#ifdef DEBUG_VELOFILL
		std::cout<<"Top Wall"<<std::endl;
		std::cout<<i<<" "<<j-1<<std::endl;
		std::cout<<i<<" "<<j<<std::endl;
	#endif
	top_send.arr[top_send.counter++]=flowField.getVelocity().getVector(i,j-1)[0];
	top_send.arr[top_send.counter++]=flowField.getVelocity().getVector(i,j-1)[1];
	top_send.arr[top_send.counter++]=flowField.getVelocity().getVector(i,j)[0];
	top_send.arr[top_send.counter++]=flowField.getVelocity().getVector(i,j)[1];

}


void VelocityBufferFillStencil::applyLeftWall   (FlowField & flowField, int i, int j, int k){
	#ifdef DEBUG_VELOFILL
		std::cout<<"Left Wall"<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k<<std::endl;
	#endif
	left_send.arr[left_send.counter++]=flowField.getVelocity().getVector(i,j,k)[0];
	left_send.arr[left_send.counter++]=flowField.getVelocity().getVector(i,j,k)[1];
}
void VelocityBufferFillStencil::applyRightWall  (FlowField & flowField, int i, int j, int k){
	#ifdef DEBUG_VELOFILL
		std::cout<<"Right Wall"<<std::endl;
		std::cout<<i-1<<" "<<j<<" k "<<k<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k<<std::endl;
	#endif
	right_send.arr[right_send.counter++]=flowField.getVelocity().getVector(i-1,j,k)[0];
	right_send.arr[right_send.counter++]=flowField.getVelocity().getVector(i-1,j,k)[1];
	right_send.arr[right_send.counter++]=flowField.getVelocity().getVector(i,j,k)[0];
	right_send.arr[right_send.counter++]=flowField.getVelocity().getVector(i,j,k)[1];
}
void VelocityBufferFillStencil::applyBottomWall (FlowField & flowField, int i, int j, int k){
	#ifdef DEBUG_VELOFILL
		std::cout<<"Bottom Wall"<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k<<std::endl;
	#endif
	//bottom_send.arr[bottom_send.counter++]=flowField.getPressure().getScalar(i,j,k);
	bottom_send.arr[bottom_send.counter++]=flowField.getVelocity().getVector(i,j,k)[0];
	bottom_send.arr[bottom_send.counter++]=flowField.getVelocity().getVector(i,j,k)[1];
}
void VelocityBufferFillStencil::applyTopWall    (FlowField & flowField, int i, int j, int k){
	#ifdef DEBUG_VELOFILL
		std::cout<<"Top Wall"<<std::endl;
		std::cout<<i<<" "<<j-1<<" "<<k<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k<<std::endl;
	#endif
	//top_send.arr[top_send.counter++]=flowField.getPressure().getScalar(i,j-1,k);
	//top_send.arr[top_send.counter++]=flowField.getPressure().getScalar(i,j,k);
	top_send.arr[top_send.counter++]=flowField.getVelocity().getVector(i,j-1,k)[0];
	top_send.arr[top_send.counter++]=flowField.getVelocity().getVector(i,j-1,k)[1];
	top_send.arr[top_send.counter++]=flowField.getVelocity().getVector(i,j,k)[0];
	top_send.arr[top_send.counter++]=flowField.getVelocity().getVector(i,j,k)[1];
}
void VelocityBufferFillStencil::applyFrontWall  (FlowField & flowField, int i, int j, int k){
	#ifdef DEBUG_VELOFILL
		std::cout<<"Front Wall"<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k-1<<std::endl;
	#endif
	front_send.arr[front_send.counter++]=flowField.getVelocity().getVector(i,j,k-1)[0];
	front_send.arr[front_send.counter++]=flowField.getVelocity().getVector(i,j,k-1)[1];
}
void VelocityBufferFillStencil::applyBackWall   (FlowField & flowField, int i, int j, int k){
	#ifdef DEBUG_VELOFILL
		std::cout<<"Back Wall"<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k-1<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k<<std::endl;
	#endif
	back_send.arr[back_send.counter++]=flowField.getVelocity().getVector(i,j,k-1)[0];
	back_send.arr[back_send.counter++]=flowField.getVelocity().getVector(i,j,k-1)[1];
	back_send.arr[back_send.counter++]=flowField.getVelocity().getVector(i,j,k)[0];
	back_send.arr[back_send.counter++]=flowField.getVelocity().getVector(i,j,k)[1];
}

/*
 * This functions needs to be triggered when all data is stored in the buffers
 */
void VelocityBufferFillStencil::resetCounter(){
	left_send.counter=0;
	right_send.counter=0;
	top_send.counter=0;
	bottom_send.counter=0;
	front_send.counter=0;
	back_send.counter=0;
}


ComBuf* VelocityBufferFillStencil::getLeft(){return &left_send;}
ComBuf* VelocityBufferFillStencil::getRight(){return &right_send;}
ComBuf* VelocityBufferFillStencil::getTop(){return &top_send;}
ComBuf* VelocityBufferFillStencil::getBottom(){return &bottom_send;}
ComBuf* VelocityBufferFillStencil::getFront(){return &front_send;}
ComBuf* VelocityBufferFillStencil::getBack(){return &back_send;}

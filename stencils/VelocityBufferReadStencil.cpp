#include <stencils/VelocityBufferReadStencil.h>

//#define DEBUG_VELOREAD

/**
 * By using this class you need to send the buffers to the class manually. Additionally you need to reset the buffer
 * counters by executing resetCounter!
 */
VelocityBufferReadStencil::VelocityBufferReadStencil(const Parameters & parameters) : BoundaryStencil<FlowField> ( parameters ){
	unsigned int value_size=(unsigned int)sizeof(FLOAT);
	int x = parameters.parallel.localSize[0];
	int y = parameters.parallel.localSize[1];
	int z = parameters.parallel.localSize[2];

	//3d
	if(parameters.geometry.dim==2){
		//todo reserve exactly the buffers which are needed
		left_recv.size=2*2*y;
		right_recv.size=2*y;
		top_recv.size=2*x;
		bottom_recv.size=2*2*x;
		back_recv.size=0;
		front_recv.size=0;


		back_recv.arr=NULL;
		front_recv.arr=NULL;

	}else{
		std::cerr<<"Not implemented yet!"<<std::endl;
		system("1");
		right_recv.size=3*y*z;
		left_recv.size=3*2*y*z;

		top_recv.size=3*x*z;
		bottom_recv.size=3*2*x*z;

		front_recv.size=3*x*y;
		back_recv.size=3*2*x*y;

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
VelocityBufferReadStencil::~VelocityBufferReadStencil(){
	free(left_recv.arr);
	free(right_recv.arr);
	free(top_recv.arr);
	free(bottom_recv.arr);
	free(front_recv.arr);
	free(back_recv.arr);
}


void VelocityBufferReadStencil::applyLeftWall(FlowField & flowField, int i, int j){
	if(j==1) return;
	#ifdef DEBUG_VELOREAD
		std::cout<<"Left Wall"<<std::endl;
		std::cout<<i-2<<" "<<j<<std::endl;
		std::cout<<i-1<<" "<<j<<std::endl;
		if(left_recv.arr==NULL) std::cerr<<"Left Wall was not initialized!"<<std::endl;
	#endif

	flowField.getVelocity().getVector(i-2,j)[0]=left_recv.arr[left_recv.counter++];
	flowField.getVelocity().getVector(i-2,j)[1]=left_recv.arr[left_recv.counter++];

	flowField.getVelocity().getVector(i-1,j)[0]=left_recv.arr[left_recv.counter++];
	flowField.getVelocity().getVector(i-1,j)[1]=left_recv.arr[left_recv.counter++];


}
void VelocityBufferReadStencil::applyRightWall(FlowField & flowField, int i, int j){
	if(j==1) return;
	#ifdef DEBUG_VELOREAD
		std::cout<<"Right Wall"<<std::endl;
		std::cout<<i+1<<" "<<j<<std::endl;
		if(right_recv.arr==NULL) std::cerr<<"Right Wall was not initialized!"<<std::endl;
	#endif
	flowField.getVelocity().getVector(i+1,j)[0]=right_recv.arr[right_recv.counter++];
	flowField.getVelocity().getVector(i+1,j)[1]=right_recv.arr[right_recv.counter++];
}
void VelocityBufferReadStencil::applyBottomWall(FlowField & flowField, int i, int j){
	#ifdef DEBUG_VELOREAD
		std::cout<<"Bottom Wall"<<std::endl;
		std::cout<<i<<" "<<j-2<<std::endl;
		std::cout<<i<<" "<<j-1<<std::endl;
		if(bottom_recv.arr==NULL) std::cerr<<"Bottom Wall was not initialized!"<<std::endl;
	#endif
	flowField.getVelocity().getVector(i,j-2)[0]=bottom_recv.arr[bottom_recv.counter++];
	flowField.getVelocity().getVector(i,j-2)[1]=bottom_recv.arr[bottom_recv.counter++];
	flowField.getVelocity().getVector(i,j-1)[0]=bottom_recv.arr[bottom_recv.counter++];
	flowField.getVelocity().getVector(i,j-1)[1]=bottom_recv.arr[bottom_recv.counter++];

}
void VelocityBufferReadStencil::applyTopWall(FlowField & flowField, int i, int j){
	#ifdef DEBUG_VELOREAD
		std::cout<<"Top Wall"<<std::endl;
		std::cout<<i<<" "<<j+1<<std::endl;
		if(top_recv.arr==NULL) std::cerr<<"Top Wall was not initialized!"<<std::endl;
	#endif
	flowField.getVelocity().getVector(i,j+1)[0]=top_recv.arr[top_recv.counter++];
	flowField.getVelocity().getVector(i,j+1)[1]=top_recv.arr[top_recv.counter++];
}






void VelocityBufferReadStencil::applyLeftWall   (FlowField & flowField, int i, int j, int k){
	std::cerr<<"VelocityBufferReadStencil not implemented"<<std::endl;
	#ifdef DEBUG_VELOREAD
		std::cout<<"Left Wall"<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k<<std::endl;
		if(left_recv.arr==NULL) std::cerr<<"Left Wall was not initialized!"<<std::endl;
	#endif

}
void VelocityBufferReadStencil::applyRightWall  (FlowField & flowField, int i, int j, int k){
	std::cerr<<"VelocityBufferReadStencil not implemented"<<std::endl;
	#ifdef DEBUG_VELOREAD
		std::cout<<"Right Wall"<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k<<std::endl;
		if(right_recv.arr==NULL) std::cerr<<"Right Wall was not initialized!"<<std::endl;
	#endif
}
void VelocityBufferReadStencil::applyBottomWall (FlowField & flowField, int i, int j, int k){
	std::cerr<<"VelocityBufferReadStencil not implemented"<<std::endl;
	#ifdef DEBUG_VELOREAD
		std::cout<<"Bottom Wall"<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k<<std::endl;
		if(bottom_recv.arr==NULL) std::cerr<<"Bottom Wall was not initialized!"<<std::endl;
	#endif
}


void VelocityBufferReadStencil::applyTopWall    (FlowField & flowField, int i, int j, int k){
	std::cerr<<"VelocityBufferReadStencil not implemented"<<std::endl;
	#ifdef DEBUG_VELOREAD
		std::cout<<"Top Wall"<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k<<std::endl;
		if(top_recv.arr==NULL) std::cerr<<"Top Wall was not initialized!"<<std::endl;
	#endif
}
void VelocityBufferReadStencil::applyFrontWall  (FlowField & flowField, int i, int j, int k){
	std::cerr<<"VelocityBufferReadStencil not implemented"<<std::endl;
	#ifdef DEBUG_VELOREAD
		std::cout<<"Front Wall"<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k<<std::endl;
		if(front_recv.arr==NULL) std::cerr<<"Front Wall was not initialized!"<<std::endl;
	#endif
}
void VelocityBufferReadStencil::applyBackWall   (FlowField & flowField, int i, int j, int k){
	std::cerr<<"VelocityBufferReadStencil not implemented"<<std::endl;
	#ifdef DEBUG_VELOREAD
		std::cout<<"Back Wall"<<std::endl;
		std::cout<<i<<" "<<j<<" "<<k<<std::endl;
		if(back_recv.arr==NULL) std::cerr<<"Back Wall was not initialized!"<<std::endl;
	#endif
}

/*
 * This functions needs to be triggered when all data is stored in the buffers
 */
void VelocityBufferReadStencil::resetCounter(){
	left_recv.counter=0;
	right_recv.counter=0;
	top_recv.counter=0;
	bottom_recv.counter=0;
	front_recv.counter=0;
	back_recv.counter=0;
}

ComBuf* VelocityBufferReadStencil::getLeft(){return &left_recv;}
ComBuf* VelocityBufferReadStencil::getRight(){return &right_recv;}
ComBuf* VelocityBufferReadStencil::getTop(){return &top_recv;}
ComBuf* VelocityBufferReadStencil::getBottom(){return &bottom_recv;}
ComBuf* VelocityBufferReadStencil::getFront(){return &front_recv;}
ComBuf* VelocityBufferReadStencil::getBack(){return &back_recv;}

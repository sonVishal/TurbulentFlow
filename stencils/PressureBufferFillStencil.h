#ifndef _PRESSURE_BUFFER_FILL_STENCIL_H_
#define _PRESSURE_BUFFER_FILL_STENCIL_H_

#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../FlowField.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "CommBuffer.h"

//#define DEBUG_PRESFILL

class PressureBufferFillStencil : public BoundaryStencil<FlowField> {
    private:
	/*
	   Rank0   Rank1
	  gggggggg gggggggg
	  gscccccg gscccccg
	  gscccccg gscccccg
	  gscccccg gscccccg
	  gscccccg gscccccg
	  gscccccg gscccccg
	  gssssssg gssssssg
	  gggggggg gggggggg
           xx-->>
            <---x

     for 3d the communication layer is on the front side
	 */
	//Buffers for each of the six directions
	ComBuf left_send;
	ComBuf right_send;
	ComBuf top_send;
	ComBuf bottom_send;
	ComBuf front_send;
	ComBuf back_send;


	void markPressureAndVelo(FlowField & flowField, int i, int j);

    public:
		PressureBufferFillStencil( const Parameters& parameters );
		~PressureBufferFillStencil();

		void resetCounter();

		/** Represents an operation in the left wall of a 2D domain.
		 *
		 * @param flowField State of the flow field
		 * @param i Index in the x direction
		 * @param j Index in the y direction
		 */
		void applyLeftWall   (FlowField & flowField, int i, int j);

		/** Represents an operation in the right wall of a 2D domain.
		 *
		 * @param flowField State of the flow field
		 * @param i Index in the x direction
		 * @param j Index in the y direction
		 */
		void applyRightWall  (FlowField & flowField, int i, int j);

		/** Represents an operation in the bottom wall of a 2D domain.
		 *
		 * @param flowField State of the flow field
		 * @param i Index in the x direction
		 * @param j Index in the y direction
		 */
		void applyBottomWall (FlowField & flowField, int i, int j);

		/** Represents an operation in the top wall of a 2D domain.
		 *
		 * @param flowField State of the flow field
		 * @param i Index in the x direction
		 * @param j Index in the y direction
		 */
		void applyTopWall    (FlowField & flowField, int i, int j);


		/** Represents an operation in the left wall of a 3D domain.
		 *
		 * @param flowField State of the flow field
		 * @param i Index in the x direction
		 * @param j Index in the y direction
		 * @param k Index in the z direction
		 */
		void applyLeftWall   (FlowField & flowField, int i, int j, int k);

		/** Represents an operation in the right wall of a 3D domain.
		 *
		 * @param flowField State of the flow field
		 * @param i Index in the x direction
		 * @param j Index in the y direction
		 * @param k Index in the z direction
		 */
		void applyRightWall  (FlowField & flowField, int i, int j, int k);

		/** Represents an operation in the bottom wall of a 3D domain.
		 *
		 * @param flowField State of the flow field
		 * @param i Index in the x direction
		 * @param j Index in the y direction
		 * @param k Index in the z direction
		 */
		void applyBottomWall (FlowField & flowField, int i, int j, int k);

		/** Represents an operation in the top wall of a 3D domain.
		 *
		 * @param flowField State of the flow field
		 * @param i Index in the x direction
		 * @param j Index in the y direction
		 * @param k Index in the z direction
		 */
		void applyTopWall    (FlowField & flowField, int i, int j, int k);

		/** Represents an operation in the front wall of a 3D domain.
		 *
		 * @param flowField State of the flow field
		 * @param i Index in the x direction
		 * @param j Index in the y direction
		 * @param k Index in the z direction
		 */
		void applyFrontWall  (FlowField & flowField, int i, int j, int k);

		/** Represents an operation in the back wall of a 3D domain.
		 *
		 * @param flowField State of the flow field
		 * @param i Index in the x direction
		 * @param j Index in the y direction
		 * @param k Index in the z direction
		 */
		void applyBackWall   (FlowField & flowField, int i, int j, int k);

		ComBuf* getLeft();
		ComBuf* getRight();
		ComBuf* getTop();
		ComBuf* getBottom();
		ComBuf* getFront();
		ComBuf* getBack();


};

#endif

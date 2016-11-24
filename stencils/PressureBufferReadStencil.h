#ifndef _PRESSURE_BUFFER_READ_STENCIL_H_
#define _PRESSURE_BUFFER_READ_STENCIL_H_

#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../FlowField.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "CommBuffer.h"


//#define DEBUG_PRESREAD

class PressureBufferReadStencil : public BoundaryStencil<FlowField> {
    private:
	//Buffer for 2d case also used in 3d
	ComBuf left_recv; //receives 2 column
	ComBuf right_recv; //receives 1 column
	ComBuf top_recv; //receives 1 row
	ComBuf bottom_recv; //receives 2 rows

	//3d
	ComBuf front_recv; //receives 2 rows
	ComBuf  back_recv; //receives 1 row

	/*
	FLOAT* left_recv; //receives 2 column
	FLOAT* right_recv; //receives 1 column
	FLOAT* top_recv; //receives 1 row
	FLOAT* bottom_recv; //receives 2 rows

	//3d
	FLOAT* front_recv; //receives 2 rows
	FLOAT* back_recv; //receives 1 row

	int left_recv_counter;
	int right_recv_counter;
	int top_recv_counter;
	int bottom_recv_counter;

	int front_recv_counter;
	int back_recv_counter;
*/

	void markPressureAndVelo(FlowField & flowField, int i, int j);

    public:
		PressureBufferReadStencil( const Parameters& parameters );
		~PressureBufferReadStencil();

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

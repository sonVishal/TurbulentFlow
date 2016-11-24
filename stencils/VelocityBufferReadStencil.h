#ifndef _VELOCITY_BUFFER_READ_STENCIL_H_
#define _VELOCITY_BUFFER_READ_STENCIL_H_

#include "../Definitions.h"
#include "../Parameters.h"
#include "../Stencil.h"
#include "../FlowField.h"
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "CommBuffer.h"


//#define DEBUG_VELOREAD

class VelocityBufferReadStencil : public BoundaryStencil<FlowField> {
    private:
	//Buffer for 2d case also used in 3d
	ComBuf left_recv; //receives 2 column
	ComBuf right_recv; //receives 1 column
	ComBuf top_recv; //receives 1 row
	ComBuf bottom_recv; //receives 2 rows

	//3d
	ComBuf front_recv; //receives 2 rows
	ComBuf  back_recv; //receives 1 row

	void markVelocityAndVelo(FlowField & flowField, int i, int j);

    public:
		VelocityBufferReadStencil( const Parameters& parameters );
		~VelocityBufferReadStencil();

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

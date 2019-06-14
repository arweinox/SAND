#pragma once


namespace renderer {

	/*enum PrimitiveType {
		CubeType,
		CylinderType,
		SphereType,
		None
	};*/

	// Can reduce precision
	struct BoundingBox2D {
		int x_min_, x_max_,
			y_min_, y_max_;
	};

	// Can reduce precision
	struct BoundingBox2Df {
		float x_min_, x_max_,
			y_min_, y_max_;
	};


	// Can reduce precision
	struct Workspace {
		float x_min_, x_max_,
			y_min_, y_max_,
			z_min_, z_max_;
	};

	struct Pose {
		float pos_[3];
		float euler_[3];
		void setPose(float x, float y, float z) {
			pos_[0] = x;
			pos_[1] = y;
			pos_[2] = z;
		}
		void setEuler(float x, float y, float z) {
			euler_[0] = x;
			euler_[1] = y;
			euler_[2] = z;
		}
	};

	// Objects in scene
	typedef int ObjectName;

	// define all object names
	const ObjectName blue_cup = 0;
	const ObjectName clorox=  1;
	const ObjectName coke = 2;
	const ObjectName detergent= 3;
	const ObjectName downy= 4;
	const ObjectName ranch= 5;
	const ObjectName red_bowl= 6;
	const ObjectName salt= 7;
	const ObjectName scotch_brite= 8;
	const ObjectName spray_bottle= 9;
	const ObjectName sugar= 10;
	const ObjectName sunscreen= 11;
	const ObjectName tide= 12;
	const ObjectName toy= 13;
	const ObjectName waterpot= 14;
	const ObjectName none = 15;

	// 14 32-bit words
	struct Object {
		//PrimitiveType type_;
		ObjectName name_;	// String name is demarcated as number corresponding to sample
		Pose pose_;			// 6-32 bit
		float dim_[3];	// 3 32-bit
		float trans_mat_[4];	//4 - 32-bit	// Most likely not needed?

		void setDim(int x, int y, int z) {
			dim_[0] = x;
			dim_[1] = y;
			dim_[2] = z;
		}
		void setTransMat(int w, int x, int y, int z) {
			trans_mat_[0] = x;
			trans_mat_[1] = y;
			trans_mat_[2] = z;
			trans_mat_[3] = w;
		}
	};

	struct Heuristic {
		Workspace workspace_;
		BoundingBox2D bbox_;
		int total_num_points_;
	};

	// Deprecated
	struct ObjectHeuristicPair {
		ObjectName object_name;
		Heuristic heuristic;
	};
}

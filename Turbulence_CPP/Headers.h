#ifndef header_hh_inluded
#define header_hh_inluded

class FVCell // cell class
{
public:
	// Location and index
	double x_loc;
	double y_loc;
	double z_loc;
	int x_index;
	int y_index;
	int z_index;

	// u-velocity component
	double u;

	// v-velocity component
	double v;

	// w-velocity component
	double w;
	
	// Pressure
	double p;

	// Cell neighbors
	std::shared_ptr< FVCell > xm1;
	std::shared_ptr< FVCell > xp1;
	std::shared_ptr< FVCell > ym1;
	std::shared_ptr< FVCell > yp1;
	std::shared_ptr< FVCell > zm1;
	std::shared_ptr< FVCell > zp1;

	// Constructor
	FVCell() :
		x_loc( 0.0 ),
		y_loc( 0.0 ),
		z_loc( 0.0 ),
		x_index( 0 ),
		y_index( 0 ),
		z_index( 0 ),

		u( 0.0 ),
		v( 0.0 ),
		w( 0.0 ),
		
		p( 0.0 ),
	
		xm1( nullptr ),
		xp1( nullptr ),
		ym1( nullptr ),
		yp1( nullptr ),
		zm1( nullptr ),
		zp1( nullptr )
	{};

	// Destructor
	~FVCell(){};

};

class Global_domain
{
public:

	int nx = 16;
	int ny = 16;
	int nz = 16;

	class Map_3D
	{
	public:
		Map_3D(int i=0, int j=0, int k=0):
			x_index( i ),
			y_index( j ),
			z_index( k )
		{}

		//Overload operator for equality in the map
		bool operator==(const Map_3D& other) const {
    	return (x_index==other.x_index &&
            y_index==other.y_index &&
            z_index==other.z_index);
  		}

		int x_index, y_index, z_index;

		// Destructor
		~Map_3D(){};
	};

	class Map_2D
	{
	public:
		Map_2D(int i=0, int j=0):
			x_index( i ),
			y_index( j )
		{}

		//Overload operator for equality in the map
		bool operator==(const Map_2D& other) const {
    	return (x_index==other.x_index &&
            y_index==other.y_index);
  		}

		int x_index, y_index;

		// Destructor
		~Map_2D(){};
	};


	// Default constructor
	Global_domain(){};
	// Destructor
	~Global_domain(){};

	//Global domain member functions
	struct hasher_key_3D
	{
    	size_t operator()(const Map_3D& map_val) const
    	{
        	return (std::hash<int>()(map_val.x_index) << 1 ^ std::hash<int>()(map_val.y_index) << 1 ^ std::hash<int>()(map_val.z_index));
    	}
	};

	struct hasher_key_2D
	{
    	size_t operator()(const Map_2D& map_val) const
    	{
        	return (std::hash<int>()(map_val.x_index) << 1 ^ std::hash<int>()(map_val.y_index));
    	}
	};

	//Virtual function here
	virtual void init_domain()=0;

};

class TGVClass : public :: Global_domain
{
public:
	double Re;
	std::unordered_map< Map_3D, std::shared_ptr< FVCell >, hasher_key_3D > FVCell_Map;

	// Default constructor
	TGVClass() :
		Re( 0.0 )
	{}

	// Member constructor
	TGVClass(
		double const Re_
	) :
		Re( Re_ )
	{}

	// Destructor
	~TGVClass(){};

	// Member Functions
	void init_domain();//Initialize this in main.cpp
};


class TGV2DClass : public :: Global_domain
{
public:
	double Re;
	std::unordered_map< Map_2D, std::shared_ptr< FVCell >, hasher_key_2D > FVCell_Map;

	// Default constructor
	TGV2DClass() :
		Re( 0.0 )
	{}

	// Member constructor
	TGV2DClass(
		double const Re_
	) :
		Re( Re_ )
	{}

	// Destructor
	~TGV2DClass(){};

	// Member Functions
	void init_domain();//Initialize this in main.cpp
};























#endif
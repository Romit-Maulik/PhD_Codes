// C++ headers
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <memory>
#include <stdio.h>
#include <sstream>
#include <string>
#include <time.h>
#include <vector>
#include<unordered_map>

// My headers
#include "Headers.h"

int main() { // Start of program
    /**
    Docs auto generated with DOxygen
    Program starts here
    */

	// Add runs here
	std::shared_ptr< TGVClass > run1( new TGVClass( 1600.0 ) );
	run1->init_domain();

	// Check for 2D connections
	std::shared_ptr< TGV2DClass > run2( new TGV2DClass( 1600.0 ) );
	run2->init_domain();
}

//Add functions defined in TGVClass here (port to different header file later)

void TGVClass::init_domain() // Initialize domain correctly for each run.
{
	//new FVCell (s) are made and pushed to unordered_map for each i,j,k
	//New map_key, new FVCell and push both into Cellvect (unordered_map)
	//https://stackoverflow.com/questions/17016175/c-unordered-map-using-a-custom-class-type-as-the-key

	for (int k=0; k<nz; k++)
	{
		for (int j=0;j<ny;j++)
		{
			for (int i=0; i<nx; i++)
			{
				std::shared_ptr< FVCell > thisCell( new FVCell );

				thisCell->x_index = i;
				thisCell->y_index = j;
				thisCell->z_index = k;

				FVCell_Map[Map_3D(i,j,k)]=thisCell;
			}
		}
	}

	//Connect neighbors in interior
	for (int k=1; k<nz-1; k++)
	{
		for (int j=1;j<ny-1;j++)
		{
			for (int i=1; i<nx-1; i++)
			{
				FVCell_Map[Map_3D(i,j,k)]->xm1 = FVCell_Map[Map_3D(i-1,j,k)];
				FVCell_Map[Map_3D(i,j,k)]->xp1 = FVCell_Map[Map_3D(i+1,j,k)];
				FVCell_Map[Map_3D(i,j,k)]->ym1 = FVCell_Map[Map_3D(i,j-1,k)];
				FVCell_Map[Map_3D(i,j,k)]->yp1 = FVCell_Map[Map_3D(i,j+1,k)];
				FVCell_Map[Map_3D(i,j,k)]->zm1 = FVCell_Map[Map_3D(i,j,k-1)];
				FVCell_Map[Map_3D(i,j,k)]->zp1 = FVCell_Map[Map_3D(i,j,k+1)];
			}
		}
	}

	//Can add string argument based if statements here - default periodic for TGV class
	//Connect neighbors in boundary - in x
	for (int k=1; k<nz-1; k++)
	{
		for (int j=1;j<ny-1;j++)
		{
			FVCell_Map[Map_3D(0,j,k)]->xm1 = FVCell_Map[Map_3D(nx-1,j,k)];
			FVCell_Map[Map_3D(0,j,k)]->xp1 = FVCell_Map[Map_3D(1,j,k)];

			FVCell_Map[Map_3D(nx-1,j,k)]->xm1 = FVCell_Map[Map_3D(nx-2,j,k)];
			FVCell_Map[Map_3D(nx-1,j,k)]->xp1 = FVCell_Map[Map_3D(0,j,k)];

		}

	}


	//Connect neighbors in boundary - in y
	for (int i=0; i<nx; i++)
	{
		for (int k=1;k<nz-1;k++)
		{
			FVCell_Map[Map_3D(i,0,k)]->ym1 = FVCell_Map[Map_3D(i,ny-1,k)];
			FVCell_Map[Map_3D(i,0,k)]->yp1 = FVCell_Map[Map_3D(i,1,k)];

			FVCell_Map[Map_3D(i,ny-1,k)]->ym1 = FVCell_Map[Map_3D(i,ny-2,k)];
			FVCell_Map[Map_3D(i,ny-1,k)]->yp1 = FVCell_Map[Map_3D(i,0,k)];
		}

	}

	//Connect neighbors in boundary - in z
	for (int j=0; j<ny; j++)
	{
		for (int i=0;i<nx;i++)
		{
			FVCell_Map[Map_3D(i,j,0)]->zm1 = FVCell_Map[Map_3D(i,j,nz-1)];
			FVCell_Map[Map_3D(i,j,0)]->zp1 = FVCell_Map[Map_3D(i,j,1)];

			FVCell_Map[Map_3D(i,j,nz-1)]->zm1 = FVCell_Map[Map_3D(i,j,nz-2)];
			FVCell_Map[Map_3D(i,j,nz-1)]->zp1 = FVCell_Map[Map_3D(i,j,0)];
		}

	}
	// All connected now

	std::cout<<"Chosen cell: "<<FVCell_Map[Map_3D(0,1,2)]->x_index<<std::endl;	
	std::cout<<"Left neighbor x index: "<<FVCell_Map[Map_3D(0,1,2)]->xm1->x_index<<std::endl;
	std::cout<<"Right neighbor x index: "<<FVCell_Map[Map_3D(0,1,2)]->xp1->x_index<<std::endl;
};



void TGV2DClass::init_domain() // Initialize domain correctly for each run.
{
	//new FVCell (s) are made and pushed to unordered_map for each i,j - This is a 2D Domain
	//New map_key, new FVCell and push both into Cellvect (unordered_map)
	//https://stackoverflow.com/questions/17016175/c-unordered-map-using-a-custom-class-type-as-the-key

	for (int j=0;j<ny;j++)
	{
		for (int i=0; i<nx; i++)
		{
			std::shared_ptr< FVCell > thisCell( new FVCell );

			thisCell->x_index = i;
			thisCell->y_index = j;

			FVCell_Map[Map_2D(i,j)]=thisCell;
		}
	}

	//Connect neighbors in interior
	for (int j=1;j<ny-1;j++)
	{
		for (int i=1; i<nx-1; i++)
		{
			FVCell_Map[Map_2D(i,j)]->xm1 = FVCell_Map[Map_2D(i-1,j)];
			FVCell_Map[Map_2D(i,j)]->xp1 = FVCell_Map[Map_2D(i+1,j)];
			FVCell_Map[Map_2D(i,j)]->ym1 = FVCell_Map[Map_2D(i,j-1)];
			FVCell_Map[Map_2D(i,j)]->yp1 = FVCell_Map[Map_2D(i,j+1)];
			//Null pointers for z-direction
		}
	}

	//Can add string argument based if statements here - default periodic for TGV class
	//Connect neighbors in boundary - in x
	for (int j=1;j<ny-1;j++)
	{
		FVCell_Map[Map_2D(0,j)]->xm1 = FVCell_Map[Map_2D(nx-1,j)];
		FVCell_Map[Map_2D(0,j)]->xp1 = FVCell_Map[Map_2D(1,j)];

		FVCell_Map[Map_2D(nx-1,j)]->xm1 = FVCell_Map[Map_2D(nx-2,j)];
		FVCell_Map[Map_2D(nx-1,j)]->xp1 = FVCell_Map[Map_2D(0,j)];

	}


	//Connect neighbors in boundary - in y
	for (int i=0; i<nx; i++)
	{
			FVCell_Map[Map_2D(i,0)]->ym1 = FVCell_Map[Map_2D(i,ny-1)];
			FVCell_Map[Map_2D(i,0)]->yp1 = FVCell_Map[Map_2D(i,1)];

			FVCell_Map[Map_2D(i,ny-1)]->ym1 = FVCell_Map[Map_2D(i,ny-2)];
			FVCell_Map[Map_2D(i,ny-1)]->yp1 = FVCell_Map[Map_2D(i,0)];
	}
	// All connected now

	std::cout<<"Chosen cell: "<<FVCell_Map[Map_2D(0,1)]->x_index<<std::endl;	
	std::cout<<"Left neighbor x index: "<<FVCell_Map[Map_2D(0,1)]->xm1->x_index<<std::endl;
	std::cout<<"Right neighbor x index: "<<FVCell_Map[Map_2D(0,1)]->xp1->x_index<<std::endl;
};
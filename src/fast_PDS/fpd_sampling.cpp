#ifndef FPDS_HPP
#define FPDS_HPP

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <cassert>
#include <random>
#include <fstream>
#include <cstring>
#include "fpd.h"

using namespace std;

#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) < (y) ? (y) : (x))

// GLOBAL FUNCTIONS
// Use C standard rand() to get uniform random number in (0,1)
inline double URand()
{
    return (double)rand() / (double)RAND_MAX;
}

// Calculate distance between 2 points.
inline double LengthSqrt(const Point2D* p1, const Point2D* p2)
{
    double a = p2->X - p1->X;
    double b = p2->Y - p1->Y;
    return sqrtf(a*a + b*b);
}

inline void PrintBool2(bool b)
{
	printf("%s", b ? "TRUE" : "FALSE");
}

PDSampling::PDSampling(const double* inp_domain_size, const char* inp_type, const int inp_grid_size, const int* inp_k)
{
	this->domain_size = *inp_domain_size;
	this->type = new char[128];
	strcpy(type, inp_type);
	this->grid_size = inp_grid_size;
	this->k = *inp_k;
	sample_points = new PointSet;
	created_points = new PointSet;
}

PDSampling::~PDSampling() {
	delete sample_points;
	PointSet::iterator it;
	for (it = created_points->begin(); it != created_points->end(); it++) {
		delete (*it);
	}
	delete created_points;
}

// generate a random point in the region of a referenced point between radius (R,2R)
Point2D* PDSampling::GenerateRandomPointInNeighborhood(const Point2D* ref_point, const double* dist) const {
	double t = URand()*M_PI*2.0;
	double r = (URand()+1)*(ref_point->R);
	double X_increment = r*cosf(t);
	double Y_increment = r*sinf(t);
	Point2D* current_point = new Point2D;
	current_point->X = X_increment + ref_point->X;  
	current_point->Y = Y_increment + ref_point->Y;
	current_point->R = *dist;
	created_points->push_back(current_point);
	return current_point;
}

// check the points in a cell if they are too close to a newly created point
bool PDSampling::CheckPointsInCell(const Grid2D* grid, const Index2D* cell_id, const Point2D* point) const {
	bool accept = true;
	if (cell_id->first >-1 && cell_id->first < this->grid_size && cell_id->second > -1 && cell_id->second < this->grid_size) {
		PointSet* points_in_cell = (*grid)[cell_id->first][cell_id->second];
		PointSet::iterator it;
		int npts = points_in_cell->size();
		if (npts > 0) {
			for (it = points_in_cell->begin(); it != points_in_cell->end(); it++) {
				assert((*it) != point);
				double this_dist = LengthSqrt(*it, point);
				if ((this_dist - MAX((*it)->R, point->R)) < 0.0f) {
					accept = false;
					break;
				}
			}
		}
	}
	return accept;	
}

bool PDSampling::CheckNeighboringCells(const Grid2D* ref_grid, const Point2D* point) const {
	bool A, B, C, D, E, F, G, H, I;
	Index2D left, right, above, below, upper_left, upper_right, lower_left, lower_right;
	Index2D ref_cell;
	double cell_size = GetCellSize();
	
	// Which cell the test point located in
	ref_cell.first = (int)floorf(point->X / cell_size);
	ref_cell.second = (int)floorf(point->Y / cell_size);

	A = B = C = D = E = F = G = H = I = false;

	left.first = ref_cell.first - 1;
	left.second = ref_cell.second;

	right.first = ref_cell.first + 1;
	right.second = ref_cell.second;

	above.first = ref_cell.first;
	above.second = ref_cell.second + 1;
	
	below.first = ref_cell.first;
	below.second = ref_cell.second - 1;

	upper_left.first = ref_cell.first - 1;
	upper_left.second = ref_cell.second + 1;

	lower_left.first = ref_cell.first - 1;
	lower_left.second = ref_cell.second - 1;
	
	upper_right.first = ref_cell.first + 1;
	upper_right.second = ref_cell.second + 1;

	lower_right.first = ref_cell.first + 1;
	lower_right.second = ref_cell.second - 1;

/*	
	printf("current "); PrintPointsInCell(ref_grid, ref_cell);
	printf("left "); 	PrintPointsInCell(ref_grid, left);
	printf("right "); 	PrintPointsInCell(ref_grid, right);
	printf("above "); 	PrintPointsInCell(ref_grid, above);
	printf("below "); 	PrintPointsInCell(ref_grid, below);
	printf("UL "); 	PrintPointsInCell(ref_grid, upper_left);
	printf("UR "); 	PrintPointsInCell(ref_grid, upper_right);
	printf("LL "); 	PrintPointsInCell(ref_grid, lower_left);
	printf("LR "); 	PrintPointsInCell(ref_grid, lower_right);			
*/

	I = CheckPointsInCell(ref_grid, &ref_cell, point);	
	//printf("Check current cell "); PrintBool2(I); printf("\n");
	A = CheckPointsInCell(ref_grid, &left, point);
	//printf("Check left cell "); PrintBool2(A); printf("\n");
	B = CheckPointsInCell(ref_grid, &right, point);
	//printf("Check right cell "); PrintBool2(B); printf("\n");
	C = CheckPointsInCell(ref_grid, &above, point);
	//printf("Check above cell "); PrintBool2(C); printf("\n");
	D = CheckPointsInCell(ref_grid, &below, point);
	//printf("Check below cell "); PrintBool2(D); printf("\n");
	E = CheckPointsInCell(ref_grid, &upper_left, point);
	//printf("Check UL cell "); PrintBool2(E); printf("\n");
	F = CheckPointsInCell(ref_grid, &upper_right, point);
	//printf("Check UR cell "); PrintBool2(F); printf("\n");
	G = CheckPointsInCell(ref_grid, &lower_left, point);
	//printf("Check LL cell "); PrintBool2(G); printf("\n");
	H = CheckPointsInCell(ref_grid, &lower_right, point);
	//printf("Check LR cell "); PrintBool2(H); printf("\n");
	//if (I == false) printf("Generate point function does not do what it means to \n");
	
	return A&B&C&D&E&F&G&H&I;
}

void PDSampling::FastPoissonDiskSampling2D() {
	sample_points->clear();
	
	int K = 0;
	if (this->k < 30) {
		K = 30;
	}
	else {
		K = this->k;
	}

	// Allocate a grid to accelerate.
	// Each cell could have one and only one sampler in the case of uniform distribution
	Grid2D Grid (boost::extents[this->grid_size][this->grid_size]);
	for (int i = 0; i<grid_size; i++)
		for (int j = 0; j<grid_size; j++) {
			Grid[i][j] = new PointSet;
		}

	std::vector<int> active_list;
	Point2D* P0 = new Point2D;
	Index2D I0;
	double dist = GetDist();
	double cell_size = GetCellSize();
	
	P0->X = URand()*this->domain_size;
	P0->Y = URand()*this->domain_size;
	P0->R = dist;
	I0.first = (int)floorf( P0->X / cell_size );
	I0.second = (int)floorf( P0->Y / cell_size );
	//printf("Idx I0: %d %d \n", I0.first, I0.second);
    
	active_list.push_back(0);
	created_points->push_back(P0);
	this->sample_points->push_back(P0);
	Grid[I0.first][I0.second]->push_back(P0);

	int done = 1;
	
	FILE * r_distr;
	r_distr = fopen("r_distr.txt","w");
	
	while (active_list.size() != 0) {
		// Initialize a sampler.
		size_t magic_idx = (size_t)floorf(URand() * active_list.size());
		size_t start_idx = active_list[magic_idx];
		//printf("Start index: %d \n", start_idx);
		Point2D* start_point = *(sample_points->begin() + start_idx);

		// whether a new point is found (i.e. accepted to the output point list)
		bool found = false;

		// generate upto K points in the neighborhood of a referenced point
		for( size_t i=0; i<K; ++i ) {
			dist = GetDist();
			Point2D* current_point = GenerateRandomPointInNeighborhood(start_point, &dist);
			//printf("Address of current point: "); cout << current_point << endl;
			//Point2D* current_point = &tmp_point;
			Index2D* current_cell = new Index2D;
					
			//cout << "Current Point X Y: " << current_point->X << " " << current_point->Y << endl;
			//cout << "R associated with this point: " << current_point->R << endl;
			
			current_cell->first = (int)floorf(current_point->X / cell_size);
			current_cell->second = (int)floorf(current_point->Y / cell_size);

			// Discard if out of domain
			if (current_point->X < 0.0f || current_point->X >= this->domain_size)
				continue;
			if (current_point->Y < 0.0f || current_point->Y >= this->domain_size)
				continue;

			if (current_cell->first < 0 || current_cell->first >= this->grid_size)
				continue;

			if (current_cell->second < 0 || current_cell->second >= this->grid_size)
				continue;

			bool accepted = CheckNeighboringCells(&Grid, current_point);
			//printf("Accept this point? "); PrintBool2(found); printf("\n");

			if (!accepted) {
				continue;
			}
			else {
				fprintf(r_distr, "%f \n", dist);
				//printf("Accept this current point \n");
				//printf("Current cell IDX: %d %d \n", current_cell->first, current_cell->second);
				Grid[current_cell->first][current_cell->second]->push_back(current_point);
				sample_points->push_back(current_point);
				active_list.push_back(done);
				found = true;
				++done;
				break;
			}
		}
	
	//if (!found) cout << "Remove point from the active list" << endl;
        // We have to remove this test sampler from active list.
		if (found == false && active_list.size() > 0) {
			//printf("Erase ID: %d \n", *(active_list.begin() + MagicIdx));
			active_list.erase( active_list.begin() + magic_idx );
		}
	//if (active_list.size() > MaxPTS) break;
	}

	//for (int i = 0; i < grid_size; i++)
	//	for (int j = 0; j < grid_size; j++) {
	//		printf("%d ", Grid[i][j]->size());
	//		if (j == grid_size - 1) printf("\n");
	//	}
		
	for (int i = 0; i<grid_size; i++)
		for (int j = 0; j<grid_size; j++)
			delete Grid[i][j];
	//delete current_point;
	//delete current_cell;
	
	//delete ref_cell;
	fflush(r_distr);
	fclose(r_distr);
	//return sample_points;
}

void PDSampling::PrintOutput(const char* output_name) const {
	FILE* output;
	PointSet::iterator vit;
	output = fopen(output_name, "w");
	if (!output) {
		printf("Unable to open output file. \n");
		exit(1);
	}
	
	printf("%d points created \n", this->sample_points->size());
	fprintf(output, "%d %d %d %d %d\n", this->sample_points->size(), 0, 0, 0, 0);
	int i = 0;
	for (vit = this->sample_points->begin(); vit != this->sample_points->end(); vit++)
	{
		++i;
		fprintf(output, "%d %f %f %f \n", i, (*vit)->X, (*vit)->Y, 0.0);
	}
	fflush(output);
	fclose(output);
}

const PointSet* PDSampling::GetPoints() const {
	return sample_points;
}

/*
void PDSUniform::PrintPointsInCell(const Grid2D grid, const Index2D idx) const {
	PointSet::iterator it;
	PointSet ps_ptr = grid[idx.first][idx.second];
	for (it = ps_ptr.begin(); it != ps_ptr.end(); it++) {
		printf("Number of points: %d \n", ps_ptr.size());
		printf("Point coordinates: %f %f \n", (*it).X, (*it).Y);
	}
}
*/

PDSUniform::PDSUniform(const double* inp_domain_size, const char* inp_type, const double* inp_min_dist, const int* inp_k):
PDSampling(inp_domain_size, inp_type, (int)ceilf((*inp_domain_size)/(*inp_min_dist)), inp_k) {
	this->min_dist = *inp_min_dist;
	this->cell_size = this->min_dist;
}

PDSUniform::~PDSUniform() {
	// Do nothing
}

// print some info
void PDSUniform::PrintInfo() const {
	printf("Geometric dimension: %d \n", geom_dim);
	printf("Domain size: %f \n", domain_size);
	printf("Type of PDS: %s \n", type);
	printf("Minimum distance: %f \n", min_dist);
	printf("Cell size: %f \n", cell_size);
	printf("Grid size: %d \n", grid_size);
	printf("Density parameter K: %d \n", k);
}

double PDSUniform::GetDist() {
	return this->min_dist;
}

double PDSUniform::GetCellSize() const {
	return this->cell_size;
}

PDSLNorm::PDSLNorm(const double* inp_domain_size, const char* inp_type, const unsigned* inp_seed, const double* inp_max_dist,
						const double* inp_min_dist, const double* inp_mean, const double* inp_std_dev, const int* inp_k):
						PDSampling(inp_domain_size, inp_type, (int)ceilf((*inp_domain_size)/(*inp_max_dist)), inp_k),
						generator(*inp_seed), distribution(*inp_mean, *inp_std_dev) {
	this->max_dist = *inp_max_dist;
	this->min_dist = *inp_min_dist;
	this->mean = *inp_mean;
	this->std_dev = *inp_std_dev;
	this->cell_size	 = this->max_dist;
}

PDSLNorm::~PDSLNorm() {
	// Do nothing
}

void PDSLNorm::PrintInfo() const {
	printf("Geometric dimension: %d \n", geom_dim);
	printf("Domain size: %f \n", domain_size);
	printf("Type of PDS: %s \n", type);
	printf("Maximum distance: %f \n", max_dist);
	printf("Minimum distance: %f \n", min_dist);
	printf("Mean: %f \n", mean);
	printf("Standard Deviation: %f \n", std_dev);
	printf("Cell size: %f \n", cell_size);
	printf("Grid size: %d \n", grid_size);
	printf("Density parameter K: %d \n", k);
}

double PDSLNorm::GetCellSize() const {
	return this->cell_size;
}

double PDSLNorm::GetDist() {
	double rand_lnorm = this->distribution(this->generator);
	while ((rand_lnorm - this->min_dist < 0.0f) || (rand_lnorm - this->max_dist > 0.0f)) {
		rand_lnorm = distribution(generator);
	}
	return rand_lnorm;
}
#endif


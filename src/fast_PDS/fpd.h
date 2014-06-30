#include <vector>
#include <map>
#include <boost/multi_array.hpp>

#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) < (y) ? (y) : (x))

// Define our basis data structure.
typedef unsigned long long Long64;
typedef std::pair<int,int> Index2D;

struct Point2D {
	double X, Y; //Coordinates X, Y
	double R; //Minimum radial distance from other points
};

typedef std::vector<Point2D*> PointSet;
typedef boost::multi_array<PointSet*, 2> Grid2D;
//typedef Grid2D::extent_gen GridDim;

// GLOBAL FUNCTIONS
// Use C standard rand() to get uniform random number in (0,1)
inline double URand();

// Calculate distance between 2 points.
inline double LengthSqrt(const Point2D* p1, const Point2D* p2);
inline void PrintBool2(bool b);

class PDSampling {
	public:
		PDSampling();
		PDSampling(const double* inp_domain_size, const char* inp_type, const int inp_grid_size, const int* inp_k);
		~PDSampling();
		
		virtual void PrintInfo() const = 0;
		const PointSet* GetPoints() const;
		virtual double GetDist() = 0;
		virtual double GetCellSize() const = 0;
		
		Point2D* GenerateRandomPointInNeighborhood(const Point2D* ref_point, const double* dist) const;
		bool CheckPointsInCell(const Grid2D* grid, const Index2D* cell_id, const Point2D* point) const;
		bool CheckNeighboringCells(const Grid2D* ref_grid, const Point2D* point) const;
		void FastPoissonDiskSampling2D();
		void PrintOutput(const char* output_name) const;
		void PrintPointsInCell(const Grid2D* grid, const Index2D* idx) const;
	
	protected:
		static const int geom_dim = 2;
		double domain_size;
		int k;
		char* type;
		int grid_size;
		PointSet* sample_points;
		PointSet* created_points;
};

class PDSUniform: public PDSampling {
	public:
		PDSUniform(const double* inp_domain_size, const char* inp_type, const double* inp_min_dist, const int* inp_k);
		~PDSUniform();
		
		double GetDist();
		double GetCellSize() const;
		void PrintInfo() const;
	
	protected:
		double min_dist;
		double cell_size;
};

class PDSLNorm: public PDSampling {
	public:
		PDSLNorm(const double* inp_domain_size, const char* inp_type, const unsigned* inp_seed, const double* inp_max_dist,
					const double* inp_min_dist, const double* inp_mean, const double* inp_std_dev, const int* inp_k);
		~PDSLNorm();
		
		double GetDist();
		double GetCellSize() const;
		void PrintInfo() const;
			
	protected:
		double max_dist;
		double min_dist;
		double mean;
		double std_dev;
		double cell_size;
		std::default_random_engine generator;
		std::lognormal_distribution<double> distribution;
};


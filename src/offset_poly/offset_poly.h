#include "clipper.hpp"
#include<vector>

using namespace ClipperLib;

#define FACTOR 1.0e+12

typedef std::pair<double, double> DoublePoint;
typedef std::vector<DoublePoint> DPolygon;

IntPoint double2int(const DoublePoint dbl_p, const double factor);
DoublePoint int2double(const IntPoint int_p, const double factor);

DoublePoint normalize(const DoublePoint v);
DoublePoint calc_normal(const DoublePoint p1, const DoublePoint p2, const double factor);
DoublePoint ComputeMidPoint(DoublePoint* p1, DoublePoint* p2);

double calc_dot_product(const DoublePoint v1, const DoublePoint v2);
double ComputeLength(DoublePoint* p1, DoublePoint* p2);
int reorder_pos(DPolygon poly, DPolygon ref_poly, const double factor);

bool PointInPolygon(DPolygon poly, DoublePoint *test_point);
bool CheckOutputOffset(DPolygon input_poly, DPolygon output_poly, const double offset_value);
bool HandleSmallEdges(DPolygon& polygon, double min_edge_length);

std::vector<int> CheckSmallEdges(DPolygon polygon, double min_edge_length);
DPolygon OffsetPolygonsWrapper(DPolygon input_poly, const double offset_value, const double factor);

std::string removeExtension(const std::string filename);

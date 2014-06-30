#include<cmath>
#include"offset_poly.h"
using namespace ClipperLib;

//The clipper library uses all large integer instead of double
//Convert double to integer
IntPoint double2int(const DoublePoint dbl_p, const double factor){
	IntPoint int_p(floor(dbl_p.first*factor),floor(dbl_p.second*factor));
	return int_p;
};

//convert integer to double
DoublePoint int2double(const IntPoint int_p, const double factor){
	DoublePoint dbl_p;
	dbl_p.first = int_p.X/factor;
	dbl_p.second = int_p.Y/factor;
	return dbl_p;
};

//Normalize a vector
DoublePoint normalize(const DoublePoint v){
	DoublePoint norm_v;
	double norm_factor = sqrt(v.first*v.first + v.second*v.second);
	norm_v.first = v.first/norm_factor;
	norm_v.second = v.second/norm_factor;
	return norm_v;
}

//calculate normal of a vector
DoublePoint calc_normal(const DoublePoint p1, const DoublePoint p2, const double factor){
	DoublePoint norm_vector;
	norm_vector.first = p1.second - p2.second;
	norm_vector.second = p2.first - p1.first;
	return normalize(norm_vector);
}

//calculate dot product
double calc_dot_product(const DoublePoint v1, const DoublePoint v2){
	return v1.first*v2.first + v1.second*v2.second;
}

//calculate distance between 2 points
double ComputeLength(DoublePoint* p1, DoublePoint* p2) {
	double dx = p2->first - p1->first;
	double dy = p2->second - p1->second;
	return(sqrtf(dx*dx + dy*dy));
}

//Return the mid point of a line segment
DoublePoint ComputeMidPoint(DoublePoint* p1, DoublePoint* p2) {
	DoublePoint mid_point;
	mid_point.first = (p1->first + p2->first)/2;
	mid_point.second = (p1->second + p2->second)/2;
	return mid_point;
}

// This function returns a value which can be used to rotate the offset
// polygon to match the original polygons.
int reorder_pos(DPolygon poly, DPolygon ref_poly, const double factor){
	std::vector<DoublePoint> ref_norm_list, norm_list;
	std::vector<DoublePoint>::iterator it, jt;
	DoublePoint p1, p2;
	DoublePoint this_norm;
	
	for (it = poly.begin(); it != poly.end(); it++){
		p1 = *it;
		if (it == (poly.end()-1)){
			p2 = *(poly.begin());
		}
		else{
			p2 = *(it+1);
		}
		this_norm = calc_normal(p1, p2, factor);
		norm_list.push_back(this_norm);
	}

	for (jt = ref_poly.begin(); jt < ref_poly.end(); jt++){
		p1 = *jt;
		if (jt == (ref_poly.end()-1)){
			p2 = *(poly.begin());
		}
		else{
			p2 = *(jt+1);
		}
		this_norm = calc_normal(p1, p2, factor);
		ref_norm_list.push_back(this_norm);
	}

	int pos = 0;
	unsigned int k = 0;
	DoublePoint target_norm = norm_list[0];
	double max_dot_product = 0;
	while (k < ref_norm_list.size()){
		++k;
		DoublePoint tmp_norm = ref_norm_list[k-1];
		double dot_product = calc_dot_product(target_norm, tmp_norm);
		if (dot_product < 0) {
			continue;
		}
		else{
			if ((max_dot_product - dot_product) < 0){
				max_dot_product = dot_product;
				pos = k;				
			}
			else continue;
		}
	}
	return pos;
}

//Check whether all the vertices of the offset polygon lie inside
//the original polygon
bool PointInPolygon(DPolygon poly, DoublePoint *test_point) {
  int   i, j=poly.size()-1 ;
  double x, y;
  bool  oddNodes=false;
  
  x = test_point->first;
  y = test_point->second;

  for (i=0; i<poly.size(); i++) {
    if (poly[i].second<y && poly[j].second>=y
    ||  poly[j].second<y && poly[i].second>=y) {
      if (poly[i].first+(y-poly[i].second)/(poly[j].second-poly[i].second)*(poly[j].first-poly[i].first)<x) {
        oddNodes=!oddNodes;
      }
    }
    j=i;
  }
  return oddNodes;
}

// Check output of offset polygons
bool CheckOutputOffset(DPolygon input_poly, DPolygon output_poly, const double offset_value){
	bool is_valid = true;
	bool point_in_poly = false;
	DPolygon::iterator vit;
	if (offset_value < 0 ){
		for (vit = output_poly.begin(); vit != output_poly.end(); vit++){
			point_in_poly = PointInPolygon(input_poly, &(*vit));
			is_valid = point_in_poly && is_valid;
		}
	}
	else{
		// Do nothing
	}
	return is_valid;
}

//Check where small edges occur
std::vector<int> CheckSmallEdges(DPolygon polygon, double min_edge_length) {
	std::vector<int> positions;
	double edge_length;
	DoublePoint p1, p2;
	for (unsigned int i = 0; i < polygon.size() - 1; i++) {
		DoublePoint p1 = polygon[i];
		DoublePoint p2 = polygon[i+1];
		edge_length = ComputeLength(&p1, &p2);
		if ((edge_length - min_edge_length) < 0) {
			positions.push_back(i);
		}
	}

	//Wrap around
	p1 = *(polygon.end()-1);
	p2 = *(polygon.begin());
	edge_length = ComputeLength(&p1, &p2);
	if ((edge_length - min_edge_length) < 0) {
		positions.push_back(polygon.size()-1);
	}
	return positions;
}

//Recursively elimiate small edges until there are no small edges 
bool HandleSmallEdges(DPolygon& polygon, double min_edge_length) {
	DPolygon poly_save = polygon;
	std::vector<int> positions;
	positions = CheckSmallEdges(polygon, min_edge_length);
	bool no_short_edges = (positions.size() == 0);
	if (no_short_edges) {
		printf("No short edges \n");
		return true;
	}
	else {
		printf("Number of short edges: %d \n", positions.size());
		printf("Eliminating short edges \n");
		DoublePoint new_point;
		for (unsigned int i = 0; i < positions.size(); i++) {
			std::vector<int> rm_pos;
			unsigned int this_pos = positions[i];
			if (this_pos < poly_save.size()) {
				new_point = ComputeMidPoint(&poly_save[this_pos], &poly_save[this_pos + 1]);
			}
			else {
				new_point = ComputeMidPoint(&poly_save[0], &poly_save[poly_save.size()-1]);
			}
			polygon.erase(polygon.begin() + this_pos);
			polygon.insert(polygon.begin() + this_pos, new_point);
		}

		unsigned int erase_cnt = 0;
		for (unsigned int i = 0; i < positions.size(); i++) {
			int this_pos = positions[i];
			int next_pos;
			if (i < positions.size()) {
				next_pos = positions[i+1];
			}
			else if (i == positions.size()) {
				continue;
			}
			if (next_pos == this_pos + 1) {
				continue;	
			}
			else {
				int rm_pos = positions[i] + 1 - erase_cnt;
				polygon.erase(polygon.begin() + rm_pos);
				++erase_cnt;
			}
		}
		for (unsigned int vertex = 0; vertex < polygon.size(); vertex++) {
			printf("Cooridnates of vertex no.%d of new polygon: %f %f \n", vertex + 1, polygon[vertex].first, polygon[vertex].second);
		}

		no_short_edges = HandleSmallEdges(polygon, min_edge_length);
	}
	return no_short_edges;
}

//Wrapper for polygon offset function in clipper library
DPolygon OffsetPolygonsWrapper(DPolygon input_poly, const double offset_value, const double factor){
	Polygons in_polys, out_polys(1);
	Polygon int_input_poly;
	DPolygon dbl_output_poly;
	DPolygon::iterator iit;
	std::vector<IntPoint>::iterator oit;

	for (iit = input_poly.begin(); iit != input_poly.end(); iit++){
		int_input_poly.push_back(double2int((*iit), factor));
	}

	in_polys.push_back(int_input_poly);
	OffsetPolygons(in_polys, out_polys, offset_value*factor, jtSquare, 2, false);

	for (oit = out_polys[0].begin(); oit != out_polys[0].end(); oit++){
		dbl_output_poly.push_back(int2double((*oit), factor));
	}

	return dbl_output_poly;
}


std::string removeExtension(const std::string filename) {
    size_t lastdot = filename.find_last_of(".");
    if (lastdot == std::string::npos) return filename;
    return filename.substr(0, lastdot); 
}

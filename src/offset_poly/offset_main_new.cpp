#include<cstdio>
#include<fstream>
#include<sstream>
#include<iostream>
#include<algorithm>
#include<string.h>
#include "offset_poly.h"
#include <iomanip>
using namespace std;

int main( int argc, char** argv)
{
	if (argc < 3) {
		cout << "Usage: PolyOffset input_file params_file" << endl;
		exit(1);
	}

	FILE *output;
	DPolygon in_poly, out_poly;
	unsigned int num_layers;
	static double min_edge_length;
	std::vector<int> small_edge_pos;
	std::vector<double> offset_values;
	std::ifstream fin;

	//Read input file (AVS format) that contains polygon information
	fin.open(argv[1]);
	if (!fin){
		printf("Unable to read input file. \n");
		exit(1);
	}
	printf("READING INPUT FILE\n");
	printf("INPUT VERTEX COORDINATES:\n");
	std::string line;
	long64 line_num = 0;
	while (std::getline(fin,line)) {
		++line_num;
		if (line_num == 1) {
			//skip first line
		}
		else {
			std::istringstream this_line;
			this_line.str(line);
			double x,y;
			int id;
			this_line >> id >> x >> y;
			printf("%.12f %.12f \n", x, y);
			DoublePoint this_point(x,y);
			in_poly.push_back(this_point);
		}
	}
	fin.close();	

	//Read parameter file
	fin.open(argv[2]);
	if (!fin){
		printf("ERROR!!! Unable to read params file. \n");
		exit(1);
	}
	line_num = 0;
	while (std::getline(fin,line)) {
		++line_num;
		std:istringstream this_line;
		this_line.str(line);
		if (line_num == 1) {
			this_line >> min_edge_length;
		}
		else if (line_num == 2) {
			this_line >> num_layers;
		}
		else {
			double h;
			this_line >> h;
			offset_values.push_back(h);
		}
	}
	fin.close();
	if (offset_values.size() != num_layers) {
		printf("ERROR!!! In params_file, the number of layers is not equal to the number of offset values \n");
		exit(1);
	}
	
	//Print out some info after reading input and param files
	std::vector<double>::iterator it;
	printf("MINIMUM EDGE LENGTH:\n %f\n", min_edge_length);
	printf("OFFSET VALUES: \n");
	for(it = offset_values.begin(); it != offset_values.end(); it++)
		printf("%.12f \n", *(it));

	std::string output_name[num_layers+1];
	for (unsigned int i = 0; i < num_layers+1; ++i){
		std::string fname_base = removeExtension(argv[1]);
		std::string fname_ext = ".inp";
		ostringstream convert;
		convert << (i+1);
		fname_base.append("_sim");
		fname_base.append(convert.str());
		fname_base.append(fname_ext);
		output_name[i] = fname_base;
	}

	//Do polygon offsetting
	try 
	{
		int i = num_layers;
		// Start processing each layer until there is none left
		while (i > 0) {
			--i;
			DPolygon::iterator vit;
			int write_flag = 0;
			int num_points_in_layer;

			// Get the offset value for this layer
			double current_offset_value = offset_values[num_layers-i-1];

			// If the offset value is zero, return the original polygon immediately and start next iterator
			if (current_offset_value == 0){
				out_poly = in_poly;
				continue;
			}
						
			// Offset the polygon, output is stored in out_poly
			out_poly = OffsetPolygonsWrapper(in_poly, current_offset_value, FACTOR);

			// Now check if the output is valid
			// If the size of the output polygon is zero, there is an external error with Clipper Library
			// Otherwise the output is not valid if the offset points are not in correct positions
			// In this case, the original polygon is returned
			if ((out_poly.size() == 0) && abs(current_offset_value)>1.0e-10) {
				printf("ERROR!!! Fail to offset polygons \n");
				exit(-10);
			}
			
			DPolygon save_out;
			
			bool is_valid = CheckOutputOffset(in_poly, out_poly, current_offset_value);
			if (!is_valid) {
				printf("The output polygon is not valid \n");
				printf("Return the original polygon, write output to file and exit \n");
				out_poly = in_poly;
				save_out = out_poly;
				write_flag = 1;
			}
			else {
				printf("============================= \n");
				for (unsigned int vertex = 0; vertex < out_poly.size(); vertex ++) {
					printf("Vertex coordinates: %f %f \n", out_poly[vertex].first, out_poly[vertex].second);
				}

				// Save the output as the input for next iteration			
				save_out = out_poly;

				// Now reorder the nodes of the output polygon so that it can be triangulated using LaGriT
				int pos = reorder_pos (out_poly, in_poly, FACTOR);
				int rot_pos;

				if (in_poly.size() == out_poly.size()){
					rot_pos = out_poly.size() - pos;
				}
				else{
					rot_pos = out_poly.size() - pos - (in_poly.size()-out_poly.size()) + 1;
					if (rot_pos < 0) rot_pos = 0;
				}
				
				//printf("Rotate positions: %d %d \n", pos, rot_pos);
				// reverse the orientation of the offset polygon
				reverse (out_poly.begin()+1, out_poly.end());

				// Collapse small edges in output polygon
				bool test = HandleSmallEdges(out_poly, min_edge_length);
				
				// Reorder offset polygon so that the vertices of the offset polygon
				// correspond to the original vertices. This is also necessary to output
				// the polygon in a way that can be triangulated by LaGriT.
				rotate (in_poly.begin(), in_poly.begin() + pos, in_poly.end());
			}
			
			// Open output file
			output = fopen(output_name[num_layers-i-1].c_str(),"w");
			if (!output) {
				printf("ERROR!!! Unable to open output file. \n");
				exit(-1);
			}
			printf("Output file: %s \n", output_name[num_layers-i-1].c_str());

			// Start writing output files in avs format
			if (write_flag == 0) {
				num_points_in_layer = in_poly.size() + out_poly.size() + 2;
				fprintf(output, "%d %d %d %d %d \n", num_points_in_layer, 0, 0, 0, 0);
				int id = 1;
				for (vit = in_poly.begin(); vit != in_poly.end(); vit ++) {
					fprintf (output, "%d %.12f %.12f %f \n", id, vit->first, vit->second, 0.0);
					++id;
				}

				for (vit = out_poly.begin(); vit != out_poly.end(); vit ++) {
					fprintf (output, "%d %.12f %.12f %f \n", id, vit->first, vit->second, 0.0);
					++id;	
				}

				//Add 2 duplicate points to connect the polygons.
				DoublePoint second_last_point = *(out_poly.begin());
				DoublePoint last_point = *(in_poly.end() - 1);
				fprintf(output, "%d %.12f %.12f %f \n", id, second_last_point.first, second_last_point.second, 0.0);
				++id;
				fprintf(output, "%d %.12f %.12f %f \n", id, last_point.first, last_point.second, 0.0);
				if (i == 0){
					fflush(output);
					fclose(output);
					output = fopen(output_name[num_layers].c_str(),"w");
					if (!output) {
						printf("ERROR!!! Unable to open last output file. \n");
						exit(-1);
					}
					printf("Output file: %s \n", output_name[num_layers].c_str());
					int sub_id = 1;
					reverse(out_poly.begin(), out_poly.end());
					fprintf(output, "%d %d %d %d %d \n", out_poly.size() , 0, 0, 0, 0);	
					for (vit = out_poly.begin(); vit != out_poly.end(); vit ++) {
						fprintf (output, "%d %.12f %.12f %f \n", sub_id, vit->first, vit->second, 0.0);
						++sub_id;
					}
					fclose(output);
				}
				out_poly.clear();
				in_poly = save_out;
			}
			else if (write_flag == 1){
				num_points_in_layer = in_poly.size();
				fprintf(output, "%d %d %d %d %d \n", num_points_in_layer, 0, 0, 0, 0);
				int id = 1;
				for (vit = in_poly.begin(); vit != in_poly.end(); vit ++) {
					fprintf (output, "%d %.12f %.12f %f \n", id, vit->first, vit->second, 0.0);
					++id;
				}
				break;
			}
			else{
				printf("ERROR!!! Unknown write flag \n");
				exit(-1);
			}			
		}
	}
	catch (...)
	{
		printf("Check exception error \n");
	}
	printf("============================= \n");
	printf("Done!!! \n");
	return 0;
}


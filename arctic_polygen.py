#!/usr/lanl/bin/python
import sys
from string import *
import os
import time
import glob
import random
from math import *

def FindMaxDist(centroid, polygon):
	dist_list = []
	#print "Centroid"
	#print centroid
	#print "Polygon"
	for i in polygon:
		#print i
		a = i[0] - centroid[0]
		b = i[1] - centroid[1]
		dist = sqrt(a*a + b*b)
		dist_list.append(dist)
	return max(dist_list)
	

def FindMaxDistList(centroid_list, polygon_list):
	if (len(centroid_list) != len(polygon_list)):
		print "ERROR!!! Number of centroids and number of polygons do not match"
		exit()
	max_dist_list = []
	for i in range(0,len(centroid_list)):
		max_dist_list.append(FindMaxDist(centroid_list[i], polygon_list[i]))
	return max_dist_list

# Function to calculate actual elevation from a polynomial. This polynomial (could be other analytic function)
# is obtained from the statistics for elevation profile (contact Cathy/Chandana/Alexei)
def MakeElevation(max_dist):
	elevation = 0.00059*(max_dist**6) - 0.0086*(max_dist**5) + 0.046*(max_dist**4) - 0.11*(max_dist**3) + 0.16*(max_dist**2) - 0.056*(max_dist) + 4.2
	return elevation

# Get the maximum distance field value for each polygon from LaGriT output	
def GetDfieldMax(file_name):
	dfield = []
	f = open(file_name,"r")
	f.readline()
	f.readline()
	f.readline()
	current_line = f.readline()
	while current_line:
		dfield.append(float(current_line.split()[0]))
		current_line = f.readline()
	return max(dfield)

# Get the ouput directory from command line argument, otherwise use default "output"
if (len(sys.argv) < 3):
	print("ERROR!!! Please give input file and output directory.")
	print("Usage: arctic_polygen.py input_file output_dir")
	exit()
else:
	input_fname = str(sys.argv[1])
	output_folder = str(sys.argv[2])

os.system("date")

current_dir = os.getcwd()
fast_pdsampling = current_dir + "/bin/fast_pdsampling "
mkdual = current_dir + "/bin/mkdual "
lagrit = current_dir + "/bin/lagrit "
offset_polygon = current_dir + "/bin/offset_polygon "
#input_fname = current_dir + "/input/input.txt"
#params_file = current_dir + "input/params.txt"

infile = open(input_fname,"r")
first_line = infile.readline().rstrip()
if (first_line != '**PolygonGenerationInput**'):
	print('Error!!! Please check input file format.')
	exit()	

current_line = infile.readline().rstrip()
if (current_line != '**SeedValue**'):
	print('Error!!! Please check input file format.')
	exit()	

seed_value = infile.readline().rstrip()
if (seed_value == 'NULL'):
	random.seed()
	seed_value = random.randint(0,4294967295)
else:
	seed_value = abs(int(seed_value))

current_line = infile.readline().rstrip()
if (current_line != '**PolygonGeneration**'):
	print('Error!!! Please check input file format.')
	exit()	

domain_size = float(infile.readline())
pdist_type = infile.readline().rstrip()
if (pdist_type == 'LogNormal'):
	in_max_dist = float(infile.readline())
	in_min_dist = float(infile.readline())
	in_mean = float(infile.readline())
	in_stdev = float(infile.readline())
	in_density_param = int(infile.readline())
	#fpd_output_name = infile.readline().rstrip()
	command = fast_pdsampling + str(seed_value) + " " + str(domain_size) + " " + pdist_type + " " + str(in_max_dist) + " " + str(in_min_dist) + " " + str(in_mean) + " " + str(in_stdev) + " " + str(in_density_param)
elif (pdist_type == 'Uniform'):
	in_min_dist = float(infile.readline())
	in_density_param = int(infile.readline())
	#fpd_output_name = infile.readline().rstrip()
	command = fast_pdsampling + str(seed_value) + " " + str(domain_size) + " " + pdist_type + " " + str(in_min_dist) + " " + str(in_density_param)
else:
	print("Error reading type of point distribution. Options are 'Uniform' and 'LogNormal'")

current_line = infile.readline().rstrip()
while (current_line != '**Model**'):
	current_line = infile.readline().rstrip()

model_type = infile.readline().rstrip()
if (model_type == 'IntermediateScale'):
	f_params = open('input/params.txt','w')
	min_edge_length = float(infile.readline())
	f_params.write(str(min_edge_length) + '\n')
	num_layers = int(infile.readline())
	f_params.write(str(num_layers) + '\n')
	i = 0
	current_line = infile.readline()
	while (i < num_layers and current_line):
		f_params.write(infile.readline())
		i = i + 1
		current_line = infile.readline()
	if (i != num_layers):
		print("ERROR!!! The number of layers does not match with the number of values you provide.")
		print("Please check input file.")
		exit()
	f_params.flush()
	f_params.close()
	refine_scale = domain_size / 100
else:
	refine_scale = domain_size / 100
infile.close()

# produce poisson disk sampling point distribution
print('****GeneratingPoints****')
print('Type of point distribution: ' + pdist_type)
print('Calling command: ')
print(command)
os.system(command)

# Writing lagrit control file to make delauney triangulation from the PD Sampling point distribution
f_la = open("tmp_process_pdpts.lgi","w")
f_la.write("cmo / create / mo / / / triplane \n")
f_la.write("read / avs / tmp_pd_pts.inp / mo \n")
f_la.write("dump / gmv / tmp_pd_pts.gmv / mo \n")
f_la.write("\n")
f_la.write("cmo / setatt / mo / imt / 1,0,0 / 1 \n")
f_la.write("cmo / setatt / mo / itp / 1,0,0 / 0 \n")
f_la.write("connect \n")
f_la.write("resetpts / itp \n")
f_la.write("\n")
#f_la.write("recon 1 \n")
f_la.write("recon 0 \n")
f_la.write("\n")
f_la.write("cmo / status / brief \n")
f_la.write("\n")
f_la.write("resetpts / itp \n")
f_la.write("dump / gmv / tmp_pts_tri.gmv / mo / ascii \n")
f_la.write("\n")
#f_la.write("cmo / create / mo1 \n")
#f_la.write("copypts / mo1 / mo \n")
f_la.write("pset / p_boundary / attribute / itp / 1,0,0 / eq / 10 \n")
f_la.write("pset / p_boundary / write / boundary_nodes / ascii \n")
#f_la.write("rmpoint / pset get p_boundary \n")
#f_la.write("rmpoint / compress \n")
#f_la.write("cmo / status / brief \n")
#f_la.write("\n")
#f_la.write("cmo / modatt / mo1 / imt / ioflag / l \n")
#f_la.write("cmo / modatt / mo1 / itp / ioflag / l \n")
#f_la.write("cmo / modatt / mo1 / icr / ioflag / l \n")
#f_la.write("cmo / modatt / mo1 / isn / ioflag / l \n")
#f_la.write("\n")
#f_la.write("dump / avs / pset.inp / mo1 \n")
#f_la.write("\n")
f_la.write("finish \n")
f_la.write("\n")
f_la.flush()
f_la.close()

#Now run the LaGriT script to get the output
command = lagrit + "<tmp_process_pdpts.lgi >tmp_process_pdpts.out"
os.system(command)

#Clean up a little bit
command = "mv tmp_process_pdpts.lgi tmp_process_pdpts.out temporary_files"
os.system(command)
os.system("rm *x3dgen")

#run mkdual to get the dual polygonal mesh
print("=================================")
print('****GeneratingPolygons****')
command = mkdual + "tmp_pts_tri.gmv"
os.system(command)
print("=================================")

#Remove outer polygons
outer_polys = set()
f = open("boundary_nodes.vertexset", "r")
f.readline()
f.readline()
this_line = f.readline()
while this_line:
	this_line = this_line.split()
	for i in this_line:
		outer_polys.add(int(i))
	this_line = f.readline()

print("Processing polygons")
print("There are " + str(len(outer_polys)) + " polygons on the boundary")

#for i in outer_polys:
#	fourdigit = str(i)
#	while len(fourdigit) < len(str(num_polys)):
#		fourdigit = '0' + fourdigit
#	os.system("rm polys/poly" + fourdigit + ".avs")
#os.system("date")

# Extract the polygons and write each polygon to a file
centroid_list = []
f_centroids = open("tmp_pd_pts.inp","r")
current_line = f_centroids.readline()
current_line = current_line.split()
for i in range(1,int(current_line[0])+1):
	if (i in outer_polys):
		f_centroids.readline()
		continue
	current_line = f_centroids.readline()
	current_line = current_line.split()
	current_line = [float(x) for x in current_line[1:]]
	centroid_list.append(current_line)

print "Read " + str(len(centroid_list)) + " centroids"

os.system("mkdir polys")

f_polys = open("tmp_pts_tri-dual.gmv", "r")
node_list = []

current_line = f_polys.readline()
while (lower(current_line[0:4]) != 'node'):
	current_line = f_polys.readline()

current_line = current_line.split()
num_nodes = int(current_line[1])
node_cnt = 0
current_line = f_polys.readline()

while (lower(current_line[0:4]) != 'cell'):
	node_cnt = node_cnt + 1
	node_list.append(current_line)
	current_line = f_polys.readline()

if (node_cnt != num_nodes):
	print "Error in reading nodes"
#else:
	#print "Read " + str(num_nodes) + " nodes"

current_line = current_line.split()
num_polys = int(current_line[1])
current_line = f_polys.readline()
poly_id = 0
name_id = 0
polygons = []

f_poly_accept = open("accepted_poly.txt","w")

while (current_line[0:3] == 'gen' or current_line[0:4] == 'quad' or current_line[0:3] == 'tri'):
	poly_id = poly_id + 1
	if (poly_id in outer_polys):
		print "Skip reading poly no." + str(poly_id) + " as it lies on the boundary"
		if (current_line[0:3] == 'gen'): 
			current_line = f_polys.readline()
		current_line = f_polys.readline()
		continue
		
	name_id += 1
	if ((name_id % 100) == 0):
		print "Writing poly no. " + str(poly_id)

	fourdigit = str(name_id)
	while len(fourdigit) < len(str(num_polys-len(outer_polys))):
		fourdigit = '0' + fourdigit

	poly_name = 'polys/poly' + fourdigit + '.avs'
	polygon = []
	if (current_line[0:3] == 'gen'):
		current_line = f_polys.readline()
		current_line = current_line.split()
		f_out = open(poly_name,"w")
		f_out.write(str(current_line[0]) + " 0 0 0 0 \n")
		for i in range(1,int(current_line[0])+1):
			node_index = int(current_line[i]) - 1
			coordinates = node_list[node_index]
			coordinates = coordinates.split()
			coordinates = [float(x) for x in coordinates[0:]]
			#print coordinates
			polygon.append(coordinates)
			f_out.write(str(i) + node_list[node_index])
		f_out.flush()
		f_out.close()
	else:
		current_line = current_line.split()
		f_out = open(poly_name,"w")
		f_out.write(str(current_line[1]) + " 0 0 0 0 \n")
		for i in range(1,int(current_line[1])+1):
			node_index = int(current_line[i+1]) - 1
			coordinates = node_list[node_index]
			coordinates = coordinates.split()
			coordinates = [float(x) for x in coordinates[0:]]
			#print coordinates
			polygon.append(coordinates)
			f_out.write(str(i) + node_list[node_index])
		f_out.flush()
		f_out.close()
	polygons.append(polygon)
	f_poly_accept.write(fourdigit + " " + str(poly_id) + "\n")
	current_line = f_polys.readline()

f_poly_accept.flush()
f_poly_accept.close()

num_valid_polys = len(polygons)

if (poly_id != num_polys):
	print "Error in reading information for polygons"
	print "Number of polygons, poly_id: " + str(num_polys) + " " + str(poly_id)
else:
	print "Done spliting " + str(num_polys) + " polygons"
	print "Total number of valid polygons: " + str(num_valid_polys)

os.system("date")

#dist_to_centroids = FindMaxDistList(centroid_list, polygons)
#print max_dist_list
#elevation_norm_const = [MakeElevation(x) for x in max_dist_list]
#max_elevation = max(elevation_norm_const)
#glob_min_dist = min(dist_to_centroids)

if (model_type == 'FineScale'):
	os.chdir(current_dir + "/polys/")
	file_list = glob.glob("poly*.avs")
	lg_merge = "la_merge_all.lgi"
	f = open(lg_merge,"w")
	mat_id = 0
	for j in file_list:
		mat_id += 1
		f.write("cmo / create / mo / / / triplane \n")
		tmp = "read / avs / " + j + " / mo \n"
		f.write(tmp)
		f.write("cmo / setatt / mo / imt / 1,0,0 / 1 \n")
		f.write("cmo / setatt / mo / itp / 1,0,0 / 0 \n")
		f.write("triangulate / counterclockwise \n")
		f.write("recon 0 \n")
		#print "Material id" + str(mat_id)
		tmp = "cmo / setatt / mo / imt1 / 1,0,0 / " + str(mat_id) + " \n"
		f.write(tmp)
		tmp = "cmo / setatt / mo / itetclr / 1,0,0 / " + str(mat_id) + " \n"
		f.write(tmp)
		#f.write("dump / gmv / " + j[0:len(j)-4] + ".gmv / mo \n")
		f.write("addmesh / merge / mo_merge / mo_merge / mo \n")
		f.write("cmo / delete / mo \n")
		f.write("\n")
		#if (i == num_files):
		#	mat_id  = mat_id + 1

	f.write("cmo / select / mo_merge \n")
	f.write("filter \n")
	f.write("rmpoint / compress \n")
	#f.write("resetpts / itp \n")
	#f.write("dump / gmv / test_full_poly.gmv / mo_all \n")
	#f.write("\n")
	f.write("dump / gmv / tmp_merge.gmv / mo_merge \n")
	f.write(" \n")
	f.write("massage / " + str(refine_scale*5) + " / 1.e-5 / 1.e-5 / strictmergelength \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("\n")
	f.write("massage / " + str(refine_scale*3) + " / 1.e-5 / 1.e-5 / strictmergelength \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("\n")
	f.write("massage / " + str(refine_scale*1.5) + " / 1.e-5 / 1.e-5 / strictmergelength \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("\n")
	f.write("massage / " + str(refine_scale) + " / 1.e-5 / 1.e-5 / strictmergelength \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("smooth / position / esug / 1,0,0; recon 0; \n")
	f.write("\n")
	f.write("dump / gmv / tmp_tri_massage.gmv / mo_merge \n")
	f.write("\n")
	f.write("finish \n")
	f.write("\n")
	f.flush()
	f.close()

	print "Merging polygons"
	os.system(lagrit + " <la_merge_all.lgi >merge_all.txt")
	print "Done merging polygons!!!"

	f = open("calc_elevation.mlgi","w")
	f.write("define / A0 / -2.4 \n")
	f.write("define / A1 / 6 \n")
	f.write("define / A2 / -6.1 \n")
	f.write("define / A3 / 2.8 \n")
	f.write("define / A4 / -0.12 \n")
	f.write("define / A5 / 0.061 \n")
	f.write("define / A6 / 0.0 \n")
	f.write(" \n")
	f.write("math / add / mo / ele_norm / 1,0,0 / mo / ele_norm / A6 \n")
	f.write(" \n")
	f.write("math / power / mo / tmp / 1,0,0 / mo / dnorm / 1 \n")
	f.write("math / multiply / mo / tmp / 1,0,0 / mo / tmp / A5 \n")
	f.write("math / add / mo / ele_norm / 1,0,0 / mo / ele_norm / mo / tmp \n")
	f.write(" \n")
	f.write("math / power / mo / tmp / 1,0,0 / mo / dnorm / 2 \n")
	f.write("math / multiply / mo / tmp / 1,0,0 / mo / tmp / A4 \n")
	f.write("math / add / mo / ele_norm / 1,0,0 / mo / ele_norm / mo / tmp \n")
	f.write(" \n")
	f.write("math / power / mo / tmp / 1,0,0 / mo / dnorm / 3 \n")
	f.write("math / multiply / mo / tmp / 1,0,0 / mo / tmp / A3 \n")
	f.write("math / add / mo / ele_norm / 1,0,0 / mo / ele_norm / mo / tmp \n")
	f.write(" \n")
	f.write("math / power / mo / tmp / 1,0,0 / mo / dnorm / 4 \n")
	f.write("math / multiply / mo / tmp / 1,0,0 / mo / tmp / A2 \n")
	f.write("math / add / mo / ele_norm / 1,0,0 / mo / ele_norm / mo / tmp \n")
	f.write(" \n")
	f.write("math / power / mo / tmp / 1,0,0 / mo / dnorm / 5 \n")
	f.write("math / multiply / mo / tmp / 1,0,0 / mo / tmp / A1 \n")
	f.write("math / add / mo / ele_norm / 1,0,0 / mo / ele_norm / mo / tmp \n")
	f.write(" \n")
	f.write("math / power / mo / tmp / 1,0,0 / mo / dnorm / 6 \n")
	f.write("math / multiply / mo / tmp / 1,0,0 / mo / tmp / A0 \n")
	f.write("math / add / mo / ele_norm / 1,0,0 / mo / ele_norm / mo / tmp \n")
	f.write("\n")
	#f.write("cmo / copyatt / mo / mo / elevation_norm / elevation \n")
	#f.write("\n")
	#f.write("math / multiply / mo / elevation / 1,0,0 / mo / elevation / NORM_CONST \n")
	#f.write(" \n")
	f.write("finish \n")
	f.write(" \n")
	f.flush()
	f.close()

	EPS = 1.e-6
	max_dist_list = []
	for i in range(1,len(polygons)+1):
		fourdigit = str(i)
		while len(fourdigit) < len(str(num_valid_polys)):
			fourdigit = '0' + fourdigit
		
		file_name_dfield = "get_dfield" + fourdigit + ".lgi"
		f_dfield_out = "dfield" + fourdigit + ".inp"
		f = open(file_name_dfield,"w")
		f.write("read / gmv / tmp_tri_massage.gmv / mo \n")
		f.write("eltset / edelete / itetclr / ne / " + str(i) + " \n")
		f.write("rmpoint / element / eltset get edelete \n")
		f.write("rmpoint / compress \n")
		f.write("\n")
		f.write("cmo / copy / mo_refine / mo \n")
		f.write("refine2d \n")
		f.write("refine2d \n")
		f.write("resetpts / itp \n")
		f.write("\n")
		f.write("dump / avs / poly_refine" + fourdigit + ".inp / mo_refine \n")
		f.write("cmo / delete / mo_refine \n")
		f.write("\n")
		f.write("read / avs / poly_refine" + fourdigit + ".inp / mo_refine \n")
		f.write("extract / surfmesh / 1,0,0 / mo_perimeter / mo_refine / external  \n")
		f.write("\n")
		f.write("compute / distance_field / mo / mo_perimeter / dfield \n")
		f.write("define / EPS / " + str(EPS) + " \n")
		f.write("pset / p_perimeter / attribute / dfield / 1,0,0 / lt / EPS \n")
		f.write("cmo / setatt / mo / dfield / pset get p_perimeter / 0.0 \n")
		f.write("dump / avs / poly" + fourdigit + ".inp / mo \n")
		f.write("cmo / modatt / mo / imt / ioflag / l \n")
		f.write("cmo / modatt / mo / itp / ioflag / l \n")
		f.write("cmo / modatt / mo / icr / ioflag / l \n")
		f.write("cmo / modatt / mo / isn / ioflag / l \n")
		f.write("\n")
		f.write("dump / avs / " + f_dfield_out + " / mo / 0 0 2 0 \n")
		f.write("\n")
		f.write("finish \n")
		f.write("\n")
		f.flush()
		f.close()
		os.system(lagrit + "<" + file_name_dfield + " >get_dfield_la_out" + fourdigit + ".txt")
		print "Output dfield file for polygon" + fourdigit
		print f_dfield_out
		this_max_dist = GetDfieldMax(f_dfield_out)
		max_dist_list.append(GetDfieldMax(f_dfield_out))
	
	#print(max_dist_list)
	elevation_norm_const = [MakeElevation(x) for x in max_dist_list]
	#max_elevation = max(elevation_norm_const)
	for i in range(1,len(polygons)+1):
		fourdigit = str(i)
		while len(fourdigit) < len(str(num_valid_polys)):
			fourdigit = '0' + fourdigit
	
		file_name_ele = "make_elevation" + fourdigit + ".lgi"
		f = open(file_name_ele,"w")
		f.write("read / avs / poly" + fourdigit + ".inp / mo \n")
		f.write("cmo / addatt / mo / dmax / vdouble / scalar / nnodes \n")
		f.write("cmo / setatt / mo / dmax / 1,0,0 / " + str(max_dist_list[i-1]) + " \n")
		f.write("cmo / addatt / mo / dnorm / vdouble / scalar / nnodes \n")
		f.write("cmo / setatt / mo / dnorm / 1,0,0 / 0.0 \n")
		f.write("cmo / addatt / mo / tmp / vdouble / scalar / nnodes \n")
		f.write("cmo / setatt / mo / tmp / 1,0,0 / 0.0 \n")
		f.write("cmo / addatt / mo / elevation / vdouble / scalar / nnodes \n")
		f.write("cmo / setatt / mo / elevation / 1,0,0 / 0.0 \n")
		f.write("cmo / addatt / mo / ele_norm / vdouble / scalar / nnodes \n")
		f.write("cmo / setatt / mo / ele_norm / 1,0,0 / 0.0 \n")
		f.write("\n")
		f.write("math / divide / mo / dnorm / 1,0,0 / mo / dfield / mo / dmax \n")
		#f.write("math / multiply / mo / dnorm / 1,0,0 / mo / dnorm / -1
		#f.write("math / add / mo / dnorm / 1,0,0 / mo / dnorm / 1
		#f.write("define / NORM_CONST / " + str(elevation_norm_const[i-1]) + " \n")
		#f.write("define / NORM_CONST / " + str(max_elevation) + " \n")
		f.write("infile calc_elevation.mlgi \n")
		f.write("\n")
		f.write("cmo / printatt / mo / ele_norm / minmax \n")
		f.write("define / ELE_MIN / 4.2 \n")
		f.write("define / ELE_MAX / " + str(elevation_norm_const[i-1]) + " \n")	
		f.write("cmo / addatt / mo / ele_max / vdouble / scalar / nnodes \n")
		f.write("cmo / setatt / mo / ele_max / 1,0,0 / ELE_MAX \n")
		f.write("cmo / addatt / mo / ele_min / vdouble / scalar / nnodes \n")
		f.write("cmo / setatt / mo / ele_min / 1,0,0 / ELE_MIN \n")
		f.write("cmo / addatt / mo / ediff / vdouble / scalar / nnodes \n")
		f.write("math / sub / mo / ediff / 1,0,0 / mo / ele_max / mo / ele_min \n")
		f.write("cmo / printatt / mo / ediff / minmax \n")
		f.write("math / multiply / mo / elevation / 1,0,0 / mo / ele_norm / mo / ediff \n")
		f.write("cmo / printatt / mo / elevation / minmax \n")
		f.write("math / add / mo / elevation / 1,0,0 / mo / elevation / mo / ele_min \n")
		f.write("cmo / printatt / mo / elevation / minmax \n")
		f.write("cmo / copyatt / mo / mo / zic_save / zic \n")
		f.write("cmo / copyatt / mo / mo / zic / elevation \n")
		f.write("\n")
		f.write("cmo / DELATT / mo / tmp \n")
		f.write("cmo / printatt / mo / ele_norm / minmax \n")
		f.write("dump / gmv / poly_surf" + fourdigit + ".gmv / mo \n")
		f.write("\n")
		f.write("cmo / modatt / mo / imt / ioflag / l \n")
		f.write("cmo / modatt / mo / itp / ioflag / l \n")
		f.write("cmo / modatt / mo / icr / ioflag / l \n")
		f.write("cmo / modatt / mo / isn / ioflag / l \n")
		f.write("cmo / modatt / mo / dmax / ioflag / l \n")
		f.write("cmo / modatt / mo / ele_norm / ioflag / l \n")
		f.write("cmo / modatt / mo / dnorm / ioflag / l \n")
		f.write("dump / avs / elevation_profile" + fourdigit + ".inp / mo / 0 0 2 0 \n")
		f.write("\n")
		f.write("finish \n")
		f.write("\n")
		f.flush()
		f.close()
		os.system(lagrit + "<" + file_name_ele + ">make_elevation_la_out" + fourdigit + ".txt")

	file_name = "final_merge.lgi"
	f = open(file_name,"w")
	for i in range(1, len(polygons)+1):
		fourdigit = str(i)
		while len(fourdigit) < len(str(num_valid_polys)):
			fourdigit = '0' + fourdigit

		f.write("read / gmv / poly_surf" + fourdigit + ".gmv / mo \n")
		f.write("addmesh / merge / mo_merge / mo_merge / mo \n")
		f.write("cmo / delete / mo \n")
		f.write("\n")

	f.write("cmo / select / mo_merge \n")
	f.write("filter \n")
	f.write("rmpoint / compress \n")
	f.write("dump / gmv / final_polygons.gmv / mo_merge \n")
	f.write("finish \n")
	f.write("\n")
	f.flush()
	f.close()

	os.system(lagrit + "<final_merge.lgi >final_merge_la_out.txt")
	os.system("date")
	
	# Clean up results
	os.chdir(current_dir)
	os.system("mkdir " + output_folder)
	os.system("mkdir " + output_folder + "/input/")
	os.system("mkdir " + output_folder + "/temporary_files/")
	os.system("mv tmp_pd_pts.gmv " + output_folder + "/pd_point_distribution.gmv")
	os.system("mv tmp_pts_tri.gmv " + output_folder + "/triangulation.gmv")
	os.system("mv tmp_pts_tri-dual.gmv " + output_folder + "/original_polygons.gmv")
	os.system("mv polys/final_polygons.gmv " + output_folder + "/final_polygons.gmv")
	os.system("cp input/input.txt " + output_folder + "/input/")
	os.system("mv polys/ boundary_nodes.vertexset accepted_poly.txt tmp_pd_pts.inp r_distr.txt temporary_files/* " + output_folder + "/temporary_files/")
	
elif (model_type == 'IntermediateScale'):
	if (num_layers > 0):
		#Write a shell script to offset polygons
		f_shell = open("offset_polygon.sh","w")
		f_shell.write("#!/bin/csh -f \n")
		f_shell.write("echo 'Offsetting Polygons' \n")
		f_shell.write("cd ./polys/ \n")
		f_shell.write("ln -s " + current_dir + "/bin/offset_polygon \n")
		f_shell.write("ln -s " + current_dir + "/input/params.txt \n")
		f_shell.write("foreach file (*.avs) \n")
		#f_shell.write("   echo $file \n")
		f_shell.write("   ./offset_polygon $file params.txt >$file.out \n")
		f_shell.write("end \n")
		f_shell.write("exit \n")
		f_shell.flush()
		f_shell.close()

		#change mode and execute the script
		print "================================="
		os.system("chmod u+rx offset_polygon.sh")
		os.system("./offset_polygon.sh")


	#change working directory to 'polys'
	os.chdir(current_dir + "/polys/")

	#Now write LaGriT control to merge each layer
	num_files = num_layers + 1
	#set_all_polys = set(range(1,num_polys+1))
	#set_outer_polys = set()
	for i in range(1,num_files + 1):
		#print "Writing lagrit control file for layer no." + str(i)
		if (num_layers > 0):
			layer_name = "*sim" + str(i) + ".inp"
		else:
			layer_name = "poly*.avs"
		file_list = glob.glob(layer_name)
		lg_merge_layer = "la_merge_layer" + str(i) + ".lgi"
		f = open(lg_merge_layer,"w")
		mat_id = i
		for j in file_list:
			f.write("cmo / create / mo / / / triplane \n")
			tmp = "read / avs / " + j + " / mo \n"
			f.write(tmp)
			#f.write("dump / gmv / " + j[0:len(j)-4] + ".gmv / mo \n")
			f.write("triangulate / counterclockwise \n")
			f.write("recon 0 \n")
			tmp = "cmo / setatt / mo / imt1 / 1,0,0 / " + str(mat_id) + " \n"
			f.write(tmp)
			tmp = "cmo / setatt / mo / itetclr / 1,0,0 / " + str(mat_id) + " \n"
			f.write(tmp)
			f.write("dump / gmv / " + j[0:len(j)-4] + "_tri.gmv / mo \n")
			f.write("addmesh / merge / mo_merge / mo_merge / mo \n")
			f.write("cmo / delete / mo \n")
			f.write("\n")
			#if (i == num_files):
			#	mat_id  = mat_id + 1

		f.write("cmo / select / mo_merge \n")
		f.write("filter \n")
		f.write("rmpoint / compress \n")
		f.write("resetpts / itp \n")
		#f.write("dump / gmv / test_full_poly.gmv / mo_all \n")
		#f.write("\n")
		f.write("dump / gmv / tmp_merge_layer" + str(i) + ".gmv / mo_merge \n")
		f.write(" \n")
		f.write("finish \n")
		f.write("\n")
		f.flush()
		f.close()
		#print "Done writing LaGriT merge file for layer no." + str(i)

	# Execute LaGriT to merge the each layer
	os.system("date")
	for i in range(1,num_files+1):
		print "Merging layer no." + str(i)
		os.system(lagrit + "<la_merge_layer" + str(i) + ".lgi >log_merge_layer" + str(i) + ".txt")
		print "Done merging layer no." + str(i)
	os.system("date")

	# Write the LaGriT control file to merge all the layers and refine
	f = open("la_merge_all.lgi","w")
	for i in range(1,num_files+1):
		tmp = "read / gmv / tmp_merge_layer" + str(i) + ".gmv / mo \n"
		f.write(tmp)
		f.write("addmesh / merge / mo_all / mo_all / mo \n")
		f.write("\n")
	f.write("filter \n")
	f.write("rmpoint / compress \n")
	f.write("resetpts / parents \n")
	f.write("\n")
	f.write("dump / gmv / proto_poly_all.gmv / mo_all \n")
	f.write("\n")
	f.write("finish \n")
	f.write("\n")
	f.flush()
	f.close()

	# Execute LaGriT to merge all layers
	print "================================="
	command = lagrit + "<la_merge_all.lgi > log_merge_all.txt"
	os.system(command)
	os.system("mv proto_poly_all.gmv ..")
	print "Final output is 'proto_poly_all.gmv'"
	os.system("date")

	# Clean up results
	os.chdir(current_dir)
	os.system("mkdir " + output_folder)
	os.system("mkdir " + output_folder + "/input/")
	os.system("mkdir " + output_folder + "/temporary_files/")
	os.system("mv tmp_pd_pts.gmv " + output_folder + "/pd_point_distribution.gmv")
	os.system("mv tmp_pts_tri.gmv " + output_folder + "/triangulation.gmv")
	os.system("mv tmp_pts_tri-dual.gmv " + output_folder + "/original_polygons.gmv")
	os.system("mv proto_poly_all.gmv " + output_folder + "/final_polygons.gmv")
	os.system("cp input/input.txt input/params.txt " + output_folder + "/input/")
	os.system("mv polys/ tmp_pd_pts.inp boundary_nodes.vertexset accepted_poly.txt r_distr.txt offset_polygon.sh temporary_files/* " + output_folder + "/temporary_files/")

		
	
	





			

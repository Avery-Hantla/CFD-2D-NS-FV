/////////////////////////////////////////////////////////
//               Function To Read In Mesh
/////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
//#include <numeric>
//#include <algorithm>

#include "../support/class_mesh.hpp"
#include "../support/struct_size.hpp"
#include "../support/struct_BC.hpp"
#include "../support/struct_inputs.hpp"

#include "../cgns-read/meshReader.hpp"

void readmesh(class_mesh* mesh, std::string grid_file, struct_size* size, struct_BC* BC, struct_inputs* inputs) {
  int num_points, num_faces, num_cells;

  std::ifstream grid_in;
  grid_in.open(grid_file); // Open Grid File
  if (grid_in.is_open()) {
    double d1, d2;
    int i1, i2;
    mesh->num_of_BC = 0;
    size->numBC1 = 0;
    size->numBC2 = 0;

    // Read the number of points, faces, and cell
    grid_in >> num_points >> num_faces >> num_cells;
    std::cout << "Reading " << grid_file << std::endl;
    std::cout << num_points << " Points \n";
    std::cout << num_faces << " Faces \n";
    std::cout << num_cells << " Cells \n"; 

    // Assign vectors size
    mesh->cell_faces.assign (6*num_cells,-101);
    mesh->cell_face_count.assign (num_cells,0);
    mesh->face2BCf.assign(num_faces, -101);

    // Read Points
    for (int idx = 0; idx < num_points; idx++) {
      grid_in >> d1 >> d2;
      mesh->x.push_back(d1);
      mesh->y.push_back(d2);
    }

    // Read Faces 
    for (int idx = 0; idx < num_faces; idx++) {
      grid_in >> i1 >> i2; 
      mesh->face_point1.push_back(i1-1);
      mesh->face_point2.push_back(i2-1);

    }

    // Read Cells 
    for (int idx = 0; idx < num_faces; idx++) {
      grid_in >> i1 >> i2; 
      mesh->face_cell1.push_back(i1-1);
      if (i2-1 == -1) {
        mesh->face_cell2.push_back(BC->BC1);
        size->numBC1++;
      } else if (i2-1 == -2) {
        mesh->face_cell2.push_back(BC->BC2);
        size->numBC2++;
      } else {
        mesh->face_cell2.push_back(i2-1);
      }

      // Assign Faces to Cell Array 
      mesh->cell_faces[(i1-1)*6+(mesh->cell_face_count[(i1-1)])] = idx;
      if (i2 > 0) {
        mesh->cell_faces[(i2-1)*6+(mesh->cell_face_count[(i2-1)])] = idx;
        mesh->cell_face_count[(i2-1)] = mesh->cell_face_count[(i2-1)]+=1;
      }
      // Count number of faces in each cell 
      mesh->cell_face_count[(i1-1)] = mesh->cell_face_count[(i1-1)]+=1;

      if (i2-1 < 0) {
        mesh->BC_faces.push_back(idx);
        mesh->face2BCf[idx] = mesh->num_of_BC;
        mesh->num_of_BC ++;
      }
    }

    // Count of many of each element is in the mesh
    int tri = 0, quad = 0, pent = 0, hex = 0;
    for (int idx = 0; idx < num_cells; idx ++ ) {
      if (mesh->cell_face_count[idx] == 3 ) {
        tri++;
      } else if (mesh->cell_face_count[idx] == 4 ) { 
        quad++;
      } else if (mesh->cell_face_count[idx] == 5 ) {
        pent++;
      } else if (mesh->cell_face_count[idx] == 6 ) {
        hex++;
      } else {
        std::cout << "ERROR: CELL DOES NOT HAVE 3-6 SIDES\n";
      }
    }
    // Display Mesh Metric
    std::cout << "Elements: \n" << "  Triangles: " << tri << std::endl;
    std::cout << "  Quads: " << quad << std::endl;
    std::cout << "  Pentagons: " << pent << std::endl;
    std::cout << "  Hexagons:" << hex << std::endl;

    grid_in.close(); // Close Grid File
  } else {
    std::cout << "ERROR: Cannot Open Grid File\n";
  }

  // Find the face coordinates
  std::vector<double> face_coords;
  int point1, point2; 
  double x1, y1, x2, y2, area, dx, dy;
  for (int idx = 0; idx < num_faces; idx++) {
    point1 = mesh->face_point1[idx];
    point2 = mesh->face_point2[idx];
    x1 = mesh->x[point1];
    y1 = mesh->y[point1];
    x2 = mesh->x[point2];
    y2 = mesh->y[point2];

    dx = x2-x1;
    dy = y2-y1;

    // Find the area of the faces
    area = sqrt( dx*dx + dy*dy );
    mesh->face_area.push_back(area);

    // Find normal vectors for faces
    mesh->face_nx.push_back(dy/area);
    mesh->face_ny.push_back(-dx/area);

    // Find mid coordinate of faces
    mesh->face_centerx.push_back(dx/2+x1);
    mesh->face_centery.push_back(dy/2+y1);
  }

  double meanx, meany, temp[12], temp1x, temp2x, temp1y, temp2y;
  double center_x, center_y;
  std::vector<double> save_cell_pointsx, save_cell_pointsy;
  save_cell_pointsx.assign(num_cells*6,-101);
  save_cell_pointsy.assign(num_cells*6,-101);
  std::vector<double> cell_pointsx, cell_pointsy;
  for (int idx = 0; idx < num_cells; idx++) {
    int vertex_count = mesh->cell_face_count[idx], found_count = 0;
    cell_pointsx.assign(vertex_count, 0);
    cell_pointsy.assign(vertex_count, 0);
    // Find points for each cell and makes sures there are no duplicates 
    for (int jdx = 0; jdx < 6; jdx ++) { // Split cell up
      bool foundtemp1 = false, foundtemp2 = false;
      if (mesh->find_cell_face(idx,jdx) == -101) {
        break;
      }
      temp1x = mesh->x[mesh->face_point1[mesh->cell_faces[idx*6 + jdx]]];
      temp1y = mesh->y[mesh->face_point1[mesh->cell_faces[idx*6 + jdx]]];

      temp2x = mesh->x[mesh->face_point2[mesh->cell_faces[idx*6 + jdx]]];
      temp2y = mesh->y[mesh->face_point2[mesh->cell_faces[idx*6 + jdx]]];

      for (int zdx = 0; zdx < vertex_count; zdx ++) {
        if ((cell_pointsx[zdx] == temp1x) && (cell_pointsy[zdx] == temp1y)) {
          foundtemp1 = true;
          break;
        }
      }
      if (foundtemp1 == false) {
        cell_pointsx[found_count] = temp1x;
        cell_pointsy[found_count] = temp1y;
        found_count ++;
      }
      for (int zdx = 0; zdx < vertex_count; zdx ++) {
        if ((cell_pointsx[zdx] == temp2x) && (cell_pointsy[zdx] == temp2y)) {
          foundtemp2 = true;
          break;
        }
      }
      if (foundtemp2 == false) {
        cell_pointsx[found_count] = temp2x;
        cell_pointsy[found_count] = temp2y;
        found_count ++;
      }
    }

    // Get mean center of cell 
    double sum_centerx = 0, sum_centery = 0;
    for (int jdx = 0; jdx < vertex_count; jdx++){
      save_cell_pointsx[idx*6+jdx] = cell_pointsx[jdx];
      save_cell_pointsy[idx*6+jdx] = cell_pointsy[jdx];
      sum_centerx += cell_pointsx[jdx];
      sum_centery += cell_pointsy[jdx];
    }
    center_x = sum_centerx/vertex_count;
    center_y = sum_centery/vertex_count;

    double vol_cell_splits[vertex_count], x1, cent_sub_x[vertex_count], cent_sub_y[vertex_count];
    double x1x, y1x, x2x, y2x, x3x, y3x;
    double sum_vol = 0, sum_centroidx = 0, sum_centroidy = 0;
    for (int jdx = 0; jdx < mesh->cell_face_count[idx]; jdx ++) { // Split cell up
      // Find volume of each cell component
      x1x = mesh->x[mesh->face_point1[mesh->cell_faces[idx*6 + jdx]]];
      y1x = mesh->y[mesh->face_point1[mesh->cell_faces[idx*6 + jdx]]];

      x2x = mesh->x[mesh->face_point2[mesh->cell_faces[idx*6 + jdx]]];
      y2x = mesh->y[mesh->face_point2[mesh->cell_faces[idx*6 + jdx]]];

      x3x = center_x;
      y3x = center_y;

      vol_cell_splits[jdx] = 0.5*std::abs( x1x*(y2x-y3x) + x2x*(y3x-y1x) + x3x*(y1x-y2x) );
      sum_vol += vol_cell_splits[jdx];

      // Find centroid of each sub elements
      cent_sub_x[jdx] = (x1x+x2x+x3x)/3;
      cent_sub_y[jdx] = (y1x+y2x+y3x)/3;

      sum_centroidx += cent_sub_x[jdx]* vol_cell_splits[jdx];
      sum_centroidy += cent_sub_y[jdx]* vol_cell_splits[jdx];
    }
    // Find volume of the cell
    mesh->cell_vol.push_back(sum_vol);

    // Find the centroid of entire cell
    mesh->cell_centerx.push_back(sum_centroidx/sum_vol);
    mesh->cell_centery.push_back(sum_centroidy/sum_vol);
  }

  // Calculate Variables for Gradient
  // Find ghoast cell center 
  mesh->BC_cell_centerx.assign(mesh->num_of_BC,-101);
  mesh->BC_cell_centery.assign(mesh->num_of_BC,-101);
  double norm_temp, rx_temp, ry_temp, face_num;
  for (int idx = 0; idx < mesh->num_of_BC; idx ++) {
    face_num = mesh->BC_faces[idx];

    rx_temp = mesh->face_centerx[face_num] - mesh->cell_centerx[mesh->face_cell1[face_num]];
    ry_temp = mesh->face_centery[face_num] - mesh->cell_centery[mesh->face_cell1[face_num]];
    norm_temp = 2.0*(rx_temp * mesh->face_nx[face_num] + ry_temp * mesh->face_ny[face_num]);

    mesh->BC_cell_centerx[idx] = mesh->cell_centerx[mesh->face_cell1[face_num]] + norm_temp*mesh->face_nx[face_num];
    mesh->BC_cell_centery[idx] = mesh->cell_centery[mesh->face_cell1[face_num]] + norm_temp*mesh->face_ny[face_num];
  }

  // Find distance from face center to cell center
  mesh->face_cell1_dx.assign(num_faces, -101);
  mesh->face_cell1_dy.assign(num_faces, -101);
  mesh->face_cell2_dx.assign(num_faces, -101);
  mesh->face_cell2_dy.assign(num_faces, -101);
  int BC_num;
  for (int idx = 0; idx < num_faces; idx ++) {
    
    mesh->face_cell1_dx[idx] = mesh->face_centerx[idx] - mesh->cell_centerx[mesh->face_cell1[idx]];
    mesh->face_cell1_dy[idx] = mesh->face_centery[idx] - mesh->cell_centery[mesh->face_cell1[idx]];

    if (mesh->face_cell2[idx] >= 0) {
      mesh->face_cell2_dx[idx] = mesh->face_centerx[idx] - mesh->cell_centerx[mesh->face_cell2[idx]];
      mesh->face_cell2_dy[idx] = mesh->face_centery[idx] - mesh->cell_centery[mesh->face_cell2[idx]];
    } else {
      BC_num = mesh->face2BCf[idx];
      mesh->face_cell2_dx[idx] = mesh->face_centerx[idx] - mesh->BC_cell_centerx[BC_num];
      mesh->face_cell2_dy[idx] = mesh->face_centery[idx] - mesh->BC_cell_centery[BC_num];
    }
  }

  // Find delta
  int cell1, cell2;
  std::vector<double> del_temp1, del_temp2, del_temp3;
  del_temp1.assign(num_cells,0);
  del_temp2.assign(num_cells,0);
  del_temp3.assign(num_cells,0);

  double xi, xj, yi, yj;
  mesh->delta.assign(num_cells, -101);
  mesh->face_dxj.assign(num_faces, -101);
  mesh->face_dyj.assign(num_faces, -101);
  for (int idx = 0; idx < num_faces; idx ++) {
    cell1 = mesh->face_cell1[idx];
    cell2 = mesh->face_cell2[idx];

    xi = mesh->cell_centerx[cell1];
    yi = mesh->cell_centery[cell1];

    if (cell2 >= 0) {
      xj = mesh->cell_centerx[cell2];
      yj = mesh->cell_centery[cell2];
    } else {
      BC_num = mesh->face2BCf[idx];
      xj = mesh->BC_cell_centerx[BC_num];
      yj = mesh->BC_cell_centery[BC_num];
    }

    del_temp1[cell1] += (xi - xj)*(xi - xj);
    del_temp2[cell1] += (yi - yj)*(yi - yj);
    del_temp3[cell1] += ((xi - xj)*(yi - yj));

    if (cell2 >= 0) {
      del_temp1[cell2] += (xj - xi)*(xj - xi);
      del_temp2[cell2] += (yj - yi)*(yj - yi);
      del_temp3[cell2] += ((xj - xi)*(yj - yi));
    }

    mesh->face_dxj[idx] = xj-xi;
    mesh->face_dyj[idx] = yj-yi;
  }

  for (int idx = 0; idx < num_cells; idx ++) {
    mesh->delta[idx] = del_temp1[idx]*del_temp2[idx] - (del_temp3[idx]*del_temp3[idx]);
  }

  // Find Ixx, Iyy, Ixy
  mesh->Ixx.assign(num_cells, -101);
  mesh->Iyy.assign(num_cells, -101);
  mesh->Ixy.assign(num_cells, -101);
  for (int idx = 0; idx < num_cells; idx ++) {
    mesh->Ixx[idx] = del_temp1[idx]/mesh->delta[idx];
    mesh->Iyy[idx] = del_temp2[idx]/mesh->delta[idx];
    mesh->Ixy[idx] = del_temp3[idx]/mesh->delta[idx];
  }

  // // Cell Centroids
  // std::ofstream temp_centroids;
  // temp_centroids.open ("mesh/cell_data_ghost.dat");
  // if (temp_centroids.is_open()) {
  //   temp_centroids << "# cell_centerx, cell_centery, x points for cell, y points for cell" << std::endl;
  //   for (int idx = 0; idx < mesh->num_of_BC; idx++) {
  //     temp_centroids << mesh->BC_cell_centerx[idx] << " " << mesh->BC_cell_centery[idx] << std::endl;   
  //   }
  // } else {
  //   std::cout << "ERROR: Cannot Save File \n";
  // }
  // temp_centroids.close();

  // Update Size Struct 
  size->num_points = num_points;
  size->num_faces = num_faces;
  size->num_cells = num_cells;

  ////////// Make Conectivity Matrix for Tecplot and CGNS Outputs  //////////
  // Make conectivity matrix
  for (int idx = 0; idx < size->num_cells; idx ++) {
    for (int jdx = 0; jdx < mesh->cell_face_count[idx]; jdx ++) {
      int face_num = mesh->find_cell_face(idx, jdx);

      int point1 = mesh->face_point1[face_num];
      int point2 = mesh->face_point2[face_num];

      mesh->connect_out.push_back(point1+1);
      mesh->connect_out.push_back(point2+1);
      mesh->connect_out.push_back(size->num_points + idx + 1);
    }
  }

  // Set number of wall elements
  if ((BC->BC1 == -2) || (BC->BC1 == -4)) {
    size->numWALL = size->numBC1;
  } else if ((BC->BC2 == -2) || (BC->BC2 == -4)) {
    size->numWALL = size->numBC2;
  }

  // Find the distances between cell centers
  mesh->face_dist.assign(num_faces, -101);
  double tempx, tempy;
  for (int idx = 0; idx < num_faces; idx ++) {
    int cell1 = mesh->face_cell1[idx];
    int cell2 = mesh->face_cell2[idx];

    if (cell2 >= 0) {
      tempx = mesh->cell_centerx[cell2] - mesh->cell_centerx[cell1];
      tempy = mesh->cell_centery[cell2] - mesh->cell_centery[cell1];
    } else {
      face_num = mesh->BC_faces[idx];
      tempx = mesh->BC_cell_centerx[face_num] - mesh->cell_centerx[cell1];
      tempy = mesh->BC_cell_centery[face_num] - mesh->cell_centery[cell1];
    }

    mesh->face_dist[idx] = std::sqrt(tempx*tempx + tempy*tempy);
  }

  mesh->face_lx.assign(num_faces, -101);
  mesh->face_ly.assign(num_faces, -101);
  mesh->face_mx.assign(num_faces, -101);
  mesh->face_my.assign(num_faces, -101);
  for (int idx = 0; idx < num_faces; idx ++) {
    // Find the lx and ly of the face
    int cell1 = mesh->face_cell1[idx];
    int cell2 = mesh->face_cell2[idx];

    if (cell2 >= 0) {
      mesh->face_lx[idx] = mesh->cell_centerx[cell2] - mesh->cell_centerx[cell1];
      mesh->face_ly[idx] = mesh->cell_centery[cell2] - mesh->cell_centery[cell1];
    } else {
      face_num = mesh->BC_faces[idx];
      mesh->face_lx[idx] = mesh->BC_cell_centerx[face_num] - mesh->cell_centerx[cell1];
      mesh->face_ly[idx] = mesh->BC_cell_centery[face_num] - mesh->cell_centery[cell1];
    }

    // Find the mx and my of the face
    mesh->face_mx[idx] = (mesh->x[mesh->face_point2[idx]] - mesh->x[mesh->face_point1[idx]]);
    mesh->face_my[idx] = (mesh->y[mesh->face_point2[idx]] - mesh->y[mesh->face_point1[idx]]);

    // Normalize
    double temp = std::sqrt(mesh->face_mx[idx]*mesh->face_mx[idx] + mesh->face_my[idx]*mesh->face_my[idx]);
    mesh->face_mx[idx] = mesh->face_mx[idx]/temp;
    mesh->face_my[idx] = mesh->face_my[idx]/temp;

    double temp2 = std::sqrt(mesh->face_lx[idx]*mesh->face_lx[idx] + mesh->face_ly[idx]*mesh->face_ly[idx]);
    mesh->face_lx[idx] = mesh->face_lx[idx]/temp2;
    mesh->face_ly[idx] = mesh->face_ly[idx]/temp2;
  }

  /////////////////////// Save metric metrics to file ///////////////////////

  // Cell Centroids
  std::ofstream output;
  output.open ("mesh/cell_data.dat");
  if (output.is_open()) {
    output << "# cell_vol, cell_centerx, cell_centery, x points for cell, y points for cell" << std::endl;
    for (int idx = 0; idx < size->num_cells; idx++) {
      output << mesh->cell_vol[idx] << " " << mesh->cell_centerx[idx] << " " << mesh->cell_centery[idx] << " ";   
      output << save_cell_pointsx[idx*6] << " " << save_cell_pointsx[idx*6+1] << " " << save_cell_pointsx[idx*6+2] << " " << save_cell_pointsx[idx*6+3] << " " << save_cell_pointsx[idx*6+4] << " " << save_cell_pointsx[idx*6+5] << " ";
      output << save_cell_pointsy[idx*6] << " " << save_cell_pointsy[idx*6+1] << " " << save_cell_pointsy[idx*6+2] << " " << save_cell_pointsy[idx*6+3] << " " << save_cell_pointsy[idx*6+4] << " " << save_cell_pointsy[idx*6+5] << std::endl;
    }
  } else {
    std::cout << "ERROR: Cannot Save File \n";
  }
  output.close();

  // Faces
  output.open ("mesh/face_data.dat");
  if (output.is_open()) {
    output << "# face point 1 x, face point 1 y, face point 2 x, face point 2y, face nx, face ny, face center x, face center y" << std::endl;
    for (int idx = 0; idx < size->num_faces; idx++) {
      output << mesh->x[mesh->face_point1[idx]] << " " <<mesh->y[mesh->face_point1[idx]] << " "; 
      output << mesh->x[mesh->face_point2[idx]] << " " <<mesh->y[mesh->face_point2[idx]] << " "; 

      output << mesh->face_nx[idx] << " " << mesh->face_ny[idx] <<  " "; 

      output << mesh->face_centerx[idx] << " " << mesh->face_centery[idx] << std::endl; 
    }
  } else {
    std::cout << "ERROR: Cannot Save File \n";
  }
  output.close();

}
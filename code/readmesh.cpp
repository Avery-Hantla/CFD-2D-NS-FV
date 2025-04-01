/////////////////////////////////////////////////////////
//               Function To Read In Mesh
/////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
//#include <numeric>
//#include <algorithm>

#include "class_mesh.hpp"
#include "struct_size.hpp"
#include "struct_BC.hpp"

void readmesh(class_mesh* mesh, std::string grid_file, struct_size* size, struct_BC* BC) {
  int num_points, num_faces, num_cells;

  std::ifstream grid_in;
  grid_in.open(grid_file); // Open Grid File
  if (grid_in.is_open()) {
    double d1, d2;
    int i1, i2;
    // Read the number of points, faces, and cell
    grid_in >> num_points >> num_faces >> num_cells;
    std::cout << "Reading " << grid_file << std::endl;
    std::cout << num_points << " Points \n";
    std::cout << num_faces << " Faces \n";
    std::cout << num_cells << " Cells \n"; 

    // Assign vectors size
    mesh->cell_faces.assign (6*num_cells,-101);
    mesh->cell_face_count.assign (num_cells,0);

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
      } else if (i2-1 == -2) {
        mesh->face_cell2.push_back(BC->BC2);
      } else {
        mesh->face_cell2.push_back(i2-1);
      }

      // Assign Faces to Cell Array 
      mesh->cell_faces[(i1-1)*6+(mesh->cell_face_count[(i1-1)])] = idx;
      if ( (i2 != -1) && (i2 != 0) ) {
        mesh->cell_faces[(i2-1)*6+(mesh->cell_face_count[(i2-1)])] = idx;
        mesh->cell_face_count[(i2-1)] = mesh->cell_face_count[(i2-1)]+=1;
      }
      // Count number of faces in each cell 
      mesh->cell_face_count[(i1-1)] = mesh->cell_face_count[(i1-1)]+=1;
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
  bool foundtemp1, foundtemp2;
  for (int idx = 0; idx < num_cells; idx++) {
    int vertex_count = mesh->cell_face_count[idx];
    double cell_pointsx[vertex_count], cell_pointsy[vertex_count];
    // Find points for each cell and makes sures there are no duplicates 
    for (int jdx = 0; jdx < 6; jdx ++) { // Split cell up
      if (mesh->find_cell_face(idx,jdx) == -101) {
        break;
      }
      temp1x = mesh->x[mesh->face_point1[mesh->cell_faces[idx*6 + jdx]]];
      temp1y = mesh->y[mesh->face_point1[mesh->cell_faces[idx*6 + jdx]]];

      temp2x = mesh->x[mesh->face_point2[mesh->cell_faces[idx*6 + jdx]]];
      temp2y = mesh->y[mesh->face_point2[mesh->cell_faces[idx*6 + jdx]]];

      for (int zdx = 0; zdx < vertex_count; zdx ++) {
        if ((cell_pointsx[zdx] = temp1x) && (cell_pointsy[zdx] = temp1y)) {
          foundtemp1 = true;
        }
      }
      if (foundtemp1 == false) {
        cell_pointsx[jdx] = temp1x;
        cell_pointsy[jdx] = temp1y;
      }
      for (int zdx = 0; zdx < vertex_count; zdx ++) {
        if ((cell_pointsx[zdx] = temp2x) && (cell_pointsy[zdx] = temp2y)) {
          foundtemp2 = true;
        }
      }
      if (foundtemp2 == false) {
        cell_pointsx[jdx] = temp2x;
        cell_pointsy[jdx] = temp2y;
      }
    }

    // Get mean center of cell 
    double sum_centerx = 0, sum_centery = 0;
    for (int jdx = 0; jdx < vertex_count; jdx++){
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

  // // Find the rx and ry vectors
  // double mag_temp, dot_temp, cell1, cell2;
  // mesh->init(num_faces);
  // for (int idx = 0; idx < num_faces; idx++) {
  //     // Cell 1
  //     cell1 = mesh->face_cell1[idx];

  //     mesh->face_cell1_rx[idx] = mesh->face_centerx[idx]-mesh->cell_centerx[cell1];
  //     mesh->face_cell1_ry[idx] = mesh->face_centery[idx]-mesh->cell_centery[cell1];

  //     mag_temp = std::sqrt(mesh->face_cell1_rx[idx] * mesh->face_cell1_rx[idx] + mesh->face_cell1_ry[idx] * mesh->face_cell1_ry[idx]);
  //     mesh->face_cell1_rx_norm[idx] = (mesh->face_cell1_rx[idx])/mag_temp; 
  //     mesh->face_cell1_ry_norm[idx] = (mesh->face_cell1_ry[idx])/mag_temp;
    
  //     dot_temp = (mesh->face_cell1_rx_norm[idx] * mesh->face_centerx[cell1] + mesh->face_cell1_ry_norm[idx] * mesh->face_centery[cell1]);

  //     // Check if normal vector goes in or out of cell
  //     if (dot_temp > 0) { // Out of face
  //       mesh->face_cell1_n_out[idx] = 1;
  //     } else if (dot_temp < 0) { // Into cell
  //       mesh->face_cell1_n_out[idx] = -1;
  //     } else {
  //       std::cout << "Error: Norm Doesnt Add Up, Cell 1: readmesh.cpp\n";
  //     }

  //     // Cell 2
  //     cell2 = mesh->face_cell2[idx];

  //     if ((cell2 == BC->freestream_patch) || (cell2 == BC->wall_patch)) {
  //       continue;
  //     }

  //     mesh->face_cell2_rx[idx] = mesh->face_centerx[idx]-mesh->cell_centerx[cell2];
  //     mesh->face_cell2_ry[idx] = mesh->face_centery[idx]-mesh->cell_centery[cell2];

  //     mag_temp = std::sqrt(mesh->face_cell2_rx[idx] * mesh->face_cell2_rx[idx] + mesh->face_cell2_ry[idx] * mesh->face_cell2_ry[idx]);
  //     mesh->face_cell2_rx_norm[idx] = (mesh->face_cell2_rx[idx])/mag_temp; 
  //     mesh->face_cell2_ry_norm[idx] = (mesh->face_cell2_ry[idx])/mag_temp;
    
  //     dot_temp = (mesh->face_cell2_rx_norm[idx] * mesh->face_centerx[cell1] + mesh->face_cell2_ry_norm[idx] * mesh->face_centery[cell1]);

  //     // Check if normal vector goes in or out of cell
  //     if (dot_temp > 0) { // Out of face
  //       mesh->face_cell2_n_out[idx] = 1;
  //     } else if (dot_temp < 0) { // Into cell
  //       mesh->face_cell2_n_out[idx] = -1;
  //     } else {
  //       std::cout << "Error: Norm Doesnt Add Up, Cell 2: readmesh.cpp\n";
  //     }
  // }

  // Find the rx and ry vectors
  int face_num;
  double mag_temp, dot_temp;
  mesh->cell_face_rx.assign(num_cells*6,-101);
  mesh->cell_face_ry.assign(num_cells*6,-101);
  mesh->cell_face_rx_norm.assign(num_cells*6,-101);
  mesh->cell_face_ry_norm.assign(num_cells*6,-101);
  mesh->cell_face_n_out.assign(num_cells*6,-101);
  for (int idx = 0; idx < num_cells; idx++) {
    for (int jdx = 0; jdx < mesh->cell_face_count[idx]; jdx++) {
      if (mesh->cell_face_count[idx] == -101) {
        break;
      }
      face_num = mesh->find_cell_face(idx,jdx);

      mesh->cell_face_rx[6*idx + jdx] = mesh->face_centerx[face_num]-mesh->cell_centerx[idx];
      mesh->cell_face_ry[6*idx + jdx] = mesh->face_centery[face_num]-mesh->cell_centery[idx];
      
      mag_temp = std::sqrt(mesh->cell_face_rx[6*idx + jdx] * mesh->cell_face_rx[6*idx + jdx] + mesh->cell_face_ry[6*idx + jdx] * mesh->cell_face_ry[6*idx + jdx]);
      mesh->cell_face_rx_norm[6*idx + jdx] = (mesh->cell_face_rx[6*idx + jdx])/mag_temp; 
      mesh->cell_face_ry_norm[6*idx + jdx] = (mesh->cell_face_ry[6*idx + jdx])/mag_temp;
    
      dot_temp = (mesh->cell_face_rx_norm[6*idx + jdx] * mesh->face_nx[face_num] + mesh->cell_face_ry_norm[6*idx + jdx] * mesh->face_ny[face_num]);
      // Check if normal vector goes in or out of cell
      if (dot_temp > 0) { // Out of face
        mesh->cell_face_n_out[6*idx+jdx] = 1;
      } else if (dot_temp < 0) { // Into cell
        mesh->cell_face_n_out[6*idx+jdx] = -1;
      } else {
        std::cout << "Error: Norm Doesnt Add Up: readmesh.cpp\n";
      }
    }
  }

  // Find the corresponding Qi face num for the face cell 1 and 2
  // double cell_num;
  // for (int idx = 0; idx < num_faces; idx ++) {
  //   cell_num = mesh->face_cells[2*idx]; // Find cell number of cell 1 of face idx
  //   for (int jdx = 0; jdx < mesh->cell_face_count[cell_num]; jdx ++) { // Search cell for face num
  //     if (mesh->find_cell_face(cell_num, jdx) == idx ) {
  //       mesh->face_cell1_Qface_num.push_back(6*cell_num+jdx);
  //       break;
  //     }
  //   }

  //   cell_num = mesh->face_cells[2*idx+1];
  //   if (cell_num == -1) {
  //     mesh->face_cell2_Qface_num.push_back(-1);
  //   } else if (cell_num == -2) {
  //     mesh->face_cell2_Qface_num.push_back(-2);
  //   } else {
  //     for (int jdx = 0; jdx < mesh->cell_face_count[cell_num]; jdx ++) {
  //       if (mesh->find_cell_face(cell_num, jdx) == idx ) {
  //         mesh->face_cell2_Qface_num.push_back(6*cell_num + jdx);
  //         break; 
  //       }
  //     }
  //   }
  // }

  // Update Size Struct 
  size->num_points = num_points;
  size->num_faces = num_faces;
  size->num_cells = num_cells;

  // Save metric metrics to file

  std::ofstream output;
  output.open ("../mesh/mesh.dat");
  if (output.is_open()) {
    for (int idx = 0; idx < size->num_cells; idx++) {
      output << mesh->cell_centerx[idx] << " " << mesh->cell_centery[idx] << std::endl;    }
  } else {
    std::cout << "ERROR: Cannot Save File \n";
  }
  output.close();
}
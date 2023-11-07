// Copyright (C) 2007-2023  CEA, EDF, OPEN CASCADE
//
// Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
// CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

//  File   : GMSHPlugin_Runner_main.cxx
//  Author : Cesar Conopoima, Open Cascade
//  Module : GMSH
//

#include "GMSHPlugin_GMSH_3D_SA.hxx"

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <chrono>

/**
 * @brief Main function
 *
 * @param argc Number of arguments
 * @param argv Arguments
 *
 * @return error code
 */
int main(int argc, char *argv[]){

  if(argc!=8||(argc==2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help")==0))){
    std::cout << "Error in number of arguments "<< argc-1<<" given expected 7" <<std::endl;
    std::cout << "Syntax:"<<std::endl;
    std::cout << "run_mesher MESHER INPUT_MESH_FILE SHAPE_FILE HYPO_FILE" << std::endl;
    std::cout << "           ELEM_ORIENT_FILE " << std::endl;
    std::cout << "           NEW_ELEMENT_FILE OUTPUT_MESH_FILE" << std::endl;
    std::cout << std::endl;
    std::cout << " Set argument to NONE to ignore them " << std::endl;
    std::cout << std::endl;
    std::cout << "Args:" << std::endl;
    std::cout << "  MESHER: mesher to use from (NETGEN3D, NETGEN2D)" << std::endl;
    std::cout << "  INPUT_MESH_FILE: MED File containing lower-dimension-elements already meshed" << std::endl;
    std::cout << "  SHAPE_FILE: STEP file containing the shape to mesh" << std::endl;
    std::cout << "  HYPO_FILE: Ascii file containint the list of parameters" << std::endl;
    std::cout << "  (optional) ELEM_ORIENT_FILE: binary file containing the list of element from INPUT_MESH_FILE associated to the shape and their orientation" << std::endl;
    std::cout << "  (optional) NEW_ELEMENT_FILE: (out) contains elements and nodes added by the meshing" << std::endl;
    std::cout << "  (optional) OUTPUT_MESH_FILE: (out) MED File containing the mesh after the run of the mesher" << std::endl;
    return 1;
  }
  std::string mesher=argv[1];
  std::string input_mesh_file=argv[2];
  std::string shape_file=argv[3];
  std::string hypo_file=argv[4];
  std::string element_orientation_file=argv[5];
  std::string new_element_file=argv[6];
  std::string output_mesh_file=argv[7];

  //std::string thing;
  //std::cin >> thing;

  if (output_mesh_file == "NONE")
    output_mesh_file = "";
  if (element_orientation_file == "NONE")
    element_orientation_file = "";
  if (new_element_file == "NONE")
    new_element_file = "";

  std::cout << mesher << "\n";
  std::cout << input_mesh_file << "\n";
  std::cout << shape_file << "\n";
  std::cout << hypo_file << "\n";
  std::cout << element_orientation_file << "\n";
  std::cout << new_element_file << "\n";
  std::cout << output_mesh_file << "\n";
  
  if (mesher=="GMSH3D"){
    GMSHPlugin_GMSH_3D_SA myplugin;
    std::cout << input_mesh_file << "\n";
    std::cout << shape_file << "\n";
    std::cout << hypo_file << "\n";
    std::cout << new_element_file << "\n";
    std::cout << output_mesh_file << "\n";
    int run = myplugin.run(input_mesh_file,
                            shape_file,
                            hypo_file,
                            element_orientation_file,
                            new_element_file,
                            output_mesh_file );
  }  
  else {
    std::cerr << "Unknown mesher:" << mesher << std::endl;
    return 1;
  }
  return 0;
}

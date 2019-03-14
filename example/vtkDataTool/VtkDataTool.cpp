/**
 * @file   VtkDataTool.h
 * @author Feng Han <fengerhan.@gmail.com>
 * @date   Fri Mar 11 13:36:59 2017
 * 
 * @brief  本程序主要为了解决 AFEPack 中二维问题的物理量的显示问题。AFEPack 可以将二维问题中
 *         的物理量输出为 easymesh 网格文件及其网格节点上的标量值。本程序用于将 easymesh网格
 *         文件，以及网格节点对应的标量值转化成 vtk 或 vtu 文件，并且可以将多个时间点下
 *         的 vtu 文件保存为 paraview 的 pvd 文件，从而实现在 paraview 下显示相关计算结果
 *         和结果动画显示。
 *         例如：
 *              EasymeshData2vtk("data/init_u0");
 *              会生成init_u0.vtk。easymesh是二维网格。如果把网格点上的物理量当作第三个个坐标，
 *              则可以生成 surface 图。调用如下函数可以生成 surface 图。
 *              EasymeshData2vtk3D("data/init_u0");
 *
 *              .vtk文件格式是ASC码文件，为了生成pvd格式的动画，需要xml格式的vtk文件。
 *               下面两个函数用于生成 xml 格式的 vtu 文件
 *
 *              EasymeshData2vtu("data/init_u0");
 *              EasymeshData2vtu3D("data/init_u0");
 * 
 *              pvd 文件的格式，请参考自行在网上搜索参考。当有了一系列 vtu 文件以及这些文件
 *              对应的时间后，调用
 *              creatpvd( time_t, file_vtu, pvdfilename);
 *              
 *              就可以生成 pvd 文件了，其中 time_t 是时间序列数组，file_vtu 是 vtu 文件
 *              名序列数组，pvdfilename 是生成的 pvd 文件名。
 *          
 *              由于要生成 xml 格式文件，我在这里使用了 tinyxml2 软件包。
 *
 *              添加了将AFEPack 二维问题中的物理量直接输出为vtu文件的函数，这样更方便了．
 *         　　　例如
 *　　　　　　　　　　　　　　　　　　　　WritevtuData(FEMFunction<double, 2>& u_h, const char * ofilename, u_int flag = 0)
 *              将输出标量函数 u_h 的函数值，当 flag = 1 时，将输出　u_h 的梯度值，目前
*　　　　　　　　　这个函数只针对２维问题，后续会添加对３维函数的处理
 */
 
#ifndef __VtkDataTool_h__
#define __VtkDataTool_h__

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>

#include <sstream>
#include "tinyxml2.h"

#include <AFEPack/DGFEMSpace.h>

#include <AFEPack/HGeometry.h>

#include "VtkDataTool.h"

using namespace tinyxml2;

int WritevtpData(RegularMesh<2>& mesh, std::vector<double> celldata, const char * ofilename)
{
		///构造 vtu 文件的 xml 模板
		XMLDocument* doc1 = new XMLDocument();
				
		XMLDeclaration * declaration = doc1->NewDeclaration();
		doc1->InsertEndChild(declaration);

/**     添加xml文件头的第二种方法
		const char* str = "<?xml version=\"1.0\" ?>";
		doc1->Parse( str );
*/
	
		XMLElement * element_VTKFile,
				   * element_PolyData,
				   * element_Piece,
		           * element_Points,
		           * element_Point_DataArray,
		           * element_CellData,
		           * element_Cells_DataArray,
		           * element_Polys,
		           * element_Polys_DataArray;
		           
		XMLNode * node_VTKFile, 
		        * node_PolyData,
		        * node_Piece,
		        * node_Points,
		        * node_Point_DataArray,
		        * node_CellData,
		        * node_Cells_DataArray,
		        * node_Polys,
		        * node_Polys_DataArray_1,
		        * node_Polys_DataArray_2;

		
		element_VTKFile = doc1->NewElement("VTKFile");
		element_VTKFile->SetAttribute( "type", "PolyData" );
		node_VTKFile = doc1->InsertEndChild( element_VTKFile );
		
		node_PolyData = node_VTKFile->InsertEndChild( doc1->NewElement( "PolyData" ) );
		element_Piece = doc1->NewElement( "Piece" );
		node_Piece = node_PolyData->InsertEndChild( element_Piece );

		element_Points = doc1->NewElement( "Points" );
		node_Points = node_Piece->InsertEndChild( element_Points );

		element_Point_DataArray = doc1->NewElement( "DataArray" );
		element_Point_DataArray->SetAttribute( "type", "Float32" );	
		element_Point_DataArray->SetAttribute( "NumberOfComponents", 3 );	
		element_Point_DataArray->SetAttribute( "format", "ascii" );	
		node_Point_DataArray = node_Points->InsertFirstChild( element_Point_DataArray );
		
		element_CellData = doc1->NewElement( "CellData" );
		element_CellData->SetAttribute( "Scalars","cell_scalars" );	
		node_CellData = node_Piece->InsertEndChild( element_CellData );

		element_Cells_DataArray = doc1->NewElement( "DataArray" );
		element_Cells_DataArray->SetAttribute( "type", "Float32" );	
		element_Cells_DataArray->SetAttribute( "Name", "cell_scalars" );	
		element_Cells_DataArray->SetAttribute( "format", "ascii" );	
		node_Cells_DataArray = node_CellData->InsertEndChild( element_Cells_DataArray );

		element_Polys = doc1->NewElement( "Polys" );
		node_Polys = node_Piece->InsertEndChild( element_Polys );

		element_Polys_DataArray = doc1->NewElement( "DataArray" );
		element_Polys_DataArray->SetAttribute( "type", "Int32" );	
		element_Polys_DataArray->SetAttribute( "Name", "connectivity" );	
		element_Polys_DataArray->SetAttribute( "format", "ascii" );	
		node_Polys_DataArray_1 = node_Polys->InsertEndChild( element_Polys_DataArray );
		
		element_Polys_DataArray = doc1->NewElement( "DataArray" );
		element_Polys_DataArray->SetAttribute( "type", "Int32" );	
		element_Polys_DataArray->SetAttribute( "Name", "offsets" );	
		element_Polys_DataArray->SetAttribute( "format", "ascii" );	
		node_Polys_DataArray_2 = node_Polys->InsertEndChild( element_Polys_DataArray );

   /// vtu 文件模板构造完成，开始写入相关数据     

 		const std::string& vtufilename(ofilename);
   		
	    std::ostringstream s1;
		std::string str1,str2;
	    s1.precision(15);
	    s1.setf(std::ios::scientific, std::ios::floatfield);

   		std::cout << "Creat vtp file..." << std::endl;	
   		std::cout << "\tWriting node data ..." << std::flush;
   
        int n_node = mesh.n_point();

        int n_element = mesh.n_geometry(2);

		element_Piece->SetAttribute( "NumberOfPoints", n_node );	
		element_Piece->SetAttribute( "NumberOfPolys", n_element );	

	    for (int i = 0; i < n_node; i++)
	    {
		      s1.str("");
		      s1  << std::endl << mesh.point(i)[0] << "  " << mesh.point(i)[1] << "  " << "0.0";
		      str1 = s1.str();
 			  node_Point_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));
	    }
 

	    double temp = *std::max_element(celldata.begin(), celldata.end());

        std::cout << " OK!" << std::endl;
   		std::cout << "\tWriting element data ..." << std::flush;

		int skip = 0;
	    for (int i = 0 ; i < n_element; ++i)
   	    {
           const GeometryBM& geo = mesh.geometry(2,i);
           switch (geo.n_vertex())
           {
              case 3: 
 		         s1.str("");
		         s1  << std::endl
		             << mesh.geometry(0, geo.vertex(0)).vertex(0) << "  "
		             << mesh.geometry(0, geo.vertex(1)).vertex(0) << "  " 
		             << mesh.geometry(0, geo.vertex(2)).vertex(0) << std::endl;
		         str1 = s1.str();
			     node_Polys_DataArray_1->InsertEndChild( doc1->NewText( str1.c_str() ));

 		         skip = skip + 3;
 		         s1.str("");
		         s1  << std::endl << skip;
		         str1 = s1.str();
			     node_Polys_DataArray_2->InsertEndChild( doc1->NewText( str1.c_str() ));
                 break;

              case 4:
 		         s1.str("");
		         s1  << std::endl
		             << mesh.geometry(0, geo.vertex(0)).vertex(0) << "  "
		             << mesh.geometry(0, geo.vertex(1)).vertex(0) << "  " 
		             << mesh.geometry(0, geo.vertex(2)).vertex(0) << "  " 
		             << mesh.geometry(0, geo.vertex(3)).vertex(0) << std::endl;
		         str1 = s1.str();
			     node_Polys_DataArray_1->InsertEndChild( doc1->NewText( str1.c_str() ));
 		         
 		         skip = skip + 4;
 		         s1.str("");
		         s1  << std::endl << skip;
		         str1 = s1.str();
			     node_Polys_DataArray_2->InsertEndChild( doc1->NewText( str1.c_str() ));
                 break;

              default:
                 Assert(false, ExcInternalError());
           }

           s1.str("");
           s1  << std::endl
	           << celldata[i] << "  "
	           << std::endl;
	       str1 = s1.str();
	       node_Cells_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));
        }

		doc1->SaveFile( (vtufilename + ".vtp").c_str() );
		std::cout << " OK!" << std::endl;
		std::cout << "\tcreat vtp file: " << (vtufilename + ".vtp").c_str() << " OK!" << std::endl;
        
	    delete doc1;
	
	    return 0;
}

int WritevtpData3D(RegularMesh<2>& mesh, std::vector<double> celldata, const char * ofilename)
{
		///构造 vtu 文件的 xml 模板
		XMLDocument* doc1 = new XMLDocument();
				
		XMLDeclaration * declaration = doc1->NewDeclaration();
		doc1->InsertEndChild(declaration);

/**     添加xml文件头的第二种方法
		const char* str = "<?xml version=\"1.0\" ?>";
		doc1->Parse( str );
*/
	
		XMLElement * element_VTKFile,
				   * element_PolyData,
				   * element_Piece,
		           * element_Points,
		           * element_Point_DataArray,
		           * element_CellData,
		           * element_Cells_DataArray,
		           * element_Polys,
		           * element_Polys_DataArray;
		           
		XMLNode * node_VTKFile, 
		        * node_PolyData,
		        * node_Piece,
		        * node_Points,
		        * node_Point_DataArray,
		        * node_CellData,
		        * node_Cells_DataArray,
		        * node_Polys,
		        * node_Polys_DataArray_1,
		        * node_Polys_DataArray_2;

		
		element_VTKFile = doc1->NewElement("VTKFile");
		element_VTKFile->SetAttribute( "type", "PolyData" );
		node_VTKFile = doc1->InsertEndChild( element_VTKFile );
		
		node_PolyData = node_VTKFile->InsertEndChild( doc1->NewElement( "PolyData" ) );
		element_Piece = doc1->NewElement( "Piece" );
		node_Piece = node_PolyData->InsertEndChild( element_Piece );

		element_Points = doc1->NewElement( "Points" );
		node_Points = node_Piece->InsertEndChild( element_Points );

		element_Point_DataArray = doc1->NewElement( "DataArray" );
		element_Point_DataArray->SetAttribute( "type", "Float32" );	
		element_Point_DataArray->SetAttribute( "NumberOfComponents", 3 );	
		element_Point_DataArray->SetAttribute( "format", "ascii" );	
		node_Point_DataArray = node_Points->InsertFirstChild( element_Point_DataArray );
		
		element_CellData = doc1->NewElement( "CellData" );
		element_CellData->SetAttribute( "Scalars","cell_scalars" );	
		node_CellData = node_Piece->InsertEndChild( element_CellData );

		element_Cells_DataArray = doc1->NewElement( "DataArray" );
		element_Cells_DataArray->SetAttribute( "type", "Float32" );	
		element_Cells_DataArray->SetAttribute( "Name", "cell_scalars" );	
		element_Cells_DataArray->SetAttribute( "format", "ascii" );	
		node_Cells_DataArray = node_CellData->InsertEndChild( element_Cells_DataArray );

		element_Polys = doc1->NewElement( "Polys" );
		node_Polys = node_Piece->InsertEndChild( element_Polys );

		element_Polys_DataArray = doc1->NewElement( "DataArray" );
		element_Polys_DataArray->SetAttribute( "type", "Int32" );	
		element_Polys_DataArray->SetAttribute( "Name", "connectivity" );	
		element_Polys_DataArray->SetAttribute( "format", "ascii" );	
		node_Polys_DataArray_1 = node_Polys->InsertEndChild( element_Polys_DataArray );
		
		element_Polys_DataArray = doc1->NewElement( "DataArray" );
		element_Polys_DataArray->SetAttribute( "type", "Int32" );	
		element_Polys_DataArray->SetAttribute( "Name", "offsets" );	
		element_Polys_DataArray->SetAttribute( "format", "ascii" );	
		node_Polys_DataArray_2 = node_Polys->InsertEndChild( element_Polys_DataArray );

   /// vtu 文件模板构造完成，开始写入相关数据     

 		const std::string& vtufilename(ofilename);
   		
	    std::ostringstream s1;
		std::string str1,str2;
	    s1.precision(15);
	    s1.setf(std::ios::scientific, std::ios::floatfield);

   		std::cout << "Creat vtp file..." << std::endl;	
   		std::cout << "\tWriting node data ..." << std::flush;
   
        int n_node = mesh.n_point();

        int n_element = mesh.n_geometry(2);


/*	    for (int i = 0; i < n_node; i++)
	    {
		      s1.str("");
		      s1  << std::endl << mesh.point(i)[0] << "  " << mesh.point(i)[1] << "  " << "0.0";
		      str1 = s1.str();
 			  node_Point_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));
	    }
*/ 
        std::cout << " OK!" << std::endl;
   		std::cout << "\tWriting element data ..." << std::flush;
		
	    double temp = *std::max_element(celldata.begin(), celldata.end());

		int skip = 0;
	    for (int i = 0 ; i < n_element; ++i)
   	    {
           const GeometryBM& geo = mesh.geometry(2,i);
           switch (geo.n_vertex())
           {
              case 3:
		     	 s1.str("");
		      	 s1  << std::endl << mesh.point(mesh.geometry(0, geo.vertex(0)).vertex(0))[0] << "  " 
		      	                  << mesh.point(mesh.geometry(0, geo.vertex(0)).vertex(0))[1] << "  " 
		      	                  << celldata[i]/temp;
		      	 str1 = s1.str();
 			  	 node_Point_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));

		     	 s1.str("");
		      	 s1  << std::endl << mesh.point(mesh.geometry(0, geo.vertex(1)).vertex(0))[0] << "  " 
		      	                  << mesh.point(mesh.geometry(0, geo.vertex(1)).vertex(0))[1] << "  " 
		      	                  << celldata[i]/temp;
		      	 str1 = s1.str();
 			  	 node_Point_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));

		     	 s1.str("");
		      	 s1  << std::endl << mesh.point(mesh.geometry(0, geo.vertex(2)).vertex(0))[0] << "  " 
		      	                  << mesh.point(mesh.geometry(0, geo.vertex(2)).vertex(0))[1] << "  " 
		      	                  << celldata[i]/temp;
		      	 str1 = s1.str();
 			  	 node_Point_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));
              
 		         s1.str("");
		         s1  << std::endl
		             << skip + 0 << "  "
		             << skip + 1 << "  " 
		             << skip + 2 << std::endl;
              	 
              	 skip = skip + 3;
		         str1 = s1.str();
			     node_Polys_DataArray_1->InsertEndChild( doc1->NewText( str1.c_str() ));

 		         s1.str("");
		         s1  << std::endl << skip;
		         str1 = s1.str();
			     node_Polys_DataArray_2->InsertEndChild( doc1->NewText( str1.c_str() ));
                 break;

              case 4:
		     	 s1.str("");
		      	 s1  << std::endl << mesh.point(mesh.geometry(0, geo.vertex(0)).vertex(0))[0] << "  " 
		      	                  << mesh.point(mesh.geometry(0, geo.vertex(0)).vertex(0))[1] << "  " 
		      	                  << celldata[i]/temp;
		      	 str1 = s1.str();
 			  	 node_Point_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));

		     	 s1.str("");
		      	 s1  << std::endl << mesh.point(mesh.geometry(0, geo.vertex(1)).vertex(0))[0] << "  " 
		      	                  << mesh.point(mesh.geometry(0, geo.vertex(1)).vertex(0))[1] << "  " 
		      	                  << celldata[i]/temp;
		      	 str1 = s1.str();
 			  	 node_Point_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));

		     	 s1.str("");
		      	 s1  << std::endl << mesh.point(mesh.geometry(0, geo.vertex(2)).vertex(0))[0] << "  " 
		      	                  << mesh.point(mesh.geometry(0, geo.vertex(2)).vertex(0))[1] << "  " 
		      	                  << celldata[i]/temp;
		      	 str1 = s1.str();
 			  	 node_Point_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));

		     	 s1.str("");
		      	 s1  << std::endl << mesh.point(mesh.geometry(0, geo.vertex(3)).vertex(0))[0] << "  " 
		      	                  << mesh.point(mesh.geometry(0, geo.vertex(3)).vertex(0))[1] << "  " 
		      	                  << celldata[i]/temp;
		      	 str1 = s1.str();
 			  	 node_Point_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));
              
 		         s1.str("");
		         s1  << std::endl
		             << skip + 0 << "  "
		             << skip + 1 << "  " 
		             << skip + 2 << "  " 
		             << skip + 3 << std::endl;
              	 
              	 skip = skip + 4;
		         str1 = s1.str();
			     node_Polys_DataArray_1->InsertEndChild( doc1->NewText( str1.c_str() ));

 		         s1.str("");
		         s1  << std::endl << skip;
		         str1 = s1.str();
			     node_Polys_DataArray_2->InsertEndChild( doc1->NewText( str1.c_str() ));
                 break;

              default:
                 Assert(false, ExcInternalError());
           }

           s1.str("");
           s1  << std::endl
	           << celldata[i] << "  "
	           << std::endl;
	       str1 = s1.str();
	       node_Cells_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));
        }
		element_Piece->SetAttribute( "NumberOfPoints", skip );	
		element_Piece->SetAttribute( "NumberOfPolys", n_element );	

		doc1->SaveFile( (vtufilename + "_3D.vtp").c_str() );
		std::cout << " OK!" << std::endl;
		std::cout << "\tcreat vtp file: " << (vtufilename + "_3D.vtp").c_str() << " OK!" << std::endl;
        
	    delete doc1;
	
	    return 0;
}

int WritevtuData(FEMFunction<double, 3>& u_h, const char * ofilename, u_int flag = 0)
{
		///构造 vtu 文件的 xml 模板
		XMLDocument* doc1 = new XMLDocument();
				
		XMLDeclaration * declaration = doc1->NewDeclaration();
		doc1->InsertEndChild(declaration);

/**     添加xml文件头的第二种方法
		const char* str = "<?xml version=\"1.0\" ?>";
		doc1->Parse( str );
*/
	
		XMLElement * element_VTKFile,
		           * element_Piece,
		           * element_Points,
		           * element_point_DataArray,
		           * element_Cells,
		           * element_Cells_DataArray_1,
		           * element_Cells_DataArray_2,
		           * element_Cells_DataArray_3,
		           * element_PointData,
		           * element_PointData_DataArray;
		           
		XMLNode * node_VTKFile, 
		        * node_UnstructuredGrid,
		        * node_piece,
		        * node_point,
		        * node_point_DataArray,
		        * node_Cells,
		        * node_Cells_DataArray_1,
		        * node_Cells_DataArray_2,
		        * node_Cells_DataArray_3,
		        * node_PointData,
		        * node_PointData_DataArray;

		
		element_VTKFile = doc1->NewElement("VTKFile");
		element_VTKFile->SetAttribute( "type", "UnstructuredGrid" );
		node_VTKFile = doc1->InsertEndChild( element_VTKFile );
		
		node_UnstructuredGrid = node_VTKFile->InsertEndChild( doc1->NewElement( "UnstructuredGrid" ) );
		element_Piece = doc1->NewElement( "Piece" );
		node_piece = node_UnstructuredGrid->InsertEndChild( element_Piece );

		element_Points = doc1->NewElement( "Points" );
		node_point = node_piece->InsertEndChild( element_Points );

		element_point_DataArray = doc1->NewElement( "DataArray" );
		element_point_DataArray->SetAttribute( "type", "Float32" );	
		element_point_DataArray->SetAttribute( "NumberOfComponents", 3 );	
		element_point_DataArray->SetAttribute( "format", "ascii" );	
		node_point_DataArray = node_point->InsertFirstChild( element_point_DataArray );
		
		element_Cells = doc1->NewElement( "Cells" );
		node_Cells = node_piece->InsertEndChild( element_Cells );

		element_Cells_DataArray_1 = doc1->NewElement( "DataArray" );
		element_Cells_DataArray_1->SetAttribute( "type", "Int32" );	
		element_Cells_DataArray_1->SetAttribute( "Name", "connectivity" );	
		element_Cells_DataArray_1->SetAttribute( "format", "ascii" );	
		node_Cells_DataArray_1 = node_Cells->InsertEndChild( element_Cells_DataArray_1 );

		element_Cells_DataArray_2 = doc1->NewElement( "DataArray" );
		element_Cells_DataArray_2->SetAttribute( "type", "Int32" );	
		element_Cells_DataArray_2->SetAttribute( "Name", "offsets" );	
		element_Cells_DataArray_2->SetAttribute( "format", "ascii" );	
		node_Cells_DataArray_2 = node_Cells->InsertEndChild( element_Cells_DataArray_2 );

		element_Cells_DataArray_3 = doc1->NewElement( "DataArray" );
		element_Cells_DataArray_3->SetAttribute( "type", "UInt8" );	
		element_Cells_DataArray_3->SetAttribute( "Name", "types" );	
		element_Cells_DataArray_3->SetAttribute( "format", "ascii" );	
		node_Cells_DataArray_3 = node_Cells->InsertEndChild( element_Cells_DataArray_3 );

		element_PointData = doc1->NewElement( "PointData" );
		if ( flag == 0)
		{
		   element_PointData->SetAttribute( "Scalars", "u_h" );
		   node_PointData = node_piece->InsertEndChild( element_PointData );
		   
		   element_PointData_DataArray = doc1->NewElement( "DataArray" );
		   element_PointData_DataArray->SetAttribute( "type", "Float32" );	
		   element_PointData_DataArray->SetAttribute( "Name", "u_h" );	
		   element_PointData_DataArray->SetAttribute( "format", "ascii" );
		   node_PointData_DataArray = node_PointData->InsertEndChild( element_PointData_DataArray );
		}
		else if ( flag == 1)
		{
		   element_PointData->SetAttribute( "Vectors", "u_h" );
		   node_PointData = node_piece->InsertEndChild( element_PointData );

		   element_PointData_DataArray = doc1->NewElement( "DataArray" );
		   element_PointData_DataArray->SetAttribute( "type", "Float32" );	
		   element_PointData_DataArray->SetAttribute( "NumberOfComponents", 3 );
		   element_PointData_DataArray->SetAttribute( "Name", "u_h" );	
		   element_PointData_DataArray->SetAttribute( "format", "ascii" );
		   node_PointData_DataArray = node_PointData->InsertEndChild( element_PointData_DataArray );
		}

   /// vtu 文件模板构造完成，开始写入相关数据     


 		const std::string& vtufilename(ofilename);
   		
	    std::ostringstream s1;
		std::string str1,str2;
	    s1.precision(15);
	    s1.setf(std::ios::scientific, std::ios::floatfield);

   		std::cout << "Creat vtu file..." << std::endl;	
   		std::cout << "\tWriting node data ..." << std::flush;
   
        const FEMSpace<double,3>& fem_space = u_h.femSpace();
        const Mesh<3,3>& mesh = fem_space.mesh();
        int n_node = mesh.n_point();
 
        int n_element = 0;
        FEMSpace<double,3>::ConstElementIterator 
        the_element = fem_space.beginElement(),
        end_element = fem_space.endElement();
        for (;the_element != end_element;++ the_element)
        {
          const GeometryBM& geo = the_element->geometry();
          switch (geo.n_vertex()) 
          {
            case 4: n_element += 1; break;
            case 5: n_element += 2; break;
            case 7: n_element += 4; break;
          }
        }    		

		element_Piece->SetAttribute( "NumberOfPoints", n_node );	
		element_Piece->SetAttribute( "NumberOfCells", n_element );	


	    for (int i = 0; i < n_node; i++)
	    {
		      s1.str("");
		      s1  << std::endl
		          << mesh.point(i)[0] << "  " 
		          << mesh.point(i)[1] << "  " 
		          << mesh.point(i)[2];
		      str1 = s1.str();
 			  node_point_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));
	    }
        std::cout << "Write node OK!" << std::endl;
   		std::cout << "\tWriting element data ..." << std::flush;

        int n_data = n_node;
        std::vector<int> count(n_data, 0);
        std::vector<double> val(n_data, 0);

		int skip = 0;
        the_element = fem_space.beginElement();
	    for (int i = 0 ;the_element != end_element; ++ the_element)
   	    {
           const GeometryBM& geo = the_element->geometry();
           //std::cout << "i = " << i << "\t" << geo.n_vertex() << std::endl;
           switch (geo.n_vertex())
           {
              case 4: 
 		         s1.str("");
		         s1  << std::endl
		             << mesh.geometry(0, geo.vertex(0)).vertex(0) << "  "
		             << mesh.geometry(0, geo.vertex(1)).vertex(0) << "  "
		             << mesh.geometry(0, geo.vertex(2)).vertex(0) << "  " 
		             << mesh.geometry(0, geo.vertex(3)).vertex(0) << std::endl;
		         str1 = s1.str();
			     node_Cells_DataArray_1->InsertEndChild( doc1->NewText( str1.c_str() ));
 		         skip = skip + 4;
 		         s1.str("");
		         s1  << std::endl << skip;
		         str1 = s1.str();
			     node_Cells_DataArray_2->InsertEndChild( doc1->NewText( str1.c_str() ));
                 node_Cells_DataArray_3->InsertEndChild( doc1->NewText( "\n10\n" ));
                 i = i + 1;
                 break;

              case 5:
 		         s1.str("");
		         s1  << std::endl
		             << mesh.geometry(0, geo.vertex(0)).vertex(0) << "  "
		             << mesh.geometry(0, geo.vertex(1)).vertex(0) << "  "
		             << mesh.geometry(0, geo.vertex(2)).vertex(0) << "  " 
		             << mesh.geometry(0, geo.vertex(4)).vertex(0) << std::endl;
		         str1 = s1.str();
			     node_Cells_DataArray_1->InsertEndChild( doc1->NewText( str1.c_str() ));
 		         skip = skip + 4;
 		         s1.str("");
		         s1  << std::endl << skip;
		         str1 = s1.str();
			     node_Cells_DataArray_2->InsertEndChild( doc1->NewText( str1.c_str() ));
                 node_Cells_DataArray_3->InsertEndChild( doc1->NewText( "\n10\n" ));
 		         
 		         s1.str("");
		         s1  << std::endl
		             << mesh.geometry(0, geo.vertex(0)).vertex(0) << "  "
		             << mesh.geometry(0, geo.vertex(2)).vertex(0) << "  "
		             << mesh.geometry(0, geo.vertex(3)).vertex(0) << "  "
		             << mesh.geometry(0, geo.vertex(4)).vertex(0) << std::endl;
		         str1 = s1.str();
			     node_Cells_DataArray_1->InsertEndChild( doc1->NewText( str1.c_str() ));
 		         skip = skip + 4;
 		         s1.str("");
		         s1  << std::endl << skip;
		         str1 = s1.str();
			     node_Cells_DataArray_2->InsertEndChild( doc1->NewText( str1.c_str() ));
                 node_Cells_DataArray_3->InsertEndChild( doc1->NewText( "\n10\n" ));
                 i = i + 1;
                 break;
                 
              case 7:
 		         s1.str("");
		         s1  << std::endl
		             << mesh.geometry(0, geo.vertex(0)).vertex(0) << "  "
		             << mesh.geometry(0, geo.vertex(1)).vertex(0) << "  "
		             << mesh.geometry(0, geo.vertex(6)).vertex(0) << "  " 
		             << mesh.geometry(0, geo.vertex(5)).vertex(0) << std::endl;
		         str1 = s1.str();
			     node_Cells_DataArray_1->InsertEndChild( doc1->NewText( str1.c_str() ));
 		         skip = skip + 4;
 		         s1.str("");
		         s1  << std::endl << skip;
		         str1 = s1.str();
			     node_Cells_DataArray_2->InsertEndChild( doc1->NewText( str1.c_str() ));
                 node_Cells_DataArray_3->InsertEndChild( doc1->NewText( "\n10\n" ));
 		         
 		         s1.str("");
		         s1  << std::endl
		             << mesh.geometry(0, geo.vertex(0)).vertex(0) << "  "
		             << mesh.geometry(0, geo.vertex(2)).vertex(0) << "  "
		             << mesh.geometry(0, geo.vertex(4)).vertex(0) << "  "
		             << mesh.geometry(0, geo.vertex(6)).vertex(0) << std::endl;
		         str1 = s1.str();
			     node_Cells_DataArray_1->InsertEndChild( doc1->NewText( str1.c_str() ));
 		         skip = skip + 4;
 		         s1.str("");
		         s1  << std::endl << skip;
		         str1 = s1.str();
			     node_Cells_DataArray_2->InsertEndChild( doc1->NewText( str1.c_str() ));
                 node_Cells_DataArray_3->InsertEndChild( doc1->NewText( "\n10\n" ));

 		         s1.str("");
		         s1  << std::endl
		             << mesh.geometry(0, geo.vertex(0)).vertex(0) << "  "
		             << mesh.geometry(0, geo.vertex(3)).vertex(0) << "  "
		             << mesh.geometry(0, geo.vertex(5)).vertex(0) << "  "
		             << mesh.geometry(0, geo.vertex(4)).vertex(0) << std::endl;
		         str1 = s1.str();
			     node_Cells_DataArray_1->InsertEndChild( doc1->NewText( str1.c_str() ));
 		         skip = skip + 4;
 		         s1.str("");
		         s1  << std::endl << skip;
		         str1 = s1.str();
			     node_Cells_DataArray_2->InsertEndChild( doc1->NewText( str1.c_str() ));
                 node_Cells_DataArray_3->InsertEndChild( doc1->NewText( "\n10\n" ));

 		         s1.str("");
		         s1  << std::endl
		             << mesh.geometry(0, geo.vertex(0)).vertex(0) << "  "
		             << mesh.geometry(0, geo.vertex(4)).vertex(0) << "  "
		             << mesh.geometry(0, geo.vertex(5)).vertex(0) << "  "
		             << mesh.geometry(0, geo.vertex(6)).vertex(0) << std::endl;
		         str1 = s1.str();
			     node_Cells_DataArray_1->InsertEndChild( doc1->NewText( str1.c_str() ));
 		         skip = skip + 4;
 		         s1.str("");
		         s1  << std::endl << skip;
		         str1 = s1.str();
			     node_Cells_DataArray_2->InsertEndChild( doc1->NewText( str1.c_str() ));
                 node_Cells_DataArray_3->InsertEndChild( doc1->NewText( "\n10\n" ));
                 i = i + 1;
                 break;
                 
              default:
                 std::cout << "error ======" << std::endl;
                 Assert(false, ExcInternalError());
           }
           for (int i = 0;i < geo.n_vertex();i ++)
           {
	          int k = mesh.geometry(0, geo.vertex(i)).vertex(0);
	          count[k] += 1;
	          val[k] += u_h.value(mesh.point(k), *the_element);
           }
        }
        std::cout << "Write element OK!" << std::endl;
        for (int i = 0;i < n_node;i ++)
        {
           Assert(count[i] > 0, ExcInternalError());
           val[i] /= count[i];
        }

   		for (int i = 0; i < n_node; i++)
   		{
 		   if( flag == 0 )
		   {
 		      s1.str("");
		      s1  << std::endl << val[i];
		      str1 = s1.str();
		      node_PointData_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));
		   }
		   else if( flag == 1 )
		   {
 		      s1.str("");
		      s1  << std::endl << val[i] << "0.0" << "0.0";
		      str1 = s1.str();
		      node_PointData_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));		   
		   }
		   else
		      Assert(false, ExcInternalError());
  		}

		//doc1->Print();
		doc1->SaveFile( (vtufilename + ".vtu").c_str() );
		std::cout << "\tcreat vtu file: " << (vtufilename + ".vtu").c_str() << " OK!" << std::endl;
        
	    delete doc1;
	    return 0;
}

int WritevtuData(FEMFunction<double, 2>& u_h, const char * ofilename, u_int flag = 0)
{
		///构造 vtu 文件的 xml 模板
		XMLDocument* doc1 = new XMLDocument();
				
		XMLDeclaration * declaration = doc1->NewDeclaration();
		doc1->InsertEndChild(declaration);

/**     添加xml文件头的第二种方法
		const char* str = "<?xml version=\"1.0\" ?>";
		doc1->Parse( str );
*/
	
		XMLElement * element_VTKFile,
		           * element_Piece,
		           * element_Points,
		           * element_point_DataArray,
		           * element_Cells,
		           * element_Cells_DataArray_1,
		           * element_Cells_DataArray_2,
		           * element_Cells_DataArray_3,
		           * element_PointData,
		           * element_PointData_DataArray;
		           
		XMLNode * node_VTKFile, 
		        * node_UnstructuredGrid,
		        * node_piece,
		        * node_point,
		        * node_point_DataArray,
		        * node_Cells,
		        * node_Cells_DataArray_1,
		        * node_Cells_DataArray_2,
		        * node_Cells_DataArray_3,
		        * node_PointData,
		        * node_PointData_DataArray;

		
		element_VTKFile = doc1->NewElement("VTKFile");
		element_VTKFile->SetAttribute( "type", "UnstructuredGrid" );
		node_VTKFile = doc1->InsertEndChild( element_VTKFile );
		
		node_UnstructuredGrid = node_VTKFile->InsertEndChild( doc1->NewElement( "UnstructuredGrid" ) );
		element_Piece = doc1->NewElement( "Piece" );
		node_piece = node_UnstructuredGrid->InsertEndChild( element_Piece );

		element_Points = doc1->NewElement( "Points" );
		node_point = node_piece->InsertEndChild( element_Points );

		element_point_DataArray = doc1->NewElement( "DataArray" );
		element_point_DataArray->SetAttribute( "type", "Float32" );	
		element_point_DataArray->SetAttribute( "NumberOfComponents", 3 );	
		element_point_DataArray->SetAttribute( "format", "ascii" );	
		node_point_DataArray = node_point->InsertFirstChild( element_point_DataArray );
		
		element_Cells = doc1->NewElement( "Cells" );
		node_Cells = node_piece->InsertEndChild( element_Cells );

		element_Cells_DataArray_1 = doc1->NewElement( "DataArray" );
		element_Cells_DataArray_1->SetAttribute( "type", "Int32" );	
		element_Cells_DataArray_1->SetAttribute( "Name", "connectivity" );	
		element_Cells_DataArray_1->SetAttribute( "format", "ascii" );	
		node_Cells_DataArray_1 = node_Cells->InsertEndChild( element_Cells_DataArray_1 );

		element_Cells_DataArray_2 = doc1->NewElement( "DataArray" );
		element_Cells_DataArray_2->SetAttribute( "type", "Int32" );	
		element_Cells_DataArray_2->SetAttribute( "Name", "offsets" );	
		element_Cells_DataArray_2->SetAttribute( "format", "ascii" );	
		node_Cells_DataArray_2 = node_Cells->InsertEndChild( element_Cells_DataArray_2 );

		element_Cells_DataArray_3 = doc1->NewElement( "DataArray" );
		element_Cells_DataArray_3->SetAttribute( "type", "UInt8" );	
		element_Cells_DataArray_3->SetAttribute( "Name", "types" );	
		element_Cells_DataArray_3->SetAttribute( "format", "ascii" );	
		node_Cells_DataArray_3 = node_Cells->InsertEndChild( element_Cells_DataArray_3 );
/*
		element_PointData = doc1->NewElement( "PointData" );
		element_PointData->SetAttribute( "Scalars", "u_h" );
		node_PointData = node_piece->InsertEndChild( element_PointData );

		element_PointData_DataArray = doc1->NewElement( "DataArray" );
		element_PointData_DataArray->SetAttribute( "type", "Float32" );	
		element_PointData_DataArray->SetAttribute( "Name", "u_h" );	
		element_PointData_DataArray->SetAttribute( "format", "ascii" );
		node_PointData->InsertEndChild( element_PointData_DataArray );
*/
		element_PointData = doc1->NewElement( "PointData" );
		if ( flag == 0)
		{
		   element_PointData->SetAttribute( "Scalars", "u_h" );
		   node_PointData = node_piece->InsertEndChild( element_PointData );
		   
		   element_PointData_DataArray = doc1->NewElement( "DataArray" );
		   element_PointData_DataArray->SetAttribute( "type", "Float32" );	
		   element_PointData_DataArray->SetAttribute( "Name", "u_h" );	
		   element_PointData_DataArray->SetAttribute( "format", "ascii" );
		   node_PointData_DataArray = node_PointData->InsertEndChild( element_PointData_DataArray );
		}
		else if ( flag == 1)
		{
		   element_PointData->SetAttribute( "Vectors", "u_h" );
		   node_PointData = node_piece->InsertEndChild( element_PointData );

		   element_PointData_DataArray = doc1->NewElement( "DataArray" );
		   element_PointData_DataArray->SetAttribute( "type", "Float32" );
		   element_PointData_DataArray->SetAttribute( "NumberOfComponents", 3 );
		   element_PointData_DataArray->SetAttribute( "Name", "u_h" );	
		   element_PointData_DataArray->SetAttribute( "format", "ascii" );
		   node_PointData_DataArray = node_PointData->InsertEndChild( element_PointData_DataArray );
		}

   /// vtu 文件模板构造完成，开始写入相关数据     


 		const std::string& vtufilename(ofilename);


   		int n_node, n_element;
   		
	    std::ostringstream s1;
		std::string str1,str2;
	    s1.precision(15);
	    s1.setf(std::ios::scientific, std::ios::floatfield);

   		std::cout << "Creat vtu file..." << std::endl;	
   		std::cout << "\tWriting node data ..." << std::flush;
   
        const FEMSpace<double,2>& fem_space = u_h.femSpace();
        const Mesh<2,2>& mesh = fem_space.mesh();
        n_node = mesh.n_point();
 
        n_element = 0;
        n_element = mesh.n_geometry(2);

		element_Piece->SetAttribute( "NumberOfPoints", n_node );	
		element_Piece->SetAttribute( "NumberOfCells", n_element );	


	    for (int i = 0; i < n_node; i++)
	    {
		      s1.str("");
		      s1  << std::endl << mesh.point(i)[0] << "  " << mesh.point(i)[1] << "  " << "0.0";
		      str1 = s1.str();
 			  node_point_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));
	    }
 
        std::cout << " OK!" << std::endl;
   		std::cout << "\tWriting element data ..." << std::flush;

        int n_data = n_node;
        std::vector<int> count(n_data, 0);
        std::vector<double> val(n_data, 0);
		std::vector<double> grad0(n_data, 0);
		std::vector<double> grad1(n_data, 0);
		
        int skip = 0;		

        FEMSpace<double,2>::ConstElementIterator 
        the_element = fem_space.beginElement(),
        end_element = fem_space.endElement();
        the_element = fem_space.beginElement();
	    for (int i = 0 ;the_element != end_element; ++ the_element)
   	    {
           const GeometryBM& geo = the_element->geometry();
           switch (geo.n_vertex())
           {
              case 3: 
 		         s1.str("");
		         s1  << std::endl
		             << mesh.geometry(0, geo.vertex(0)).vertex(0) << "  "
		             << mesh.geometry(0, geo.vertex(1)).vertex(0) << "  " 
		             << mesh.geometry(0, geo.vertex(2)).vertex(0) << std::endl;
		         str1 = s1.str();
			     node_Cells_DataArray_1->InsertEndChild( doc1->NewText( str1.c_str() ));
				 skip = skip + 3;
 		         s1.str("");
		         s1  << std::endl << skip;
		         str1 = s1.str();
			     node_Cells_DataArray_2->InsertEndChild( doc1->NewText( str1.c_str() ));
                 node_Cells_DataArray_3->InsertEndChild( doc1->NewText( "\n5\n" ));
                 i = i + 1;
                 break;

              case 4:
 		         s1.str("");
		         s1  << std::endl
		             << mesh.geometry(0, geo.vertex(0)).vertex(0) << "  "
		             << mesh.geometry(0, geo.vertex(1)).vertex(0) << "  " 
		             << mesh.geometry(0, geo.vertex(2)).vertex(0) << "  " 
		             << mesh.geometry(0, geo.vertex(3)).vertex(0) << std::endl;
		         str1 = s1.str();
			     node_Cells_DataArray_1->InsertEndChild( doc1->NewText( str1.c_str() ));
 		         skip = skip + 4;
 		         s1.str("");
		         s1  << std::endl << skip;
		         str1 = s1.str();
			     node_Cells_DataArray_2->InsertEndChild( doc1->NewText( str1.c_str() ));
                 node_Cells_DataArray_3->InsertEndChild( doc1->NewText( "\n9\n" ));
                 i = i + 1;
                 break;

              default:
                 Assert(false, ExcInternalError());
           }
		   if ( flag == 0 )
		   {
	           for (int i = 0;i < geo.n_vertex();i ++)
    	       {
		          int k = mesh.geometry(0, geo.vertex(i)).vertex(0);
		          count[k] += 1;
		          val[k] += u_h.value(mesh.point(k), *the_element);
    	       }
 		   }
		   else if ( flag == 1 )
		   {
	           for (int i = 0;i < geo.n_vertex();i ++)
    	       {
		          int k = mesh.geometry(0, geo.vertex(i)).vertex(0);
		          count[k] += 1;
		          grad0[k] += u_h.gradient(mesh.point(k), *the_element)[0];
		          grad1[k] += u_h.gradient(mesh.point(k), *the_element)[1];
		       }
		   }
       }
        std::cout << " OK!" << std::endl;
        for (int i = 0;i < n_node;i ++)
        {
           Assert(count[i] > 0, ExcInternalError());
           val[i] /= count[i];
           grad0[i] /= count[i];
           grad1[i] /= count[i];
        }

   		for (int i = 0; i < n_node; i++)
   		{
		   if ( flag == 0 )
		   {
 		      s1.str("");
		      s1  << std::endl << val[i];
		      str1 = s1.str();
		      node_PointData_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));
		   }
		   else if ( flag == 1 )
		   {
 		      s1.str("");
		      s1  << std::endl << grad0[i] << "  " << grad1[i] << "  0.0";
		      str1 = s1.str();
		      node_PointData_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));		   
		   }
   		}

		//doc1->Print();
		doc1->SaveFile( (vtufilename + ".vtu").c_str() );
		std::cout << "\tcreat vtu file: " << (vtufilename + ".vtu").c_str() << " OK!" << std::endl;
        
	    delete doc1;
	
	    return 0;
}

//　按矢量分量输出矢量
int WritevtuData(FEMFunction<double, 2>& v0, FEMFunction<double, 2>& v1, const char * ofilename)
{
		///构造 vtu 文件的 xml 模板
		XMLDocument* doc1 = new XMLDocument();
				
		XMLDeclaration * declaration = doc1->NewDeclaration();
		doc1->InsertEndChild(declaration);

/**     添加xml文件头的第二种方法
		const char* str = "<?xml version=\"1.0\" ?>";
		doc1->Parse( str );
*/
	
		XMLElement * element_VTKFile,
		           * element_Piece,
		           * element_Points,
		           * element_point_DataArray,
		           * element_Cells,
		           * element_Cells_DataArray_1,
		           * element_Cells_DataArray_2,
		           * element_Cells_DataArray_3,
		           * element_PointData,
		           * element_PointData_DataArray;
		           
		XMLNode * node_VTKFile, 
		        * node_UnstructuredGrid,
		        * node_piece,
		        * node_point,
		        * node_point_DataArray,
		        * node_Cells,
		        * node_Cells_DataArray_1,
		        * node_Cells_DataArray_2,
		        * node_Cells_DataArray_3,
		        * node_PointData,
		        * node_PointData_DataArray;

		
		element_VTKFile = doc1->NewElement("VTKFile");
		element_VTKFile->SetAttribute( "type", "UnstructuredGrid" );
		node_VTKFile = doc1->InsertEndChild( element_VTKFile );
		
		node_UnstructuredGrid = node_VTKFile->InsertEndChild( doc1->NewElement( "UnstructuredGrid" ) );
		element_Piece = doc1->NewElement( "Piece" );
		node_piece = node_UnstructuredGrid->InsertEndChild( element_Piece );

		element_Points = doc1->NewElement( "Points" );
		node_point = node_piece->InsertEndChild( element_Points );

		element_point_DataArray = doc1->NewElement( "DataArray" );
		element_point_DataArray->SetAttribute( "type", "Float32" );	
		element_point_DataArray->SetAttribute( "NumberOfComponents", 3 );	
		element_point_DataArray->SetAttribute( "format", "ascii" );	
		node_point_DataArray = node_point->InsertFirstChild( element_point_DataArray );
		
		element_Cells = doc1->NewElement( "Cells" );
		node_Cells = node_piece->InsertEndChild( element_Cells );

		element_Cells_DataArray_1 = doc1->NewElement( "DataArray" );
		element_Cells_DataArray_1->SetAttribute( "type", "Int32" );	
		element_Cells_DataArray_1->SetAttribute( "Name", "connectivity" );	
		element_Cells_DataArray_1->SetAttribute( "format", "ascii" );	
		node_Cells_DataArray_1 = node_Cells->InsertEndChild( element_Cells_DataArray_1 );

		element_Cells_DataArray_2 = doc1->NewElement( "DataArray" );
		element_Cells_DataArray_2->SetAttribute( "type", "Int32" );	
		element_Cells_DataArray_2->SetAttribute( "Name", "offsets" );	
		element_Cells_DataArray_2->SetAttribute( "format", "ascii" );	
		node_Cells_DataArray_2 = node_Cells->InsertEndChild( element_Cells_DataArray_2 );

		element_Cells_DataArray_3 = doc1->NewElement( "DataArray" );
		element_Cells_DataArray_3->SetAttribute( "type", "UInt8" );	
		element_Cells_DataArray_3->SetAttribute( "Name", "types" );	
		element_Cells_DataArray_3->SetAttribute( "format", "ascii" );	
		node_Cells_DataArray_3 = node_Cells->InsertEndChild( element_Cells_DataArray_3 );
/*
		element_PointData = doc1->NewElement( "PointData" );
		element_PointData->SetAttribute( "Scalars", "u_h" );
		node_PointData = node_piece->InsertEndChild( element_PointData );

		element_PointData_DataArray = doc1->NewElement( "DataArray" );
		element_PointData_DataArray->SetAttribute( "type", "Float32" );	
		element_PointData_DataArray->SetAttribute( "Name", "u_h" );	
		element_PointData_DataArray->SetAttribute( "format", "ascii" );
		node_PointData->InsertEndChild( element_PointData_DataArray );
*/
		element_PointData = doc1->NewElement( "PointData" );
		{
		   element_PointData->SetAttribute( "Vectors", "Vector_v" );
		   node_PointData = node_piece->InsertEndChild( element_PointData );

		   element_PointData_DataArray = doc1->NewElement( "DataArray" );
		   element_PointData_DataArray->SetAttribute( "type", "Float32" );
		   element_PointData_DataArray->SetAttribute( "NumberOfComponents", 3 );
		   element_PointData_DataArray->SetAttribute( "Name", "u_h" );	
		   element_PointData_DataArray->SetAttribute( "format", "ascii" );
		   node_PointData_DataArray = node_PointData->InsertEndChild( element_PointData_DataArray );
		}

   /// vtu 文件模板构造完成，开始写入相关数据     


 		const std::string& vtufilename(ofilename);


   		int n_node, n_element;
   		
	    std::ostringstream s1;
		std::string str1,str2;
	    s1.precision(15);
	    s1.setf(std::ios::scientific, std::ios::floatfield);

   		std::cout << "Creat vtu file..." << std::endl;	
   		std::cout << "\tWriting node data ..." << std::flush;
   
        const FEMSpace<double,2>& fem_space = v0.femSpace();
        const Mesh<2,2>& mesh = fem_space.mesh();
        n_node = mesh.n_point();
 
        n_element = 0;
        n_element = mesh.n_geometry(2);

		element_Piece->SetAttribute( "NumberOfPoints", n_node );	
		element_Piece->SetAttribute( "NumberOfCells", n_element );	


	    for (int i = 0; i < n_node; i++)
	    {
		      s1.str("");
		      s1  << std::endl << mesh.point(i)[0] << "  " << mesh.point(i)[1] << "  " << "0.0";
		      str1 = s1.str();
 			  node_point_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));
	    }
 
        std::cout << " OK!" << std::endl;
   		std::cout << "\tWriting element data ..." << std::flush;

        int n_data = n_node;
        std::vector<int> count(n_data, 0);
        std::vector<double> val0(n_data, 0);
        std::vector<double> val1(n_data, 0);

        int skip = 0;		

        FEMSpace<double,2>::ConstElementIterator 
        the_element = fem_space.beginElement(),
        end_element = fem_space.endElement();
        the_element = fem_space.beginElement();
	    for (int i = 0 ;the_element != end_element; ++ the_element)
   	    {
           const GeometryBM& geo = the_element->geometry();
           switch (geo.n_vertex())
           {
              case 3: 
 		         s1.str("");
		         s1  << std::endl
		             << mesh.geometry(0, geo.vertex(0)).vertex(0) << "  "
		             << mesh.geometry(0, geo.vertex(1)).vertex(0) << "  " 
		             << mesh.geometry(0, geo.vertex(2)).vertex(0) << std::endl;
		         str1 = s1.str();
			     node_Cells_DataArray_1->InsertEndChild( doc1->NewText( str1.c_str() ));
				 skip = skip + 3;
 		         s1.str("");
		         s1  << std::endl << skip;
		         str1 = s1.str();
			     node_Cells_DataArray_2->InsertEndChild( doc1->NewText( str1.c_str() ));
                 node_Cells_DataArray_3->InsertEndChild( doc1->NewText( "\n5\n" ));
                 i = i + 1;
                 break;

              case 4:
 		         s1.str("");
		         s1  << std::endl
		             << mesh.geometry(0, geo.vertex(0)).vertex(0) << "  "
		             << mesh.geometry(0, geo.vertex(1)).vertex(0) << "  " 
		             << mesh.geometry(0, geo.vertex(2)).vertex(0) << "  " 
		             << mesh.geometry(0, geo.vertex(3)).vertex(0) << std::endl;
		         str1 = s1.str();
			     node_Cells_DataArray_1->InsertEndChild( doc1->NewText( str1.c_str() ));
 		         skip = skip + 4;
 		         s1.str("");
		         s1  << std::endl << skip;
		         str1 = s1.str();
			     node_Cells_DataArray_2->InsertEndChild( doc1->NewText( str1.c_str() ));
                 node_Cells_DataArray_3->InsertEndChild( doc1->NewText( "\n9\n" ));
                 i = i + 1;
                 break;

              default:
                 Assert(false, ExcInternalError());
           }
           for (int i = 0;i < geo.n_vertex();i ++)
           {
	          int k = mesh.geometry(0, geo.vertex(i)).vertex(0);
	          count[k] += 1;
	          val0[k] += v0.value(mesh.point(k), *the_element);
	          val1[k] += v1.value(mesh.point(k), *the_element);
           }
        }
        std::cout << " OK!" << std::endl;
        
        for (int i = 0;i < n_node;i ++)
        {
           Assert(count[i] > 0, ExcInternalError());
           val0[i] /= count[i];
           val1[i] /= count[i];
        }

   		for (int i = 0; i < n_node; i++)
   		{
		   {
 		      s1.str("");
		      s1  << std::endl << val0[i] << "  " << val1[i] << "  0.0";
		      str1 = s1.str();
		      node_PointData_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));		   
		   }
   		}

		//doc1->Print();
		doc1->SaveFile( (vtufilename + ".vtu").c_str() );
		std::cout << "\tcreat vtu file: " << (vtufilename + ".vtu").c_str() << " OK!" << std::endl;
        
	    delete doc1;
	
	    return 0;
}

// 直接将解写成vtu 文件.　直接将每个三角形单元的三个点输出为节点，将三个点上的函数值作为ｚ轴坐标，从而将解输出为三维曲面

int WritevtuData3D(FEMFunction<double, 2>& u_h, const char * ofilename, u_int flag = 0)
{
		///构造 vtu 文件的 xml 模板
		XMLDocument* doc1 = new XMLDocument();
				
		XMLDeclaration * declaration = doc1->NewDeclaration();
		doc1->InsertEndChild(declaration);

/**     添加xml文件头的第二种方法
		const char* str = "<?xml version=\"1.0\" ?>";
		doc1->Parse( str );
*/
	
		XMLElement * element_VTKFile,
		           * element_Piece,
		           * element_Points,
		           * element_point_DataArray,
		           * element_Cells,
		           * element_Cells_DataArray_1,
		           * element_Cells_DataArray_2,
		           * element_Cells_DataArray_3,
		           * element_PointData,
		           * element_PointData_DataArray;
		           
		XMLNode * node_VTKFile, 
		        * node_UnstructuredGrid,
		        * node_piece,
		        * node_point,
		        * node_point_DataArray,
		        * node_Cells,
		        * node_Cells_DataArray_1,
		        * node_Cells_DataArray_2,
		        * node_Cells_DataArray_3,
		        * node_PointData,
		        * node_PointData_DataArray;

		
		element_VTKFile = doc1->NewElement("VTKFile");
		element_VTKFile->SetAttribute( "type", "UnstructuredGrid" );
		node_VTKFile = doc1->InsertEndChild( element_VTKFile );
		
		node_UnstructuredGrid = node_VTKFile->InsertEndChild( doc1->NewElement( "UnstructuredGrid" ) );
		element_Piece = doc1->NewElement( "Piece" );
		node_piece = node_UnstructuredGrid->InsertEndChild( element_Piece );

		element_Points = doc1->NewElement( "Points" );
		node_point = node_piece->InsertEndChild( element_Points );

		element_point_DataArray = doc1->NewElement( "DataArray" );
		element_point_DataArray->SetAttribute( "type", "Float32" );	
		element_point_DataArray->SetAttribute( "NumberOfComponents", 3 );	
		element_point_DataArray->SetAttribute( "format", "ascii" );	
		node_point_DataArray = node_point->InsertFirstChild( element_point_DataArray );
		
		element_Cells = doc1->NewElement( "Cells" );
		node_Cells = node_piece->InsertEndChild( element_Cells );

		element_Cells_DataArray_1 = doc1->NewElement( "DataArray" );
		element_Cells_DataArray_1->SetAttribute( "type", "Int32" );	
		element_Cells_DataArray_1->SetAttribute( "Name", "connectivity" );	
		element_Cells_DataArray_1->SetAttribute( "format", "ascii" );	
		node_Cells_DataArray_1 = node_Cells->InsertEndChild( element_Cells_DataArray_1 );

		element_Cells_DataArray_2 = doc1->NewElement( "DataArray" );
		element_Cells_DataArray_2->SetAttribute( "type", "Int32" );	
		element_Cells_DataArray_2->SetAttribute( "Name", "offsets" );	
		element_Cells_DataArray_2->SetAttribute( "format", "ascii" );	
		node_Cells_DataArray_2 = node_Cells->InsertEndChild( element_Cells_DataArray_2 );

		element_Cells_DataArray_3 = doc1->NewElement( "DataArray" );
		element_Cells_DataArray_3->SetAttribute( "type", "UInt8" );	
		element_Cells_DataArray_3->SetAttribute( "Name", "types" );	
		element_Cells_DataArray_3->SetAttribute( "format", "ascii" );	
		node_Cells_DataArray_3 = node_Cells->InsertEndChild( element_Cells_DataArray_3 );
/*
		element_PointData = doc1->NewElement( "PointData" );
		element_PointData->SetAttribute( "Scalars", "u_h" );
		node_PointData = node_piece->InsertEndChild( element_PointData );

		element_PointData_DataArray = doc1->NewElement( "DataArray" );
		element_PointData_DataArray->SetAttribute( "type", "Float32" );	
		element_PointData_DataArray->SetAttribute( "Name", "u_h" );	
		element_PointData_DataArray->SetAttribute( "format", "ascii" );
		node_PointData->InsertEndChild( element_PointData_DataArray );
*/
		element_PointData = doc1->NewElement( "PointData" );
		if ( flag == 0)
		{
		   element_PointData->SetAttribute( "Scalars", "u_h" );
		   node_PointData = node_piece->InsertEndChild( element_PointData );
		   
		   element_PointData_DataArray = doc1->NewElement( "DataArray" );
		   element_PointData_DataArray->SetAttribute( "type", "Float32" );	
		   element_PointData_DataArray->SetAttribute( "Name", "u_h" );	
		   element_PointData_DataArray->SetAttribute( "format", "ascii" );
		   node_PointData_DataArray = node_PointData->InsertEndChild( element_PointData_DataArray );
		}
		else if ( flag == 1)
		{
		   element_PointData->SetAttribute( "Vectors", "u_h" );
		   node_PointData = node_piece->InsertEndChild( element_PointData );

		   element_PointData_DataArray = doc1->NewElement( "DataArray" );
		   element_PointData_DataArray->SetAttribute( "type", "Float32" );
		   element_PointData_DataArray->SetAttribute( "NumberOfComponents", 3 );
		   element_PointData_DataArray->SetAttribute( "Name", "u_h" );	
		   element_PointData_DataArray->SetAttribute( "format", "ascii" );
		   node_PointData_DataArray = node_PointData->InsertEndChild( element_PointData_DataArray );
		}

   /// vtu 文件模板构造完成，开始写入相关数据     


 		const std::string& vtufilename(ofilename);


   		int n_node, n_element;

        int n_data = n_node;
        double val_1, val_2, val_3, max_val = 0;
        double pointx,pointy;
        int kk = 0;
   		
	    std::ostringstream s1;
		std::string str1,str2;
	    s1.precision(15);
	    s1.setf(std::ios::scientific, std::ios::floatfield);

   		std::cout << "Creat vtu file..." << std::endl;	
   
        const FEMSpace<double,2>& fem_space = u_h.femSpace();
        const Mesh<2,2>& mesh = fem_space.mesh();
        
        n_node = 0;
        n_element = 0;
        FEMSpace<double,2>::ConstElementIterator 
        the_element = fem_space.beginElement(),
        end_element = fem_space.endElement();
        for (;the_element != end_element;++ the_element)
        {
          const GeometryBM& geo = the_element->geometry();
          switch (geo.n_vertex()) 
          {
            case 3: 
            	n_element += 1;
            	n_node += 3;
          		for (int k = 0; k < 3; k++)
            	{
               		kk = mesh.geometry(0, geo.vertex(k)).vertex(0);
               		val_1 = u_h.value(mesh.point(kk),*the_element);
               		if (max_val < fabs(val_1))
                  		max_val = fabs(val_1);
            	}
            	break;
            case 4: 
            	n_element += 1; 
            	n_node += 4; 
          		for (int k = 0; k < 4; k++)
            	{
               		kk = mesh.geometry(0, geo.vertex(k)).vertex(0);
               		val_1 = u_h.value(mesh.point(kk),*the_element);
               		if (max_val < fabs(val_1))
                  		max_val = fabs(val_1);
            	}
            	break;
          }
        }    		
		max_val = 1.0;
		element_Piece->SetAttribute( "NumberOfPoints", n_node );	
		element_Piece->SetAttribute( "NumberOfCells", n_element );	

   		std::cout << "\tWriting element data ..." << std::flush;


        int skip = 0;		
        the_element = fem_space.beginElement();
	    for (;the_element != end_element; ++ the_element)
   	    {
           const GeometryBM& geo = the_element->geometry();
           switch (geo.n_vertex())
           {
              case 3: 
                 for (int k = 0; k <3; k++)
                 {
                    kk = mesh.geometry(0, geo.vertex(k)).vertex(0);
                    pointx = mesh.point(kk)[0];
                    pointy = mesh.point(kk)[1];
                    
                    if ( flag ==0 )
                       val_1 = u_h.value(mesh.point(kk),*the_element);
                    else if ( flag == 1 )
                    {
                       val_1 = u_h.gradient(mesh.point(kk),*the_element)[0];
                       val_2 = u_h.gradient(mesh.point(kk),*the_element)[1];  
                    }                                               
 		            s1.str("");
		            s1  << std::endl << pointx << "  " << pointy << "  " << val_1/max_val;
		            str1 = s1.str();
 			        node_point_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));

		            if ( flag == 0 )
		            {
 		               s1.str("");
		               s1  << std::endl << val_1;
		               str1 = s1.str();
		               node_PointData_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));
		            }
		            else if ( flag == 1 )
		            {
 		               s1.str("");
		               s1  << std::endl << val_1 << "  " << val_2;
		               str1 = s1.str();
		               node_PointData_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));		   
		            }
                 }
		         s1.str("");
		         s1  << std::endl
		             << skip + 0 << "  "
		             << skip + 1 << "  " 
		             << skip + 2 << std::endl;
		         str1 = s1.str();
			     node_Cells_DataArray_1->InsertEndChild( doc1->NewText( str1.c_str() ));

				 skip = skip + 3;
 		         s1.str("");
		         s1  << std::endl << skip;
		         str1 = s1.str();
			     node_Cells_DataArray_2->InsertEndChild( doc1->NewText( str1.c_str() ));
                 node_Cells_DataArray_3->InsertEndChild( doc1->NewText( "\n5\n" ));
                 break;

              case 4:
                 for (int k = 0; k < 4; k++)
                 {
                    kk = mesh.geometry(0, geo.vertex(k)).vertex(0);
                    pointx = mesh.point(kk)[0];
                    pointy = mesh.point(kk)[1];
                    
                    if ( flag ==0 )
                       val_1 = u_h.value(mesh.point(kk),*the_element);
                    else if ( flag == 1 )
                    {
                       val_1 = u_h.gradient(mesh.point(kk),*the_element)[0];
                       val_2 = u_h.gradient(mesh.point(kk),*the_element)[1];  
                    }                                               
 		            s1.str("");
		            s1  << std::endl << pointx << "  " << pointy << "  " << val_1/max_val;
		            str1 = s1.str();
 			        node_point_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));

		            if ( flag == 0 )
		            {
 		               s1.str("");
		               s1  << std::endl << val_1;
		               str1 = s1.str();
		               node_PointData_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));
		            }
		            else if ( flag == 1 )
		            {
 		               s1.str("");
		               s1  << std::endl << val_1 << "  " << val_2;
		               str1 = s1.str();
		               node_PointData_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));		   
		            }
                 }
		         s1.str("");
		         s1  << std::endl
		             << skip + 0 << "  "
		             << skip + 1 << "  " 
		             << skip + 2 << "  "
		             << skip + 3 << std::endl;
		         str1 = s1.str();
			     node_Cells_DataArray_1->InsertEndChild( doc1->NewText( str1.c_str() ));

				 skip = skip + 4;
 		         s1.str("");
		         s1  << std::endl << skip;
		         str1 = s1.str();
			     node_Cells_DataArray_2->InsertEndChild( doc1->NewText( str1.c_str() ));
                 node_Cells_DataArray_3->InsertEndChild( doc1->NewText( "\n9\n" ));
                 break;
 

              default:
                 Assert(false, ExcInternalError());
           }
        }
        std::cout << " OK!" << std::endl;

		//doc1->Print();
		doc1->SaveFile( (vtufilename + "_3D.vtu").c_str() );
		std::cout << "\tcreat vtu file: " << (vtufilename + "_3D.vtu").c_str() << " OK!" << std::endl;
        
	    delete doc1;
	
	    return 0;
}
// 将并行计算得到个各个分区的 vtu 文件合成为一个 pvtu 文件，方便显示。
// 各个分区的文件名前段相同，后面加0,1,2,3...等一个数字表示不同分区。
// 例如：data_0.vtu,data_1.vtu。
// piece_filename 表示文件名前段，n_rank 表示分区的个数。
// filenamepvtu 表示生成的 pvtu 文件，例如 data_all.pvtu。

int creatpvtu( std::string piece_filename, 
               u_int n_rank,
			   std::string filenamepvtu)
{
   std::ostringstream s1;
   std::string str1,str2;


	XMLDocument* doc = new XMLDocument();
			
	XMLDeclaration * declaration = doc->NewDeclaration();
	doc->InsertEndChild(declaration);

		XMLElement * element_VTKFile,
				   * element_PUnstructuredGrid,
		           * element_PPiece,
		           * element_PPoints,
		           * element_PPoint_DataArray,
		           * element_PPointData,
		           * element_PPointData_DataArray;
		           
		XMLNode * node_VTKFile, 
		        * node_PUnstructuredGrid,
		        * node_PPoints,
		        * node_PPointData;

		element_VTKFile = doc->NewElement("VTKFile");
		element_VTKFile->SetAttribute( "type", "PUnstructuredGrid" );
		node_VTKFile = doc->InsertEndChild( element_VTKFile );
		
		  element_PUnstructuredGrid = doc->NewElement( "PUnstructuredGrid" );
		  element_PUnstructuredGrid->SetAttribute( "GhostLevel", "0" );
		  node_PUnstructuredGrid = node_VTKFile->InsertEndChild( element_PUnstructuredGrid );

		     element_PPointData = doc->NewElement( "PPointData" );
		     element_PPointData->SetAttribute( "Scalars", "u_h" );
		     node_PPointData = node_PUnstructuredGrid->InsertEndChild( element_PPointData );

		        element_PPointData_DataArray = doc->NewElement( "PDataArray" );
		        element_PPointData_DataArray->SetAttribute( "type", "Float32" );	
		        element_PPointData_DataArray->SetAttribute( "Name", "u_h" );	
		        element_PPointData_DataArray->SetAttribute( "format", "ascii" );
		        node_PPointData->InsertEndChild( element_PPointData_DataArray );

 		     element_PPoints = doc->NewElement( "PPoints" );
		     node_PPoints = node_PUnstructuredGrid->InsertEndChild( element_PPoints );

		        element_PPoint_DataArray = doc->NewElement( "PDataArray" );
		        element_PPoint_DataArray->SetAttribute( "type", "Float32" );	
		        element_PPoint_DataArray->SetAttribute( "NumberOfComponents", 3 );	
		        element_PPoint_DataArray->SetAttribute( "format", "ascii" );	
		        node_PPoints->InsertFirstChild( element_PPoint_DataArray );


	        for( u_int i = 0; i < n_rank; i++ )
	        {
               s1.str("");
               s1 << piece_filename.c_str() << i << ".vtu";
               str1 = s1.str();

		       element_PPiece = doc->NewElement( "Piece" );
		       element_PPiece->SetAttribute( "Source", str1.c_str() );	
		       node_PUnstructuredGrid->InsertEndChild( element_PPiece );
	        }
	doc->SaveFile( (filenamepvtu).c_str() );
	delete doc;

	return 0;
}

int creatpvd( std::vector< double > time_t, 
			  std::vector< std::string > filename,
			  std::string filenamepvd)
{
	if( time_t.capacity() != filename.capacity() )
	{
		std::cout << "number of time step is not equal to number of data files"  << std::endl;
		return 1;
	}
	XMLDocument* doc = new XMLDocument();
			
	XMLDeclaration * declaration = doc->NewDeclaration();
	doc->InsertEndChild(declaration);

		XMLElement * element_VTKFile,
		           * element_Collection,
		           * element_DataSet;
		           
		XMLNode * node_VTKFile, 
		        * node_Collection;

		element_VTKFile = doc->NewElement("VTKFile");
		element_VTKFile->SetAttribute( "type", "Collection" );
		node_VTKFile = doc->InsertEndChild( element_VTKFile );

		element_Collection = doc->NewElement( "Collection" );
		node_Collection = node_VTKFile->InsertEndChild( element_Collection );

	for( u_int i = 0; i < time_t.capacity(); i++ )
	{
		element_DataSet = doc->NewElement( "DataSet" );
		element_DataSet->SetAttribute( "timestep", time_t.at(i) );	
		element_DataSet->SetAttribute( "group", "" );	
		element_DataSet->SetAttribute( "part", 0 );	
		element_DataSet->SetAttribute( "file", (filename.at(i)).c_str() );	
		node_Collection->InsertEndChild( element_DataSet );
	}
//	doc->Print();
	doc->SaveFile( (filenamepvd + ".pvd").c_str() );
	std::cout << "\t creat pvd file: " << (filenamepvd  + ".pvd").c_str() << " OK!" << std::endl;
	delete doc;

	return 0;
}



int EasymeshData2vtu(const char * file1)
{

		///构造 vtu 文件的 xml 模板
		XMLDocument* doc1 = new XMLDocument();
				
		XMLDeclaration * declaration = doc1->NewDeclaration();
		doc1->InsertEndChild(declaration);

/**     添加xml文件头的第二种方法
		const char* str = "<?xml version=\"1.0\" ?>";
		doc1->Parse( str );
*/
	
		XMLElement * element_VTKFile,
		           * element_Piece,
		           * element_Points,
		           * element_point_DataArray,
		           * element_Cells,
		           * element_Cells_DataArray_1,
		           * element_Cells_DataArray_2,
		           * element_Cells_DataArray_3,
		           * element_PointData,
		           * element_PointData_DataArray;
		           
		XMLNode * node_VTKFile, 
		        * node_UnstructuredGrid,
		        * node_piece,
		        * node_point,
		        * node_point_DataArray,
		        * node_Cells,
		        * node_Cells_DataArray_1,
		        * node_Cells_DataArray_2,
		        * node_Cells_DataArray_3,
		        * node_PointData,
		        * node_PointData_DataArray;

		
		element_VTKFile = doc1->NewElement("VTKFile");
		element_VTKFile->SetAttribute( "type", "UnstructuredGrid" );
		node_VTKFile = doc1->InsertEndChild( element_VTKFile );
		
		node_UnstructuredGrid = node_VTKFile->InsertEndChild( doc1->NewElement( "UnstructuredGrid" ) );
		element_Piece = doc1->NewElement( "Piece" );
//		element_Piece->SetAttribute( "NumberOfPoints", a );	
//		element_Piece->SetAttribute( "NumberOfCells", b );	
		node_piece = node_UnstructuredGrid->InsertEndChild( element_Piece );

		element_Points = doc1->NewElement( "Points" );
		node_point = node_piece->InsertEndChild( element_Points );

		element_point_DataArray = doc1->NewElement( "DataArray" );
		element_point_DataArray->SetAttribute( "type", "Float32" );	
		element_point_DataArray->SetAttribute( "NumberOfComponents", 3 );	
		element_point_DataArray->SetAttribute( "format", "ascii" );	
		node_point_DataArray = node_point->InsertFirstChild( element_point_DataArray );
		
		element_Cells = doc1->NewElement( "Cells" );
		node_Cells = node_piece->InsertEndChild( element_Cells );

		element_Cells_DataArray_1 = doc1->NewElement( "DataArray" );
		element_Cells_DataArray_1->SetAttribute( "type", "Int32" );	
		element_Cells_DataArray_1->SetAttribute( "Name", "connectivity" );	
		element_Cells_DataArray_1->SetAttribute( "format", "ascii" );	
		node_Cells_DataArray_1 = node_Cells->InsertEndChild( element_Cells_DataArray_1 );

		element_Cells_DataArray_2 = doc1->NewElement( "DataArray" );
		element_Cells_DataArray_2->SetAttribute( "type", "Int32" );	
		element_Cells_DataArray_2->SetAttribute( "Name", "offsets" );	
		element_Cells_DataArray_2->SetAttribute( "format", "ascii" );	
		node_Cells_DataArray_2 = node_Cells->InsertEndChild( element_Cells_DataArray_2 );

		element_Cells_DataArray_3 = doc1->NewElement( "DataArray" );
		element_Cells_DataArray_3->SetAttribute( "type", "UInt8" );	
		element_Cells_DataArray_3->SetAttribute( "Name", "types" );	
		element_Cells_DataArray_3->SetAttribute( "format", "ascii" );	
		node_Cells_DataArray_3 = node_Cells->InsertEndChild( element_Cells_DataArray_3 );

		element_PointData = doc1->NewElement( "PointData" );
		element_PointData->SetAttribute( "Scalars", "u_h" );
		node_PointData = node_piece->InsertEndChild( element_PointData );

		element_PointData_DataArray = doc1->NewElement( "DataArray" );
		element_PointData_DataArray->SetAttribute( "type", "Float32" );	
		element_PointData_DataArray->SetAttribute( "Name", "u_h" );	
		element_PointData_DataArray->SetAttribute( "format", "ascii" );
		node_PointData_DataArray = node_PointData->InsertEndChild( element_PointData_DataArray );

   /// vtu 文件模板构造完成，开始写入相关数据，首先从计算结果的 easymesh 文件中读取网格信息     


 		const std::string& easymeshfilename(file1);

   		int i, j, p1, p2, p3, mark;
   		u_int n_node, n_side, n_element;
   		double x, y, u;
   		char text[64];

	    std::ostringstream s1;
		   std::string str1,str2;
	   s1.precision(15);
	   s1.setf(std::ios::scientific, std::ios::floatfield);
		//    s1.str("");
		//    s1 << "data/u_" << 10.2 << "_0_S";
		//    str1 = s1.str();

   		std::cout << "Reading easymesh data file to creat vtu file..." << std::endl;	
   		std::cout << "\treading node data ..." << std::flush;
   
   		std::ifstream is((easymeshfilename + ".n").c_str());
   		if( !is ) std::cout << "Error" << std::endl;
   		
   		is >> n_node >> n_element >> n_side;
   		is.getline(text, 64);

		element_Piece->SetAttribute( "NumberOfPoints", n_node );	
		element_Piece->SetAttribute( "NumberOfCells", n_element );	


	    for (u_int i = 0; i < n_node; i++)
	    {
	       is >> j
	          >> x 
	          >> y
    	      >> mark;
 		      s1.str("");
		      s1  << std::endl << x << "  " << y << "  " << "0.0";
		      str1 = s1.str();
 			  node_point_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));
	    }
   		is.close();
 
        std::cout << " OK!" << std::endl;
   		std::cout << "\treading element data ..." << std::flush;

   		is.open((easymeshfilename + ".e").c_str());
   		is >> i;
	    is >> i;
	    is >> i;
 	    is.getline(text, 64);

	    for (u_int i = 0; i < n_element; i++)
   	    {
      	   is >> j
         	  >> p1 
              >> p2
              >> p3
              >> j >> j >> j
              >> j >> j >> j;

 		      s1.str("");
		      s1  << std::endl << p1 << "  " << p2 << "  " << p3 << std::endl;
		      str1 = s1.str();
			  node_Cells_DataArray_1->InsertEndChild( doc1->NewText( str1.c_str() ));

 		      s1.str("");
		      s1  << std::endl << (i+1)*3;
		      str1 = s1.str();
			  node_Cells_DataArray_2->InsertEndChild( doc1->NewText( str1.c_str() ));

			  node_Cells_DataArray_3->InsertEndChild( doc1->NewText( "\n5\n" ));
        }
	    is.close();
   

	    is.open((easymeshfilename + ".dat").c_str());
   		for (u_int i = 0; i < n_node; i++)
   		{
      		 is >> u;
      		// os << u << std::endl;
 		     s1.str("");
		     s1  << std::endl << u;
		     str1 = s1.str();
			 node_PointData->InsertEndChild( doc1->NewText( str1.c_str() ));
   		}
   		is.close();
 	    std::cout << " OK!" << std::endl;

		//doc1->Print();
		doc1->SaveFile( (easymeshfilename + ".vtu").c_str() );
		std::cout << "\tcreat vtu file: " << (easymeshfilename + ".vtu").c_str() << " OK!" << std::endl;

	delete doc1;


	return 0;
}

int EasymeshData2vtu3D(const char * file1)
{

		///构造 vtu 文件的 xml 模板
		XMLDocument* doc1 = new XMLDocument();
				
		XMLDeclaration * declaration = doc1->NewDeclaration();
		doc1->InsertEndChild(declaration);

/**     添加xml文件头的第二种方法
		const char* str = "<?xml version=\"1.0\" ?>";
		doc1->Parse( str );
*/
	
		XMLElement * element_VTKFile,
		           * element_Piece,
		           * element_Points,
		           * element_point_DataArray,
		           * element_Cells,
		           * element_Cells_DataArray_1,
		           * element_Cells_DataArray_2,
		           * element_Cells_DataArray_3,
		           * element_PointData,
		           * element_PointData_DataArray;
		           
		XMLNode * node_VTKFile, 
		        * node_UnstructuredGrid,
		        * node_piece,
		        * node_point,
		        * node_point_DataArray,
		        * node_Cells,
		        * node_Cells_DataArray_1,
		        * node_Cells_DataArray_2,
		        * node_Cells_DataArray_3,
		        * node_PointData,
		        * node_PointData_DataArray;

		
		element_VTKFile = doc1->NewElement("VTKFile");
		element_VTKFile->SetAttribute( "type", "UnstructuredGrid" );
		node_VTKFile = doc1->InsertEndChild( element_VTKFile );
		
		node_UnstructuredGrid = node_VTKFile->InsertEndChild( doc1->NewElement( "UnstructuredGrid" ) );
		element_Piece = doc1->NewElement( "Piece" );
//		element_Piece->SetAttribute( "NumberOfPoints", a );	
//		element_Piece->SetAttribute( "NumberOfCells", b );	
		node_piece = node_UnstructuredGrid->InsertEndChild( element_Piece );

		element_Points = doc1->NewElement( "Points" );
		node_point = node_piece->InsertEndChild( element_Points );

		element_point_DataArray = doc1->NewElement( "DataArray" );
		element_point_DataArray->SetAttribute( "type", "Float32" );	
		element_point_DataArray->SetAttribute( "NumberOfComponents", 3 );	
		element_point_DataArray->SetAttribute( "format", "ascii" );	
		node_point_DataArray = node_point->InsertFirstChild( element_point_DataArray );
		
		element_Cells = doc1->NewElement( "Cells" );
		node_Cells = node_piece->InsertEndChild( element_Cells );

		element_Cells_DataArray_1 = doc1->NewElement( "DataArray" );
		element_Cells_DataArray_1->SetAttribute( "type", "Int32" );	
		element_Cells_DataArray_1->SetAttribute( "Name", "connectivity" );	
		element_Cells_DataArray_1->SetAttribute( "format", "ascii" );	
		node_Cells_DataArray_1 = node_Cells->InsertEndChild( element_Cells_DataArray_1 );

		element_Cells_DataArray_2 = doc1->NewElement( "DataArray" );
		element_Cells_DataArray_2->SetAttribute( "type", "Int32" );	
		element_Cells_DataArray_2->SetAttribute( "Name", "offsets" );	
		element_Cells_DataArray_2->SetAttribute( "format", "ascii" );	
		node_Cells_DataArray_2 = node_Cells->InsertEndChild( element_Cells_DataArray_2 );

		element_Cells_DataArray_3 = doc1->NewElement( "DataArray" );
		element_Cells_DataArray_3->SetAttribute( "type", "UInt8" );	
		element_Cells_DataArray_3->SetAttribute( "Name", "types" );	
		element_Cells_DataArray_3->SetAttribute( "format", "ascii" );	
		node_Cells_DataArray_3 = node_Cells->InsertEndChild( element_Cells_DataArray_3 );

		element_PointData = doc1->NewElement( "PointData" );
		element_PointData->SetAttribute( "Scalars", "u_h" );
		node_PointData = node_piece->InsertEndChild( element_PointData );

		element_PointData_DataArray = doc1->NewElement( "DataArray" );
		element_PointData_DataArray->SetAttribute( "type", "Float32" );	
		element_PointData_DataArray->SetAttribute( "Name", "u_h" );	
		element_PointData_DataArray->SetAttribute( "format", "ascii" );
		node_PointData_DataArray = node_PointData->InsertEndChild( element_PointData_DataArray );

   /// vtu 文件模板构造完成，开始写入相关数据，首先从计算结果的 easymesh 文件中读取网格信息     


 		const std::string& easymeshfilename(file1);

   		int i, j, p1, p2, p3, mark;
   		u_int n_node, n_side, n_element;
   		double x, y, u, max_x_y, max_u;
   		max_x_y = 1.0;
   		max_u = 1.0;
   		
   		char text[64];

	    std::ostringstream s1;
		   std::string str1,str2;
	   s1.precision(15);
	   s1.setf(std::ios::scientific, std::ios::floatfield);
		//    s1.str("");
		//    s1 << "data/u_" << 10.2 << "_0_S";
		//    str1 = s1.str();

   		std::cout << "Reading easymesh data file to creat vtu file..." << std::endl;	
   		std::cout << "\treading node data ..." << std::flush;
   
   		std::ifstream is((easymeshfilename + ".n").c_str());
   		std::ifstream is1((easymeshfilename + ".dat").c_str());
   		if( !is || ! is1 ) std::cout << "Error" << std::endl;

   		is >> n_node >> n_element >> n_side;
   		is.getline(text, 64);

	    for (u_int i = 0; i < n_node; i++)
	    {
	       is >> j
	          >> x 
	          >> y
    	      >> mark;
    	   is1 >> u;
    	   
    	   if ( max_x_y <= fabs(x)) max_x_y = x;
    	   if ( max_x_y <= fabs(y)) max_x_y = y;
      	   if ( max_u <= fabs(u)) max_u = u;
  	   
 	    }
   		is.close();
   		is1.close();

        is.open((easymeshfilename + ".n").c_str());
        is1.open((easymeshfilename + ".dat").c_str());
   		if( !is || ! is1 ) std::cout << "Error" << std::endl;
   		
   		is >> n_node >> n_element >> n_side;
   		is.getline(text, 64);

		element_Piece->SetAttribute( "NumberOfPoints", n_node );	
		element_Piece->SetAttribute( "NumberOfCells", n_element );	


	    for (u_int i = 0; i < n_node; i++)
	    {
	       is >> j
	          >> x 
	          >> y
    	      >> mark;
    	   is1 >> u;
 		   s1.str("");
		   s1  << std::endl << x << "  " << y << "  " << max_x_y*u/max_u;
		   str1 = s1.str();
 		   node_point_DataArray->InsertEndChild( doc1->NewText( str1.c_str() ));
	    }
   		is.close();
   		is1.close();
 
        std::cout << " OK!" << std::endl;
   		std::cout << "\treading element data ..." << std::flush;

   		is.open((easymeshfilename + ".e").c_str());
   		is >> i;
	    is >> i;
	    is >> i;
 	    is.getline(text, 64);

	    for (u_int i = 0; i < n_element; i++)
   	    {
      	   is >> j
         	  >> p1 
              >> p2
              >> p3
              >> j >> j >> j
              >> j >> j >> j;

 		      s1.str("");
		      s1  << std::endl << p1 << "  " << p2 << "  " << p3 << std::endl;
		      str1 = s1.str();
			  node_Cells_DataArray_1->InsertEndChild( doc1->NewText( str1.c_str() ));

 		      s1.str("");
		      s1  << std::endl << (i+1)*3;
		      str1 = s1.str();
			  node_Cells_DataArray_2->InsertEndChild( doc1->NewText( str1.c_str() ));

			  node_Cells_DataArray_3->InsertEndChild( doc1->NewText( "\n5\n" ));
        }
	    is.close();
   

	    is.open((easymeshfilename + ".dat").c_str());
   		for (u_int i = 0; i < n_node; i++)
   		{
      		 is >> u;
      		// os << u << std::endl;
 		     s1.str("");
		     s1  << std::endl << u;
		     str1 = s1.str();
			 node_PointData->InsertEndChild( doc1->NewText( str1.c_str() ));
   		}
   		is.close();
 	    std::cout << " OK!" << std::endl;

		//doc1->Print();
		doc1->SaveFile( (easymeshfilename + "_3D.vtu").c_str() );
		std::cout << "\tcreat vtu file: " << (easymeshfilename + "_3D.vtu").c_str() << " OK!" << std::endl;

	delete doc1;


	return 0;
}

int EasymeshData2vtk(const char * file1)
{
   const std::string& easymeshfilename(file1);
//   const std::string& vtkfilename(file2);
   
 //  easymeshfilename = argv[1];
  // vtkfilename = argv[2];
   


   int i, j, p1, p2, p3, mark;
   int n_node, n_side, n_element;
   double x, y, u;
   char text[64];

   std::ofstream os((easymeshfilename + ".vtk").c_str());
   os.precision(15);
   os.setf(std::ios::scientific, std::ios::floatfield);
   os << "# vtk DataFile Version 3.0" << std::endl;
   os << "2D scalar data" << std::endl;
   os << "ASCII" << std::endl;
   os << std::endl;
   os << "DATASET UNSTRUCTURED_GRID" << std::endl;

   std::cout << "Reading easymesh data file to creat vtk file..." << std::endl;	
   std::cout << "\treading node data ..." << std::flush;
   
   std::ifstream is((easymeshfilename + ".n").c_str());
   is >> n_node >> n_element >> n_side;
   is.getline(text, 64);

   os << "POINTS  " << n_node << "  float" << std::endl;
   
   for (i = 0;i < n_node;i ++)
   {
      is >> j
         >> x 
         >> y
         >> mark;
      os << x << "  " << y << "  " << 0.0 << std::endl;
   }
   is.close();
 
   std::cout << " OK!" << std::endl;
   std::cout << "\treading element data ..." << std::flush;

   is.open((easymeshfilename + ".e").c_str());
   is >> i;
//   Assert(i == n_element, ExcMeshData("in element file: element number error"));
   is >> i;
//   Assert(i == n_node, ExcMeshData("in element file: node number error"));
   is >> i;
//   Assert(i == n_side, ExcMeshData("in element file: side number error"));
   is.getline(text, 64);

   int temp = 0;
   temp = 4 * n_element;
   os << std::endl;
   os << "CELLS  " << n_element << "  " << temp << std::endl;
   
   for (i = 0;i < n_element;i ++)
   {
      is >> j
         >> p1 
         >> p2
         >> p3
         >> j >> j >> j
         >> j >> j >> j;
      os << 3 << "  " << p1 << "  " << p2 << "  " << p3 << std::endl;
   }
   os << std::endl;
   os << "CELL_TYPES " << n_element << std::endl;
   for (i = 0;i < n_element;i ++)
   {
      os << 5 << std::endl;
   }
   is.close();
   
   os << std::endl;
   os << "POINT_DATA  " << n_node << std::endl;
   os << "SCALARS u_h float 1" << std::endl;
   os << "LOOKUP_TABLE default" << std::endl;

   is.open((easymeshfilename + ".dat").c_str());
   for (i = 0;i < n_node;i ++)
   {
      is >> u;
      os << u << std::endl;
   }
   is.close();
   os << std::endl;
   
   os.close(); 

   std::cout << " OK!" << std::endl;
   std::cout << "\tcreat vtk file: " << (easymeshfilename + ".vtk").c_str() << " OK!" << std::endl;
   return 0;
}

int EasymeshData2vtk3D(const char * file1)
{
   const std::string& easymeshfilename(file1);
//   const std::string& vtkfilename(file2);
   
 //  easymeshfilename = argv[1];
  // vtkfilename = argv[2];
   


   int i, j, p1, p2, p3, mark;
   int n_node, n_side, n_element;
   double x, y, u;
   char text[64];

   std::ofstream os((easymeshfilename + "_3D.vtk").c_str());
   os.precision(15);
   os.setf(std::ios::scientific, std::ios::floatfield);
   os << "# vtk DataFile Version 3.0" << std::endl;
   os << "2D scalar data" << std::endl;
   os << "ASCII" << std::endl;
   os << std::endl;
   os << "DATASET UNSTRUCTURED_GRID" << std::endl;

   std::cout << "Reading easymesh data file to creat vtk file..." << std::endl;	
   std::cout << "\treading node data ..." << std::flush;
   
   std::ifstream is((easymeshfilename + ".n").c_str());
   is >> n_node >> n_element >> n_side;
   is.getline(text, 64);

   std::ifstream is1((easymeshfilename + ".dat").c_str());

   os << "POINTS  " << n_node << "  float" << std::endl;
   
   for (i = 0;i < n_node;i ++)
   {
      is >> j
         >> x 
         >> y
         >> mark;
      is1 >> u;
      os << x << "  " << y << "  " << u << std::endl;
   }
   is.close();
   is1.close();
 
   std::cout << " OK!" << std::endl;
   std::cout << "\treading element data ..." << std::flush;

   is.open((easymeshfilename + ".e").c_str());
   is >> i;
//   Assert(i == n_element, ExcMeshData("in element file: element number error"));
   is >> i;
//   Assert(i == n_node, ExcMeshData("in element file: node number error"));
   is >> i;
//   Assert(i == n_side, ExcMeshData("in element file: side number error"));
   is.getline(text, 64);

   int temp = 0;
   temp = 4 * n_element;
   os << std::endl;
   os << "CELLS  " << n_element << "  " << temp << std::endl;
   
   for (i = 0;i < n_element;i ++)
   {
      is >> j
         >> p1 
         >> p2
         >> p3
         >> j >> j >> j
         >> j >> j >> j;
      os << 3 << "  " << p1 << "  " << p2 << "  " << p3 << std::endl;
   }
   os << std::endl;
   os << "CELL_TYPES " << n_element << std::endl;
   for (i = 0;i < n_element;i ++)
   {
      os << 5 << std::endl;
   }
   is.close();
   
   os << std::endl;
   os << "POINT_DATA  " << n_node << std::endl;
   os << "SCALARS u_h float 1" << std::endl;
   os << "LOOKUP_TABLE default" << std::endl;

   is.open((easymeshfilename + ".dat").c_str());
   for (i = 0;i < n_node;i ++)
   {
      is >> u;
      os << u << std::endl;
   }
   is.close();
   os << std::endl;
   
   os.close(); 

   std::cout << " OK!" << std::endl;
   std::cout << "\tcreat vtk file: " << (easymeshfilename + "_3D.vtk").c_str() << " OK!" << std::endl;
   return 0;
}

#endif // __VtkDataTool_h__



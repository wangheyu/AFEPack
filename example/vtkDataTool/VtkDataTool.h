/**
 * @file   VtkDataTool.h
 * @author Feng Han <fengerhan.@gmail.com>
 * @date   Fri Mar 11 13:36:59 2017
 * 
 * @brief  本程序主要为了解决 AFEPack 中二维问题物理量的显示问题。AFEPack 可以将二维问题中
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
 *              pvd 文件的格式，请自行在网上搜索参考。当有了一系列 vtu 文件以及这些文件
 *              对应的时间后，调用
 *              creatpvd( time_t, file_vtu, pvdfilename);
 *              
 *              就可以生成 pvd 文件了，其中 time_t 是时间序列数组，file_vtu 是 vtu 文件
 *              名序列数组，pvdfilename 是生成的 pvd 文件名。
 *          
 *              由于要生成 xml 格式的vtu文件，我在这里使用了 tinyxml2 软件包。
 *
 *              添加了将AFEPack 二维问题中的物理量直接输出为vtu文件的函数，这样更方便了．
 *         　　　例如
 *　　　　　　　　　　　　　　　　　　　　WritevtuData(FEMFunction<double, 2>& u_h, const char * ofilename, u_int flag = 0)
 *              将输出标量函数　u_h　的函数值，当 flag = 1 时，将输出　u_h　的梯度值，目前
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

using namespace tinyxml2;

/// 按片输出结果，每个网格作为一个cell，对应一个值，用于输出误差指示子的结果
int WritevtpData(RegularMesh<2>& mesh, std::vector<double> celldata, const char * ofilename);

/// 按片输出结果，每个网格作为一个cell，对应一个值，根据该值的大小，输出到三维空间中，用于输出误差指示子的结果
int WritevtpData3D(RegularMesh<2>& mesh, std::vector<double> celldata, const char * ofilename);

int WritevtuData(FEMFunction<double, 3>& u_h, const char * ofilename, u_int flag = 0);

int WritevtuData(FEMFunction<double, 2>& u_h, const char * ofilename, u_int flag = 0);

///　按矢量分量输出矢量
int WritevtuData(FEMFunction<double, 2>& v0, FEMFunction<double, 2>& v1, const char * ofilename);


/// 直接将解写成vtu 文件.　直接将每个三角形单元的三个点输出为节点，将三个点上的函数值作为ｚ轴坐标，从而将解输出为三维曲面
int WritevtuData3D(FEMFunction<double, 2>& u_h, const char * ofilename, u_int flag = 0);

/* 将并行计算得到个各个分区的 vtu 文件合成为一个 pvtu 文件，方便显示。
*  各个分区的文件名前段相同，后面加0,1,2,3...等一个数字表示不同分区。
*  例如：data_0.vtu,data_1.vtu。
*  piece_filename 表示文件名前段，n_rank 表示分区的个数。
*  filenamepvtu 表示生成的 pvtu 文件，例如 data_all.pvtu。
*/
int creatpvtu( std::string piece_filename, 
               u_int n_rank,
			   std::string filenamepvtu);

int creatpvd( std::vector< double > time_t, 
			  std::vector< std::string > filename,
			  std::string filenamepvd);



int EasymeshData2vtu(const char * file1);

int EasymeshData2vtu3D(const char * file1);

int EasymeshData2vtk(const char * file1);

int EasymeshData2vtk3D(const char * file1);

#endif // __VtkDataTool_h__



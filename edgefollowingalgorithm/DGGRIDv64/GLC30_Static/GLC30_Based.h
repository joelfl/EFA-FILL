#pragma once
#ifndef _GLC30_BASED_H_
#define _GLC30_BASED_H_


#include "gdal.h" 
#include "gdal_priv.h" 
#include <ogr_spatialref.h>
#include <gdalwarper.h>
#include <ogr_api.h>
#include <ogrsf_frmts.h>
#include <dggrid.h>
#include <gridgen.h>
#include <clipper_region.h>
#include <cmath>
#include <climits>

#include "clipper.hpp"
#include "DgIVec2D.h"
#include "DgInputStream.h"
#include "DgInAIGenFile.h"
#include "DgInShapefile.h"
#include "DgOutShapefile.h"
#include "DgInShapefileAtt.h"
#include "DgOutLocFile.h"
#include "DgIDGG.h"
#include "DgBoundedIDGG.h"
#include "DgGeoSphRF.h"
#include "DgBoundedRF2D.h"
#include "DgParamList.h"
#include "DgProjGnomonicRF.h"
#include "DgGeoProjConverter.h"
#include "DgRandom.h"
#include "DgCell.h"
#include "DgLocList.h"
#include "DgDmdD4Grid2D.h"
#include "DgDmdD4Grid2DS.h"
#include "DgTriGrid2D.h"
#include "DgOutRandPtsText.h"
#include "DgProjISEA.h"
#include "DgHexSF.h"

using namespace std;
typedef int Status;
#define ERR_OK 0
#define ERR_NO 1
#define AREA_Half 371.275969/2.0
#define AREA 371.275969
#define DIST 30
#define NLEVEL 18
//quam=5时，三角形公共边旋转到x轴的旋转参数
#define _AX 581832.459600000
#define _AY 4201628.67820000
#define _BX 581850.358400000
#define _BY 4201610.61650000
#define COSt 0.703896386240366
#define SINt -0.710302666078168

#define _HEIGHT 29 //三角形的高度
#define _TYPE_ 30 //地类类别10 30 40 50 60 90 100 
#define random(x) (rand() % x)

typedef struct stackNode {
	DgQ2DICoord COOR;
	stackNode * next;
}SeqStack;
SeqStack * init_SeqStack(SeqStack * top);
int is_Empty(SeqStack * top);
SeqStack * push_Stack(SeqStack * top, DgQ2DICoord coor);
SeqStack * pop_Stack(SeqStack * top, DgQ2DICoord &coor);
SeqStack * top_Stack(SeqStack * top, DgQ2DICoord &coor);

//定义节点
typedef struct queueNode {
	DgQ2DICoord COOR;
	struct queueNode * next;
}queueNode;

//定义队列（保存队首和队尾指针）
typedef struct queue_link {
	queueNode * front;
	queueNode * rear;
}que;

typedef struct {
	int numpoint;
	Vec2D *pointlist;
}polygon;

typedef struct {
	int polyID;//多边形索引
	int Class;//多边形属性
	int BDRaster_count;//边界栅格总数
	vector<DgQ2DICoord> QIJ;//多边形边界栅格
	vector<double> Ratio;//边界栅格与多边形相交面积，与QIJ对应
}Poly_BDRasters;

vector<Vec2D> GridVercoord(DgPolygon cell, OGRCoordinateTransformation *coordTransInv);
int getImgInfo(char *szInFile, GDALDataset **poDataset, int *nbands, double **geoTrans, int *width, int *height, GDALDataType *gdt, const char** projRef, GDALRasterBand *** poBand, OGRCoordinateTransformation **pocoordTrans, OGRCoordinateTransformation **pocoordTransInv);
unsigned short * getImgData(GDALDataset *poDataset, int nbands, int width, int height, GDALDataType *gdt);
void getValByPolygonIntesect(vector<Vec2D> gridVerCor,
	unsigned short *imgData, double *geoTrans, int width, int height, int *val);
void genGrid_(GridGenParam& dp, char *SHP_path, string polyID);
void DGGRID_GLC30(DgGridPList &plist, char *SHP_path, int **GLC30, string polyID, char * txt_path, char * outshp_path, char * outtxt_Path);
void setParam_(DgGridPList *PLIST, int level);
int GetClass(int gridval);
void GetPath(string *shppath, string *txtPath, string in);

void Grid_located(GridGenParam& dp, unsigned short * imgData, int Width, int Height, double *geoTrans,
	OGRCoordinateTransformation *coordTrans, OGRCoordinateTransformation *coordTransInv, int **GLC30);
void Grid_Voronoi(GridGenParam& dp, unsigned short * imgData, int Width, int Height, double *geoTrans,
	OGRCoordinateTransformation *coordTrans, OGRCoordinateTransformation *coordTransInv, int **GLC30);
vector<DgQ2DICoord> Edge_Proximity(DgLocation& Loc, int Quam, int i, int j, DgPolygon verts, int maxI, int maxJ);
void LL2Coord(DgQ2DICoord coor, double lon, double lat,
	DgGeoSphRF geoRF, DgIDGG dgg, int nDensify);
void QIJ2LL(DgQ2DICoord coor, DgLocation& Loc, DgPolygon verts,
	DgGeoSphRF geoRF, DgIDGG dgg, int nDensify);
bool is_Same_Side(vector<Vec2D> gridVerCor, Vec2D A, Vec2D B);
bool is_In_On_Polygon(vector<Vec2D> gridVerCor, Vec2D A);
double clipper_intersection(vector<Vec2D> grid, vector<Vec2D> BD);
double Angle(double pix, double piy, double pjx, double pjy);
bool PolygonOrti(vector<Vec2D> gridVerCor);
vector<DgQ2DICoord> BDT_EdgeADJ(DgQ2DICoord cood_i, bool Orti, double Angle);
vector<DgQ2DICoord> BDT_AllADJ(DgQ2DICoord cood_i, bool Orti);
vector<DgQ2DICoord> Edge_Proximity(DgQ2DICoord cood_i, bool UP, int maxI, int maxJ);
double Point2LineDist(double x1, double y1, double p1x, double p1y, double p2x, double p2y);
double PointsDist(double p1x, double p1y, double p2x, double p2y);
bool Point_internal(double x1, double y1, double p1x, double p1y, double p2x, double p2y);
void Delete_same(vector<DgQ2DICoord> *coor, vector<DgQ2DICoord> COOD);
double clipper_intersection_area(double x1, double y1, double x2, double y2,
	double x3, double y3, OGRLineString *boundary);
bool Line_Grid_intersection(vector<Vec2D> gridVerCor, double x1, double y1, double x2, double y2,
	double cost, double sint, double AB);
bool Line_Line_intersection(double Ax, double Ay, double Bx, double By,
	double Px, double Py, double Qx, double Qy, double AB, double cost, double sint);
Vec2D Grid_center(vector<Vec2D> gridVerCor);
OGRGeometry* ConvertPolygonToPolyline(OGRGeometry* polygon);
OGRGeometry* ConvertPolylineToPolygon(OGRGeometry* polyline);
void SplitString(string line, int **str2i);
void SplitString_OURSGRID(string line, int **STR2I);
int BDOri(double xi, double xj, double yi, double yj);
double clipper_intersection_area2(double v1_x, double v1_y, double v2_x, double v2_y,
	double v3_x, double v3_y, Paths subject);
bool PointInPolygon(double x, double y, Paths subject, int pathID);
void BDTracing(GridGenParam& dp, char *SHP_path, string polyID, string outshp_path, string outtxt_path);
void getEdgeGridSHPwithType(GridGenParam& dp, char *SHP_path, char *txt_path, char *outshp_path);
void ShpSnyderFwd(GridGenParam& dp, char *SHP_path, char *outshp_path);


#endif // !_GLC30_BASED_H_

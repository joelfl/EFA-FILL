#pragma once
#ifndef _GLC30_BASED_H_
#define _GLC30_BASED_H_


#include "gdal.h" 
#include "gdal_priv.h" 
#include <ogr_spatialref.h>
#include <gdalwarper.h>
#include <ogr_api.h>
#include <ogrsf_frmts.h>
#include <omp.h>
#include <dggrid.h>
#include <gridgen.h>
#include <clipper_region.h>
//#include <iostream>
//#include <istream>
//#include <string>
#include <fstream>
#include <sstream> 
//#include <string>

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
#include "DgHexSF.h"

#define Grid_Area 0.38e06//level=13
		//6.0805e06//level=11
		//389.1492e06//level=8
#define Area_Threshold 0.68*Grid_Area
#define UTMN 48
#define NLEVEL 18

class Val_ {
public:
	int val;
	bool used = false;
	int index[19] = { 0 };
};
class QuadVals_ {
public:
	bool isUsed;
	DgIVec2D offset;     // offset of min (i, j)
	DgIVec2D upperRight; // (maxi, maxj) relative to offset
	int numI;
	int numJ;
	Val_** vals;
};

class BD_GRID {
public:
	double area;
	DgQ2DICoord coor;
};

class Pt_GRID
{
public:
	DgQ2DICoord coor;
	int ID;
	double ai;
	vector<int> polyID;
	vector<double> area;
};
class Pt_topo
{
public:
	int ID;
	double lon, lat;
	vector<int> polyID;
};

class BoundaryG
{
public:
	int PolyID;
	DgQ2DICoord coor;
	double area;
	bool IS = 0;//标记coor是否属于当前PolyID.
	double weight_node = 1;
};

class BoundarGs
{
public:
	vector<BoundaryG> bdg;
};

class Counting
{
public:
	int PolyID;
	int current_count;
	int actual_count;
};
class tmpClass1
{
public:
	int index;
	int deltaN;
	int currentN;
	int actualN;
};
class Location
{
public:
	int pos1;
	int pos2;
	int DistIndexPOS;
};

using namespace std;
vector<Vec2D> GridVercoord(DgPolygon cell, OGRCoordinateTransformation *coordTransInv);
int getImgInfo(char *szInFile, GDALDataset **poDataset, int *nbands, double **geoTrans, int *width, int *height, GDALDataType *gdt, const char** projRef, GDALRasterBand *** poBand, OGRCoordinateTransformation **pocoordTrans, OGRCoordinateTransformation **pocoordTransInv);
unsigned char * getImgData(GDALDataset *poDataset, int nbands, int width, int height, GDALDataType *gdt);
void getValByPolygonIntesect(vector<Vec2D> gridVerCor,
	unsigned char *imgData, double *geoTrans, int width, int height, int *val);
void genGrid_(GridGenParam& dp,
	OGRSpatialReference oSRS,
	OGRCoordinateTransformation *coordTrans,
	OGRCoordinateTransformation *coordTransInv,
	OGRGeometry* poGeometry,
	char *outpath,
	double PolyGon_Area, int type, int polyID, int &count,
	vector<BoundaryG> &bd_grid);
void binValsPartial4RandomPt_(BinValsParam& dp,
	OGRCoordinateTransformation *coordTrans, OGRCoordinateTransformation *coordTransInv,
	char * PtTXT, char * PolygonSHPPath, vector<Pt_GRID> &Pts_Grid);
void binValsPartial4RandomPt_Res(BinValsParam& dp, unsigned char * imgData, int Width, int Height, double *geoTrans,
	OGRCoordinateTransformation *coordTrans, OGRCoordinateTransformation *coordTransInv, 
	unsigned char * imgData1, int Width1, int Height1, double *geoTrans1, string RandomPtpath, string outpath);
void DGGRID_GLC30(DgGridPList &plist, char * SHPPath, char * MBRPath, char * txt_path, char * OutPath, int *NIGN);
void setParam_(DgGridPList *PLIST, int level);
int GetClass(int gridval);
void GetPath(string *shppath, string *txtPath, string in);
OGRLineString *SHPBoundary(char *shppath);
double clipper_intersection_area(double v1_x, double v1_y, double v2_x, double v2_y,
	double v3_x, double v3_y, OGRLineString *poly);
double clipper_intersection_area1(double v1_x, double v1_y, double v2_x, double v2_y,
	double v3_x, double v3_y, double v4_x, double v4_y, Paths subject);//OGRLineString *poly
double evalCell_(GridGenParam& dp, const DgIDGG& dgg, const DgContCartRF& cc1,
	const DgDiscRF2D& grid, OGRGeometry* poGeometry, DgQuadClipRegion& clipRegion,
	const DgIVec2D& add2D);
int *Sort_bdgrid(int total, vector<BD_GRID> bd_grid);
double BD_Grid_AREA(OGRGeometry* poGeometry, vector<Vec2D> gridVerCor);
void PolygonPTGrid(DgGridPList &plist, char * SHPPath, char * PTtxt, vector<Pt_GRID> &Pts_grid);
bool IS_same_ptgrid(vector<Pt_GRID> Pts_grid, DgQ2DICoord coor);
bool IS_same_pttopo(vector<Pt_topo> Pts_topo, vector<string> iter);
vector<string> Gen_PtTopo(char * Ptpath);
void BDGTraver(vector<BoundaryG> &bdg, int PolyID, int T);
void NodeBoundaryG(vector<BoundaryG> bd_grid, BoundarGs *&bdg);
void NodeReDistr(BoundarGs *&bdg, vector<Counting> Counter);
void frquency(vector<tmpClass1> data, int &count, int &val);
bool isused(DgQ2DICoord coor, vector<BoundaryG> bdg);
void BDGWight(vector<BoundaryG> &bdg);
bool Line_Line_intersection(double Ax, double Ay, double Bx, double By, double Px, double Py, double Qx, double Qy);
bool Line_BD_intersect(double Ax, double Ay, double Bx, double By, OGRGeometry *poGeometry);
void Coor2SHP(DgGridPList &plist, vector<BoundaryG> bdg, OGRSpatialReference oSRS, vector<Counting> counter);
void Coor2SHP1(DgGridPList &plist, BoundarGs *bdg, OGRSpatialReference oSRS, vector<Counting> counter);
int *sortdouble(int n, float* data);
BoundarGs * GetMissedGrid(BoundarGs *bdg);
int FindActualCounter(int polyid, vector<Counting> counter);
int MinDistID(double *dist, int size);
double CoorMBRDist(double cx, double cy, OGREnvelope *poEnvelope);
int Find_CounterPolyID(vector<Counting> counter, int polyid);
vector<BoundaryG> Find_PolyIDBDG(BoundarGs *bdg, int polyid);
int * Dist2MBR(DgGridPList &plist, vector<BoundaryG> bdgi, char *MBRShpPath);
void Find_PolyIDCoor(BoundarGs *bdg, int polyID, DgQ2DICoord coor, int &i, int &j);
bool Same_LOC(vector<Location> LOC, int pos1, int pos2);
int PolyGonChange(int polyid, DgQ2DICoord coor, BoundarGs *bdg, vector<Counting> counter);
int MaxPolyID1(BoundarGs *bdg, vector<Counting> counter, int PolyID, vector<int> &usePolyID, DgQ2DICoord coor);
int CandidatePolygon(DgGridPList &plist, int polyID, char *MBRShpPath, BoundarGs *&bdg, vector<Counting> counter, vector<int> & usePolyID, DgQ2DICoord &coord);
int MaxCountPolyID(vector<BoundaryG> bdg,vector<Counting> counter, int i, int j, int &pos);
void SpaceDistConst(DgGridPList &plist, char *MBRShpPath, BoundarGs *&bdg, vector<Counting> counter);
Vec2D Grid_center(vector<Vec2D> gridVerCor);
void GetRotate(double Ax, double Ay, double Bx, double By, double sint, double cost);
bool Line_Line_intersection1(double Ax, double Ay, double sint, double cost, double Px, double Py, double Qx, double Qy, Vec2D *point);
vector<Vec2D> GetIntersectPoint(Vec2D A, Vec2D B, OGRGeometry *poGeometry);
float *evaluation(vector<Vec2D> POINT);
int *evaluation_(vector<Vec3D> QIJ);
Vec2D FindSamePoint(vector<Vec2D> P1, vector<Vec2D> P2);
Paths PreDefineSubject(OGRLineString *poly);
void ScanlineFill(GridGenParam& dp, char *SHP_path, char *SHP_out_path, double x0, double y0, double x3, double y3);
int ScanlineFill_(GridGenParam& dp, char *SHP_path, char *txt_path, char *SHP_out_path, double x0, double y0, double x3, double y3);
void FillPolygong4Paper2(GridGenParam& dp, char *SHP_path, char *txt_path, char *SHP_out_path, char *NGIN_path);
#endif // !_GLC30_BASED_H_

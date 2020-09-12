#include "GLC30_Based.h"


int main()
{
	clock_t start, finish;
	double duration = 0;
	start = clock();
	GDALAllRegister();

	int GLC30[12] = { 0 };

	int argc = 1;
	char* token = "dggrid_operation";
	char* remainder = "GENERATE_GRID";
	{ 
		ofstream outfile;
		int *tmp;
		DgGridPList plist;
		setParam_(&plist, NLEVEL);
		int k = 0;
		string filename;
		string line;
		string tifPath, shpPath, txtPath, outshpPath;
		shpPath = "..\\..\\..\\data\\shp\\region_proj.shp";
		plist.setParam("clip_region_files", (char *)(shpPath.data()));
		shpPath = "..\\..\\..\\data\\shp\\region.shp";
		txtPath = "..\\..\\..\\data\\txt.txt";
		outshpPath = "..\\..\\..\\data\\result\\edgegrid\\" + to_string(_TYPE_) + "_";
		string outtxtPath = "..\\..\\..\\data\\result\\edgegrid_txt\\" + to_string(_TYPE_) + "_";
		start = clock();
		DGGRID_GLC30(plist, (char*)(shpPath.data()), &tmp, line, (char *)(txtPath.data()), (char *)(outshpPath.data()), (char *)(outtxtPath.data()));
		finish = clock();
		duration = finish - start;
	}
	GDALDestroyDriverManager();
	system("pause");
}
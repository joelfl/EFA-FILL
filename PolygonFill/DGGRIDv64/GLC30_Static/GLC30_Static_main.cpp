#include "GLC30_Based.h"

int main()
{
	clock_t start, finish;
	double duration = 0;
	start = clock();
	float startTime = omp_get_wtime();
	GDALAllRegister();
	unsigned int GLC30[12] = { 0 };
	int argc = 1;
	char* token = "dggrid_operation";
	char* remainder = "GENERATE_GRID";
	ofstream outfile;
	int *tmp;
	DgGridPList plist;
	setParam_(&plist, NLEVEL);
	int NGIN = 0;
	int k = 0;
	ifstream in("D:\\data\\paper2_Shandong\\PropertyforFill\\filename.txt");//以TIF中的文件为标准
	string filename;
	string line;
	string tifPath, shpMBRPath, shpPath, shpPathProj, txtPath, outPath, PTtxt;
	PTtxt = "D:\\data\\SAMPLE\\ptout.txt";
	vector<Pt_GRID> Pts_grid;
	vector<BoundaryG> bd_grid;
	vector<string> FileLine;
	shpPath = "..\\..\\..\\Data\\Shp\\region_proj.shp";
	txtPath = "..\\..\\..\\Data\\result\\edgegrid_txt\\";
	outPath = "..\\..\\..\\Data\\result\\interiorgrid\\ ";
	cout << line << endl;
	plist.setParam("clip_region_files", (char *)(shpPath.data()));
	DGGRID_GLC30(plist, (char*)(shpPath.data()), (char*)(shpPath.data()),
		(char*)(txtPath.data()), (char*)(outPath.data()), &NGIN);
	GDALDestroyDriverManager();
	double endTime = omp_get_wtime();
	double time_tmp;
	time_tmp = endTime - startTime;
	system("pause");
}
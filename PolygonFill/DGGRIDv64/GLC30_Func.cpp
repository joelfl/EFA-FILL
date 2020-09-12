#include "GLC30_Based.h"

void GetPath(string *shppath, string *txtPath, string in)
{
	int pos = 0, pos1 = 0, pos2 = 0;
	string smybol1 = "2000_";
	string smybol2 = "2010_";
	pos = in.find(".tif");
	in.erase(pos, 4);
	*txtPath = in;
	pos1 = in.find(smybol1);
	pos2 = in.find(smybol2);
	pos = pos1 < 0 ? pos2 : pos1;
	in.erase(pos, 5);
	*shppath = in;
}

vector<Vec2D> GridVercoord(DgPolygon cell, OGRCoordinateTransformation *coordTransInv)
{
	const DgGeoSphRF* geoRF = dynamic_cast<const DgGeoSphRF*>(&cell.rf());
	vector<Vec2D> gridVerCor;
	int count = cell.size();
	for (int i = 0; i < count; i++)//此处即可得到cell的顶点坐标和中心点坐标
								   //	                                 
								   //直接使用verts时，所得点区域时不闭合的，此处要注意
	{
		const DgGeoCoord& v = *geoRF->getAddress(cell[i]);
		double x = v.lat();
		double y = v.lon();
		/*x = x * M_PI_180;
		y = y * M_PI_180;*/
		x = x * M_180_PI;
		y = y * M_180_PI;
		coordTransInv->Transform(1, &y, &x);//(lon,lat)
		Vec2D tmp;
		tmp.x = y;
		tmp.y = x;
		gridVerCor.push_back(tmp);
	}
	return gridVerCor;
}


int getImgInfo(char *szInFile, GDALDataset **poDataset, int *nbands, double **geoTrans, int *width, int *height, GDALDataType *gdt, const char** projRef, GDALRasterBand *** poBand, OGRCoordinateTransformation **pocoordTrans, OGRCoordinateTransformation **pocoordTransInv)
{

	GDALDataset *poDatasetTmp = *poDataset;
	poDatasetTmp = (GDALDataset*)GDALOpen(szInFile, GA_ReadOnly);

	//poDatasetTmp = (GDALDataset*)GDALOpen(szInFile.c_str(), GA_ReadOnly);

	int widthTmp = *width, heightTmp = *height, nbandsTmp = *nbands;
	widthTmp = poDatasetTmp->GetRasterXSize();
	heightTmp = poDatasetTmp->GetRasterYSize();
	nbandsTmp = poDatasetTmp->GetRasterCount();

	GDALDataType gdtTmp = *gdt;
	gdtTmp = poDatasetTmp->GetRasterBand(1)->GetRasterDataType();

	double *geoTransTmp = *geoTrans;
	geoTransTmp = new double[6];
	poDatasetTmp->GetGeoTransform(geoTransTmp);//获取地理坐标信息，地理坐标信息是一个含6个double型数据的数组，
	const char* projRefTmp = *projRef;
	projRefTmp = poDatasetTmp->GetProjectionRef();  //获取投影信息

	GDALRasterBand ** poBandTmp = *poBand;
	poBandTmp = new GDALRasterBand *[nbandsTmp];
	if (poBand == NULL)
	{
		cout << "GDALRasterBand ** poBand = new GDALRasterBand *[nBands]; failed!" << endl;
	}
	for (int i = 0; i < nbandsTmp; i++)
	{
		poBandTmp[i] = poDatasetTmp->GetRasterBand(i + 1);
	}
	OGRSpatialReference oSRS = OGRSpatialReference(poDatasetTmp->GetProjectionRef());//直接将tif的Dataset中的project信息添加到spatialreference中即可；
	OGRSpatialReference oTRS;
	oTRS.SetWellKnownGeogCS("WGS84");
	OGRCoordinateTransformation *coordTrans = OGRCreateCoordinateTransformation(&oSRS, &oTRS);//from plane to sphere
	OGRCoordinateTransformation *coordTransInv = OGRCreateCoordinateTransformation(&oTRS, &oSRS);//from sphere to plane

	*poDataset = poDatasetTmp;
	*nbands = nbandsTmp;
	*geoTrans = geoTransTmp;
	*width = widthTmp;
	*height = heightTmp;
	*gdt = gdtTmp;
	*projRef = projRefTmp;
	*poBand = poBandTmp;
	*pocoordTrans = coordTrans;
	*pocoordTransInv = coordTransInv;

	//释放内存
	//free(geoTransTmp);
	//GDALClose(poDatasetTmp);

	return 0;
}

unsigned char * getImgData(GDALDataset *poDataset, int nbands, int width, int height, GDALDataType *gdt)
{
	//GDT_Byte: 8bit正整型(C++中对应unsigned char) MODIS数据
	//GDT_UInt16 : 16bit正整型 (C++中对应 unsigned short)
	unsigned char *imgBuf = new unsigned char[width*height];
	if (nbands == 1)
	{
		unsigned char *imgBuf = new unsigned char[width*height];
		GDALRasterBand *poBand = poDataset->GetRasterBand(1);
		int re = poBand->RasterIO(GF_Read, 0, 0, width, height, imgBuf, width, height, GDT_Byte, 0, 0);
		return imgBuf;
	}
	else
	{
		return imgBuf;
	}
}

int GetClass(int gridval)
{
	if (gridval == 254)
		return 17;
	if (gridval == 255)
		return 18;
	return gridval;
	/*if (gridval == 0)
		return 11;
	else if (gridval == 255)
	{
		return 10;
	}
	else if (gridval == 100)
	{
		return 9;
	}
	else if (gridval == 90)
	{
		return 8;
	}
	else if (gridval == 80)
	{
		return 7;
	}
	else if (gridval == 70)
	{
		return 6;
	}
	else if (gridval == 60)
	{
		return 5;
	}
	else if (gridval == 50)
	{
		return 4;
	}
	else if (gridval == 40)
	{
		return 3;
	}
	else if (gridval == 30)
	{
		return 2;
	}
	else if (gridval == 20)
	{
		return 1;
	}
	else if (gridval == 10)
	{
		return 0;
	}
	return 11;*/
}

int sig(double d) {
	return(d>1E-8) - (d<-1E-8);
}
struct Point__ {
	double x, y; Point__() {}
	Point__(double x, double y) :x(x), y(y) {}
	bool operator==(const Point__&p)const {
		return sig(x - p.x) == 0 && sig(y - p.y) == 0;
	}
};
double cross(Point__ o, Point__ a, Point__ b) {
	return(a.x - o.x)*(b.y - o.y) - (b.x - o.x)*(a.y - o.y);
}
double cross_(double o_x, double o_y, double a_x, double a_y, double b_x, double b_y) {
	return(a_x - o_x)*(b_y - o_y) - (b_x - o_x)*(a_y - o_y);
}

double area(Point__* ps, int n) {
	ps[n] = ps[0];
	double res = 0;
	for (int i = 0; i<n; i++) {
		res += ps[i].x*ps[i + 1].y - ps[i].y*ps[i + 1].x;
	}
	return res / 2.0;
}
double area_(double* ps_x, double *ps_y, int n) {
	ps_x[n] = ps_x[0];
	ps_y[n] = ps_y[0];

	double res = 0;
	for (int i = 0; i<n; i++) {
		res += ps_x[i]* ps_y[i + 1] - ps_y[i]*ps_x[i + 1];
	}
	return res / 2.0;
}

int lineCross(Point__ a, Point__ b, Point__ c, Point__ d, Point__&p) {
	double s1, s2;
	s1 = cross(a, b, c);
	s2 = cross(a, b, d);
	if (sig(s1) == 0 && sig(s2) == 0) return 2;
	if (sig(s2 - s1) == 0) return 0;
	p.x = (c.x*s2 - d.x*s1) / (s2 - s1);
	p.y = (c.y*s2 - d.y*s1) / (s2 - s1);
	return 1;
}

int lineCross_(double a_x, double a_y, double b_x, double b_y, 
	double c_x, double c_y, double d_x, double d_y,
	double *p_x, double *p_y) {
	double s1, s2;
	s1 = cross_(a_x, a_y, b_x, b_y, c_x, c_y);
	s2 = cross_(a_x, a_y, b_x, b_y, d_x, d_y);
	if (sig(s1) == 0 && sig(s2) == 0) return 2;
	if (sig(s2 - s1) == 0) return 0;

	*p_x = (c_x*s2 - d_x*s1) / (s2 - s1);
	*p_y = (c_y*s2 - d_y*s1) / (s2 - s1);
	return 1;
}

//多边形切割
//用直线ab切割多边形p，切割后的在向量(a,b)的左侧，并原地保存切割结果
//如果退化为一个点，也会返回去,此时n为1
void polygon_cut(Point__*p, int&n, Point__ a, Point__ b) {
	static Point__ pp[510];
	int m = 0; p[n] = p[0];
	for (int i = 0; i<n; i++) {
		if (sig(cross(a, b, p[i]))>0) pp[m++] = p[i];
		if (sig(cross(a, b, p[i])) != sig(cross(a, b, p[i + 1])))
			lineCross(a, b, p[i], p[i + 1], pp[m++]);
	}
	n = 0;
	for (int i = 0; i<m; i++)
		if (!i || !(pp[i] == pp[i - 1]))
			p[n++] = pp[i];
	while (n>1 && p[n - 1] == p[0])n--;
}
void polygon_cut_(double *p_x, double *p_y, int *n, double a_x, double a_y, double b_x, double b_y) {
	static double pp_x[10];
	static double pp_y[10];

	int m = 0; 
	int tmp_n = *n;
	p_x[tmp_n] = p_x[0];
	p_y[tmp_n] = p_y[0];

	for (int i = 0; i<tmp_n; i++) {
		if (sig(cross_(a_x, a_y, b_x, b_y, p_x[i], p_y[i])) > 0)
		{
			pp_x[m] = p_x[i];
			pp_y[m] = p_y[i];
			m++;
		}
		if (sig(cross_(a_x, a_y, b_x, b_y, p_x[i], p_y[i])) != sig(cross_(a_x, a_y, b_x, b_y, p_x[i + 1], p_y[i + 1])))

		{
			lineCross_(a_x, a_y, b_x, b_y, p_x[i], p_y[i], p_x[i + 1], p_y[i + 1], &pp_x[m], &pp_y[m]);
			m++;
		}
		}
	tmp_n = 0;
	for (int i = 0; i<m; i++)
		if (!i || !(pp_x[i] == pp_x[i - 1] && pp_y[i] == pp_y[i - 1]))
		{
			p_x[tmp_n] = pp_x[i];
			p_y[tmp_n] = pp_y[i];
			tmp_n++;
		}
	while (tmp_n>1 && p_x[tmp_n - 1] == p_x[0] && p_y[tmp_n - 1] == p_y[0]) tmp_n--;
	*n = tmp_n;
}

//---------------华丽的分隔线-----------------//
//返回三角形oab和三角形ocd的有向交面积,o是原点//
double intersectArea(Point__ a, Point__ b, Point__ c, Point__ d) {
	Point__ o(0, 0);
	int s1 = sig(cross(o, a, b));
	int s2 = sig(cross(o, c, d));
	if (s1 == 0 || s2 == 0)return 0.0;//退化，面积为0
	if (s1 == -1) swap(a, b);
	if (s2 == -1) swap(c, d);
	Point__ p[10] = { o,a,b };
	int n = 3;
	polygon_cut(p, n, o, c);
	polygon_cut(p, n, c, d);
	polygon_cut(p, n, d, o);
	double res = fabs(area(p, n));
	if (s1*s2 == -1) res = -res; return res;
}
//求两多边形的交面积
double intersectArea(Point__*ps1, int n1, Point__*ps2, int n2) {

	if (area(ps1, n1)<0) reverse(ps1, ps1 + n1);
	if (area(ps2, n2)<0) reverse(ps2, ps2 + n2);//若面积小于0，则将点顺序倒置
	ps1[n1] = ps1[0];
	ps2[n2] = ps2[0];
	double res = 0;
	for (int i = 0; i<n1; i++) {
		for (int j = 0; j<n2; j++) {
			res += intersectArea(ps1[i], ps1[i + 1], ps2[j], ps2[j + 1]);
		}
	}
	return res;//assumeresispositive!
}

double intersectArea_(double a_x, double a_y, double b_x, double b_y, 
	double c_x, double c_y, double d_x, double d_y) {
	double res = 0;
	double ox = 0, oy = 0;
	int s1 = sig(cross_(ox, oy, a_x, a_y, b_x, b_y));
	int s2 = sig(cross_(ox, oy, c_x, c_y, d_x, d_y));
	if (s1 == 0 || s2 == 0)return 0.0;//退化，面积为0
	if (s1 == -1)
	{
		swap(a_x, b_x);
		swap(a_y, b_y);
	}
	if (s2 == -1)
	{
		swap(c_x, d_x);
		swap(c_y, d_y);
	}
	double p_x[10] = { ox,a_x,b_x };
	double p_y[10]= { oy,a_y,b_y };
	int n = 3;
	polygon_cut_(p_x, p_y, &n, ox, oy, c_x, c_y);
	polygon_cut_(p_x, p_y, &n, c_x, c_y, d_x, d_y);
	polygon_cut_(p_x, p_y, &n, d_x, d_y, ox, oy);
	res = fabs(area_(p_x, p_y, n));
	if (s1*s2 == -1) 
		res = -res; 
	return res;
}
double intersectArea__(double *ps1_x,  double *ps1_y, int n1, double *ps2_x, double *ps2_y, int n2) {

	if (area_(ps1_x, ps1_y, n1)<0)
	{
		reverse(ps1_x, ps1_x+ n1);
		reverse(ps1_y, ps1_y + n1);
		//double tmp_[4];
		//for (int i = 0; i < n1; i++)
		//{
		//	tmp_[i] = ps1_x[n1 - i - 1];
		//}
		//reverse(ps1_x, ps1_x + n1);
	}
	double tmp_[4];
	if (area_(ps2_x, ps2_y, n1)<0)
	{
		reverse(ps2_x, ps2_x + n2);
		reverse(ps2_y, ps2_y + n2);
	}
	ps1_x[n1] = ps1_x[0];
	ps1_y[n1] = ps1_y[0];

	ps2_x[n2] = ps2_x[0];
	ps2_y[n2] = ps2_y[0];

	double res = 0;
	for (int i = 0; i<n1; i++) {
		for (int j = 0; j<n2; j++) {
			res += intersectArea_(ps1_x[i], ps1_y[i], ps1_x[i + 1], 
				ps1_y[i+1], ps2_x[j], ps2_y[j], ps2_x[j + 1], ps2_y[j + 1]);
		}
	}
	return res;//assumeresispositive!
}


void getValByPolygonIntesect(vector<Vec2D> gridVerCor,
	unsigned char *imgData, double *geoTrans, int width, int height, int *val)
{
	//定义变量
	int k;
	int count, total_size;
	int col, row;
	int center_col, center_row;
	int glc_class = 11, img_k;
	double g0, g1, g2, g3, g4, g5;
	double x, y, xi, yi, xj, yj;
	double center_x, center_y;
	double ai, area;
	double gird_area[19] = { 0.0 };
	Point__ ps1[5], ps2[5];

	//平面坐标转行列号，根据tif geoTrans
	//vector<Vec2D> IJ = GridCor2IJ(gridVerCor, geoTrans);
	//精确定位IJ
	//遍历IJ，确定IJ的9邻域，确定每个pixel的地理坐标，用gridVerCor判断
	vector<Vec3D> IJ;
	Vec2D ij;
	Vec3D tmp;
	count = 0;
	total_size = (int)gridVerCor.size();
	g0 = geoTrans[0];
	g1 = geoTrans[1];
	g2 = geoTrans[2];//0
	g3 = geoTrans[3];
	g4 = geoTrans[4];//0
	g5 = geoTrans[5];
	//遍历gridVerCor，获得中心点坐标
	x = 0;
	y = 0;
	for (k = 0; k < total_size; k++)
	{
		x += gridVerCor[k].x;
		y += gridVerCor[k].y;
	}
	center_x = x / total_size;
	center_y = y / total_size;
	//根据中心点坐标计算行列号
	center_col = (int)((center_x - g0) / g1);
	center_row = (int)((center_y - g3) / g5);
	//根据中心点行列号遍历九邻域
	for (k = 0; k < 9; k++)
	{
		if (k == 0)
		{
			col = center_col;
			row = center_row;
		}
		else if (k == 1)
		{
			col = center_col - 1;
			row = center_row - 1;
		}
		else if (k == 2)
		{
			col = center_col - 1;
			row = center_row;
		}
		else if (k == 3)
		{
			col = center_col - 1;
			row = center_row + 1;
		}
		else if (k == 4)
		{
			col = center_col;
			row = center_row - 1;
		}
		else if (k == 5)
		{
			col = center_col;
			row = center_row + 1;
		}
		else if (k == 6)
		{
			col = center_col + 1;
			row = center_row - 1;
		}
		else if (k == 7)
		{
			col = center_col + 1;
			row = center_row;
		}
		else if (k == 8)
		{
			col = center_col + 1;
			row = center_row + 1;
		}
		xi = g0 + col * g1 + row * g2;//像素的左上角坐标
		yi = g3 + col * g4 + row * g5;
		xj = g0 + (col + 1) * g1 + (row + 1) * g2;//像素的右下角坐标
		yj = g3 + (col + 1) * g4 + (row + 1) * g5;

		//引入clipper.lib
		//ISEA4T
		if (total_size == 3)
			ai = clipper_intersection_area(gridVerCor[0].x, gridVerCor[0].y,
				gridVerCor[1].x, gridVerCor[1].y,
				gridVerCor[2].x, gridVerCor[2].y,
				xi, yi, xj, yj);
		////ISEA4D
		/*if (total_size == 4)
			ai = clipper_intersection_area_D(gridVerCor[0].x, gridVerCor[0].y,
				gridVerCor[1].x, gridVerCor[1].y,
				gridVerCor[2].x, gridVerCor[2].y,
				gridVerCor[3].x, gridVerCor[3].y,
				xi, yi, xj, yj);*/
		
		////ISEA4H
		//if (total_size == 6)
		//	ai = clipper_intersection_area_H6(gridVerCor[0].x, gridVerCor[0].y,
		//		gridVerCor[1].x, gridVerCor[1].y,
		//		gridVerCor[2].x, gridVerCor[2].y,
		//		gridVerCor[3].x, gridVerCor[3].y,
		//		gridVerCor[4].x, gridVerCor[4].y,
		//		gridVerCor[5].x, gridVerCor[5].y,
		//		xi, yi, xj, yj);
		//if (total_size == 5)
		//	ai = clipper_intersection_area_H5(gridVerCor[0].x, gridVerCor[0].y,
		//		gridVerCor[1].x, gridVerCor[1].y,
		//		gridVerCor[2].x, gridVerCor[2].y,
		//		gridVerCor[3].x, gridVerCor[3].y,
		//		gridVerCor[4].x, gridVerCor[4].y,
		//		xi, yi, xj, yj);
		//越界判断
		if (col < 0)
			col = 0;
		if (row < 0)
			row = 0;
		if (col >= width)
			col = width - 1;
		if (row >= height)
			row = height - 1;
		img_k = row * width + col;
		*val = imgData[img_k];
		glc_class = GetClass(*val);
		double ps1_x[5], ps1_y[5];
		double ps2_x[5], ps2_y[5];

		/*ps1_x[0] = gridVerCor[0].x;
		ps1_y[0] = gridVerCor[0].y;
		ps1_x[1] = gridVerCor[1].x;
		ps1_y[1] = gridVerCor[1].y;
		ps1_x[2] = gridVerCor[2].x;
		ps1_y[2] = gridVerCor[2].y;
		ps1_x[3] = gridVerCor[3].x;
		ps1_y[3] = gridVerCor[3].y;
		ps1_x[4] = gridVerCor[0].x;
		ps1_y[4] = gridVerCor[0].y;

		ps2_x[0] = xi;
		ps2_y[0] = yi;
		ps2_x[1] = xj;
		ps2_y[1] = yi;
		ps2_x[2] = xj;
		ps2_y[2] = yj;
		ps2_x[3] = xi;
		ps2_y[3] = yj;
		ps2_x[4] = xi;
		ps2_y[4] = yi;*/
		//ai = intersectArea__(ps1_x, ps1_y, 4, ps2_x, ps2_y, 4);
		gird_area[glc_class] += ai;

		//if (abs(ai - area_new_) > 1)
		//{
		//	cout << "un equal" << endl;
		//}
	}
	area = -10;
	glc_class = -1;
	for (k = 0; k < 19; k++)
	{
		ai = gird_area[k];
		if (ai >= area)
		{
			area = ai;
			glc_class = k;
		}
	}
	if (area == 0)
		glc_class = 18;

	*val = glc_class;
	//vector<Vec3D>().swap(IJ);

	//return;
}

void genGrid_(GridGenParam& dp, unsigned char * imgData, int Width, int Height, double *geoTrans,
	OGRCoordinateTransformation *coordTrans, OGRCoordinateTransformation *coordTransInv, int **GLC30)
{
	int GLC[19] = { 0 };
	int gridval = 0;
	int cols[4] = { 0, Width, 0, Width };
	int rows[4] = { 0, 0,  Height, Height };
	/////////
	//string shpoutpath = "D:\\data\\SAMPLE\\N03_05_2000LC030_41.shp";
	//GDALDriver *poDriver;
	//poDriver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
	////创建shp文件
	//GDALDataset *poDS = poDriver->Create(shpoutpath.data(), 0, 0, 0, GDT_Unknown, NULL);
	////设置Dst的空间参考系
	//OGRSpatialReference DstSPF;
	//DstSPF.SetProjCS("UTM 45(WGS84) in northern hemisphere.");
	//DstSPF.SetWellKnownGeogCS("WGS84");
	//DstSPF.SetUTM(45, TRUE);
	////创建图层文件，一般为1个图层，图层中添加空间参考与几何类型
	//OGRLayer *poLayer = poDS->CreateLayer("DGGRID_T", &DstSPF, wkbPolygon, NULL);
	//OGRFeature *poFeature;
	//poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
	/////////
	double upper_x, upper_y, down_x, down_y;
	double g0 = geoTrans[0],
		g1 = geoTrans[1],
		g2 = geoTrans[2],//0
		g3 = geoTrans[3],
		g4 = geoTrans[4],//0
		g5 = geoTrans[5];
	for (int i = 0; i < 4; i++)
	{
		double xi = g0 + cols[i] * g1 + rows[i] * g2;
		double yi = g3 + cols[i] * g4 + rows[i] * g5;
		if (i == 0)
		{
			upper_x = xi;
			upper_y = yi;
		}
		if (i == 3)
		{
			down_x = xi;
			down_y = yi;
		}
	}

	////// create the reference frames ////////

	DgRFNetwork net0;
	DgGeoSphRF geoRF(net0, dp.datum, dp.earthRadius);
	DgIDGG dgg(geoRF, dp.vert0, dp.azimuthDegs, dp.aperture, dp.actualRes,
		"DDG", dp.gridTopo, dp.projType, dp.isMixed43, dp.numAp4,
		dp.isSuperfund, dp.sfRes, dp.precision);

	//cout << "Res " << dgg.outputRes() << " " << dgg.gridStats();

	// set-up to convert to degrees
	DgGeoSphDegRF deg(geoRF, geoRF.name() + "Deg");

	// create output files that rely on having the RF's created

	dp.nCellsOutputToFile = 0;
	dp.nOutputFile = 1;

	string cellOutFileName = dp.cellOutFileName;
	string ptOutFileName = dp.ptOutFileName;
	string randPtsOutFileName = dp.randPtsOutFileName;
	if (dp.maxCellsPerFile)
	{
		cellOutFileName += "_1";
		ptOutFileName += "_1";
		randPtsOutFileName += "_1";
	}

	dp.cellOut = DgOutLocFile::makeOutLocFile(dp.cellOutType, cellOutFileName,
		deg, false, dp.precision, dp.shapefileIdLen,
		dp.kmlColor, dp.kmlWidth, dp.kmlName, dp.kmlDescription);

	dp.cellOutShp = NULL;
	if (dp.outCellAttributes)
	{
		dp.cellOutShp = static_cast<DgOutShapefile*>(dp.cellOut);
		dp.cellOutShp->setDefIntAttribute(dp.shapefileDefaultInt);
		dp.cellOutShp->setDefDblAttribute(dp.shapefileDefaultDouble);
		dp.cellOutShp->setDefStrAttribute(dp.shapefileDefaultString);
	}

	dp.ptOut = DgOutLocFile::makeOutLocFile(dp.pointOutType, ptOutFileName,
		deg, true, dp.precision, dp.shapefileIdLen,
		dp.kmlColor, dp.kmlWidth, dp.kmlName, dp.kmlDescription);

	dp.ptOutShp = NULL;
	if (dp.outPointAttributes)
	{
		dp.ptOutShp = static_cast<DgOutShapefile*>(dp.ptOut);
		dp.ptOutShp->setDefIntAttribute(dp.shapefileDefaultInt);
		dp.ptOutShp->setDefDblAttribute(dp.shapefileDefaultDouble);
		dp.ptOutShp->setDefStrAttribute(dp.shapefileDefaultString);
	}

	dp.randPtsOut = NULL;
	if (dp.doRandPts)
	{
		if (dp.curGrid == 1 || !dp.concatPtOut)
		{
			if (!dp.randPtsOutType.compare("TEXT"))
				dp.randPtsOut = new DgOutRandPtsText(deg, randPtsOutFileName,
					dp.precision);
			else
				dp.randPtsOut = DgOutLocFile::makeOutLocFile(dp.randPtsOutType,
					randPtsOutFileName, deg, true, dp.precision, dp.shapefileIdLen,
					dp.kmlColor, dp.kmlWidth, dp.kmlName, dp.kmlDescription);
		}
	}

	////// do whole earth grid if applicable /////

	if (dp.seqToPoly)
	{
		dp.nCellsAccepted = 0;
		dp.nCellsTested = 0;

		set<unsigned long int> seqnums; //To ensure each cell is printed once

										// read-in the sequence numbers
		for (int i = 0; i < dp.regionFiles.size(); i++)
		{
			DgInputStream fin(dp.regionFiles[i].c_str(), "", DgBase::Fatal);
			unsigned long int seqnum;
			const int maxLine = 1000;
			char buff[maxLine];

			while (1) {
				dp.nCellsTested++;

				fin.getline(buff, maxLine);
				if (fin.eof()) break;

				unsigned long int sNum;
				if (sscanf(buff, "%ld", &sNum) != 1)
					::report("doTransform(): invalid SEQNUM " + string(buff), DgBase::Fatal);

				seqnums.insert(sNum);
			}

			fin.close();
		}

		// generate the cells
		for (set<unsigned long int>::iterator i = seqnums.begin(); i != seqnums.end(); i++) {

			DgLocation* loc = static_cast<const DgIDGG&>(dgg).bndRF().locFromSeqNum(*i);
			if (!dgg.bndRF().validLocation(*loc)) {
				std::cerr << "genGrid(): SEQNUM " << (*i) << " is not a valid location" << std::endl;
				::report("genGrid(): Invalid SEQNUM found.", DgBase::Fatal);
			}

			dp.nCellsAccepted++;
			outputStatus(dp);

			DgPolygon verts(dgg);
			dgg.setVertices(*loc, verts, dp.nDensify);

			//outputCellAdd2D(dp, dgg, *loc, verts, deg);

			delete loc;
		}
	}
	else if (dp.wholeEarth)
	{
		dp.nCellsAccepted = 0;
		dp.nCellsTested = 0;
		if (!dp.isSuperfund)
		{
			DgLocation* addLoc = new DgLocation(dgg.bndRF().first());
			ofstream out("tmp_grid.txt");
			while (1)
			{
				dp.nCellsAccepted++;
				dp.nCellsTested++;
				outputStatus(dp);

				DgPolygon verts(dgg);
				dgg.setVertices(*addLoc, verts, dp.nDensify);
				const DgAddress<DgQ2DICoord> *add= (DgAddress<DgQ2DICoord>*)((addLoc)->address());
				//DgLocation tLoc((DgLocation)(*addLoc));
				/*DgQ2DICoord add(tLoc.address());*/
				int Quam = add->address().quadNum();
				// tmp->quadNum();
				//outputCellAdd2D(dp, dgg, *addLoc, verts, deg);
				const DgGeoSphRF* geoRF = dynamic_cast<const DgGeoSphRF*>(&verts.rf());
				//out << "Quam=" <<<< endl;
				for (int i = 0; i < 3; i++)
				{
					const DgGeoCoord& v = *geoRF->getAddress(verts[i]);
					double x = v.lat();
					double y = v.lon();
					x = x * M_180_PI;
					y = y * M_180_PI;
					out << x << "	" << y << endl;

				}
				
				dgg.bndRF().incrementLocation(*addLoc);
				if (!dgg.bndRF().validLocation(*addLoc)) break;
			}
			delete addLoc;
		}
		else // dp.isSuperfund
		{
			for (int q = 0; q < 12; q++)
			{
				DgHexSF baseTile(0, 0, 0, 0, true, q);
				baseTile.setType('P');
				baseTile.depthFirstTraversal(dp, dgg, deg, 2);
			}
		}
	}
	else // use clip regions
	{
		DgQuadClipRegion clipRegions[12]; // clip regions for each quad
		set<DgIVec2D> overageSet[12];     // overage sets
		map<DgIVec2D, set<DgDBFfield> > overageFields[12]; // associated fields

		try {
			createClipRegions(dp, dgg, clipRegions, overageSet, overageFields);
		}
		catch (ClipperLib::clipperException& e) {
			cerr << "ERROR: a clipping polygon vertex exceeds the range for the clipping library.\n";
			report("Try reducing the value of parameter clipper_scale_factor and/or breaking-up large clipping polygons.", DgBase::Fatal);
		}

		if (dp.buildShapeFileAttributes)
		{
			if (dp.outCellAttributes)
				dp.cellOutShp->addFields(dp.allFields);

			if (dp.outPointAttributes)
				dp.ptOutShp->addFields(dp.allFields);
		}

		//// now process the cells by quad ////

		const DgContCartRF& cc1 = dgg.ccFrame();
		const DgDiscRF2D& grid = dgg.grid2D();

		cout << "\n";
		for (int q = 0; q < 12; q++)
		{
			if (overageSet[q].empty() && !clipRegions[q].isQuadUsed())
			{
				cout << string("* No intersections in quad ")
					<< dgg::util::to_string(q) << "." << endl;
				continue;
			}

			cout << string("* Testing quad ") << dgg::util::to_string(q)
				<< "... " << endl;

			if (dp.megaVerbose)
				cout << "Generating: " << q << " " << clipRegions[q].offset()
				<< " " << clipRegions[q].upperRight() << endl;

			DgIVec2D lLeft;
			DgIVec2D uRight;

			if (clipRegions[q].isQuadUsed())
			{
				lLeft = clipRegions[q].offset();
				uRight = clipRegions[q].upperRight();
			}

			// assume dp.isSuperfund
			if (dp.isSuperfund)
			{
				DgEvalData ed(dp, dgg, cc1, grid, clipRegions[q], overageSet[q],
					overageFields[q], deg, lLeft, uRight);

				DgHexSF baseTile(0, 0, 0, 0, true, q);
				baseTile.setType('P');
				baseTile.depthFirstTraversal(dp, dgg, deg, 2, &ed);
			}
			else // !dp.isSuperfund
			{
				DgBoundedRF2D b1(grid, DgIVec2D(0, 0), (uRight - lLeft));
				DgIVec2D tCoord = lLeft; // where are we on the grid?
				while (!overageSet[q].empty() || clipRegions[q].isQuadUsed())
				{
					DgIVec2D coord = tCoord;
					bool accepted = false;

					// first check if there are cells on the overage set

					if (!overageSet[q].empty())
					{
						if (clipRegions[q].isQuadUsed())
						{
							set<DgIVec2D>::iterator it = overageSet[q].find(tCoord);
							if (it != overageSet[q].end()) // found tCoord
							{
								accepted = true;
								overageSet[q].erase(it);
								if (dp.megaVerbose) cout << "found OVERAGE coord " << coord << endl;

								tCoord -= lLeft;
								tCoord = b1.incrementAddress(tCoord);
								if (tCoord == b1.invalidAdd())
									clipRegions[q].setIsQuadUsed(false);

								tCoord += lLeft;
							}
							else
							{
								set<DgIVec2D>::iterator it = overageSet[q].begin();
								if (*it < tCoord)
								{
									accepted = true;
									coord = *it;
									overageSet[q].erase(it);
									if (dp.megaVerbose) cout << "processing OVERAGE " << coord << endl;
								}
								else
								{
									tCoord -= lLeft;
									tCoord = b1.incrementAddress(tCoord);
									if (tCoord == b1.invalidAdd())
										clipRegions[q].setIsQuadUsed(false);

									tCoord += lLeft;
								}
							}
						}
						else
						{
							set<DgIVec2D>::iterator it = overageSet[q].begin();
							coord = *it;
							overageSet[q].erase(it);
							accepted = true;
							if (dp.megaVerbose) cout << "processing OVERAGE " << coord << endl;
						}
					}
					else if (clipRegions[q].isQuadUsed())
					{
						tCoord -= lLeft;
						tCoord = b1.incrementAddress(tCoord);
						if (tCoord == b1.invalidAdd())
							clipRegions[q].setIsQuadUsed(false);

						tCoord += lLeft;
					}

					// skip subfrequency cells as appropriate if doing classII
					// (this should all be done using the seqNum methods, would be
					// much cleaner)

					if (!dgg.isClassI())
						if ((coord.j() + coord.i()) % 3) continue;

					outputStatus(dp);

					if (!accepted)
						accepted = evalCell(dp, dgg, cc1, grid, clipRegions[q], coord);//注意此函数，提供了求交和求差等操作

					if (!accepted) continue;
					//if (accepted) continue;

					// if we're here we have a good one

					dp.nCellsAccepted++;
					//cout << "XX " << q << " " << coord << endl;

					DgLocation* addLoc = dgg.makeLocation(DgQ2DICoord(q, coord));
					DgPolygon verts(dgg);
					dgg.setVertices(*addLoc, verts, dp.nDensify);
					//////
					//根据tif信息，将verts投影至平面
					const DgGeoSphRF* geoRF = dynamic_cast<const DgGeoSphRF*>(&verts.rf());
					vector<Vec2D> gridVerCor = GridVercoord(verts, coordTransInv);
					////判断gridVerCor是否在tif范围内
					//if (!judgeIn(gridVerCor, upper_x, upper_y, down_x, down_y))
					//{
					//	gridVerCor.clear();
					//	continue;
					//}
					//三角形与pixel求交确定val
					getValByPolygonIntesect(gridVerCor, imgData, geoTrans, Width, Height, &gridval);
					////////
					GLC[gridval] = GLC[gridval] + 1;//统计地类个数

					//////
					vector<Vec2D>().swap(gridVerCor);
					//outputCellAdd2D(dp, dgg, *addLoc, verts, deg);
					//////

					//OGRLinearRing ring;
					//for (int i = 0; i < verts.size(); i++)//此处即可得到cell的顶点坐标和中心点坐标
					//							   //	                                 
					//							   //直接使用verts时，所得点区域时不闭合的，此处要注意
					//{
					//	const DgGeoCoord& v = *geoRF->getAddress(verts[i]);
					//	double x = v.lat();
					//	double y = v.lon();
					//	x = x * M_180_PI;
					//	y = y * M_180_PI;
					//	ring.addPoint(x, y);
					//}
					////OGRLinearRing ring;
					//for (int i = 0; i < gridVerCor.size(); i++)
					//{
					//	ring.addPoint(gridVerCor[i].x, gridVerCor[i].y);
					//}
					//	ring.closeRings();
					//	//环加入到polygon中
					//	OGRPolygon polygon;
					//	polygon.addRing(&ring);
					//	//
					//	poFeature->SetGeometry(&polygon);
					//	poLayer->CreateFeature(poFeature);
					// check for special cases 
					if (q == 0 || q == 11) break; // only one cell
				}
			} // else !dp.isSuperfund

			  //cout << "...quad " << q << " complete." << endl;
		}

	} // end if wholeEarth else

	  // close the output files
	delete dp.cellOut;
	dp.cellOut = NULL;
	delete dp.ptOut;
	dp.ptOut = NULL;

	if (dp.numGrids == 1 || !dp.concatPtOut)
	{
		delete dp.randPtsOut;
		dp.randPtsOut = NULL;
	}
	*GLC30 = GLC;
} // void genGrid_

void binValsPartial(BinValsParam& dp, unsigned char * imgData, int Width, int Height, double *geoTrans,
	OGRCoordinateTransformation *coordTrans, OGRCoordinateTransformation *coordTransInv, int **GLC30)
{
	int GLC[19] = { 0 };
	////// create the reference frames ////////
	DgRFNetwork net0;
	DgGeoSphRF geoRF(net0, dp.datum, dp.earthRadius);
	DgIDGG dgg(geoRF, dp.vert0, dp.azimuthDegs, dp.aperture, dp.actualRes,
		"DDG", dp.gridTopo, dp.projType, dp.isMixed43, dp.numAp4,
		dp.isSuperfund, dp.sfRes, dp.precision);

	cout << "Res " << dgg.outputRes() << " " << dgg.gridStats() << endl;

	// set-up to convert to degrees
	DgGeoSphDegRF deg(geoRF, geoRF.name() + "Deg");
	// set-up the output reference frame

	/*dp.outSeqNum = false;
	const DgRFBase* pOutRF = NULL;
	if (dp.outAddType == "PROJTRI") pOutRF = &dgg.projTriRF();
	else if (dp.outAddType == "VERTEX2DD") pOutRF = &dgg.vertexRF();
	else if (dp.outAddType == "Q2DD") pOutRF = &dgg.q2ddRF();
	else if (dp.outAddType == "INTERLEAVE") pOutRF = &dgg.intRF();
	else if (dp.outAddType == "PLANE") pOutRF = &dgg.planeRF();
	else if (dp.outAddType == "Q2DI") pOutRF = &dgg;
	else if (dp.outAddType == "SEQNUM")
	{
		dp.outSeqNum = true;
		pOutRF = &dgg;
	}
	else
	{
		::report("binValsPartial(): invalid output_address_type " +
			dp.outAddType, DgBase::Fatal);
	}

	const DgRFBase& outRF = *pOutRF;*/

	// create a place to store the values by quad

	//QuadVals_ qvals[12];
	QuadVals_ *qvals = (QuadVals_*)malloc(12 * sizeof(QuadVals_));
	memset(qvals, 0, 12 * sizeof(QuadVals_));
	for (int q = 0; q < 12; q++)
	{
		qvals[q].isUsed = false;
		qvals[q].offset = DgIVec2D(dgg.maxI() + 1, dgg.maxJ() + 1);
		qvals[q].upperRight = DgIVec2D(-1, -1);
		qvals[q].numI = 0;
		qvals[q].numJ = 0;
		qvals[q].vals = 0;
	}

	// now make a first pass through the input files and determine what
	// cells are represented

	cout << "determing quad bounds..." << endl;

	const int maxLine = 100;
	char buff[maxLine];
	double lon, lat, val;
	//遍历img
	double g0 = geoTrans[0],
		g1 = geoTrans[1],
		g2 = geoTrans[2],//0
		g3 = geoTrans[3],
		g4 = geoTrans[4],//0
		g5 = geoTrans[5];
	for (int rows = 0; rows < Height; rows++)
	{
		for (int cols = 0; cols < Width; cols++)
		{
			// 计算(rows,cols)的中心点坐标
			double x = g0 + cols * g1 + rows * g2;
			double y = g3 + cols * g4 + rows * g5;
			//（x,y）转换为经纬度
			coordTrans->Transform(1, &x, &y);
			DgLocation* tloc = geoRF.makeLocation(DgGeoCoord(x, y, false));
			dgg.convert(tloc);
			int q = dgg.getAddress(*tloc)->quadNum();
			const DgIVec2D& coord = dgg.getAddress(*tloc)->coord();
			QuadVals_& qv = qvals[q];

			qv.isUsed = true;
			if (coord.i() < qv.offset.i()) qv.offset.setI(coord.i());
			if (coord.i() > qv.upperRight.i()) qv.upperRight.setI(coord.i());
			if (coord.j() < qv.offset.j()) qv.offset.setJ(coord.j());
			if (coord.j() > qv.upperRight.j()) qv.upperRight.setJ(coord.j());

			delete tloc;
		}
	}
	// now initialize the vals storage in the quads which are used
	for (int q = 0; q < 12; q++)
	{
		QuadVals_& qv = qvals[q];
		if (!qv.isUsed) continue;

		qv.upperRight = qv.upperRight - qv.offset; // make relative

		qv.numI = qv.upperRight.i() + 1;
		qv.numJ = qv.upperRight.j() + 1;
		//qv.vals = new Val_*[qv.numI];
		qv.vals = (Val_ **)malloc(qv.numI * sizeof(Val_ *));
		memset(qv.vals, 0, qv.numI * sizeof(Val_ *));
		for (int i = 0; i < qv.numI; i++)
		{
			//qv.vals[i] = new Val_[qv.numJ];
			qv.vals[i] = (Val_ *)malloc(qv.numJ * sizeof(Val_ ));
			memset(qv.vals[i], 0, qv.numJ * sizeof(Val_ ));
			/*for (int j = 0; j < qv.numJ; j++)
			{
				qv.vals[i][j].val = false;
				qv.vals[i][j].index[19] = {0};
			}*/
		}
	}
	// now process the points in each input file
	cout << "binning values..." << endl;
	for (int rows = 0; rows < Height; rows++)
	{
		for (int cols = 0; cols < Width; cols++)
		{
			int index = imgData[rows * Width + cols];
			// 计算(rows,cols)的中心点坐标
			double x = g0 + cols * g1 + rows * g2;
			double y = g3 + cols * g4 + rows * g5;
			//（x,y）转换为经纬度
			coordTrans->Transform(1, &x, &y);
			//确定qij
			DgLocation* tloc = geoRF.makeLocation(DgGeoCoord(x, y, false));
			dgg.convert(tloc);
			int q = dgg.getAddress(*tloc)->quadNum();
			QuadVals_& qv = qvals[q];
			DgIVec2D coord = dgg.getAddress(*tloc)->coord() - qv.offset;
			delete tloc;
			if (index == 254)index = 17;
			if (index == 255) index = 18;
			qv.vals[coord.i()][coord.j()].used = true;
			qv.vals[coord.i()][coord.j()].index[index]++;
		}
	}
	
	///// 找到在Qij中落入栅格像元最多的属性值 /////

	for (int q = 0; q < 12; q++)
	{
		QuadVals_& qv = qvals[q];
		if (!qv.isUsed) continue;

		for (int i = 0; i < qv.numI; i++)
		{
			for (int j = 0; j < qv.numJ; j++)
			{
				int index = -1;
				int max_num = -1;
				if (qv.vals[i][j].used)
				{
					for (int n_index = 0; n_index < 19; n_index++)
					{
						int num = qv.vals[i][j].index[n_index];
						if (num > max_num)
						{
							max_num = num;
							index = n_index;
						}
					}
					qv.vals[i][j].val = index;
					GLC[index]++;
				}
			}
		}
	}

	///// output the values /////

	//if (dp.outputAllCells)
	//{
	//	for (unsigned long int i = 0; i < dgg.bndRF().size(); i++)
	//	{
	//		unsigned long int sNum = i + 1;
	//		DgLocation* tloc = dgg.bndRF().locFromSeqNum(sNum);

	//		double val = 0.0;

	//		// check to see if there is a value for this cell

	//		int q = dgg.getAddress(*tloc)->quadNum();
	//		QuadVals_& qv = qvals[q];
	//		if (qv.isUsed)
	//		{
	//			DgIVec2D coord = dgg.getAddress(*tloc)->coord() - qv.offset;
	//			if (coord.i() >= 0 && coord.j() >= 0 &&
	//				coord.i() <= qv.upperRight.i() &&
	//				coord.j() <= qv.upperRight.j())
	//				val = qv.vals[coord.i()][coord.j()].val;
	//		}

	//		// output the value

	//		if (dp.outSeqNum)
	//			*dp.outFile << sNum << dp.outputDelimiter << val << endl;
	//		else
	//		{
	//			outRF.convert(tloc);
	//			*dp.outFile << tloc->asString(dp.outputDelimiter)
	//				<< dp.outputDelimiter << val << endl;
	//		}

	//		delete tloc;
	//	}
	//}
	//else
	//{
	//	for (int q = 0; q < 12; q++)
	//	{
	//		QuadVals_& qv = qvals[q];
	//		if (!qv.isUsed) continue;

	//		for (int i = 0; i < qv.numI; i++)
	//		{
	//			for (int j = 0; j < qv.numJ; j++)
	//			{
	//				double val = qv.vals[i][j].val;
	//				if (val == 0.0) continue;

	//				DgIVec2D coord(qv.offset.i() + i, qv.offset.j() + j);

	//				DgLocation* tloc = dgg.makeLocation(DgQ2DICoord(q, coord));

	//				if (dp.outSeqNum)
	//				{
	//					unsigned long int sNum = dgg.bndRF().seqNum(*tloc);
	//					*dp.outFile << sNum << dp.outputDelimiter << val << endl;
	//				}
	//				else
	//				{
	//					outRF.convert(tloc);
	//					*dp.outFile << tloc->asString(dp.outputDelimiter)
	//						<< dp.outputDelimiter << val << endl;
	//				}

	//				delete tloc;
	//			}
	//		}
	//	}
	//}

	///// clean-up /////

	for (int q = 0; q < 12; q++)
	{
		QuadVals_& qv = qvals[q];
		if (!qv.isUsed) continue;

		for (int i = 0; i < qv.numI; i++)
		{
			/*delete qv.vals[i];
			qv.vals[i] = NULL;*/
			free( qv.vals[i] );
		}

		/*delete qv.vals;
		qv.vals = NULL;*/
		free(qv.vals);
	}
	*GLC30 = GLC;
}
  //申明 设置DGGRID参数
void setParam_(DgGridPList *plist, int level)//, char * shp_path, char * tif_path
{
	//DgGridPList plist;
	char* token = "dggrid_operation";
	char* remainder = "GENERATE_GRID";
	plist->setParam(token, remainder);
	token = "dggs_type";
	remainder = "ISEA4T";
	//remainder = "ISEA4D";
	//remainder = "ISEA4H";

	plist->setParam(token, remainder);
	//token = "precision";
	//remainder = "6";
	//plist->setParam(token, remainder);
	//token = "rng_type";
	//remainder = "RAND";
	//plist->setParam(token, remainder);
	//token = "verbosity";
	//remainder = "0";
	//plist->setParam(token, remainder);
	////
	//token = "dggs_topology";
	//remainder = "TRIANGLE";
	//plist->setParam(token, remainder);
	////
	//token = "dggs_proj";
	//remainder = "ISEA";
	//plist->setParam(token, remainder);
	////
	//token = "proj_datum";
	//remainder = "WGS84_AUTHALIC_SPHERE";
	//plist->setParam(token, remainder);
	////
	//token = "dggs_orient_specify_type";
	//remainder = "SPECIFIED";
	//plist->setParam(token, remainder);
	////
	token = "dggs_vert0_lat";
	remainder = "90";
	plist->setParam(token, remainder);
	//
	token = "dggs_vert0_lon";
	remainder = "0";
	plist->setParam(token, remainder);
	//
	token = "dggs_vert0_azimuth";
	remainder = "0.0";
	plist->setParam(token, remainder);
	////
	//token = "dggs_num_placements";
	//remainder = "1";
	//plist->setParam(token, remainder);
	////
	//token = "dggs_res_specify_type";
	//remainder = "SPECIFIED";
	//plist->setParam(token, remainder);
	////
	token = "dggs_res_spec";
	char  nlevel[10];
	sprintf(nlevel, "%d", level);
	plist->setParam(token, nlevel);
	//当数据分辨率与格网分辨率相差不大时，即格网单元面积<=2*栅格单元面积
	token = "clip_subset_type";
	remainder = "SHAPEFILE";
	//remainder = "WHOLE_EARTH";
	plist->setParam(token, remainder);
	//当数据分辨率远小于格网分辨率时，可能出现多个栅格单元落入同一格网中
	/*token = "bin_coverage";
	remainder = "PARTIAL";
	plist->setParam(token, remainder);*/
	//
	token = "densification";
	remainder = "0";
	plist->setParam(token, remainder);
	////
	//token = "max_cells_per_output_file";
	//remainder = "0";
	//plist->setParam(token, remainder);
	//
	token = "cell_output_type";
	remainder = "KML";
	plist->setParam(token, remainder);
	//
	token = "cell_output_file_name";
	remainder = "D:\\a";
	plist->setParam(token, remainder);
	//
	token = "shapefile_id_field_length";
	remainder = "5";
	plist->setParam(token, remainder);
}
//实现函数DGGRID_GLC30
void DGGRID_GLC30(DgGridPList &plist, char *tif_path, int **GLC30)
{
	MainParam* pdp = new GridGenParam(plist);
	orientGrid(static_cast<GridGenParam&>(*pdp), plist);
	//根据tif_path获得tif的基本信息
	double *adfGeoTransform;
	GDALDataType gdt;
	GDALDataset *poDataset;//存放tif数据
	int Width, Height, nBands;
	GDALRasterBand ** poBand;
	const char* projRef;
	OGRCoordinateTransformation *coordTrans, *coordTransInv;
	int ret = getImgInfo(tif_path, &poDataset, &nBands, &adfGeoTransform, &Width, &Height, &gdt, &projRef, &poBand, &coordTrans, &coordTransInv);
	unsigned char *imgData = getImgData(poDataset, nBands, Width, Height, &gdt);

	//根据clipregion.shp确定cell

	genGrid_(static_cast<GridGenParam&>(*pdp), imgData, Width, Height, adfGeoTransform,
		coordTrans, coordTransInv, GLC30);
	/*binValsPartial(static_cast<BinValsParam&>(*pdp), imgData, Width, Height, adfGeoTransform,
		coordTrans, coordTransInv, GLC30);*/

	//projRef = NULL;
	//free(adfGeoTransform);
	//free(imgData);
	//delete pdp;
	//delete adfGeoTransform;
	//GDALClose(poDataset);
	//GDALDestroyDriverManager();
}


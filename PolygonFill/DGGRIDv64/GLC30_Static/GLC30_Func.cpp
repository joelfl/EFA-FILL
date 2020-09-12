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
	if (gridval == 0)
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
	return 11;
}

int sig(double d) {
	return(d > 1E-8) - (d < -1E-8);
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
	for (int i = 0; i < n; i++) {
		res += ps[i].x*ps[i + 1].y - ps[i].y*ps[i + 1].x;
	}
	return res / 2.0;
}
double area_(double* ps_x, double *ps_y, int n) {
	ps_x[n] = ps_x[0];
	ps_y[n] = ps_y[0];

	double res = 0;
	for (int i = 0; i < n; i++) {
		res += ps_x[i] * ps_y[i + 1] - ps_y[i] * ps_x[i + 1];
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
	for (int i = 0; i < n; i++) {
		if (sig(cross(a, b, p[i])) > 0) pp[m++] = p[i];
		if (sig(cross(a, b, p[i])) != sig(cross(a, b, p[i + 1])))
			lineCross(a, b, p[i], p[i + 1], pp[m++]);
	}
	n = 0;
	for (int i = 0; i < m; i++)
		if (!i || !(pp[i] == pp[i - 1]))
			p[n++] = pp[i];
	while (n > 1 && p[n - 1] == p[0])n--;
}
void polygon_cut_(double *p_x, double *p_y, int *n, double a_x, double a_y, double b_x, double b_y) {
	static double pp_x[10];
	static double pp_y[10];

	int m = 0;
	int tmp_n = *n;
	p_x[tmp_n] = p_x[0];
	p_y[tmp_n] = p_y[0];

	for (int i = 0; i < tmp_n; i++) {
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
	for (int i = 0; i < m; i++)
		if (!i || !(pp_x[i] == pp_x[i - 1] && pp_y[i] == pp_y[i - 1]))
		{
			p_x[tmp_n] = pp_x[i];
			p_y[tmp_n] = pp_y[i];
			tmp_n++;
		}
	while (tmp_n > 1 && p_x[tmp_n - 1] == p_x[0] && p_y[tmp_n - 1] == p_y[0]) tmp_n--;
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

	if (area(ps1, n1) < 0) reverse(ps1, ps1 + n1);
	if (area(ps2, n2) < 0) reverse(ps2, ps2 + n2);//若面积小于0，则将点顺序倒置
	ps1[n1] = ps1[0];
	ps2[n2] = ps2[0];
	double res = 0;
	for (int i = 0; i < n1; i++) {
		for (int j = 0; j < n2; j++) {
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
	double p_y[10] = { oy,a_y,b_y };
	int n = 3;
	polygon_cut_(p_x, p_y, &n, ox, oy, c_x, c_y);
	polygon_cut_(p_x, p_y, &n, c_x, c_y, d_x, d_y);
	polygon_cut_(p_x, p_y, &n, d_x, d_y, ox, oy);
	res = fabs(area_(p_x, p_y, n));
	if (s1*s2 == -1)
		res = -res;
	return res;
}
double intersectArea__(double *ps1_x, double *ps1_y, int n1, double *ps2_x, double *ps2_y, int n2) {

	if (area_(ps1_x, ps1_y, n1) < 0)
	{
		reverse(ps1_x, ps1_x + n1);
		reverse(ps1_y, ps1_y + n1);
		//double tmp_[4];
		//for (int i = 0; i < n1; i++)
		//{
		//	tmp_[i] = ps1_x[n1 - i - 1];
		//}
		//reverse(ps1_x, ps1_x + n1);
	}
	double tmp_[4];
	if (area_(ps2_x, ps2_y, n1) < 0)
	{
		reverse(ps2_x, ps2_x + n2);
		reverse(ps2_y, ps2_y + n2);
	}
	ps1_x[n1] = ps1_x[0];
	ps1_y[n1] = ps1_y[0];

	ps2_x[n2] = ps2_x[0];
	ps2_y[n2] = ps2_y[0];

	double res = 0;
	for (int i = 0; i < n1; i++) {
		for (int j = 0; j < n2; j++) {
			res += intersectArea_(ps1_x[i], ps1_y[i], ps1_x[i + 1],
				ps1_y[i + 1], ps2_x[j], ps2_y[j], ps2_x[j + 1], ps2_y[j + 1]);
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
	double gird_area[12] = { 0.0 };
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
	//for (k = 0; k < 9; k++)
	//{
	//	if (k == 0)
	//	{
	//		col = center_col;
	//		row = center_row;
	//	}
	//	else if (k == 1)
	//	{
	//		col = center_col - 1;
	//		row = center_row - 1;
	//	}
	//	else if (k == 2)
	//	{
	//		col = center_col - 1;
	//		row = center_row;
	//	}
	//	else if (k == 3)
	//	{
	//		col = center_col - 1;
	//		row = center_row + 1;
	//	}
	//	else if (k == 4)
	//	{
	//		col = center_col;
	//		row = center_row - 1;
	//	}
	//	else if (k == 5)
	//	{
	//		col = center_col;
	//		row = center_row + 1;
	//	}
	//	else if (k == 6)
	//	{
	//		col = center_col + 1;
	//		row = center_row - 1;
	//	}
	//	else if (k == 7)
	//	{
	//		col = center_col + 1;
	//		row = center_row;
	//	}
	//	else if (k == 8)
	//	{
	//		col = center_col + 1;
	//		row = center_row + 1;
	//	}
	//	xi = g0 + col * g1 + row * g2;//像素的左上角坐标
	//	yi = g3 + col * g4 + row * g5;
	//	xj = g0 + (col + 1) * g1 + (row + 1) * g2;//像素的右下角坐标
	//	yj = g3 + (col + 1) * g4 + (row + 1) * g5;

	//	//引入clipper.lib
	//	//ISEA4T
	//	if (total_size == 3)
	//		ai = clipper_intersection_area(gridVerCor[0].x, gridVerCor[0].y,
	//			gridVerCor[1].x, gridVerCor[1].y,
	//			gridVerCor[2].x, gridVerCor[2].y,
	//			xi, yi, xj, yj);
	//	////ISEA4D
	//	/*if (total_size == 4)
	//		ai = clipper_intersection_area_D(gridVerCor[0].x, gridVerCor[0].y,
	//			gridVerCor[1].x, gridVerCor[1].y,
	//			gridVerCor[2].x, gridVerCor[2].y,
	//			gridVerCor[3].x, gridVerCor[3].y,
	//			xi, yi, xj, yj);*/
	//	
	//	////ISEA4H
	//	//if (total_size == 6)
	//	//	ai = clipper_intersection_area_H6(gridVerCor[0].x, gridVerCor[0].y,
	//	//		gridVerCor[1].x, gridVerCor[1].y,
	//	//		gridVerCor[2].x, gridVerCor[2].y,
	//	//		gridVerCor[3].x, gridVerCor[3].y,
	//	//		gridVerCor[4].x, gridVerCor[4].y,
	//	//		gridVerCor[5].x, gridVerCor[5].y,
	//	//		xi, yi, xj, yj);
	//	//if (total_size == 5)
	//	//	ai = clipper_intersection_area_H5(gridVerCor[0].x, gridVerCor[0].y,
	//	//		gridVerCor[1].x, gridVerCor[1].y,
	//	//		gridVerCor[2].x, gridVerCor[2].y,
	//	//		gridVerCor[3].x, gridVerCor[3].y,
	//	//		gridVerCor[4].x, gridVerCor[4].y,
	//	//		xi, yi, xj, yj);
	//	//越界判断
	//	if (col < 0)
	//		col = 0;
	//	if (row < 0)
	//		row = 0;
	//	if (col >= width)
	//		col = width - 1;
	//	if (row >= height)
	//		row = height - 1;
	//	img_k = row * width + col;
	//	*val = imgData[img_k];
	//	glc_class = GetClass(*val);
	//	double ps1_x[5], ps1_y[5];
	//	double ps2_x[5], ps2_y[5];

	//	/*ps1_x[0] = gridVerCor[0].x;
	//	ps1_y[0] = gridVerCor[0].y;
	//	ps1_x[1] = gridVerCor[1].x;
	//	ps1_y[1] = gridVerCor[1].y;
	//	ps1_x[2] = gridVerCor[2].x;
	//	ps1_y[2] = gridVerCor[2].y;
	//	ps1_x[3] = gridVerCor[3].x;
	//	ps1_y[3] = gridVerCor[3].y;
	//	ps1_x[4] = gridVerCor[0].x;
	//	ps1_y[4] = gridVerCor[0].y;

	//	ps2_x[0] = xi;
	//	ps2_y[0] = yi;
	//	ps2_x[1] = xj;
	//	ps2_y[1] = yi;
	//	ps2_x[2] = xj;
	//	ps2_y[2] = yj;
	//	ps2_x[3] = xi;
	//	ps2_y[3] = yj;
	//	ps2_x[4] = xi;
	//	ps2_y[4] = yi;*/
	//	//ai = intersectArea__(ps1_x, ps1_y, 4, ps2_x, ps2_y, 4);
	//	gird_area[glc_class] += ai;

	//	//if (abs(ai - area_new_) > 1)
	//	//{
	//	//	cout << "un equal" << endl;
	//	//}
	//}
	//area = -10;
	//glc_class = -1;
	//for (k = 0; k < 12; k++)
	//{
	//	ai = gird_area[k];
	//	if (ai >= area)
	//	{
	//		area = ai;
	//		glc_class = k;
	//	}
	//}
	//if (area == 0)
	//	glc_class = 11;
	//越界处理
	if (center_col >= width || center_row >= height || center_col < 0 || center_row < 0)
		*val = 0;
	else
	{
		img_k = center_row * width + center_col;
		*val = imgData[img_k];
		glc_class = GetClass(*val);
		*val = glc_class;
	}
	vector<Vec3D>().swap(IJ);
	//return;
}

double clipper_intersection_area(double v1_x, double v1_y, double v2_x, double v2_y,
	double v3_x, double v3_y, OGRLineString *poly)
{
	int SCALE_ = 1000;
	Paths subject, clip, solution;
	ClipType ct = ctIntersection;
	FillRule fr = frNonZero;
	Clipper clipper;
	/************************************************************************/
	ct = ctIntersection;//ctUnion;//
	fr = frEvenOdd;//frNonZero;//
	/************************************************************************/
	subject.resize(1);
	int count = poly->getNumPoints();
	subject[0].resize(count + 1);
	for (int i = 0; i < count; i++)
	{
		subject[0][i].x = (int64_t)(poly->getX(i) * SCALE_);
		subject[0][i].y = (int64_t)(poly->getY(i) * SCALE_);
	}
	subject[0][count].x = (int64_t)(poly->getX(0) * SCALE_);
	subject[0][count].y = (int64_t)(poly->getY(0) * SCALE_);

	clip.resize(1);
	clip[0].resize(4);
	clip[0][0].x = (int64_t)(v1_x * SCALE_);		clip[0][0].y = (int64_t)(v1_y * SCALE_);
	clip[0][1].x = (int64_t)(v2_x * SCALE_);		clip[0][1].y = (int64_t)(v2_y * SCALE_);
	clip[0][2].x = (int64_t)(v3_x * SCALE_);		clip[0][2].y = (int64_t)(v3_y * SCALE_);
	clip[0][3].x = (int64_t)(v1_x * SCALE_);		clip[0][3].y = (int64_t)(v1_y * SCALE_);

	clipper.Clear();
	clipper.AddPaths(subject, ptSubject);
	clipper.AddPaths(clip, ptClip);

	clipper.Execute(ct, solution, fr);//solution没有值 定义可能冲突了

									  /*std::cout << "Result Polygon Size:" << solution.size() << std::endl;
									  for (int i = 0; i<solution.size(); i++)
									  {
									  std::cout << "Polygon: " << i + 1 << std::endl;
									  std::cout << "Point Number: " << solution[i].size() << std::endl;
									  std::cout << "x\ty" << std::endl;
									  for (int j = 0; j<solution[i].size(); j++)
									  {
									  std::cout << solution[i][j].x << "\t" << solution[i][j].y << std::endl;
									  }
									  }*/
									  //cout << solution.size() << endl;
	if (solution.size())
	{
		double re_ai = 0;
		for (int si = 0; si < solution.size(); si++)
			re_ai += Area(solution[si]);
		return (re_ai / SCALE_ / SCALE_);
		//return  ((Area(solution[0])+ Area(solution[1])) / SCALE_ / SCALE_);//此处不对，
	}
	else
	{
		return 0;//表示三角形和pixel不相交
	}
}

double clipper_intersection_area1(double v1_x, double v1_y, double v2_x, double v2_y,
	double v3_x, double v3_y, double v4_x, double v4_y, Paths subject)//OGRLineString *poly
{
	int SCALE_ = 1000;
	Paths  clip, solution;
	ClipType ct = ctIntersection;
	FillRule fr = frNonZero;
	Clipper clipper;
	/************************************************************************/
	ct = ctIntersection;//ctUnion;//
	fr = frEvenOdd;//frNonZero;//
				   /************************************************************************/
	/*subject.resize(1);
	int count = poly->getNumPoints();
	subject[0].resize(count + 1);
	for (int i = 0; i < count; i++)
	{
		subject[0][i].x = (int64_t)(poly->getX(i) * SCALE_);
		subject[0][i].y = (int64_t)(poly->getY(i) * SCALE_);
	}
	subject[0][count].x = (int64_t)(poly->getX(0) * SCALE_);
	subject[0][count].y = (int64_t)(poly->getY(0) * SCALE_);*/

	clip.resize(1);
	clip[0].resize(5);
	clip[0][0].x = (int64_t)(v1_x * SCALE_);		clip[0][0].y = (int64_t)(v1_y * SCALE_);
	clip[0][1].x = (int64_t)(v2_x * SCALE_);		clip[0][1].y = (int64_t)(v2_y * SCALE_);
	clip[0][2].x = (int64_t)(v3_x * SCALE_);		clip[0][2].y = (int64_t)(v3_y * SCALE_);
	clip[0][3].x = (int64_t)(v4_x * SCALE_);		clip[0][3].y = (int64_t)(v4_y * SCALE_);
	clip[0][4].x = (int64_t)(v1_x * SCALE_);		clip[0][4].y = (int64_t)(v1_y * SCALE_);

	clipper.Clear();
	clipper.AddPaths(subject, ptSubject);
	clipper.AddPaths(clip, ptClip);

	clipper.Execute(ct, solution, fr);//solution没有值 定义可能冲突了

									  /*std::cout << "Result Polygon Size:" << solution.size() << std::endl;
									  for (int i = 0; i<solution.size(); i++)
									  {
									  std::cout << "Polygon: " << i + 1 << std::endl;
									  std::cout << "Point Number: " << solution[i].size() << std::endl;
									  std::cout << "x\ty" << std::endl;
									  for (int j = 0; j<solution[i].size(); j++)
									  {
									  std::cout << solution[i][j].x << "\t" << solution[i][j].y << std::endl;
									  }
									  }*/
									  //cout << solution.size() << endl;
	if (solution.size())
	{
		double re_ai = 0;
		for (int si = 0; si < solution.size(); si++)
			re_ai += Area(solution[si]);
		return (re_ai / SCALE_ / SCALE_);
		//return  ((Area(solution[0])+ Area(solution[1])) / SCALE_ / SCALE_);//此处不对，
	}
	else
	{
		return 0;//表示三角形和pixel不相交
	}
}
double clipper_intersection_area2(double v1_x, double v1_y, double v2_x, double v2_y,
	double v3_x, double v3_y, Paths subject)//OGRLineString *poly
{
	int SCALE_ = 1000;
	Paths  clip, solution;
	ClipType ct = ctIntersection;
	FillRule fr = frNonZero;
	Clipper clipper;
	/************************************************************************/
	ct = ctIntersection;//ctUnion;//
	fr = frEvenOdd;//frNonZero;//
				   /************************************************************************/
				   /*subject.resize(1);
				   int count = poly->getNumPoints();
				   subject[0].resize(count + 1);
				   for (int i = 0; i < count; i++)
				   {
				   subject[0][i].x = (int64_t)(poly->getX(i) * SCALE_);
				   subject[0][i].y = (int64_t)(poly->getY(i) * SCALE_);
				   }
				   subject[0][count].x = (int64_t)(poly->getX(0) * SCALE_);
				   subject[0][count].y = (int64_t)(poly->getY(0) * SCALE_);*/

	clip.resize(1);
	clip[0].resize(4);
	clip[0][0].x = (int64_t)(v1_x * SCALE_);		clip[0][0].y = (int64_t)(v1_y * SCALE_);
	clip[0][1].x = (int64_t)(v2_x * SCALE_);		clip[0][1].y = (int64_t)(v2_y * SCALE_);
	clip[0][2].x = (int64_t)(v3_x * SCALE_);		clip[0][2].y = (int64_t)(v3_y * SCALE_);
	clip[0][3].x = (int64_t)(v1_x * SCALE_);		clip[0][3].y = (int64_t)(v1_y * SCALE_);

	clipper.Clear();
	clipper.AddPaths(subject, ptSubject);
	clipper.AddPaths(clip, ptClip);

	clipper.Execute(ct, solution, fr);//solution没有值 定义可能冲突了

									  /*std::cout << "Result Polygon Size:" << solution.size() << std::endl;
									  for (int i = 0; i<solution.size(); i++)
									  {
									  std::cout << "Polygon: " << i + 1 << std::endl;
									  std::cout << "Point Number: " << solution[i].size() << std::endl;
									  std::cout << "x\ty" << std::endl;
									  for (int j = 0; j<solution[i].size(); j++)
									  {
									  std::cout << solution[i][j].x << "\t" << solution[i][j].y << std::endl;
									  }
									  }*/
									  //cout << solution.size() << endl;
	if (solution.size())
	{
		double re_ai = 0;
		for (int si = 0; si < solution.size(); si++)
			re_ai += Area(solution[si]);
		return (re_ai / SCALE_ / SCALE_);
		//return  ((Area(solution[0])+ Area(solution[1])) / SCALE_ / SCALE_);//此处不对，
	}
	else
	{
		return 0;//表示三角形和pixel不相交
	}
}

double evalCell_(GridGenParam& dp, const DgIDGG& dgg, const DgContCartRF& cc1,
	const DgDiscRF2D& grid, vector<Vec2D>gridVerCor, OGRCoordinateTransformation *coordTransInv,
	DgQuadClipRegion& clipRegion,
	const DgIVec2D& add2D)
{
	if (!dgg.isClassI() && (add2D.j() + add2D.i()) % 3)
		return false;

	dp.nCellsTested++;

	if (dp.megaVerbose)
		cout << "Testing #" << dp.nCellsTested << ": " << add2D << endl;

	bool accepted = false;
	double ai = 0.0;
	// start by checking the points
	set<DgIVec2D>::iterator it = clipRegion.points().find(add2D);
	if (it != clipRegion.points().end())
	{
		accepted = true;
		clipRegion.points().erase(it);

		if (dp.buildShapeFileAttributes)
		{
			// add the fields for this point

			map<DgIVec2D, set<DgDBFfield> >::iterator itFields =
				clipRegion.ptFields().find(add2D);
			const set<DgDBFfield>& fields = itFields->second;
			for (set<DgDBFfield>::iterator it = fields.begin();
				it != fields.end(); it++)
				dp.curFields.insert(*it);

			clipRegion.ptFields().erase(itFields);
		}
		else // only need one intersection
			return accepted;
	}

	// generate the boundary

	DgLocation* loc = grid.makeLocation(add2D);

	DgPolygon verts;
	grid.setVertices(add2D, verts);

	cc1.convert(loc);
	DgDVec2D cp = *cc1.getAddress(*loc);
	delete loc;

	cc1.convert(verts);

	// discard cells that don't meet poly-intersect clipping criteria if
	// applicable

	if (dp.doPolyIntersect)
	{
		bool failure = true;

		// check against bounding box
		bool okminx = false;
		bool okmaxx = false;
		bool okminy = false;
		bool okmaxy = false;
		for (int i = 0; i < verts.size(); i++) {
			const DgDVec2D& p0 = *cc1.getAddress((verts)[i]);
			if (!okminx && p0.x() > clipRegion.minx()) okminx = true;
			if (!okminy && p0.y() > clipRegion.miny()) okminy = true;
			if (!okmaxx && p0.x() < clipRegion.maxx()) okmaxx = true;
			if (!okmaxy && p0.y() < clipRegion.maxy()) okmaxy = true;
		}

		if (!(okminx && okmaxx) || !(okminy && okmaxy)) {
			accepted = false;
		}
		else {
			ClipperLib::Paths cellPoly(1);

			for (int i = 0; i < verts.size(); i++)
			{
				//将投影后坐标写进来vector<Vec2D>gridVerCor, 
				cellPoly[0] << ClipperLib::IntPoint(dp.clipperFactor*gridVerCor[i].x, dp.clipperFactor*gridVerCor[i].y);
				//cellPoly[0] << ClipperLib::IntPoint(dp.clipperFactor*cc1.getAddress((verts)[i])->x(), dp.clipperFactor*cc1.getAddress((verts)[i])->y());
			}
			cout << cellPoly[0];
			//vector<ClipperLib::Paths> clip;

			for (unsigned int i = 0; i < clipRegion.clpPolys().size(); i++) //表示一共有多少条path
			{
				ClipperLib::Paths clip_i_path;
				vector<double>X, Y;
				vector<ClipperLib::Paths> clip;
				ClipperLib::Paths clip_i = clipRegion.clpPolys()[i];
				double x0 = clip_i[0][0].X;
				double y0 = clip_i[0][0].Y;
				int j = 1;
				X.push_back(x0);
				Y.push_back(y0);
				cout << x0 << "	" << y0 << endl;
				while (true)
				{
					cout << clip_i[0][j].X << "	" << clip_i[0][j].Y << endl;
					if (clip_i[0][j].X == NULL)
						break;
					X.push_back(clip_i[0][j].X);
					Y.push_back(clip_i[0][j].Y);
					j++;
				}
				clip_i_path.resize(j);
				for (j = 0; j < X.size(); j++)
				{
					x0 = X[0];
					y0 = Y[0];
					coordTransInv->Transform(1, &x0, &y0);
					clip_i_path[i][j].X = (int64_t)(dp.clipperFactor * x0);
					clip_i_path[i][j].Y = (int64_t)(dp.clipperFactor * y0);
				}
				clipRegion.clpPolys()[i] = clip_i_path;
				//clipRegion.clpPolys()[i][0].size();
				ClipperLib::Clipper c;
				c.AddPaths(cellPoly, ClipperLib::ptSubject, true);
				c.AddPaths(clipRegion.clpPolys()[i], ClipperLib::ptClip, true);
				ClipperLib::Paths solution;
				c.Execute(ClipperLib::ctIntersection, solution, ClipperLib::pftNonZero, ClipperLib::pftNonZero);

				if (solution.size() != 0) {
					//accepted = true;
					failure = false;
					ai = Area(solution[0]);
					if (dp.buildShapeFileAttributes) {
						// add the fields for this polygon
						const set<DgDBFfield>& fields = clipRegion.polyFields()[i];
						for (set<DgDBFfield>::iterator it = fields.begin();
							it != fields.end(); it++)
							dp.curFields.insert(*it);
					}
					else { // only need one intersection
						goto EVALCELL_FINISH;
					}
				}
			}
		}

		// If we are here, we did not fail:
		failure = false;

		goto EVALCELL_FINISH;

	EVALCELL_FINISH:

		if (failure)
			throw "Out of memory in evalCell()";
	}

	return ai;

} // bool evalCell

  //快排

void Sort(float a[], int n, int id[], int m)//在网上找到函数 https://zhidao.baidu.com/question/587033720.html
{
	if (m > 1)
	{
		int i = 0;
		int j = m - 1;
		int tmp = id[i];
		while (i < j)
		{
			while (j > i && a[id[j]] > a[tmp])
				--j;
			if (j > i)
				id[i++] = id[j];  //只改变索引顺序
			while (j > i && a[id[i]] < a[tmp])
				++i;
			if (j > i)
				id[j--] = id[i];  //只改变索引顺序
		}
		id[i] = tmp;
		Sort(a, n, id, i);
		Sort(a, n, id + i + 1, m - i - 1);
	}
}
int *Sort_bdgrid(int total, vector<BoundaryG> bd_grid)
{
	int LEN = bd_grid.size();
	float *a = new float[LEN];
	int *id = new int[LEN];
	//将bd_grid中的ai值赋值给a
	for (int i = 0; i < LEN; ++i)
	{
		id[i] = i;
		a[i] = bd_grid[i].area * bd_grid[i].weight_node;//
	}
	Sort(a, LEN, id, LEN);
	delete[]a;
	return id;
}
int *sortdouble(int LEN, float* data)
{
	int *id = new int[LEN];
	for (int i = 0; i < LEN; ++i)
	{
		id[i] = i;
	}
	Sort(data, LEN, id, LEN);
	return id;
}

int *sortint(int LEN, int* data)
{
	int *id = new int[LEN];
	for (int i = 0; i < LEN; ++i)
	{
		id[i] = i;
	}
	Sort((float*)data, LEN, id, LEN);
	return id;
}

bool IS_same_ptgrid(vector<Pt_GRID> Pts_grid, DgQ2DICoord coor)
{
	for (int i = 0; i < Pts_grid.size(); i++)
	{
		DgQ2DICoord coori = Pts_grid[i].coor;
		if (coor == coori) return true;
	}
	return false;
}

//void genGrid_(GridGenParam& dp,
//	OGRSpatialReference oSRS,
//	OGRCoordinateTransformation *coordTrans,
//	OGRCoordinateTransformation *coordTransInv,
//	OGRGeometry* poGeometry,
//	char *outpath,
//	double PolyGon_Area, int type, int PolyID, int &count,
//	vector<BoundaryG> &bd_grid)
//{
//	int GridCount = 0;
//	bool ai_accept = false;
//	double ai = 0.0;
//	vector<BD_GRID> tmp_bd_grid;
//	/////////
//	//输出路径控制
//	GDALDriver *poDriver;
//	poDriver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
//	//创建shp文件
//	GDALDataset *poDS = poDriver->Create(outpath, 0, 0, 0, GDT_Unknown, NULL);
//	//创建图层文件，一般为1个图层，图层中添加空间参考与几何类型
//	OGRLayer *poLayer = poDS->CreateLayer("DGGRID_T", &oSRS, wkbPolygon, NULL);
//	// 创建属性字段
//	OGRFieldDefn oFieldId("Type", OFTInteger);//多边形属性
//	oFieldId.SetWidth(5);
//	OGRFieldDefn firstField("area", OFTReal);//面积
//	firstField.SetWidth(100);
//	OGRFieldDefn secondField("Quam", OFTReal);
//	secondField.SetPrecision(2);
//	OGRFieldDefn thirdField("I", OFTReal);
//	thirdField.SetPrecision(5);
//	OGRFieldDefn forthField("J", OFTReal);
//	forthField.SetPrecision(5);
//	poLayer->CreateField(&oFieldId);
//	poLayer->CreateField(&firstField);
//	poLayer->CreateField(&secondField);
//	poLayer->CreateField(&thirdField);
//	poLayer->CreateField(&forthField);
//
//	OGRFeature *poFeature;
//	poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
//
//	////// create the reference frames ////////
//	DgRFNetwork net0;
//	DgGeoSphRF geoRF(net0, dp.datum, dp.earthRadius);
//	DgIDGG dgg(geoRF, dp.vert0, dp.azimuthDegs, dp.aperture, dp.actualRes,
//		"DDG", dp.gridTopo, dp.projType, dp.isMixed43, dp.numAp4,
//		dp.isSuperfund, dp.sfRes, dp.precision);
//
//	//cout << "Res " << dgg.outputRes() << " " << dgg.gridStats();
//
//	// set-up to convert to degrees
//	DgGeoSphDegRF deg(geoRF, geoRF.name() + "Deg");
//
//	// create output files that rely on having the RF's created
//
//	dp.nCellsOutputToFile = 0;
//	dp.nOutputFile = 1;
//
//	string cellOutFileName = dp.cellOutFileName;
//	string ptOutFileName = dp.ptOutFileName;
//	string randPtsOutFileName = dp.randPtsOutFileName;
//	if (dp.maxCellsPerFile)
//	{
//		cellOutFileName += "_1";
//		ptOutFileName += "_1";
//		randPtsOutFileName += "_1";
//	}
//
//	dp.cellOut = DgOutLocFile::makeOutLocFile(dp.cellOutType, cellOutFileName,
//		deg, false, dp.precision, dp.shapefileIdLen,
//		dp.kmlColor, dp.kmlWidth, dp.kmlName, dp.kmlDescription);
//
//	dp.cellOutShp = NULL;
//	if (dp.outCellAttributes)
//	{
//		dp.cellOutShp = static_cast<DgOutShapefile*>(dp.cellOut);
//		dp.cellOutShp->setDefIntAttribute(dp.shapefileDefaultInt);
//		dp.cellOutShp->setDefDblAttribute(dp.shapefileDefaultDouble);
//		dp.cellOutShp->setDefStrAttribute(dp.shapefileDefaultString);
//	}
//
//	dp.ptOut = DgOutLocFile::makeOutLocFile(dp.pointOutType, ptOutFileName,
//		deg, true, dp.precision, dp.shapefileIdLen,
//		dp.kmlColor, dp.kmlWidth, dp.kmlName, dp.kmlDescription);
//
//	dp.ptOutShp = NULL;
//	if (dp.outPointAttributes)
//	{
//		dp.ptOutShp = static_cast<DgOutShapefile*>(dp.ptOut);
//		dp.ptOutShp->setDefIntAttribute(dp.shapefileDefaultInt);
//		dp.ptOutShp->setDefDblAttribute(dp.shapefileDefaultDouble);
//		dp.ptOutShp->setDefStrAttribute(dp.shapefileDefaultString);
//	}
//
//	dp.randPtsOut = NULL;
//	if (dp.doRandPts)
//	{
//		if (dp.curGrid == 1 || !dp.concatPtOut)
//		{
//			if (!dp.randPtsOutType.compare("TEXT"))
//				dp.randPtsOut = new DgOutRandPtsText(deg, randPtsOutFileName,
//					dp.precision);
//			else
//				dp.randPtsOut = DgOutLocFile::makeOutLocFile(dp.randPtsOutType,
//					randPtsOutFileName, deg, true, dp.precision, dp.shapefileIdLen,
//					dp.kmlColor, dp.kmlWidth, dp.kmlName, dp.kmlDescription);
//		}
//	}
//
//	////// do whole earth grid if applicable /////
//
//	if (dp.seqToPoly)
//	{
//		dp.nCellsAccepted = 0;
//		dp.nCellsTested = 0;
//
//		set<unsigned long int> seqnums; //To ensure each cell is printed once
//
//										// read-in the sequence numbers
//		for (int i = 0; i < dp.regionFiles.size(); i++)
//		{
//			DgInputStream fin(dp.regionFiles[i].c_str(), "", DgBase::Fatal);
//			unsigned long int seqnum;
//			const int maxLine = 1000;
//			char buff[maxLine];
//
//			while (1) {
//				dp.nCellsTested++;
//
//				fin.getline(buff, maxLine);
//				if (fin.eof()) break;
//
//				unsigned long int sNum;
//				if (sscanf(buff, "%ld", &sNum) != 1)
//					::report("doTransform(): invalid SEQNUM " + string(buff), DgBase::Fatal);
//
//				seqnums.insert(sNum);
//			}
//
//			fin.close();
//		}
//
//		// generate the cells
//		for (set<unsigned long int>::iterator i = seqnums.begin(); i != seqnums.end(); i++) {
//
//			DgLocation* loc = static_cast<const DgIDGG&>(dgg).bndRF().locFromSeqNum(*i);
//			if (!dgg.bndRF().validLocation(*loc)) {
//				std::cerr << "genGrid(): SEQNUM " << (*i) << " is not a valid location" << std::endl;
//				::report("genGrid(): Invalid SEQNUM found.", DgBase::Fatal);
//			}
//
//			dp.nCellsAccepted++;
//			outputStatus(dp);
//
//			DgPolygon verts(dgg);
//			dgg.setVertices(*loc, verts, dp.nDensify);
//
//			//outputCellAdd2D(dp, dgg, *loc, verts, deg);
//
//			delete loc;
//		}
//	}
//	else if (dp.wholeEarth)
//	{
//		dp.nCellsAccepted = 0;
//		dp.nCellsTested = 0;
//		if (!dp.isSuperfund)
//		{
//			DgLocation* addLoc = new DgLocation(dgg.bndRF().first());
//			ofstream out("tmp_grid.txt");
//			while (1)
//			{
//				dp.nCellsAccepted++;
//				dp.nCellsTested++;
//				outputStatus(dp);
//
//				DgPolygon verts(dgg);
//				dgg.setVertices(*addLoc, verts, dp.nDensify);
//				const DgAddress<DgQ2DICoord> *add = (DgAddress<DgQ2DICoord>*)((addLoc)->address());
//				//DgLocation tLoc((DgLocation)(*addLoc));
//				/*DgQ2DICoord add(tLoc.address());*/
//				int Quam = add->address().quadNum();
//				// tmp->quadNum();
//				//outputCellAdd2D(dp, dgg, *addLoc, verts, deg);
//				const DgGeoSphRF* geoRF = dynamic_cast<const DgGeoSphRF*>(&verts.rf());
//				//out << "Quam=" <<<< endl;
//				for (int i = 0; i < 3; i++)
//				{
//					const DgGeoCoord& v = *geoRF->getAddress(verts[i]);
//					double x = v.lat();
//					double y = v.lon();
//					x = x * M_180_PI;
//					y = y * M_180_PI;
//					out << x << "	" << y << endl;
//
//				}
//
//				dgg.bndRF().incrementLocation(*addLoc);
//				if (!dgg.bndRF().validLocation(*addLoc)) break;
//			}
//			delete addLoc;
//		}
//		else // dp.isSuperfund
//		{
//			for (int q = 0; q < 12; q++)
//			{
//				DgHexSF baseTile(0, 0, 0, 0, true, q);
//				baseTile.setType('P');
//				baseTile.depthFirstTraversal(dp, dgg, deg, 2);
//			}
//		}
//	}
//	else // use clip regions
//	{
//		DgQuadClipRegion clipRegions[12]; // clip regions for each quad
//		set<DgIVec2D> overageSet[12];     // overage sets
//		map<DgIVec2D, set<DgDBFfield> > overageFields[12]; // associated fields
//
//		try {
//			createClipRegions(dp, dgg, clipRegions, overageSet, overageFields);
//		}
//		catch (ClipperLib::clipperException& e) {
//			cerr << "
//				: a clipping polygon vertex exceeds the range for the clipping library.\n";
//			report("Try reducing the value of parameter clipper_scale_factor and/or breaking-up large clipping polygons.", DgBase::Fatal);
//		}
//
//		if (dp.buildShapeFileAttributes)
//		{
//			if (dp.outCellAttributes)
//				dp.cellOutShp->addFields(dp.allFields);
//
//			if (dp.outPointAttributes)
//				dp.ptOutShp->addFields(dp.allFields);
//		}
//
//		//// now process the cells by quad ////
//
//		const DgContCartRF& cc1 = dgg.ccFrame();
//		const DgDiscRF2D& grid = dgg.grid2D();
//		//首先将结点对应格网写入PolyID的shp中
//		//for (int i = 0; i < Pts_grid.size(); i++)
//		//{
//		//	int pt_polyid = Pts_grid[i].ID;
//		//	if (pt_polyid == PolyID)
//		//	{
//		//		DgQ2DICoord tmp_coor = Pts_grid[i].coor;
//		//		ai = Pts_grid[i].ai;
//		//		//QIJ转经纬度
//		//		DgLocation* addLoc = dgg.makeLocation(DgQ2DICoord(tmp_coor.quadNum(),tmp_coor.coord()));
//		//		DgPolygon verts(dgg);
//		//		dgg.setVertices(*addLoc, verts, dp.nDensify);
//		//		//////
//		//		//根据tif信息，将verts投影至平面
//		//		const DgGeoSphRF* geoRF = dynamic_cast<const DgGeoSphRF*>(&verts.rf());
//		//		vector<Vec2D> gridVerCor = GridVercoord(verts, coordTransInv);
//		//		//
//		//		OGRLinearRing ring;
//		//		for (int i = 0; i < gridVerCor.size(); i++)
//		//			ring.addPoint(gridVerCor[i].x, gridVerCor[i].y);
//		//		ring.closeRings();
//		//		//环加入到polygon中
//		//		OGRPolygon polygon;
//		//		polygon.addRing(&ring);
//		//		//
//		//		poFeature->SetGeometry(&polygon);
//		//		poFeature->SetField(0, type);
//		//		poFeature->SetField(1, ai);
//		//		poLayer->CreateFeature(poFeature);
//		//		GridCount++;
//		//	}
//		//}
//		//将GridCount个数与shp最优格网数对比
//		//if (GridCount >= PolyGon_Area / Grid_Area) return;
//		//其次，按照相交面积对shp进行栅格化
//		cout << "\n";
//		for (int q = 0; q < 12; q++)
//		{
//			if (overageSet[q].empty() && !clipRegions[q].isQuadUsed())
//			{
//				cout << string("* No intersections in quad ")
//					<< dgg::util::to_string(q) << "." << endl;
//				continue;
//			}
//
//			cout << string("* Testing quad ") << dgg::util::to_string(q)
//				<< "... " << endl;
//
//			if (dp.megaVerbose)
//				cout << "Generating: " << q << " " << clipRegions[q].offset()
//				<< " " << clipRegions[q].upperRight() << endl;
//
//			DgIVec2D lLeft;
//			DgIVec2D uRight;
//
//			if (clipRegions[q].isQuadUsed())
//			{
//				lLeft = clipRegions[q].offset();
//				uRight = clipRegions[q].upperRight();
//			}
//
//			// assume dp.isSuperfund
//			if (dp.isSuperfund)
//			{
//				DgEvalData ed(dp, dgg, cc1, grid, clipRegions[q], overageSet[q],
//					overageFields[q], deg, lLeft, uRight);
//
//				DgHexSF baseTile(0, 0, 0, 0, true, q);
//				baseTile.setType('P');
//				baseTile.depthFirstTraversal(dp, dgg, deg, 2, &ed);
//			}
//			else // !dp.isSuperfund
//			{
//				DgBoundedRF2D b1(grid, DgIVec2D(0, 0), (uRight - lLeft));
//				DgIVec2D tCoord = lLeft; // where are we on the grid?
//				while (!overageSet[q].empty() || clipRegions[q].isQuadUsed())
//				{
//					DgIVec2D coord = tCoord;
//					bool accepted = false;
//
//					// first check if there are cells on the overage set
//
//					if (!overageSet[q].empty())
//					{
//						if (clipRegions[q].isQuadUsed())
//						{
//							set<DgIVec2D>::iterator it = overageSet[q].find(tCoord);
//							if (it != overageSet[q].end()) // found tCoord
//							{
//								accepted = true;
//								overageSet[q].erase(it);
//								if (dp.megaVerbose) cout << "found OVERAGE coord " << coord << endl;
//
//								tCoord -= lLeft;
//								tCoord = b1.incrementAddress(tCoord);
//								if (tCoord == b1.invalidAdd())
//									clipRegions[q].setIsQuadUsed(false);
//
//								tCoord += lLeft;
//							}
//							else
//							{
//								set<DgIVec2D>::iterator it = overageSet[q].begin();
//								if (*it < tCoord)
//								{
//									accepted = true;
//									coord = *it;
//									overageSet[q].erase(it);
//									if (dp.megaVerbose) cout << "processing OVERAGE " << coord << endl;
//								}
//								else
//								{
//									tCoord -= lLeft;
//									tCoord = b1.incrementAddress(tCoord);
//									if (tCoord == b1.invalidAdd())
//										clipRegions[q].setIsQuadUsed(false);
//
//									tCoord += lLeft;
//								}
//							}
//						}
//						else
//						{
//							set<DgIVec2D>::iterator it = overageSet[q].begin();
//							coord = *it;
//							overageSet[q].erase(it);
//							accepted = true;
//							if (dp.megaVerbose) cout << "processing OVERAGE " << coord << endl;
//						}
//					}
//					else if (clipRegions[q].isQuadUsed())
//					{
//						tCoord -= lLeft;
//						tCoord = b1.incrementAddress(tCoord);
//						if (tCoord == b1.invalidAdd())
//							clipRegions[q].setIsQuadUsed(false);
//
//						tCoord += lLeft;
//					}
//
//					// skip subfrequency cells as appropriate if doing classII
//					// (this should all be done using the seqNum methods, would be
//					// much cleaner)
//
//					if (!dgg.isClassI())
//						if ((coord.j() + coord.i()) % 3) continue;
//
//					outputStatus(dp);
//
//					if (!accepted)
//					{
//						accepted = evalCell(dp, dgg, cc1, grid, clipRegions[q], coord);
//					}//注意此函数，提供了求交和求差等操作
//
//					if (!accepted) continue;
//
//					// if we're here we have a good one
//
//					dp.nCellsAccepted++;
//					//cout << "XX " << q << " " << coord << endl;
//					////首先判断是否与结点格网重复
//					//bool is_same_ptgrid = IS_same_ptgrid(Pts_grid, DgQ2DICoord(q, coord));
//					//if (is_same_ptgrid) continue;
//					//
//					//计算两个多边形的相交面积
//					DgLocation* addLoc = dgg.makeLocation(DgQ2DICoord(q, coord));
//					DgPolygon verts(dgg);
//					dgg.setVertices(*addLoc, verts, dp.nDensify);
//					//////
//					//根据tif信息，将verts投影至平面
//					const DgGeoSphRF* geoRF = dynamic_cast<const DgGeoSphRF*>(&verts.rf());
//					vector<Vec2D> gridVerCor = GridVercoord(verts, coordTransInv);
//					OGRwkbGeometryType pGeoType = poGeometry->getGeometryType();
//					/*if (coord.i() == 190268 && coord.j() == 515120)
//					{
//						int aaaa = 0;
//						cout << aaaa;
//					}*/
//					//在此处重新修改grid与boundary相交的判断情况
//					//verts是否相交判断
//					int np_1 = 0;
//					bool intersect = 0;
//					for (int np = 0; np < gridVerCor.size(); np++)
//					{
//						double Ax = gridVerCor[np].x;
//						double Ay = gridVerCor[np].y;
//						if (np < 2) np_1 = np + 1;
//						if (np == 2)np_1 = 0;
//						double Bx = gridVerCor[np_1].x;
//						double By = gridVerCor[np_1].y;
//						intersect = Line_BD_intersect(Ax, Ay, Bx, By, poGeometry);
//						if (intersect == 1) break;
//					}
//					//计算面积
//					{
//						if (pGeoType == wkbMultiLineString)
//						{
//							OGRMultiLineString *poMultiLineString = (OGRMultiLineString*)poGeometry;
//							int nGeoCount = poMultiLineString->getNumGeometries();//当前polygon由多少条闭合polyline组成
//							OGRGeometry *poLineGeometry;
//
//							//当iLine==0时，表示外环
//							double *AI = new double[nGeoCount];
//							for (int iLine = 0; iLine < nGeoCount; iLine++)
//							{
//								poLineGeometry = poMultiLineString->getGeometryRef(iLine);
//								OGRLineString *boundary = (OGRLineString*)poLineGeometry;
//								//cout << boundary << endl;
//								ai = clipper_intersection_area(gridVerCor[0].x, gridVerCor[0].y,
//									gridVerCor[1].x, gridVerCor[1].y, gridVerCor[2].x, gridVerCor[2].y,
//									boundary);
//								AI[iLine] = ai;
//								//cout << ai << endl;
//							}
//							if (AI[0] < Area_Threshold)
//							{
//								ai = AI[0];
//							}
//							else
//							{
//								double sumai = 0;
//								for (int iLine = 1; iLine < nGeoCount; iLine++)
//								{
//									sumai += abs(AI[iLine]);
//								}
//								ai = AI[0] - sumai;
//							}
//							delete[] AI;
//						}
//						else
//						{
//							OGRLineString *boundary = (OGRLineString*)poGeometry;
//							ai = clipper_intersection_area(gridVerCor[0].x, gridVerCor[0].y,
//								gridVerCor[1].x, gridVerCor[1].y, gridVerCor[2].x, gridVerCor[2].y,
//								boundary);
//							//cout << ai << endl;
//						}
//					}
//					//
//					if (ai == 0) continue;
//					if(intersect == 0)//不相交
//					{
//						OGRLinearRing ring;
//						for (int i = 0; i < gridVerCor.size(); i++)
//							ring.addPoint(gridVerCor[i].x, gridVerCor[i].y);
//						ring.closeRings();
//						//环加入到polygon中
//						OGRPolygon polygon;
//						polygon.addRing(&ring);
//						//
//						poFeature->SetGeometry(&polygon);
//						poFeature->SetField(0, type);
//						poFeature->SetField(1, ai);
//						poFeature->SetField(2, q);
//						poFeature->SetField(3, coord.i());
//						poFeature->SetField(4, coord.j());
//						poLayer->CreateFeature(poFeature);
//						GridCount++;
//					}
//					else if(intersect == 1)
//					{
//						BoundaryG tmp;
//						tmp.PolyID = PolyID;
//						tmp.area = ai;
//						tmp.coor = DgQ2DICoord(q, coord);
//						bd_grid.push_back(tmp);
//					}
//						
//					
//					//gridVerCor.swap(vector<Vec2D>());
//
//					// check for special cases 
//					if (q == 0 || q == 11) break; // only one cell
//				}
//			} // else !dp.isSuperfund
//
//			  //cout << "...quad " << q << " complete." << endl;
//		}
//
//	} // end if wholeEarth else
//	// close the output files
//	delete dp.cellOut;
//	dp.cellOut = NULL;
//	delete dp.ptOut;
//	dp.ptOut = NULL;
//
//	if (dp.numGrids == 1 || !dp.concatPtOut)
//	{
//		delete dp.randPtsOut;
//		dp.randPtsOut = NULL;
//	}
//	count = GridCount;
//} // void genGrid_

double BD_Grid_AREA(OGRGeometry* poGeometry, vector<Vec2D> gridVerCor)
{
	OGRwkbGeometryType pGeoType = poGeometry->getGeometryType();
	double ai = 0;
	if (pGeoType == wkbMultiLineString)
	{
		OGRMultiLineString *poMultiLineString = (OGRMultiLineString*)poGeometry;
		int nGeoCount = poMultiLineString->getNumGeometries();//当前polygon由多少条闭合polyline组成
		OGRGeometry *poLineGeometry;
		//当iLine==0时，表示外环
		double *AI = new double[nGeoCount];
		for (int iLine = 0; iLine < nGeoCount; iLine++)
		{
			poLineGeometry = poMultiLineString->getGeometryRef(iLine);
			OGRLineString *boundary = (OGRLineString*)poLineGeometry;
			//cout << boundary << endl;
			ai = clipper_intersection_area(gridVerCor[0].x, gridVerCor[0].y,
				gridVerCor[1].x, gridVerCor[1].y, gridVerCor[2].x, gridVerCor[2].y,
				boundary);
			AI[iLine] = ai;
			//cout << ai << endl;
		}
		if (AI[0] < Area_Threshold)
		{
			ai = AI[0];
		}
		else
		{
			for (int iLine = 1; iLine < nGeoCount; iLine++)
			{
				if (abs(AI[iLine]) >= Area_Threshold)
				{
					ai = AI[0] - abs(AI[iLine]);
					break;
				}
				else ai = AI[0] - abs(AI[iLine]);
			}
		}
	}
	else
	{
		OGRLineString *boundary = (OGRLineString*)poGeometry;
		ai = clipper_intersection_area(gridVerCor[0].x, gridVerCor[0].y,
			gridVerCor[1].x, gridVerCor[1].y, gridVerCor[2].x, gridVerCor[2].y,
			boundary);
		//cout << ai << endl;
	}
	return ai;
}

void binValsPartial4RandomPt_Res(BinValsParam& dp, unsigned char * imgData, int Width, int Height, double *geoTrans,
	OGRCoordinateTransformation *coordTrans, OGRCoordinateTransformation *coordTransInv,
	unsigned char * imgData1, int Width1, int Height1, double *geoTrans1, string RandomPtpath, string outpath)
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

	double lon, lat, val;
	double x, y;
	//遍历img
	double g0 = geoTrans[0],
		g1 = geoTrans[1],
		g2 = geoTrans[2],//0
		g3 = geoTrans[3],
		g4 = geoTrans[4],//0
		g5 = geoTrans[5];

	double g01 = geoTrans1[0],
		g11 = geoTrans1[1],
		g21 = geoTrans1[2],//0
		g31 = geoTrans1[3],
		g41 = geoTrans1[4],//0
		g51 = geoTrans1[5];

	////open RandomPt file
	ofstream outfile(outpath.data());

	const int maxLine = 100;
	char line[maxLine];
	DgInputStream inFile(RandomPtpath, "", DgBase::Fatal);
	int gridval_grid, gridval_raster;
	while (inFile.getline(line, maxLine))
	{
		sscanf(line, "%lf	%lf", &lon, &lat);
		//确定栅格属性 利用重采样后的tiff 下标为1的
		x = lon;
		y = lat;
		coordTransInv->Transform(1, &x, &y);//坐标变换重采样前后是一致的，直接用
		int col = (int)((x - g01) / g11);
		int row = (int)((y - g31) / g51);
		//越界处理
		if (col >= Width1 || row >= Height1 || col < 0 || row < 0)
			gridval_raster = 0;
		else
		{
			gridval_raster = imgData1[row * Width1 + col];
			gridval_raster = GetClass(gridval_raster);
		}
		// 确定格网属性
		DgLocation* tloc = geoRF.makeLocation(DgGeoCoord(lon, lat, false));
		dgg.convert(tloc);
		int Q = dgg.getAddress(*tloc)->quadNum();
		const DgIVec2D& COORD = dgg.getAddress(*tloc)->coord();
		for (int rows = row - 16; rows < row + 16; rows++)
		{
			for (int cols = col - 16; cols < col + 16; cols++)
			{
				// 计算(rows,cols)的中心点坐标
				x = g0 + cols * g1 + rows * g2;
				y = g3 + cols * g4 + rows * g5;
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
		//for (int q = 0; q < 12; q++)
		{
			QuadVals_& qv = qvals[Q];
			if (!qv.isUsed) continue;

			qv.upperRight = qv.upperRight - qv.offset; // make relative

			qv.numI = qv.upperRight.i() + 1;
			qv.numJ = qv.upperRight.j() + 1;
			//qv.vals = new Val_*[qv.numI];
			qv.vals = (Val_ **)malloc(qv.numI * sizeof(Val_ *));
			memset(qv.vals, 0, qv.numI * sizeof(Val_ *));
			for (int i = 0; i < qv.numI; i++)
			{
				qv.vals[i] = (Val_ *)malloc(qv.numJ * sizeof(Val_));
				memset(qv.vals[i], 0, qv.numJ * sizeof(Val_));
			}
		}
		// now process the points in each input file
		//cout << "binning values..." << endl;
		for (int rows = row - 16; rows < row + 16; rows++)
		{
			for (int cols = col - 16; cols < col + 16; cols++)
			{
				if (rows < 0 || rows >= Height || cols < 0 || cols >= Width) continue;
				int index = imgData[rows * Width + cols];
				// 计算(rows,cols)的中心点坐标
				x = g0 + cols * g1 + rows * g2;
				y = g3 + cols * g4 + rows * g5;
				//（x,y）转换为经纬度
				coordTrans->Transform(1, &x, &y);
				//确定qij
				DgLocation* tloc = geoRF.makeLocation(DgGeoCoord(x, y, false));
				dgg.convert(tloc);
				int q = dgg.getAddress(*tloc)->quadNum();
				if (q != Q) continue;
				QuadVals_& qv = qvals[q];
				DgIVec2D coord = dgg.getAddress(*tloc)->coord();
				if (coord.i() != COORD.i() || coord.j() != COORD.j()) continue;
				coord = coord - qv.offset;
				delete tloc;
				index = GetClass(index);
				qv.vals[coord.i()][coord.j()].used = true;
				qv.vals[coord.i()][coord.j()].index[index]++;
			}
		}
		///// 找到在Qij中落入栅格像元最多的属性值 /////
		QuadVals_& qv = qvals[Q];
		int max_num = -1, index = -1;
		DgIVec2D coord = COORD - qv.offset;
		for (int n_index = 0; n_index < 12; n_index++)
		{
			int num = qv.vals[coord.i()][coord.j()].index[n_index];
			if (num > max_num)
			{
				max_num = num;
				index = n_index;
			}
		}
		outfile << setprecision(15) << lon << "	" << lat << "	" << index << "	" << gridval_raster << endl;
		for (int q = 0; q < 12; q++)
		{
			QuadVals_& qv = qvals[q];
			if (!qv.isUsed) continue;
			for (int i = 0; i < qv.numI; i++)
			{
				free(qv.vals[i]);
			}
			free(qv.vals);
		}
	}
}
//申明 设置DGGRID参数

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
	remainder = "50.7635";
	plist->setParam(token, remainder);
	//
	token = "dggs_vert0_lon";
	remainder = "154.2628";
	plist->setParam(token, remainder);
	//
	token = "dggs_vert0_azimuth";
	remainder = "226.3062";
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
	//
	token = "clip_subset_type";
	remainder = "SHAPEFILE";
	////remainder = "WHOLE_EARTH";
	plist->setParam(token, remainder);
	//
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
OGRGeometry* ConvertPolygonToPolyline(OGRGeometry* polygon)
{
	// 线生成
	OGRwkbGeometryType sourceGeometryType = polygon->getGeometryType();
	sourceGeometryType = wkbFlatten(sourceGeometryType);

	OGRwkbGeometryType targetGeometryType;
	switch (sourceGeometryType)
	{
	case OGRwkbGeometryType::wkbPolygon:
	{
		OGRPolygon* pOGRPolygon = (OGRPolygon*)polygon;
		int innerCount = pOGRPolygon->getNumInteriorRings();
		if (innerCount == 0)
		{
			targetGeometryType = OGRwkbGeometryType::wkbLineString;
			OGRLineString* pOGRLineString = (OGRLineString*)OGRGeometryFactory::createGeometry(targetGeometryType);

			OGRLinearRing* pOGRLinearRing = pOGRPolygon->getExteriorRing();
			int pointCount = pOGRLinearRing->getNumPoints();
			double x = 0; double y = 0;
			for (int i = 0; i < pointCount; i++)
			{
				x = pOGRLinearRing->getX(i);
				y = pOGRLinearRing->getY(i);
				pOGRLineString->addPoint(x, y);
			}

			return pOGRLineString;
		}
		else
		{
			targetGeometryType = OGRwkbGeometryType::wkbMultiLineString;
			OGRMultiLineString* pOGRMultiLineString = (OGRMultiLineString*)OGRGeometryFactory::createGeometry(targetGeometryType);

			// 添加外环
			OGRLineString ogrLineString;
			OGRLinearRing* pOGRLinearRing = pOGRPolygon->getExteriorRing();
			int pointCount = pOGRLinearRing->getNumPoints();
			double x = 0; double y = 0;
			for (int i = 0; i < pointCount; i++)
			{
				x = pOGRLinearRing->getX(i);
				y = pOGRLinearRing->getY(i);
				ogrLineString.addPoint(x, y);
			}
			pOGRMultiLineString->addGeometry(&ogrLineString);

			for (int i = 0; i < innerCount; i++)
			{
				// 添加内环
				OGRLineString ogrLineString0;
				OGRLinearRing* pOGRLinearRing0 = pOGRPolygon->getInteriorRing(i);
				int pointCount = pOGRLinearRing0->getNumPoints();
				double x = 0; double y = 0;
				for (int i = 0; i < pointCount; i++)
				{
					x = pOGRLinearRing0->getX(i);
					y = pOGRLinearRing0->getY(i);
					ogrLineString0.addPoint(x, y);
				}
				pOGRMultiLineString->addGeometry(&ogrLineString0);
			}

			return pOGRMultiLineString;
		}
	}
	case OGRwkbGeometryType::wkbMultiPolygon:
	{
		targetGeometryType = OGRwkbGeometryType::wkbMultiLineString;
		OGRMultiLineString* pOGRMultiLineString = (OGRMultiLineString*)OGRGeometryFactory::createGeometry(targetGeometryType);

		OGRGeometryCollection* pOGRPolygons = (OGRGeometryCollection*)polygon;
		int geometryCount = pOGRPolygons->getNumGeometries();

		for (int i = 0; i < geometryCount; i++)
		{
			OGRGeometry* pOGRGeo = ConvertPolygonToPolyline(pOGRPolygons->getGeometryRef(i));
			pOGRMultiLineString->addGeometry(pOGRGeo);
		}

		return pOGRMultiLineString;
	}
	default:
		return NULL;
	}

	return NULL;
}

/*
* @brief ConvertPolylineToPolygon        转换线为面
* @param[in] OGRGeometry* polygon        要转换的面
* @return OGRGeometry*            　　　　转换成功后的线
* @author
* @date
* @note 2015年11月04日 小八创建；
*/
OGRGeometry* ConvertPolylineToPolygon(OGRGeometry* polyline)
{
	// 线生成
	OGRwkbGeometryType sourceGeometryType = polyline->getGeometryType();
	sourceGeometryType = wkbFlatten(sourceGeometryType);

	OGRwkbGeometryType targetGeometryType;
	switch (sourceGeometryType)
	{
	case OGRwkbGeometryType::wkbLineString:
	{
		OGRLineString* pOGRLineString = (OGRLineString*)polyline;
		targetGeometryType = OGRwkbGeometryType::wkbPolygon;

		OGRPolygon* pOGRPolygon = (OGRPolygon*)OGRGeometryFactory::createGeometry(targetGeometryType);

		OGRLinearRing pOGRLinearRing;
		int pointCount = pOGRLineString->getNumPoints();
		double x = 0; double y = 0;
		for (int i = 0; i < pointCount; i++)
		{
			x = pOGRLineString->getX(i);
			y = pOGRLineString->getY(i);
			pOGRLinearRing.addPoint(x, y);
		}
		pOGRLinearRing.closeRings();
		pOGRPolygon->addRing(&pOGRLinearRing);
		return pOGRPolygon;
	}
	case OGRwkbGeometryType::wkbMultiLineString:
	{
		targetGeometryType = OGRwkbGeometryType::wkbMultiPolygon;
		OGRMultiPolygon* pOGRMultiPolygon = (OGRMultiPolygon*)OGRGeometryFactory::createGeometry(targetGeometryType);

		OGRGeometryCollection* pOGRPolylines = (OGRGeometryCollection*)polyline;
		int geometryCount = pOGRPolylines->getNumGeometries();

		for (int i = 0; i < geometryCount; i++)
		{
			OGRGeometry* pOGRGeo = ConvertPolylineToPolygon(pOGRPolylines->getGeometryRef(i));
			pOGRMultiPolygon->addGeometry(pOGRGeo);
		}

		return pOGRMultiPolygon;
	}
	default:
		return NULL;
	}

	return NULL;
}


//OGRLineString *SHPBoundary(char *shppath)
//{
//	GDALDataset   *poDS;
//	poDS = (GDALDataset*)GDALOpenEx(shppath, GDAL_OF_VECTOR, NULL, NULL, NULL);
//	OGRLayer  *poLayer;
//	poLayer = poDS->GetLayer(0); //读取层
//	poLayer->ResetReading();
//	OGRFeature *poFeature = poLayer->GetNextFeature();
//	OGRGeometry* poGeometry = ConvertPolygonToPolyline(poFeature->GetGeometryRef());
//	OGRwkbGeometryType pGeoType = poGeometry->getGeometryType();
//	OGRPolygon *rdPolygon = (OGRPolygon *)poGeometry;
//	OGRLineString *boundary;
//	if (pGeoType == wkbMultiLineString)
//	{
//		OGRMultiLineString *poMultiLineString = (OGRMultiLineString*)poGeometry;
//		int nGeoCount = poMultiLineString->getNumGeometries();
//		OGRGeometry *poLineGeometry;
//		for (int iLine = 0; iLine < nGeoCount; iLine++)
//		{
//			poLineGeometry = poMultiLineString->getGeometryRef(iLine);
//			OGRLineString* poLineString = (OGRLineString*)poLineGeometry;
//		}
//	}
//	else
//	{
//		boundary = (OGRLineString*)poGeometry;
//	}
//	return boundary;
//}

//确定边界上结点格网与不同多边形相交面积
void PolygonPTGrid(DgGridPList &plist, char * SHPPath, char * PTtxt, vector<Pt_GRID> &Pts_grid)
{
	MainParam* pdp = new GridGenParam(plist);
	orientGrid(static_cast<GridGenParam&>(*pdp), plist);
	//根据shppath确定boundary
	GDALDataset *poDS = (GDALDataset*)GDALOpenEx(SHPPath, GDAL_OF_VECTOR, NULL, NULL, NULL);
	OGRLayer  *poLayer = poDS->GetLayer(0); //读取层
	OGRFeature *poFeature = poLayer->GetNextFeature();
	OGRGeometry* poGeometry = ConvertPolygonToPolyline(poFeature->GetGeometryRef());
	//获得投影信息
	OGRCoordinateTransformation *coordTrans, *coordTransInv;
	OGRSpatialReference oSRS, oTRS;
	oTRS.SetWellKnownGeogCS("WGS84");
	oSRS = *poLayer->GetSpatialRef();
	/*oSRS.SetProjCS("UTM UTMN (wgs84) in northern hemisphere");
	oSRS.SetWellKnownGeogCS("WGS84");
	oSRS.SetUTM(UTMN, TRUE);*/
	coordTrans = OGRCreateCoordinateTransformation(&oSRS, &oTRS);//from plane to sphere
	coordTransInv = OGRCreateCoordinateTransformation(&oTRS, &oSRS);//from sphere to plane
	//根据结点坐标确定格网属性
	vector<Pt_GRID> tmp_Pts_Grid;
	binValsPartial4RandomPt_(static_cast<BinValsParam&>(*pdp),
		coordTrans, coordTransInv,
		PTtxt, SHPPath, tmp_Pts_Grid);
	Pts_grid = tmp_Pts_Grid;
}

Vec2D Grid_center(vector<Vec2D> gridVerCor)
{
	Vec2D c;
	double sum_x = 0.0, sum_y = 0.0;
	for (int j = 0; j < 3; j++)
	{
		sum_x += gridVerCor[j].x;
		sum_y += gridVerCor[j].y;
	}
	//计算coor中心坐标
	sum_x = sum_x / 3.0;
	sum_y = sum_y / 3.0;
	c.x = sum_x;
	c.y = sum_y;
	return c;
}
void GetRotate(double Ax, double Ay, double Bx, double By, double sint, double cost)
{
	double AB = sqrtf(Bx * Bx + By * By);
	sint = By / AB;
	cost = Bx / AB;
}
bool Line_Line_intersection1(double Ax, double Ay, double sint, double cost, double Px, double Py, double Qx, double Qy, Vec2D *point)
{//与AB所在的直线相交，而非AB线段，并返回交点坐标
	Vec2D Point;
	Point.x = 0;
	Point.y = 0;
	//1 按A点坐标，整体平移至原点
	
	Px = Px - Ax;
	Qx = Qx - Ax;
	Py = Py - Ay;
	Qy = Qy - Ay;

	//2 按AB与x轴正向夹角，整体旋转
	double PxR = Px * cost + Py * sint;
	double PyR = -Px * sint + Py * cost;
	double QxR = Qx * cost + Qy * sint;
	double QyR = -Qx * sint + Qy * cost;
	//3 判断PQ与AB的位置
	if (PyR * QyR > 0) return 0;//PQ位于同侧 不相交
	//4. 若相交 计算交点坐标
	double k = (PxR - QxR) / (PyR - QyR);
	double InterSectX = k * (QxR / k - QyR);
	Point.x = InterSectX * cost + Ax;
	Point.y = InterSectX * sint + Ay;
	*point = Point;
	return 1;//相交
}
vector<Vec2D> GetIntersectPoint(Vec2D A, Vec2D B, OGRGeometry *poGeometry)
{
	//定义变量
	int i_1;
	double Px, Py, Qx, Qy;
	double Ax = A.x, Ay = A.y;
	double Bx = B.x - Ax, By = B.y - Ay;
	Vec2D point;
	vector<Vec2D> POINT;
	OGRwkbGeometryType pGeoType = poGeometry->getGeometryType();
	//计算AB的平移旋转矩阵
	double AB = sqrtf(Bx * Bx + By * By);
	double sint = By / AB;
	double cost = Bx / AB;

	//暂时处理不带环的情况
	if (pGeoType == wkbMultiLineString)
	{ 
	}
	else
	{
		OGRLinearRing *boundary = (OGRLinearRing*)poGeometry;
		int pointcount = boundary->getNumPoints();
		for (int i = 0; i < pointcount - 1; i++)
		{
			i_1 = i + 1;
			Px = boundary->getX(i/* + 1144*/);
			Py = boundary->getY(i /*+ 1144*/);
			Qx = boundary->getX(i_1 /*+ 1144*/);
			Qy = boundary->getY(i_1 /*+ 1144*/);
			if (Line_Line_intersection1(Ax, Ay, sint, cost, Px, Py, Qx, Qy, &point) == 1)
				POINT.push_back(point);
		}
	}
	return POINT;
}
void exchange(double * a, double* b) {
	double temp = *a;
	*a = *b;
	*b = temp;
}
void print_arr(double *a, int size) //打印函数 
{
	cout << "打印数组：";
	for (int i = 0; i<size; i++)  //打印数组 
	{
		cout << a[i] << " ";
	}
	cout << endl << endl;
}
/*序列划分函数*/
int partition(double a[], int p, int r) {
	double key = a[r];//取最后一个
	int i = p - 1;
	for (int j = p; j < r; j++)
	{
		if (a[j] <= key)
		{
			i++;
			//i一直代表小于key元素的最后一个索引，当发现有比key小的a[j]时候，i+1 后交换     
			exchange(&a[i], &a[j]);
		}
	}
	exchange(&a[i + 1], &a[r]);//将key切换到中间来，左边是小于key的，右边是大于key的值。
	return i + 1;
}

void quickSort(double a[], int p, int r) {
	int position = 0;
	if (p<r)
	{
		position = partition(a, p, r);//返回划分元素的最终位置
		quickSort(a, p, position - 1);//划分左边递归
		quickSort(a, position + 1, r);//划分右边递归
	}
}
float *evaluation(vector<Vec2D> POINT)
{
	int NP = POINT.size();
	float *PY = new float[NP];
	for (int i = 0; i < NP; i++)
	{
		PY[i] = POINT[i].y;
	}
	return PY;
}
int *evaluation_(vector<Vec3D> QIJ)
{
	int NP = QIJ.size();
	int *PY = new int[NP];
	for (int i = 0; i < NP; i++)
	{
		PY[i] = QIJ[i].z;
	}
	return PY;
}
Vec2D FindSamePoint(vector<Vec2D> P1, vector<Vec2D> P2)
{
	for (int k = 0; k < 3; k++)
	{
		double x = P1[k].x;
		double y = P1[k].y;
		//cout << setprecision(10) << x << "	" << y << endl;
		for (int j = 0; j < 3; j++)
		{
			//cout << setprecision(10) << P2[j].x << "	" << P2[j].y << endl;
			if (abs(x - P2[j].x) <= 10e-4 && abs(y - P2[j].y) <= 10e-4)
			{
				//cout << setprecision(10) << P2[j].x << "	" << P2[j].y << endl;
				//cout<< endl;
				return P2[j];
			}
		}
	}
}
Paths PreDefineSubject(OGRLineString *poly)
{
	int SCALE_ = 1000;
	Paths subject;
	subject.resize(1);
	int count = poly->getNumPoints();
	subject[0].resize(count + 1);
	for (int i = 0; i < count; i++)
	{
		subject[0][i].x = (int64_t)(poly->getX(i) * SCALE_);
		subject[0][i].y = (int64_t)(poly->getY(i) * SCALE_);
	}
	subject[0][count].x = (int64_t)(poly->getX(0) * SCALE_);
	subject[0][count].y = (int64_t)(poly->getY(0) * SCALE_);
	return subject;
}

int *GetQIJFromstring(string s)
{
	int *QIJ = new int[4];
	int len, pos, tmp, i = 0;
	char t[256];
	t[0] = 0;
	while (1)
	{
		len = s.size();
		pos = s.find_last_of(',');
		string str1 = s.substr(pos + 1, len);
		strncpy_s(t, (char*)s.data(), pos);
		if (pos == -1) break;
		t[pos] = 0;
		s = t;
		//if (i == 1) continue;
		sscanf_s((char*)str1.data(), "%d", &tmp);
		QIJ[3 - i] = tmp;
		//cout << tmp << endl;
		i++;
		if (i >= 3) break;
	}
	sscanf_s(t, "%d", &tmp);
	QIJ[0] = tmp;
	return QIJ;
}
void FillPolygong4Paper2(GridGenParam& dp, char *filename, char *SHP_path, 
							char *txt_path,char *SHP_out_path, char *NGIN_path)
{
	//define the parameters
	GDALDataset   *poDS, *poDS1;
	OGRLayer  *poLayer, *poLayer1;
	OGRFeature *poFeature, *poFeature1;
	OGRGeometry *poGeometry;
	OGRPolygon *BDPolygon;
	OGRLinearRing  grid_ring;
	OGRSpatialReference DstSPF;

	////// create the reference frames ////////
	DgRFNetwork net0;
	DgGeoSphRF geoRF(net0, dp.datum, dp.earthRadius);
	DgIDGG dgg(geoRF, dp.vert0, dp.azimuthDegs, dp.aperture, dp.actualRes,
		"DDG", dp.gridTopo, dp.projType, dp.isMixed43, dp.numAp4,
		dp.isSuperfund, dp.sfRes, dp.precision);

	// set-up to convert to degrees
	DgGeoSphDegRF deg(geoRF, geoRF.name() + "Deg");
	DgLocation* addLoc;
	DgPolygon verts(dgg);
	DgQ2DICoord coormax, coormin, coorI, coorI1, coorI2;
	int I, Imin, Imax;
	int J;
	int polyid;
	const DgAddress<DgQ2DICoord> *add;
	DgLocation ADDLOC;
	vector<DgQ2DICoord> coor;
	vector<Vec2D> gridVerCor1, gridVerCor2, gridVerCor3;
	GDALDriver *poDriver;
	poDriver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
	poDS = (GDALDataset*)GDALOpenEx(SHP_path, GDAL_OF_VECTOR, NULL, NULL, NULL);
	if (poDS == NULL)
	{
		printf("Open failed.\n%s");
		exit;
	}
	poLayer = poDS->GetLayer(0); 
	poLayer->ResetReading();
	DstSPF = *(poLayer->GetSpatialRef());
	OGRSpatialReference oTRS;
	oTRS.SetWellKnownGeogCS("WGS84");
	OGRCoordinateTransformation *coordTrans = OGRCreateCoordinateTransformation(&DstSPF, &oTRS);//from plane to sphere
	OGRCoordinateTransformation *coordTransInv = OGRCreateCoordinateTransformation(&oTRS, &DstSPF);//from sphere to plane
	int polycount = 0;
	ifstream in(filename);
	string line;
	ofstream out(NGIN_path);
	int ID = 0;
	while (getline(in, line))
	{
		sscanf((char*)(line.data()), "%d", &ID);
		//cout << ID << endl;
		//if (ID != 325071) continue;
		while ((poFeature = poLayer->GetNextFeature()) != NULL)
		{
			int polyID = poFeature->GetFieldAsDouble(0);
			if (polyID != ID) continue;
			string GINTXTPath = SHP_out_path + to_string(polyID) + "_" + to_string(NLEVEL) + "_IN.txt";
			string GINShpPath = SHP_out_path + to_string(polyID) + "_" + to_string(NLEVEL) + "_IN.shp";
			string edgetxtPath = txt_path + to_string(polyID) + "_" + to_string(NLEVEL) + ".txt";
			poGeometry = ConvertPolygonToPolyline(poFeature->GetGeometryRef());
			OGRwkbGeometryType pGeoType = poGeometry->getGeometryType();
			vector<Paths> SUBJECTs;
			if (pGeoType == wkbMultiLineString)
			{
				OGRMultiLineString *poMultiLineString = (OGRMultiLineString*)poGeometry;
				int nGeoCount = poMultiLineString->getNumGeometries();//get the count of line of polyon
				OGRGeometry *poLineGeometry;
				for (int iLine = 0; iLine < nGeoCount; iLine++)
				{
					poLineGeometry = poMultiLineString->getGeometryRef(iLine);
					Paths subject = PreDefineSubject((OGRLinearRing*)poLineGeometry);
					SUBJECTs.push_back(subject);
				}
			}
			else
			{
				Paths subject = PreDefineSubject((OGRLinearRing*)poGeometry);
				SUBJECTs.push_back(subject);
			}
			int Nsubjects = SUBJECTs.size();
			{
				OGRFieldDefn oFieldId("Type", OFTInteger);
				oFieldId.SetWidth(5);
				OGRFieldDefn firstField("polyID", OFTInteger);
				firstField.SetWidth(255);
				OGRFieldDefn secondField("Quam", OFTInteger);
				secondField.SetWidth(2);
				OGRFieldDefn thirdField("I", OFTInteger);
				thirdField.SetWidth(255);
				OGRFieldDefn forthField("J", OFTInteger);
				forthField.SetWidth(255);
				poDS1 = poDriver->Create((char*)(GINShpPath.data()), 0, 0, 0, GDT_Unknown, NULL);
				poLayer1 = poDS1->CreateLayer("DGGRID_T", &oTRS, wkbPolygon, NULL);
				poLayer1->CreateField(&oFieldId);
				poLayer1->CreateField(&firstField);
				poLayer1->CreateField(&secondField);
				poLayer1->CreateField(&thirdField);
				poLayer1->CreateField(&forthField);
			}
			poFeature1 = OGRFeature::CreateFeature(poLayer1->GetLayerDefn());

			vector<Vec3D> QIJ;
			string line, line1;
			char name[256], name1[256];
			ifstream in(edgetxtPath);

			int i = 0, TXT_count = 0, Top_left_count = 0;;
			int NGIN = 0;
			Vec3D qij;
			while (in.getline(name, 256)) TXT_count++;
			ifstream in1(edgetxtPath);
			int  Q = 0, _I_ = 0, _J_ = 0, PolyID = 0, Type = 0;
			int originQIJ_I = 0;
			while (i <= TXT_count)
			{
				if (i < TXT_count)
				{
					in1.getline(name, 256);
					sscanf(name, "%d,%d,%d,%d", &polyid, &Q, &_I_, &_J_);
					qij.x = Q; qij.y = _I_; qij.z = _J_;
				}
				i++;
				if (i == 1)
				{
					QIJ.push_back(qij);
					originQIJ_I = QIJ[0].y;
				}
				//cout << name << endl;
				//cout << QIJ[0][1] << endl;
				//cout <<i<<"	"<< qij[2] << "	" << qij[3]<< endl;
				if (qij.y == originQIJ_I && i <= TXT_count)
				{
					if (i != 1)
						QIJ.push_back(qij);
					/*for (int j = 0; j < QIJ.size(); j++)
					{
						cout << QIJ[j].x << "	" << QIJ[j].y << "	" << QIJ[j].z << endl;
					}*/
					continue;
				}
				else
				{
					
					int NP = QIJ.size();
					/*for (int j = 0; j < NP; j++)
					{
					cout << QIJ[j][0] << "	" << QIJ[j][1] << "	" << QIJ[j][2] << endl;
					}*/
					
					//sort by J
					int *PY = evaluation_(QIJ);
					int *index = sortint(NP, PY);

					/*for (int j = 0; j < NP; j++)
					{
						cout << QIJ[j].x << "	" << QIJ[j].y << "	" << QIJ[j].z << endl;
					}*/
					//get the interval
					int quam = QIJ[index[0]].x;
					int i1 = QIJ[index[0]].y;
					int N_Continue = 0;
					vector<Vec3D> QIJ_INTERVAL;
					Vec3D tmp_interval;

					tmp_interval.x = quam;
					tmp_interval.y = i1;
					for (int j = 0; j < NP - 1; j++)
					{
						int j1 = QIJ[index[j]].z;
						int j2 = QIJ[index[j + 1]].z;
						if (abs(j1 - j2) <= 1)
						{
							continue;
						}

						int j12 = (j1 + j2) / 2;
						DgLocation* addLoc = dgg.makeLocation(DgQ2DICoord(quam, DgIVec2D(i1, j12)));
						DgPolygon verts(dgg);
						dgg.setVertices(*addLoc, verts, dp.nDensify);
						gridVerCor1 = GridVercoord(verts, coordTransInv);
						double a0 = clipper_intersection_area2(gridVerCor1[0].x, gridVerCor1[0].y, gridVerCor1[1].x,
							gridVerCor1[1].y, gridVerCor1[2].x, gridVerCor1[2].y, SUBJECTs[0]);
						if (Nsubjects == 1)
						{
							if (a0)
							{
								tmp_interval.z = j1 + 1;
								QIJ_INTERVAL.push_back(tmp_interval);
								tmp_interval.z = j2 - 1;
								QIJ_INTERVAL.push_back(tmp_interval);
								NGIN += j2 - j1 - 1;
								N_Continue = 0;
								for (int jk = j + 1; jk < NP - 1; jk++)
								{
									if (abs(QIJ[index[jk]].z - QIJ[index[jk + 1]].z) <= 1)
										N_Continue++;
									else break;
								}
								N_Continue++;
								Top_left_count += N_Continue;
							}
						}
						else
						{
							bool intersect = 1;
							for (int ss = 1; ss < Nsubjects; ss++)
							{
								double ai = clipper_intersection_area2(gridVerCor1[0].x, gridVerCor1[0].y, gridVerCor1[1].x,
									gridVerCor1[1].y, gridVerCor1[2].x, gridVerCor1[2].y, SUBJECTs[ss]);
								if (ai != 0)
								{
									intersect = 0;
									break;
								}
							}
							if (intersect == 1 && a0)
							{
								tmp_interval.z = j1 + 1;
								QIJ_INTERVAL.push_back(tmp_interval);
								tmp_interval.z = j2 - 1;
								QIJ_INTERVAL.push_back(tmp_interval);
								NGIN += j2 - j1 - 1;
								N_Continue = 0;
								for (int jk = j + 1; jk < NP - 1; jk++)
								{
									if (abs(QIJ[index[jk]].z - QIJ[index[jk + 1]].z) <= 1)
										N_Continue++;
									else break;
								}
								N_Continue++;
								Top_left_count += N_Continue;
							}
						}

					}
					bool LastGrid = 0;
					quam = QIJ[index[NP - 1]].x;
					i1 = QIJ[index[NP - 1]].y;
					int j1 = QIJ[index[NP - 1]].z;
					if ((j1 + 3) > dgg.maxJ()) j1 = dgg.maxJ();
					else
						j1 = j1 + 3;
					DgLocation* addLoc = dgg.makeLocation(DgQ2DICoord(quam, DgIVec2D(i1, j1)));
					DgPolygon verts(dgg);
					dgg.setVertices(*addLoc, verts, dp.nDensify);
					gridVerCor1 = GridVercoord(verts, coordTransInv);
					double a0 = clipper_intersection_area2(gridVerCor1[0].x, gridVerCor1[0].y, gridVerCor1[1].x,
						gridVerCor1[1].y, gridVerCor1[2].x, gridVerCor1[2].y, SUBJECTs[0]);
					if (Nsubjects == 1)
					{
						if (a0)
						{
							LastGrid = 1;
						}
					}
					else
					{
						bool intersect = 1;
						for (int ss = 1; ss < Nsubjects; ss++)
						{
							double ai = clipper_intersection_area2(gridVerCor1[0].x, gridVerCor1[0].y, gridVerCor1[1].x,
								gridVerCor1[1].y, gridVerCor1[2].x, gridVerCor1[2].y, SUBJECTs[ss]);
							if (ai != 0)
							{
								intersect = 0;
								break;
							}
						}
						if (intersect == 1 && a0)
						{
							LastGrid = 1;
						}
					}
					if (LastGrid)
					{
						int sJ = QIJ[index[NP - 1]].z;
						int sI = QIJ[index[NP - 1]].y;
						int sQ = QIJ[index[NP - 1]].x;
						int sJ1 = 2 * (dgg.maxI() - sI) + 1;
						if (sJ - 1 < dgg.maxJ())
						{
							tmp_interval.x = sQ;
							tmp_interval.y = sI;
							tmp_interval.z = sJ + 1;
							QIJ_INTERVAL.push_back(tmp_interval);
							tmp_interval.z = dgg.maxJ();
							QIJ_INTERVAL.push_back(tmp_interval);
							NGIN += dgg.maxJ() - sJ;
						}
						if (sJ - 1 == dgg.maxJ()) NGIN++;
					}
					//cout << NGIN << endl;
					for (int i = 0; i < QIJ_INTERVAL.size(); i = i + 2)
					{
						//
						Vec3D qij1 = QIJ_INTERVAL[i];
						Vec3D qij2 = QIJ_INTERVAL[i + 1];
						//cout << qij1.x << "	" << qij1.y << "	" << qij1.z << "	";
						//cout << qij2.x << "	" << qij2.y << "	" << qij2.z << "	" << endl;

						if (qij1.y == qij2.y)
						{
							int quam = qij1.x;
							int j1 = qij1.z;
							int j2 = qij2.z;
							int I = qij1.y;
							while (j1 <= j2)
							{
								grid_ring.empty();
								DgLocation* addLoc = dgg.makeLocation(DgQ2DICoord(quam, DgIVec2D(I, j1)));
								DgPolygon verts(dgg);
								dgg.setVertices(*addLoc, verts, dp.nDensify);
								const DgGeoSphRF* geoRF = dynamic_cast<const DgGeoSphRF*>(&verts.rf());
								for (int k = 0; k < 3; k++)
								{
									const DgGeoCoord& v = *geoRF->getAddress(verts[k]);
									grid_ring.addPoint(v.lon()* M_180_PI, v.lat() * M_180_PI);
								}
								grid_ring.closeRings();
								OGRPolygon polygon;
								polygon.addRing(&grid_ring);
								poFeature1->SetGeometry(&polygon);
								poFeature1->SetField(1, polyid);
								poFeature1->SetField(2, quam);
								poFeature1->SetField(3, I);
								poFeature1->SetField(4, j1);
								poLayer1->CreateFeature(poFeature1);
								j1++;
							}
						}
						if (qij1.z == qij2.z)
						{
							int quam = qij1.x;
							int i1 = qij1.y;
							int i2 = qij2.y;
							int J = qij1.z;
							while (i1 <= i2)
							{
								grid_ring.empty();
								DgLocation* addLoc = dgg.makeLocation(DgQ2DICoord(quam, DgIVec2D(i1, J)));
								DgPolygon verts(dgg);
								dgg.setVertices(*addLoc, verts, dp.nDensify);
								const DgGeoSphRF* geoRF = dynamic_cast<const DgGeoSphRF*>(&verts.rf());
								for (int k = 0; k < 3; k++)
								{
									const DgGeoCoord& v = *geoRF->getAddress(verts[k]);
									grid_ring.addPoint(v.lon()* M_180_PI, v.lat() * M_180_PI);
								}
								grid_ring.closeRings();
								OGRPolygon polygon;
								polygon.addRing(&grid_ring);
								poFeature1->SetGeometry(&polygon);
								poFeature1->SetField(1, polyid);
								poFeature1->SetField(2, quam);
								poFeature1->SetField(3, i1);
								poFeature1->SetField(4, J);
								poLayer1->CreateFeature(poFeature1);
								i1++;
							}
						}
					}
					vector<Vec3D>().swap(QIJ_INTERVAL);
					vector<Vec3D>().swap(QIJ);
					QIJ.push_back(qij);
					originQIJ_I = qij.y;
				}
			}
			cout << polyID << "," << NGIN << endl;
			out << polyID << "," << NGIN << endl;
			vector<Vec3D>().swap(QIJ);
			OGRFeature::DestroyFeature(poFeature);
			OGRFeature::DestroyFeature(poFeature1);
			GDALClose(poDS1);
			break;
		}
	}
	out.close();
	GDALClose(poDS);
	delete dp.cellOut;
	dp.cellOut = NULL;
	delete dp.ptOut;
	dp.ptOut = NULL;
	if (dp.numGrids == 1 || !dp.concatPtOut)
	{
		delete dp.randPtsOut;
		dp.randPtsOut = NULL;
	}
}

void ScanlineFill(GridGenParam& dp, char *SHP_path, char *SHP_out_path, double x0, double y0, double x3, double y3)
{
	//定义参数
	GDALDataset   *poDS, *poDS1;
	OGRLayer  *poLayer, *poLayer1;
	OGRFeature *poFeature, *poFeature1;
	OGRGeometry *poGeometry;
	OGRwkbGeometryType pGeoType;
	OGRPolygon *BDPolygon;
	OGRLinearRing  grid_ring;
	OGRSpatialReference DstSPF;

	////// create the reference frames ////////
	DgRFNetwork net0;
	DgGeoSphRF geoRF(net0, dp.datum, dp.earthRadius);
	DgIDGG dgg(geoRF, dp.vert0, dp.azimuthDegs, dp.aperture, dp.actualRes,
		"DDG", dp.gridTopo, dp.projType, dp.isMixed43, dp.numAp4,
		dp.isSuperfund, dp.sfRes, dp.precision);

	// set-up to convert to degrees
	DgGeoSphDegRF deg(geoRF, geoRF.name() + "Deg");
	DgLocation* addLoc;
	DgPolygon verts(dgg);
	DgQ2DICoord coormax, coormin, coorI, coorI1, coorI2;
	int I, Imin, Imax;
	int J;
	const DgAddress<DgQ2DICoord> *add;
	DgLocation ADDLOC;
	vector<DgQ2DICoord> coor;
	vector<Vec2D> gridVerCor1, gridVerCor2, gridVerCor3;
	//定义GDAL环境
	GDALDriver *poDriver;
	poDriver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
	
																								   
	 //读取SHP_path的几何轮廓
	poDS = (GDALDataset*)GDALOpenEx(SHP_path, GDAL_OF_VECTOR, NULL, NULL, NULL);

	if (poDS == NULL)
	{
		printf("Open failed.\n%s");
		exit;
	}
	poLayer = poDS->GetLayer(0); //读取层
	poLayer->ResetReading();
	DstSPF = *(poLayer->GetSpatialRef());
	OGRSpatialReference oTRS;
	oTRS.SetWellKnownGeogCS("WGS84");
	OGRCoordinateTransformation *coordTrans = OGRCreateCoordinateTransformation(&DstSPF, &oTRS);//from plane to sphere
	OGRCoordinateTransformation *coordTransInv = OGRCreateCoordinateTransformation(&oTRS, &DstSPF);//from sphere to plane
	poFeature = poLayer->GetNextFeature();
	poGeometry = ConvertPolygonToPolyline(poFeature->GetGeometryRef());


	//创建输出文件
	// 创建属性字段
	{
		OGRFieldDefn oFieldId("Type", OFTInteger);//多边形属性
		oFieldId.SetWidth(5);
		OGRFieldDefn firstField("Quam", OFTInteger);
		firstField.SetWidth(2);
		OGRFieldDefn secondField("I", OFTInteger);
		secondField.SetWidth(255);
		OGRFieldDefn thirdField("J", OFTInteger);
		thirdField.SetWidth(255);

		//创建shp文件
		poDS1 = poDriver->Create(SHP_out_path, 0, 0, 0, GDT_Unknown, NULL);

		//创建图层文件，一般为1个图层，图层中添加空间参考与几何类型
		poLayer1 = poDS1->CreateLayer("DGGRID_T", &oTRS, wkbPolygon, NULL);
		poLayer1->CreateField(&oFieldId);
		poLayer1->CreateField(&firstField);
		poLayer1->CreateField(&secondField);
		poLayer1->CreateField(&thirdField);
	}
	poFeature1 = OGRFeature::CreateFeature(poLayer1->GetLayerDefn());

	//1. 根据(x0,y0)(x3,y3)确定Imin和Imax
	addLoc = geoRF.makeLocation(DgGeoCoord(x0, y0, false));
	dgg.setVertices(*addLoc, verts, dp.nDensify);
	dgg.convert(addLoc);
	ADDLOC = *addLoc;
	add = (DgAddress<DgQ2DICoord>*)((ADDLOC).address());
	coormax = DgQ2DICoord(add->address().quadNum(),DgIVec2D(add->address().coord().i(),add->address().coord().j()));

	addLoc = geoRF.makeLocation(DgGeoCoord(x3, y3, false));
	dgg.setVertices(*addLoc, verts, dp.nDensify);
	dgg.convert(addLoc);
	ADDLOC = *addLoc;
	add = (DgAddress<DgQ2DICoord>*)((ADDLOC).address());
	coormin = DgQ2DICoord(add->address().quadNum(), DgIVec2D(add->address().coord().i(), add->address().coord().j()));
	coorI = coormin;
	Imax = coormax.coord().i();
	Imin = coormin.coord().i();
	//I = Imin;
	I = 3809;
	J = add->address().coord().j();
	//提前将polygon坐标点提取
	Paths subject = PreDefineSubject((OGRLinearRing*)poGeometry);

	while (1)
	{
		cout << I << endl;
		coorI = DgQ2DICoord(add->address().quadNum(), DgIVec2D(I, J));
		coorI1 = DgQ2DICoord(add->address().quadNum(), DgIVec2D(I, J + 2));
		coorI2 = DgQ2DICoord(add->address().quadNum(), DgIVec2D(I, J + 4));
		DgQ2DICoord coorI3 = DgQ2DICoord(add->address().quadNum(), DgIVec2D(I, J + 6));
		I++;
		if (I > 3810) break;
		//1. 找到(I, J)&(I, J+2)的中点
		((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coorI);
		dgg.setVertices(ADDLOC, verts, dp.nDensify);
		gridVerCor1 = GridVercoord(verts, coordTransInv);
		((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coorI1);
		dgg.setVertices(ADDLOC, verts, dp.nDensify);
		gridVerCor2 = GridVercoord(verts, coordTransInv);
		((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coorI2);
		dgg.setVertices(ADDLOC, verts, dp.nDensify);
		gridVerCor3 = GridVercoord(verts, coordTransInv);

		Vec2D P1 = FindSamePoint(gridVerCor1, gridVerCor2);
		Vec2D P2 = FindSamePoint(gridVerCor2, gridVerCor3);

		////

		//((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coorI3);
		//dgg.setVertices(ADDLOC, verts, dp.nDensify);
		//vector<Vec2D> gridVerCor4 = GridVercoord(verts, coordTransInv);
		//Vec2D P3 = FindSamePoint(gridVerCor3, gridVerCor4);

		////


		//2. 计算交点
		vector<Vec2D> POINT = GetIntersectPoint(P1, P2, poGeometry);
		
		int NP = POINT.size();
		if (NP < 2)
		{
			I++;
			continue;
		}
		//3. POINT中坐标按y值从小到大排序
		/*cout << setprecision(5) << P1.x << "	" << P1.y << endl;
		cout << setprecision(5) << P2.x << "	" << P2.y << endl;
		cout << endl;*/
		float *PY = evaluation(POINT);
		int *index = sortdouble(NP, PY);
		/*for (int i = 0; i < NP; i++)
		{
			cout << setprecision(5)<< POINT[i].x << "	" << POINT[i].y << endl;
		}
		cout << endl;
		for (int i = 0; i < NP; i++) 
		{
			cout << setprecision(5) << POINT[index[i]].x << "	" << POINT[index[i]].y << endl;
		}*/
		//4. 确定区间
		vector<Vec2D> INTERVAL;
		for (int i = 0; i < NP - 1; i++)
		{
			Vec2D P1 = POINT[index[i]];
			Vec2D P2 = POINT[index[i + 1]];
			Vec2D midP;
			midP.x = (P1.x + P2.x) / 2.0;
			midP.y = (P1.y + P2.y) / 2.0;
			double midPx1 = midP.x - 1, midPx2 = midP.x + 1, midPy1 = midP.y - 1, midPy2 = midP.y + 1;
			//计算以midP为中心的半径为0.5的正方形的相交情况
			if (clipper_intersection_area1(midPx1, midPy1, midPx2, midPy1, midPx2, midPy2, midPx1, midPy2, subject))//仅考虑无环情况
			{
				INTERVAL.push_back(P1);
				INTERVAL.push_back(P2);
			}
		}
		/*for (int i = 0; i < INTERVAL.size(); i++)
		{
			cout << INTERVAL[i].x << " " << INTERVAL[i].y << endl;
		}
		cout << endl;*/
		//5. 写入SHP文件中输出
		for (int i = 0; i < INTERVAL.size(); i = i + 2)
		{
			double px1 = INTERVAL[i].x;
			double py1 = INTERVAL[i].y;
			double px2 = INTERVAL[i + 1].x;
			double py2 = INTERVAL[i + 1].y;

			//投影至经纬度
			coordTrans->Transform(1, &px1, &py1);
			coordTrans->Transform(1, &px2, &py2);
			//根据经纬度确定格网 coor_i
			addLoc = geoRF.makeLocation(DgGeoCoord(px1, py1, false));
			dgg.setVertices(*addLoc, verts, dp.nDensify);
			dgg.convert(addLoc);
			ADDLOC = *addLoc;
			add = (DgAddress<DgQ2DICoord>*)((ADDLOC).address());
			int Ii = add->address().coord().i();
			int Ji = add->address().coord().j() + 1;
			int Qi = add->address().quadNum();

			//根据经纬度确定格网 coor_i
			addLoc = geoRF.makeLocation(DgGeoCoord(px2, py2, false));
			dgg.setVertices(*addLoc, verts, dp.nDensify);
			dgg.convert(addLoc);
			ADDLOC = *addLoc;
			add = (DgAddress<DgQ2DICoord>*)((ADDLOC).address());
			int Ii_1 = add->address().coord().i();
			int Ji_1 = add->address().coord().j();
			int Qi_1 = add->address().quadNum();
			while (Ji <= Ji_1)
			{
				grid_ring.empty();
				DgLocation* addLoc = dgg.makeLocation(DgQ2DICoord(Qi, DgIVec2D(Ii, Ji)));
				DgPolygon verts(dgg);
				dgg.setVertices(*addLoc, verts, dp.nDensify);
				const DgGeoSphRF* geoRF = dynamic_cast<const DgGeoSphRF*>(&verts.rf());
				for (int k = 0; k < 3; k++)
				{
					const DgGeoCoord& v = *geoRF->getAddress(verts[k]);
					grid_ring.addPoint( v.lon()* M_180_PI, v.lat() * M_180_PI);
				}
				grid_ring.closeRings();
				//环加入到paddRingolygon中
				OGRPolygon polygon;
				polygon.addRing(&grid_ring);
				poFeature1->SetGeometry(&polygon);
				poFeature1->SetField(1, Qi);
				poFeature1->SetField(2, Ii);
				poFeature1->SetField(3, Ji);
				poLayer1->CreateFeature(poFeature1);
				Ji++;
			}
		}
		vector<Vec2D>().swap(INTERVAL);
		vector<Vec2D>().swap(POINT);
	}
}
//实现函数DGGRID_GLC30
void DGGRID_GLC30(DgGridPList &plist, char * SHPPath, char * MBRPath, char * txt_path, char * OutPath, int *NGIN )
{
	MainParam* pdp = new GridGenParam(plist);
	orientGrid(static_cast<GridGenParam&>(*pdp), plist);
	char *filename = "..\\..\\..\\Data\\result\\filename.txt";
	char *SHP_path = "..\\..\\..\\Data\\Shp\\region_proj.shp";
	//the coordinate must be sorted by I value descending order and the same grid must be deleted before running
	char *_txt_path = "..\\..\\..\\Data\\result\\edgegrid_txt\\30_";
	char *SHP_out_path = "..\\..\\..\\Data\\result\\interiorgrid\\30_ ";
	char *NGIN_path = "..\\..\\..\\Data\\result\\interiorgrid_txt\\30_NGIN.txt";
	FillPolygong4Paper2(static_cast<GridGenParam&>(*pdp),
		filename, SHP_path, _txt_path, SHP_out_path, NGIN_path);
}

bool IS_same_pttopo(vector<string> Ptinfo, string str_ori)
{
	auto iter = Ptinfo.begin();
	for (; iter != Ptinfo.end();)
	{
		string str = *iter;
		if (str == str_ori) return true;
	}
	return false;
}

vector<string> Gen_PtTopo(char * Ptpath)
{
	double lon, lat, PolyID, d;
	int a = 0, b = 0, ptstart = 0, ptend = 0;
	int PtID = 0, ptIDend = 0, ptID1st = 0;
	const int maxLine = 100;
	char line1[maxLine];
	DgInputStream inFile1(Ptpath, "", DgBase::Fatal);
	//inFile1.getline(line1, maxLine);
	vector<string> PTinfo;
	while (inFile1.getline(line1, maxLine))
	{
		string tmpString = line1;
		PTinfo.push_back(tmpString);
	}
	//int lines = 0;
	//auto iter1 = PTinfo.begin(), iter2 = PTinfo.begin(), iter3 = PTinfo.begin();
	//string str1, str2, str3;
	//for (; iter1 != PTinfo.end();)
	//{
	//	if (*iter1 == "	")
	//	{
	//		iter1 = PTinfo.erase(iter1);
	//		continue;
	//	}
	//	str1 = *iter1;
	//	sscanf((char*)(str1.data()), "%d,%d,%lf,%lf,%lf",
	//		&a, &ptID1st, &PolyID, &lon, &lat);
	//	//cout << str1 << endl;
	//	//循环检查是否存在与PolyID/lon/lat相同的字符串，若存在，则删除
	//	iter2 = iter1 + 1;
	//	double lon2, lat2, PolyID2;
	//	for (; iter2 != PTinfo.end();)
	//	{
	//		string str2 = *iter2;
	//		if (str2 == "	")
	//		{
	//			iter2++;
	//			continue;
	//		}
	//		//cout << str2 << endl;
	//		sscanf((char*)(str2.data()), "%d,%d,%lf,%lf,%lf",
	//			&a, &ptID1st, &PolyID2, &lon2, &lat2);
	//		if (PolyID == PolyID2 && lon == lon2 && lat == lat2)
	//			*iter2 = "	";
	//		iter2++;
	//	}
	//	iter1++;
	//}
	////测试结果
	//ofstream stream("D:\\data\\SAMPLE\\ptout.txt");
	//for (iter1 = PTinfo.begin(); iter1 != PTinfo.end(); iter1++)
	//{
	//	stream << *iter1 << endl;
	//}
	//stream.close();
	return PTinfo;
}

void binValsPartial4RandomPt_(BinValsParam& dp,
	OGRCoordinateTransformation *coordTrans, OGRCoordinateTransformation *coordTransInv,
	char * PtTXT, char * PolygonSHPPath, vector<Pt_GRID> &Pts_Grid)
{
	////// create the reference frames ////////
	DgRFNetwork net0;
	DgGeoSphRF geoRF(net0, dp.datum, dp.earthRadius);
	DgIDGG dgg(geoRF, dp.vert0, dp.azimuthDegs, dp.aperture, dp.actualRes,
		"DDG", dp.gridTopo, dp.projType, dp.isMixed43, dp.numAp4,
		dp.isSuperfund, dp.sfRes, dp.precision);

	//cout << "Res " << dgg.outputRes() << " " << dgg.gridStats() << endl;

	// set-up to convert to degrees
	DgGeoSphDegRF deg(geoRF, geoRF.name() + "Deg");

	// create a place to store the values by quad
	// now make a first pass through the input files and determine what
	// cells are represented

	cout << "determing quad bounds..." << endl;

	int i, j;

	GDALDataset *poDS = (GDALDataset*)GDALOpenEx(PolygonSHPPath, GDAL_OF_VECTOR, NULL, NULL, NULL);
	//生成pttopo，检查输入的pt文本中是否存在重复的点

	int a = 0, PtID = 0, PtID2 = 0;
	double PolyID = 0, lon = 0, lat = 0;
	double PolyID2 = 0;
	vector<string> PTinfo = Gen_PtTopo(PtTXT);
	//cout << PTinfo.size() << endl;
	for (i = 0; i < PTinfo.size();)
	{
		string tmpString = PTinfo[i];
		//cout << tmpString << endl;

		sscanf((char*)(tmpString.data()), "%d,%d,%lf,%lf,%lf", &a, &PtID, &PolyID, &lon, &lat);
		//LL转QIJ
		// 确定格网属性
		DgLocation* tloc = geoRF.makeLocation(DgGeoCoord(lon, lat, false));
		dgg.convert(tloc);
		DgPolygon verts(dgg);
		dgg.setVertices(*tloc, verts, 0);
		//根据tif信息，将verts投影至平面
		const DgGeoSphRF* geoRF = dynamic_cast<const DgGeoSphRF*>(&verts.rf());
		vector<Vec2D> gridVerCor = GridVercoord(verts, coordTransInv);
		Pt_GRID pt_grid;
		dgg.convert(tloc);
		pt_grid.coor = DgQ2DICoord(dgg.getAddress(*tloc)->quadNum(), dgg.getAddress(*tloc)->coord());
		double ai_max = -1;
		for (j = i; j < PTinfo.size(); j++)
		{
			string tmpString = PTinfo[j];
			sscanf((char*)(tmpString.data()), "%d,%d,%lf,%lf,%lf", &a, &PtID2, &PolyID2, &lon, &lat);
			//cout << tmpString << endl;
			if (PtID != PtID2)
			{
				//pt_grid.ID = -1;
				break;
			}
			else
			{//根据PolyID计算格网相交面积
			 //读取polygon信息
				OGRLayer  *poLayer = poDS->GetLayer(0); //读取层
				OGRFeature *poFeature;/* = poLayer->GetNextFeature();*/
				poLayer->ResetReading();
				while ((poFeature = poLayer->GetNextFeature()) != NULL)
				{
					int oriPolyID = poFeature->GetFieldAsDouble(0);
					if (oriPolyID != PolyID2) continue;
					else
					{
						//提取oriPolyID所对应SHP
						OGRGeometry* poGeometry = ConvertPolygonToPolyline(poFeature->GetGeometryRef());
						double ai = BD_Grid_AREA(poGeometry, gridVerCor);
						if (ai_max < ai)
						{
							ai_max = ai;
							pt_grid.ID = oriPolyID;
							pt_grid.ai = ai_max;
						}
						gridVerCor.swap(vector<Vec2D>());
						/*pt_grid.area.push_back(ai);
						pt_grid.polyID.push_back(PolyID + 1);*/
						break;
					}
				}
			}
		}
		i = j;
		if (pt_grid.ID != -1)
		{
			Pts_Grid.push_back(pt_grid);
			//cout << pt_grid.ID<< "	"<<pt_grid.coor << "	" << pt_grid.ai << endl;
		}
	}
	GDALClose(poDS);
	PTinfo.swap(vector<string>());
}
	
void BDGWight(vector<BoundaryG> &bdg)
{//穷举，统计与coori相同的个数
	for (auto iter = bdg.begin(); iter != bdg.end(); iter++)
	{
		int n = 0;
		DgQ2DICoord coori = (*iter).coor;
		//if (coori.coord().i() == 190360 && coori.coord().j() == 515292)
		//{
		//	int aaaa = 0;
		//	cout << aaaa;
		//}
		for (auto iter1 = bdg.begin(); iter1 != bdg.end(); iter1++)
		{
			DgQ2DICoord coorj = (*iter1).coor;
			if (coori == coorj) n++;
		}
		(*iter).weight_node = n;
	}
}

void NodeBoundaryG(vector<BoundaryG> bd_grid, BoundarGs *&bdg)
{
	//分类结点bd_grid
	BoundarGs* BDM = (BoundarGs*)malloc(10 * sizeof(BoundarGs));
	memset(BDM, 0, 10 * sizeof(BoundarGs));
	vector<BoundaryG> BDmTmp;
	BoundaryG tmp;
	/*tmp.area = -1;
	tmp.coor = DgQ2DICoord(0, (0, 0));
	tmp.PolyID = -1;*/
	//BDmTmp.push_back(tmp);
	//for (int i = 0; i < 10; i++)
	//{//等价于先初始化10个
	//	//BDM[i].bdg.push_back(tmp);
	//}
	//BDmTmp.swap(vector<BoundaryG>());
	int N = bd_grid.size();
	for (int i = 0; i < N; i++)
	{
		BoundaryG BGi = bd_grid[i];
		//if (BGi.coor.coord.i() == 190348 && BGi.coor.coord.j() == 515295)
		//{
		//	int aaa = 0;
		//}
		if (BGi.PolyID == -1) continue;
		DgQ2DICoord coori = BGi.coor;
		BDmTmp.push_back(BGi);
		for (int j = i + 1; j < N; j++)
		{
			BoundaryG BGj = bd_grid[j];
			DgQ2DICoord coorj = BGj.coor;
			if (coori == coorj)
			{
				BDmTmp.push_back(BGj);
				bd_grid[j].PolyID = -1;
			}
		}
		int tmp_total = BDmTmp.size();
		for (int k = 0; k < tmp_total; k++)
		{
			BDM[tmp_total].bdg.push_back(BDmTmp[k]);
		}
		BDmTmp.swap(vector<BoundaryG>());

	}
	ofstream  out;
	out.open("D:\\data\\SAMPLE\\node.txt");
	for (int i= 0; i < 10; i++)
	{
		vector<BoundaryG> t = BDM[i].bdg;
		//cout << i << "	" << t.size() << endl;
		for (auto iter1 = t.begin(); iter1 != t.end(); iter1++)
		{
			//cout<< (*iter1).PolyID << "	" << (*iter1).coor << "	" << (*iter1).area << "	" << (*iter1).IS << endl;
			out << (*iter1).PolyID <<"	"<< (*iter1).coor << "	" << (*iter1).area << "	" << (*iter1).IS << endl;
		}
		t.swap(vector<BoundaryG>());
	}
	out.close();
	bdg = BDM;
}

bool isused(DgQ2DICoord coor, vector<BoundaryG> bdg)
{
	for (auto iter = bdg.begin(); iter != bdg.end(); iter++)
	{
		DgQ2DICoord coori = (*iter).coor;
		if (coori == coor)
		{
			if ((*iter).IS) 
				return true;
		}
	}
	return false;
}

void BDGTraver(vector<BoundaryG> &bdg, int PolyID, int T)
{
	vector<BoundaryG> VT;
	BoundaryG tmp;
	for (auto iter = bdg.begin(); iter != bdg.end();)
	{
		int pid = (*iter).PolyID;
		if (pid != PolyID)
		{
			iter++;
			continue;
		}
		VT.push_back(*iter);
		iter = bdg.erase(iter);
	}
	//面积排序
	int N = VT.size();
	int *Index = Sort_bdgrid(N, VT);
	//for (int i = 0; i < N; i++)
	//{
	//	cout << VT[Index[i]].coor<<"	" << VT[Index[i]].area << VT[Index[i]].weight_node<<"	" << endl;
	//}

	////取前T个，将其bool令为1
	//cout << "取后T个" << endl;
	int tmpsum = 0;
	for (int i = N - 1; i>=0 && tmpsum<T; i--)
	{
		int pos = Index[i];
		DgQ2DICoord coori = VT[pos].coor;
		if (isused(coori, bdg))//当前coori已使用
			continue;
		tmpsum++;
		//cout << VT[Index[i]].coor << "	" << VT[Index[i]].area << VT[Index[i]].weight_node << "	" << endl;

		VT[pos].IS = 1;
	}
	//cout << PolyID<< "	" << tmpsum << "	" << endl;

	//for (int i = 0; i < T; i++)
	//{
	//	int pos = 0;
	//	int key = N - i - 1;
	//	bool used = 0;
	//	while (key != -1)
	//	{
	//		pos = Index[key];
	//		DgQ2DICoord coori = VT[pos].coor;
	//		if (isused(coori, bdg))//当前coori已使用
	//		{
	//			key--;
	//		}
	//		else
	//		{
	//			
	//			VT[pos].IS = 1;
	//			break;
	//		}
	//	}
	//}
	tmpsum = 0;
	for (auto iter = VT.begin(); iter != VT.end(); iter++)
	{
		bdg.push_back(*iter);
		if((*iter).IS) tmpsum++;
	}
	//cout << tmpsum << "	" << endl;
	VT.swap(vector<BoundaryG>());
}

void NodeReDistr(BoundarGs *&bdg, vector<Counting> Counter)
{
	for (int i = 2; i < 10; i++)
	{
		vector<BoundaryG> BDGi = bdg[i].bdg;
		if (BDGi.size() == 0) continue;
		//统计BDGi中bool==1的个数；若个数==sie数，则continue；反之，对等于0的项重分配
		int sum = 0;
		for (auto iter = BDGi.begin(); iter != BDGi.end(); iter++)
			sum += (*iter).IS;
		if (sum == BDGi.size()) continue;
		//每i个bool相交，若值不等于1，则将面积排序，面积最大值bool为1
		for (int j = 0; j < BDGi.size(); j = j + i)
		{
			sum = 0;
			int id = 0, delta_id = 0,sameT = -1;
			int  MaxN = -1, maxDelatN = -1, sameDeltaN = -1;
			for (int k = 0; k < i; k++)
				sum += BDGi[j + k].IS;
			if (sum == i) continue;
			vector<Counting> tmpCount;
			vector<tmpClass1>tmpClass;
			tmpClass1 tmpClass_;
			int tmp = 0;
			int tmpNsum = 0;
			bool is_same = 1;
			for (int k = 0; k < i; k++)
			{
				int polyID = BDGi[j + k].PolyID;
				for (int l = 0; l < Counter.size(); l++)
				{
					Counting cl = Counter[l];
					if (cl.PolyID != polyID) continue;
					else
					{
						int deltaN = abs(cl.current_count - cl.actual_count);
						tmpNsum += deltaN;
						if (MaxN < cl.actual_count)
						{
							MaxN = cl.actual_count;
							id = j + k;
						}
						if (maxDelatN < deltaN)
						{
							maxDelatN = deltaN;
							delta_id = j + k;
						}
						if (deltaN != 0)
						{
							tmpClass_.index = j + k;
							tmpClass_.deltaN = deltaN;
							tmpClass_.currentN = cl.current_count;
							tmpClass_.actualN = cl.actual_count;
							tmpClass.push_back(tmpClass_);
						}
					}
				}
			}
			if(tmpNsum == 0)	BDGi[id].IS = 1;
			//统计deltacount中值相同值出现的次数count，最大次数对应的val
			int count, val;
			frquency(tmpClass, count, val);
			if (count == 1)//每个值仅出现一次
				id = delta_id;
			else
			{
				if (val < maxDelatN) id = delta_id;
				if (val == maxDelatN)
				{
					MaxN = -1;
					for (auto iter = tmpClass.begin(); iter != tmpClass.end(); iter++)
					{
						if ((*iter).deltaN != maxDelatN) continue;
						if (MaxN < (*iter).actualN)
						{
							MaxN = (*iter).actualN;
							id = (*iter).index;
						}
					}
				}
			}
			BDGi[id].IS = 1;
		}
		bdg[i].bdg = BDGi;
		BDGi.swap(vector<BoundaryG>());
	}
	//
	ofstream  out;
	out.open("D:\\data\\SAMPLE\\node.txt");
	for (int i = 0; i < 10; i++)
	{
		vector<BoundaryG> t = bdg[i].bdg;
		//cout << i << "	" << t.size() << endl;
		for (auto iter1 = t.begin(); iter1 != t.end(); iter1++)
		{
			//cout << (*iter1).PolyID << "	" << (*iter1).coor << "	" << (*iter1).area << "	" << (*iter1).IS << endl;
			out << (*iter1).PolyID << "	" << (*iter1).coor << "	" << (*iter1).area << "	" << (*iter1).IS << endl;
		}
	}
	out.close();

}
void frquency(vector<tmpClass1> data, int &count, int &val)
{
	count = 0;
	val = 0;
	int *array = new int[1000];
	int N = data.size();
	int maxCount = -1;
	for (int i = 0; i < N; i++)
	{
		int pos = data[i].deltaN;
		array[pos]++;
	}
	for (int i = 0; i < 1000; i++)
	{
		int tmpD = array[i];
		if (tmpD) 
		{
			if (maxCount < array[i])
			{
				maxCount = array[i];
				count = maxCount;
				val = i;
			}
		}
	}
}
void Coor2SHP(DgGridPList &plist, vector<BoundaryG> bdg, OGRSpatialReference oSRS, vector<Counting> counter)
{
	MainParam* pdp = new GridGenParam(plist);
	orientGrid(static_cast<GridGenParam&>(*pdp), plist);
	DgRFNetwork net0;
	GridGenParam dp = static_cast<GridGenParam&>(*pdp);
	DgGeoSphRF geoRF(net0, dp.datum, dp.earthRadius);
	DgIDGG dgg(geoRF, dp.vert0, dp.azimuthDegs, dp.aperture, dp.actualRes,
		"DDG", dp.gridTopo, dp.projType, dp.isMixed43, dp.numAp4,
		dp.isSuperfund, dp.sfRes, dp.precision);
	DgGeoSphDegRF deg(geoRF, geoRF.name() + "Deg");
	cout << "determing quad bounds..." << endl;

	GDALDriver *poDriver;
	poDriver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
	// 创建属性字段
	OGRFieldDefn oFieldId("Type", OFTInteger);//多边形属性
	oFieldId.SetWidth(5);
	OGRFieldDefn firstField("area", OFTReal);//面积
	firstField.SetWidth(100);
	OGRFieldDefn secondField("Quam", OFTReal);
	secondField.SetPrecision(2);
	OGRFieldDefn thirdField("I", OFTReal);
	thirdField.SetPrecision(5);
	OGRFieldDefn forthField("J", OFTReal);
	forthField.SetPrecision(5);
	OGRFieldDefn fifthField("wight", OFTReal);
	fifthField.SetPrecision(5);
	//获得投影信息
	OGRCoordinateTransformation *coordTrans, *coordTransInv;
	OGRSpatialReference oTRS;
	oTRS.SetWellKnownGeogCS("WGS84");
	coordTrans = OGRCreateCoordinateTransformation(&oSRS, &oTRS);//from plane to sphere
	coordTransInv = OGRCreateCoordinateTransformation(&oTRS, &oSRS);//from sphere to plane

	for (auto iter = counter.begin(); iter != counter.end(); iter++)
	{
		int polyID = (*iter).PolyID;
		//int polyID = 111;

		string path = "D:\\data\\SAMPLE\\OutSHP\\" + to_string(polyID) + "_" + to_string(NLEVEL) + "_EDGE.shp";
		GDALDataset *poDS = poDriver->Create((char *)(path.data()), 0, 0, 0, GDT_Unknown, NULL);
		OGRLayer *poLayer = poDS->CreateLayer("DGGRID_T", &oSRS, wkbPolygon, NULL);
		poLayer->CreateField(&oFieldId);
		poLayer->CreateField(&firstField);
		poLayer->CreateField(&secondField);
		poLayer->CreateField(&thirdField);
		poLayer->CreateField(&forthField);
		poLayer->CreateField(&fifthField);
		OGRFeature *poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
		//double x[23] = { 190343, 190393, 190406, 190404, 190395, 190375, 190375, 190385, 190368, 190368, 190368, 190351, 190357, 190353, 190382, 190380, 190398, 190358, 190358, 190347, 190348, 190360, 190359 };
		//double y[23] = { 515314, 515376, 515370, 515362, 515329, 515369, 515370, 515377, 515306, 515307, 515308, 515281, 515296, 515298, 515278, 515270, 515340, 515288, 515289, 515294, 515295, 515292, 515372 };
		double x[2] = { 190347,190348 };
		double y[2] = { 515294,515295 };
		int sumtmp = 0;
		for (auto iter1 = bdg.begin(); iter1 != bdg.end(); iter1++)
		//for(int i = 0; i< 2;i++)
		{
			int ID = (*iter1).PolyID;
			double ai = (*iter1).area;
			if (ID != polyID) continue;
			if (!(*iter1).IS) continue;
			//if ((*iter1).weight_node < 2) continue;
			sumtmp++;
			DgQ2DICoord coor = (*iter1).coor;
			//DgQ2DICoord coor(1, DgIVec2D(x[i], y[i]));
			DgLocation* tloc = dgg.makeLocation(coor);
			dgg.convert(tloc);
			DgPolygon verts(dgg);
			dgg.setVertices(*tloc, verts, 0);
			//根据tif信息，将verts投影至平面
			const DgGeoSphRF* geoRF = dynamic_cast<const DgGeoSphRF*>(&verts.rf());
			vector<Vec2D> gridVerCor = GridVercoord(verts, coordTransInv);
			OGRLinearRing ring;
			for (int i = 0; i < gridVerCor.size(); i++)
				ring.addPoint(gridVerCor[i].x, gridVerCor[i].y);
			ring.closeRings();
			//环加入到polygon中
			OGRPolygon polygon;
			polygon.addRing(&ring);
			//
			poFeature->SetGeometry(&polygon);
			poFeature->SetField(1, ai);
			poFeature->SetField(2, coor.quadNum());
			poFeature->SetField(3, coor.coord().i());
			poFeature->SetField(4, coor.coord().j());
			poFeature->SetField(5, (*iter1).weight_node);

			poLayer->CreateFeature(poFeature);
			gridVerCor.swap(vector<Vec2D>());
		}
		//cout << sumtmp << endl;
	}
}
void Coor2SHP1(DgGridPList &plist, BoundarGs *bdg, OGRSpatialReference oSRS,  vector<Counting> counter)
{
	MainParam* pdp = new GridGenParam(plist);
	orientGrid(static_cast<GridGenParam&>(*pdp), plist);
	DgRFNetwork net0;
	GridGenParam dp = static_cast<GridGenParam&>(*pdp);
	DgGeoSphRF geoRF(net0, dp.datum, dp.earthRadius);
	DgIDGG dgg(geoRF, dp.vert0, dp.azimuthDegs, dp.aperture, dp.actualRes,
		"DDG", dp.gridTopo, dp.projType, dp.isMixed43, dp.numAp4,
		dp.isSuperfund, dp.sfRes, dp.precision);
	DgGeoSphDegRF deg(geoRF, geoRF.name() + "Deg");
	cout << "determing quad bounds..." << endl;

	GDALDriver *poDriver;
	poDriver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
	// 创建属性字段
	OGRFieldDefn oFieldId("Type", OFTInteger);//多边形属性
	oFieldId.SetWidth(5);
	OGRFieldDefn firstField("area", OFTReal);//面积
	firstField.SetWidth(100);
	OGRFieldDefn secondField("Quam", OFTReal);
	secondField.SetPrecision(2);
	OGRFieldDefn thirdField("I", OFTReal);
	thirdField.SetPrecision(5);
	OGRFieldDefn forthField("J", OFTReal);
	forthField.SetPrecision(5);
	OGRFieldDefn fifthField("wight", OFTReal);
	fifthField.SetPrecision(5);
	//获得投影信息
	OGRCoordinateTransformation *coordTrans, *coordTransInv;
	OGRSpatialReference oTRS;
	oTRS.SetWellKnownGeogCS("WGS84");
	coordTrans = OGRCreateCoordinateTransformation(&oSRS, &oTRS);//from plane to sphere
	coordTransInv = OGRCreateCoordinateTransformation(&oTRS, &oSRS);//from sphere to plane

	for (auto iter = counter.begin(); iter != counter.end(); iter++)
	{
		int polyID = (*iter).PolyID;

		string path = "D:\\data\\SAMPLE\\OutSHP\\" + to_string(polyID) + "_" + to_string(NLEVEL) + "_EDGE.shp";
		GDALDataset *poDS = poDriver->Create((char *)(path.data()), 0, 0, 0, GDT_Unknown, NULL);
		char * wktChar;
		oSRS.exportToWkt(&wktChar);
		cout << wktChar << endl;
		OGRLayer *poLayer = poDS->CreateLayer("DGGRID_T", &oSRS, wkbPolygon, NULL);

		poLayer->CreateField(&oFieldId);
		poLayer->CreateField(&firstField);
		poLayer->CreateField(&secondField);
		poLayer->CreateField(&thirdField);
		poLayer->CreateField(&forthField);
		poLayer->CreateField(&fifthField);

		OGRFeature *poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
		//double x[23] = { 190343, 190393, 190406, 190404, 190395, 190375, 190375, 190385, 190368, 190368, 190368, 190351, 190357, 190353, 190382, 190380, 190398, 190358, 190358, 190347, 190348, 190360, 190359 };
		//double y[23] = { 515314, 515376, 515370, 515362, 515329, 515369, 515370, 515377, 515306, 515307, 515308, 515281, 515296, 515298, 515278, 515270, 515340, 515288, 515289, 515294, 515295, 515292, 515372 };
		double x[8] = { 190298, 190297, 190296, 190295, 190297, 190296, 190295, 190294 };
		double y[8] = { 515165, 515165, 515165, 515165, 515164, 515164, 515164, 515164 };
		int sumtmp = 0;
		for (int i = 1; i < 10; i++)
		{
			vector<BoundaryG>bdgi = bdg[i].bdg;
			for (auto iter1 = bdgi.begin(); iter1 != bdgi.end(); iter1++)
			//for (int i = 0; i < 8; i++)
			{
				int ID = (*iter1).PolyID;
				double ai = (*iter1).area;
				if (ID != polyID) continue;
				if (!(*iter1).IS) continue;
				//if ((*iter1).weight_node < 2) continue;
				sumtmp++;
				DgQ2DICoord coor = (*iter1).coor;
				//DgQ2DICoord coor(1, DgIVec2D(x[i], y[i]));
				DgLocation* tloc = dgg.makeLocation(coor);
				dgg.convert(tloc);
				DgPolygon verts(dgg);
				dgg.setVertices(*tloc, verts, 0);
				//根据tif信息，将verts投影至平面
				const DgGeoSphRF* geoRF = dynamic_cast<const DgGeoSphRF*>(&verts.rf());
				vector<Vec2D> gridVerCor = GridVercoord(verts, coordTransInv);
				OGRLinearRing ring;
				for (int i = 0; i < gridVerCor.size(); i++)
					ring.addPoint(gridVerCor[i].x, gridVerCor[i].y);
				ring.closeRings();
				//环加入到polygon中
				OGRPolygon polygon;
				polygon.addRing(&ring);
				//
				poFeature->SetGeometry(&polygon);
				poFeature->SetField(1, ai);
				poFeature->SetField(2, coor.quadNum());
				poFeature->SetField(3, coor.coord().i());
				poFeature->SetField(4, coor.coord().j());
				poFeature->SetField(5, (*iter1).weight_node);

				poLayer->CreateFeature(poFeature);
				//gridVerCor.swap(vector<Vec2D>());
			}
			//bdgi.swap(vector<BoundaryG>());
		}
		//cout << sumtmp << endl;
	}
}
BoundarGs * GetMissedGrid(BoundarGs *bdg)
{
	vector<BoundaryG> MissedGrid;
	BoundarGs *Missedbdg = (BoundarGs *)malloc(10 * sizeof(BoundarGs));
	memset(Missedbdg, 0, 10 * sizeof(BoundarGs));
	for (int i = 2; i < 10; i++)
	{
		vector<BoundaryG> bdgi = bdg[i].bdg;
		int total = bdgi.size();
		for(int j = 0; j < total; j = j + i)
		{
			int tmpIS = 0;
			for (int k = 0; k < i; k++)
			{
				BoundaryG tmpbdgi = bdgi[j + k];
				tmpIS += tmpbdgi.IS;
			}
			if (tmpIS) continue;
			for (int k = 0; k < i; k++)
			{
				BoundaryG tmpbdgi = bdgi[j + k];
				MissedGrid.push_back(tmpbdgi);
			}
		}
		Missedbdg[i].bdg = MissedGrid;
		MissedGrid.swap(vector<BoundaryG>());
	}
	return Missedbdg;
}
void SpaceDistConst(DgGridPList &plist, char *MBRShpPath, BoundarGs *&bdg, vector<Counting> counter)
{
	//定义变量
	int i = 0, j = 0, k = 0;
	int issum = 0;
	double areasum = 0;
	int polyid = 0;
	vector<BoundaryG> bdgi;
	BoundaryG tmpbdgjk;
	DgQ2DICoord coor;
	for (i = 1; i < 10; i++)//这里的起始是否存在问题？？
		                        //没有问题，因为此处仅处理node>1的情况，即i从2开始
	{
		bdgi = bdg[i].bdg;
		for (j = 0; j < bdgi.size(); j = j+i)
		{
			issum = 0;
			areasum = 0;
			for (k = 0; k < i; k++)
			{
				tmpbdgjk = bdgi[j + k];
				issum += tmpbdgjk.IS;
				areasum += tmpbdgjk.area;
			}
			if (issum) continue;//说明已为格网分配属性
			if (areasum < Grid_Area) continue;//说明格网处于MBR的边界上
			//上面是确定在BoundarGs中的空缺格网
			//找到面积最大的格网，并返回ID
			polyid = MaxCountPolyID(bdgi, counter, i, j, k);//在第i中node形式中，j到j+i格网中，返回对应polygon面积最大的polyid
			cout << polyid << "	" << i << "	" << j << endl;
			//20191130 在此处进行修改
			vector<Location> LOC;
			vector<BoundaryG> polyID_bdg;
			int *Index;
			int DistIndexPOS = 0;
			int NLOC = 0;
			while (1)
			{
				if (polyid == 0) break;
				polyID_bdg = Find_PolyIDBDG(bdg, polyid);
				Index = Dist2MBR(plist, polyID_bdg, MBRShpPath);
				int pos = DistIndexPOS;
				int polyid1 = 0, polyid2 = 0;
				int pos11 = 0, pos12 = 0;
				int pos21 = 0, pos22 = 0;
				DgQ2DICoord coor1;
				polyid1 = polyid;
				while (pos < polyID_bdg.size())
				{
					int ID = Index[pos];
					coor1 = polyID_bdg[ID].coor;
					if (polyID_bdg[ID].weight_node < 2)
					{
						polyid2 = 0;
						break;
					}
					polyid2 = PolyGonChange(polyid1, coor1, bdg, counter);
					Find_PolyIDCoor(bdg, polyid2, coor1, pos21, pos22);
					bool same_loc = Same_LOC(LOC, pos21, pos22);
					cout << polyid1 << "	" << polyid2 << "	" << coor1 << endl;
					if (!same_loc)
					{
						polyid = polyid2;
						DistIndexPOS = pos;
						break;
					}
					else
					{
						pos++;
						continue;
					}
				}
				if (pos == polyID_bdg.size())
				{
					auto iter = LOC.end()-1;
					DistIndexPOS = (*iter).DistIndexPOS;
					polyid = bdg[(*iter).pos1].bdg[(*iter).pos2].PolyID;
					LOC.erase(iter);
					LOC.erase(LOC.end());
					continue;
				}
				Find_PolyIDCoor(bdg, polyid1, coor1, pos11, pos12);
				LOC.push_back({ pos11, pos12, DistIndexPOS + 1});
				if (polyid2 != 0)
				{
					Find_PolyIDCoor(bdg, polyid2, coor1, pos21, pos22);
					LOC.push_back({ pos21, pos22, 0 });
				}
				polyid = polyid2;
				DistIndexPOS = 0;
			}
			bdg[i].bdg[j + k].IS = 1;
			for (k = 0; k < LOC.size(); k++)//此处应为奇数
			{
				int pos1 = LOC[k].pos1;
				int pos2 = LOC[k].pos2;
				if (k % 2 == 0) bdg[pos1].bdg[pos2].IS = 0;
				if (k % 2 == 1) bdg[pos1].bdg[pos2].IS = 1;
			}
			LOC.swap(vector<Location>());
			polyID_bdg.swap(vector<BoundaryG>());
		}
		bdgi.swap(vector<BoundaryG>());
	}
}
int MaxCountPolyID(vector<BoundaryG> bdg, vector<Counting> counter, int i, int j, int &pos)
{
	int pp = 0;
	int maxcount = -1;
	int ID = -1;
	for (int k = 0; k < i; k++)
	{
		BoundaryG tmpbdg = bdg[j + k];
		int polyid = tmpbdg.PolyID;
		int pp = Find_CounterPolyID(counter, polyid);
		if (maxcount < counter[pp].actual_count)
		{
			maxcount = counter[pp].actual_count;
			ID = polyid;
			pos = k;
		}
	}
	return ID;
}
int CandidatePolygon(DgGridPList &plist, int polyID, char *MBRShpPath, BoundarGs *&bdg, vector<Counting> counter, vector<int> &usePolyID, DgQ2DICoord &coord)
{
	int pos1 = 0, pos2 = 0;
	vector<BoundaryG> PolyIDbdg;
	//根据shppath确定boundary
	GDALDataset *poDS = (GDALDataset*)GDALOpenEx(MBRShpPath, GDAL_OF_VECTOR, NULL, NULL, NULL);
	OGRLayer  *poLayer = poDS->GetLayer(0); //读取层
	OGRFeature *poFeature = poLayer->GetNextFeature();
	OGRGeometry* poGeometry = poFeature->GetGeometryRef();
	OGREnvelope *poEnvelope = new OGREnvelope;
	poGeometry->getEnvelope(poEnvelope);

	MainParam* pdp = new GridGenParam(plist);
	orientGrid(static_cast<GridGenParam&>(*pdp), plist);
	GridGenParam dp = static_cast<GridGenParam&>(*pdp);
	////// create the reference frames ////////
	DgRFNetwork net0;
	DgGeoSphRF geoRF(net0, dp.datum, dp.earthRadius);
	DgIDGG dgg(geoRF, dp.vert0, dp.azimuthDegs, dp.aperture, dp.actualRes,
		"DDG", dp.gridTopo, dp.projType, dp.isMixed43, dp.numAp4,
		dp.isSuperfund, dp.sfRes, dp.precision);
	DgGeoSphDegRF deg(geoRF, geoRF.name() + "Deg");
	OGRCoordinateTransformation *coordTrans, *coordTransInv;
	OGRSpatialReference oSRS, oTRS;
	oTRS.SetWellKnownGeogCS("WGS84");
	oSRS = *poLayer->GetSpatialRef();
	/*oSRS.SetProjCS("UTM UTMN (wgs84) in northern hemisphere");
	oSRS.SetWellKnownGeogCS("WGS84");
	oSRS.SetUTM(UTMN, TRUE);*/
	coordTrans = OGRCreateCoordinateTransformation(&oSRS, &oTRS);//from plane to sphere
	coordTransInv = OGRCreateCoordinateTransformation(&oTRS, &oSRS);//from sphere to plane
	double minDist = 10000000;
	//在bdg中找到对应polyID的全部Grid
	for (int i = 1; i < 10; i++)//首先确定相邻polygon，所以此处也是从2开始取的。
	{
		vector<BoundaryG> bdgi = bdg[i].bdg;
		for (auto iter = bdgi.begin(); iter != bdgi.end(); iter++)
		{
			BoundaryG tmpbdgi = *iter;
			if (polyID != tmpbdgi.PolyID) continue;
			if (!tmpbdgi.IS) continue;
			PolyIDbdg.push_back(tmpbdgi);
		}
		bdgi.swap(vector<BoundaryG>());
	}
	//计算到MBR的距离
	int totalsize = PolyIDbdg.size();
	float *Dist = new float[totalsize];
	for (int i = 0; i < totalsize; i++)
	{
		BoundaryG Polybdg = PolyIDbdg[i];
		DgQ2DICoord coor = Polybdg.coor;
		//中心点坐标
		DgLocation* addLoc = dgg.makeLocation(coor);
		DgPolygon verts(dgg);
		dgg.setVertices(*addLoc, verts, dp.nDensify);
		vector<Vec2D> gridVerCor = GridVercoord(verts, coordTransInv);
		double cx = 0.0, cy = 0.0;
		for (int j = 0; j < gridVerCor.size(); j++)
		{
			cx += gridVerCor[j].x;
			cy += gridVerCor[j].y;
		}
		cx = cx / gridVerCor.size();
		cy = cy / gridVerCor.size();
		Dist[i] = CoorMBRDist(cx, cy, poEnvelope);
		gridVerCor.swap(vector<Vec2D>());
	}
	int *index = sortdouble(totalsize, Dist);
	int pos = 0;
	int candidatePID;
	DgQ2DICoord coor_mindist;
	//for (int mm = 0; mm < totalsize; ++mm)
	//{
	//	cout << Dist[index[mm]] << "	" << endl;
	//}
	while (1)
	{
		for (int mm = 0; mm < totalsize; ++mm)
		{
			
			cout << Dist[index[mm]] << "	" << PolyIDbdg[index[mm]].coor<< endl;
		}
		int coorID = index[pos];
		coor_mindist = PolyIDbdg[coorID].coor;
		if (PolyIDbdg[coorID].weight_node < 2)
		{
			candidatePID = 0;
			break;
		}
		//若未找到合适的候选polygon
		candidatePID = MaxPolyID1(bdg, counter, polyID, usePolyID, coor_mindist);
		if (candidatePID == -1)
		{
			pos++;
			if (pos >= totalsize)
				int aaa = 0;
			continue;
		}
		else break;
	}
	cout << polyID << "	" << coor_mindist <<"	0" <<endl;
	Find_PolyIDCoor(bdg, polyID, coor_mindist, pos1, pos2);
	bdg[pos1].bdg[pos2].IS = 0;
	usePolyID.push_back(candidatePID);
	for (int m = 0; m < usePolyID.size(); m++)
		cout << usePolyID[m] << endl;
	PolyIDbdg.swap(vector<BoundaryG>());
	delete[] Dist;
	coord = coor_mindist;
	return candidatePID;
}
int MinDistID(double *dist, int size)
{
	int ID = 0;
	double mindist = 10000000;
	for (int i = 0; i < size; i++)
	{
		if (mindist > dist[i])
		{
			mindist = dist[i];
			ID = i;
		}
	}
	return ID;
}
double CoorMBRDist(double cx, double cy, OGREnvelope *poEnvelope)
{
	double Xmin = poEnvelope->MinX;
	double Xmax = poEnvelope->MaxX;
	double Ymin = poEnvelope->MinY;
	double Ymax = poEnvelope->MaxY;
	double dist[4] = {0};
	dist[0] = abs(cx - Xmin);
	dist[1] = abs(cx - Xmax);
	dist[2] = abs(cy - Ymin);
	dist[3] = abs(cy - Ymax);
	double Dist = 1000000;
	for (int i = 0; i < 4; i++)
	{
		if (Dist > dist[i])
			Dist = dist[i];
	}
	return Dist;
}
int MaxPolyID1(BoundarGs *bdg, vector<Counting> counter, int PolyID, vector<int> &usePolyID, DgQ2DICoord coor)
{
	bool ok = 0;
	int ID = 0;
	double SumArea = 0;
	for (int i = 1; i < 10; i++)
	{
		int maxC = -1;
		vector<BoundaryG> bdgi = bdg[i].bdg;
		vector<BoundaryG>tmp;
		for (int j = 0; j < bdgi.size(); j = j + i)
		{

			if (bdgi[j].coor == coor)
			{
				SumArea = 0;
				for (int k = 0; k < i; k++)
				{
					int pid = bdgi[j + k].PolyID;
					SumArea += bdgi[j + k].area;
					if (pid == PolyID)
						continue;
					tmp.push_back(bdgi[j + k]);
				}
				if (SumArea < Grid_Area) return 0;
				//重复性检测
				for (auto iter = tmp.begin(); iter != tmp.end();)
				{
					int pid = (*iter).PolyID;
					bool issame = false;
					for (int m = 0; m < usePolyID.size(); m++)
					{
						if (pid == usePolyID[m])
						{
							iter = tmp.erase(iter);
							issame = true;
							break;
						}
					}
					if (!issame) iter++;
				}
				if (tmp.size() > 1)
				{
					for (auto iter = tmp.begin(); iter != tmp.end(); iter++)
					{
						int pid = (*iter).PolyID;
						int ac = FindActualCounter(pid, counter);
						if (maxC < ac)
						{
							maxC = ac;
							ID = pid;
						}
					}
				}
				else if (tmp.size() == 1)
					return(tmp[0].PolyID);
				else
					return -1;
				return ID;
			}
		}
		bdgi.swap(vector<BoundaryG>());
		tmp.swap(vector<BoundaryG>());
	}
	return ID;
}
void Find_PolyIDCoor(BoundarGs *bdg, int polyID, DgQ2DICoord coor, int &pos1, int &pos2)
{
	for (int i = 0; i < 10; i++)
	{
		vector<BoundaryG>bdgi = bdg[i].bdg;
		for (int j = 0; j < bdgi.size(); j++)
		{
			BoundaryG tmpbdg = bdgi[j];
			if (tmpbdg.PolyID == polyID && tmpbdg.coor == coor)
			{
				pos1 = i;
				pos2 = j;
				return;
			}
		}
		bdgi.swap(vector<BoundaryG>());
	}
}
int FindActualCounter(int polyid, vector<Counting> counter)
{
	for (int i = 0; i < counter.size(); i++)
	{
		int ID = counter[i].PolyID;
		if (ID == polyid)
			return counter[i].actual_count;
	}
}
int Find_CounterPolyID(vector<Counting> counter, int polyid)
{
	for (int i=0;i<counter.size();i++)
	{
		int ID = counter[i].PolyID;
		if (ID == polyid) return i;
	}
}
vector<BoundaryG> Find_PolyIDBDG(BoundarGs *bdg, int polyid)
{
	vector<BoundaryG> PolyIDbdg;
	//在bdg中找到对应polyID的全部Grid
	for (int i = 1; i < 10; i++)//首先确定相邻polygon，所以此处也是从2开始取的。
	{
		vector<BoundaryG> bdgi = bdg[i].bdg;
		for (auto iter = bdgi.begin(); iter != bdgi.end(); iter++)
		{
			BoundaryG tmpbdgi = *iter;
			if (polyid != tmpbdgi.PolyID) continue;
			if (!tmpbdgi.IS) continue;
			PolyIDbdg.push_back(tmpbdgi);
		}
		bdgi.swap(vector<BoundaryG>());
	}
	return PolyIDbdg;
}
int * Dist2MBR(DgGridPList &plist, vector<BoundaryG> bdgi, char *MBRShpPath)
{
	//根据shppath确定boundary
	GDALDataset *poDS = (GDALDataset*)GDALOpenEx(MBRShpPath, GDAL_OF_VECTOR, NULL, NULL, NULL);
	OGRLayer  *poLayer = poDS->GetLayer(0); //读取层
	OGRFeature *poFeature = poLayer->GetNextFeature();
	OGRGeometry* poGeometry = poFeature->GetGeometryRef();
	OGREnvelope *poEnvelope = new OGREnvelope;
	poGeometry->getEnvelope(poEnvelope);
	OGRCoordinateTransformation *coordTrans, *coordTransInv;
	OGRSpatialReference oSRS, oTRS;
	oTRS.SetWellKnownGeogCS("WGS84");
	oSRS = *poLayer->GetSpatialRef();
	/*oSRS.SetProjCS("UTM UTMN (wgs84) in northern hemisphere");
	oSRS.SetWellKnownGeogCS("WGS84");
	oSRS.SetUTM(UTMN, TRUE);*/
	coordTrans = OGRCreateCoordinateTransformation(&oSRS, &oTRS);//from plane to sphere
	coordTransInv = OGRCreateCoordinateTransformation(&oTRS, &oSRS);//from sphere to plane
	//Dggrid的环境设定
	MainParam* pdp = new GridGenParam(plist);
	orientGrid(static_cast<GridGenParam&>(*pdp), plist);
	GridGenParam dp = static_cast<GridGenParam&>(*pdp);
	////// create the reference frames ////////
	DgRFNetwork net0;
	DgGeoSphRF geoRF(net0, dp.datum, dp.earthRadius);
	DgIDGG dgg(geoRF, dp.vert0, dp.azimuthDegs, dp.aperture, dp.actualRes,
		"DDG", dp.gridTopo, dp.projType, dp.isMixed43, dp.numAp4,
		dp.isSuperfund, dp.sfRes, dp.precision);
	DgGeoSphDegRF deg(geoRF, geoRF.name() + "Deg");
	
	//计算到MBR的距离
	int totalsize = bdgi.size();
	float *Dist = new float[totalsize];
	for (int i = 0; i < totalsize; i++)
	{
		BoundaryG Polybdg = bdgi[i];
		DgQ2DICoord coor = Polybdg.coor;
		//中心点坐标
		DgLocation* addLoc = dgg.makeLocation(coor);
		DgPolygon verts(dgg);
		dgg.setVertices(*addLoc, verts, dp.nDensify);
		vector<Vec2D> gridVerCor = GridVercoord(verts, coordTransInv);
		double cx = 0.0, cy = 0.0;
		for (int j = 0; j < gridVerCor.size(); j++)
		{
			cx += gridVerCor[j].x;
			cy += gridVerCor[j].y;
		}
		cx = cx / gridVerCor.size();
		cy = cy / gridVerCor.size();
		Dist[i] = CoorMBRDist(cx, cy, poEnvelope);
		gridVerCor.swap(vector<Vec2D>());
	}
	int *index = sortdouble(totalsize, Dist);
	return index;
}
bool Same_LOC(vector<Location> LOC, int pos1, int pos2)
{
	for (auto iter = LOC.begin(); iter != LOC.end(); iter++)
	{
		int i = (*iter).pos1;
		int j = (*iter).pos2;
		if (i == pos1 && j == pos2) return 1;
	}
	return 0;
}
int PolyGonChange(int polyid, DgQ2DICoord coor, BoundarGs *bdg, vector<Counting> counter)
{
	int ID = 0;
	int maxCount = -1;
	double SumArea = 0;
	vector<BoundaryG>tmp;
	for (int i = 1; i < 10; i++)//首先确定相邻polygon，所以此处也是从2开始取的。
	{
		vector<BoundaryG> bdgi = bdg[i].bdg;
		for (int j = 0; j < bdgi.size(); j = j + i)
		{
			//cout << bdgi[j].coor << endl;
			if (bdgi[j].coor != coor) 
				continue;
			SumArea = 0;
			for (int k = 0; k < i; k++)
			{
				int pid = bdgi[j + k].PolyID;
				SumArea += bdgi[j + k].area;
				if (pid == polyid)
					continue;
				tmp.push_back(bdgi[j + k]);
			}
		}
		bdgi.swap(vector<BoundaryG>());
	}
	if (SumArea < Grid_Area)
	{
		ID = 0;
	}
	if (tmp.size() == 1)
	{
		ID = tmp[0].PolyID;
	}
	for (auto iter = tmp.begin(); iter != tmp.end(); iter++)
	{
		int pid = (*iter).PolyID;
		int ac = FindActualCounter(pid, counter);
		if (maxCount < ac)
		{
			maxCount = ac;
			ID = pid;
		}
	}
	tmp.swap(vector<BoundaryG>());
	return ID;
}
bool Line_Line_intersection(double Ax, double Ay, double Bx, double By, double Px, double Py, double Qx, double Qy)
{
	//1 按A点坐标，整体平移至原点
	Bx = Bx - Ax;
	Px = Px - Ax;
	Qx = Qx - Ax;
	By = By - Ay;
	Py = Py - Ay;
	Qy = Qy - Ay;
	Ax = Ax - Ax;
	Ay = Ay - Ay;

	//2 按AB与x轴正向夹角，整体旋转
	double AB = sqrtf(Bx * Bx + By * By);
	double sint = By / AB;
	double cost = Bx / AB;
	double BxR = Bx * cost + By * sint;
	double ByR = -Bx * sint + By * cost;
	double PxR = Px * cost + Py * sint;
	double PyR = -Px * sint + Py * cost;
	double QxR = Qx * cost + Qy * sint;
	double QyR = -Qx * sint + Qy * cost;
	//3 判断PQ与AB的位置
	if (PyR * QyR > 0) return 0;//PQ位于同侧
								//
	double InterSectX = PxR + (QxR - PxR) *((0 - PyR) / (QyR - PyR));
	if ((InterSectX - 0) < 1e-5) return 0;
	if ((InterSectX - AB) > 1e-4) return 0;
	return 1;
}
bool Line_BD_intersect(double Ax, double Ay, double Bx, double By, OGRGeometry *poGeometry)
{
	bool non_intersect = false;
	//计算AB旋转至X轴正向的参数
	double Bx_ = Bx - Ax;
	double By_ = By - Ay;
	double AB = sqrtf(Bx_ * Bx_ + By_ * By_);
	double sint = By_ / AB;
	double cost = Bx_ / AB;
	double BxR = Bx_ * cost + By_ * sint;
	double ByR = -Bx_ * sint + By_ * cost;
	OGRwkbGeometryType pGeoType = poGeometry->getGeometryType();
	if (pGeoType == wkbMultiLineString)
	{

	}
	else
	{
		OGRLinearRing *boundary = (OGRLinearRing *)poGeometry;
		int pointcount = boundary->getNumPoints();
		double Px = 0, Py = 0, Qx = 0, Qy = 0;
		for (int np = 0; np < pointcount - 1; np++)
		{
			Px = boundary->getX(np);
			Py = boundary->getY(np);
			Qx = boundary->getX(np + 1);
			Qy = boundary->getY(np + 1);
			bool inter = Line_Line_intersection(Ax, Ay, Bx, By, Px, Py, Qx, Qy);
			if (inter == 0)
			{
				continue;
			}
			else return 1;
		}
		return 0;
	}
}
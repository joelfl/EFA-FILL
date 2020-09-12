#include "GLC30_Based.h"
SeqStack * init_SeqStack(SeqStack * top){
	top = NULL;
	return top;
}
int is_Empty(SeqStack * top) {
	if (top == NULL)return 1;
	else return 0;
}
SeqStack * push_Stack(SeqStack * top, DgQ2DICoord coor) {
	SeqStack * New;
	New = (SeqStack *)malloc(sizeof(SeqStack));
	New->COOR = coor;
	New->next = top;
	top = New;
	return top;
}
SeqStack * pop_Stack(SeqStack * top, DgQ2DICoord &coor) {
	SeqStack * p = NULL;
	coor = top->COOR;
	p = top;
	top = top->next;
	free(p);
	return top;
}
SeqStack * top_Stack(SeqStack * top, DgQ2DICoord &coor) {
	if (!is_Empty(top)) {
		coor = top->COOR;
		return top;
	}
}

//初始化队列
que * InitQueue()
{
	que * q = (que*)malloc(sizeof(que));
	q->front = q->rear = NULL;
	return q;
}
//判断队列是否为空
int EmptyQueue(que * q)
{
	return q->front == NULL;
}
//入队
void InsertQueue(que *q, DgQ2DICoord data)
{
	queueNode * n = (queueNode *)malloc(sizeof(queueNode));
	if (n == NULL)//内存分配失败
		return;
	n->COOR = data;
	n->next = NULL;//尾插法，插入元素指向空
	if (q->rear == NULL)
	{
		q->front = n;
		q->rear = n;
	}
	else {
		q->rear->next = n;//让n成为当前的尾部节点下一节点
		q->rear = n;//尾部指针指向n
	}
}

//出队(删除队首元素)
void DeleteQueue(que *q, DgQ2DICoord &data)
{
	queueNode * n = q->front;
	data = n->COOR;
	if (q->front == NULL)//判断队列是否为空
		return;
	if (q->front == q->rear)//是否只有一个元素
	{
		q->front = NULL;
		q->rear = NULL;
	}
	else {
		q->front = q->front->next;
		free(n);
	}
}

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

unsigned short * getImgData(GDALDataset *poDataset, int nbands, int width, int height, GDALDataType *gdt)
{
	unsigned short *imgBuf = new unsigned short[width*height];
	if (nbands == 1)
	{
		unsigned short *imgBuf = new unsigned short[width*height];
		GDALRasterBand *poBand = poDataset->GetRasterBand(1);
		int re = poBand->RasterIO(GF_Read, 0, 0, width, height, imgBuf, width, height, GDT_UInt16, 0, 0);
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
	unsigned short *imgData, double *geoTrans, int width, int height, int *val)
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
	for (k = 0; k < 12; k++)
	{
		ai = gird_area[k];
		if (ai >= area)
		{
			area = ai;
			glc_class = k;
		}
	}
	if (area == 0)
		glc_class = 11;

	*val = glc_class;
	//vector<Vec3D>().swap(IJ);

	//return;
}

void genGrid_(GridGenParam& dp, char *SHP_path, string polyID)
{
	GDALDataset   *poDS, *poDS1;
	OGRLayer  *poLayer, *poLayer1;
	OGRFeature *poFeature, *poFeature1;
	Poly_BDRasters * PBDRs, *PBDRs_Line1, *PBDRs_Line2;
	OGRGeometry *poGeometry;
	OGRwkbGeometryType pGeoType;
	OGRPolygon *BDPolygon;
	OGRLinearRing  grid_ring;
	OGRSpatialReference DstSPF;
	vector<Vec2D> gridVerCor, gridVerCor_tmp;
	vector<DgQ2DICoord> Vctmp;

	//定义GDAL环境
	GDALDriver *poDriver;
	poDriver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
	/*DstSPF.SetProjCS("UTM 12(WGS84) in northern hemisphere.");
	DstSPF.SetWellKnownGeogCS("WGS84");
	DstSPF.SetUTM(12, TRUE);*/
	//读取shp文件
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

																								   //遍历PBDRs，将结果写成SHP
																								   //string shpoutpath = "D:\\data\\SAMPLE\\edge.shp";
	string shpoutpath = "D:\\data\\SAMPLE\\GlobeLand30\\SHP\\OutSHP\\" + polyID + "_" + to_string(NLEVEL) + "_ALL.shp";
	// 创建属性字段
	{
		OGRFieldDefn oFieldId("Type", OFTInteger);//多边形属性
		oFieldId.SetWidth(5);
		OGRFieldDefn firstField("Quam", OFTInteger);//编号，为浮点类型
		firstField.SetWidth(2);
		//fifthField.SetPrecision(2);
		OGRFieldDefn secondField("I", OFTInteger);//编号，为浮点类型
		secondField.SetWidth(255);
		//sixthField.SetPrecision(10);
		OGRFieldDefn thirdField("J", OFTInteger);//编号，为浮点类型
		thirdField.SetWidth(255);
		//sixthField.SetPrecision(10);
		//OGRFieldDefn seventhField("50", OFTReal);//编号，为浮点类型
		//										 //seventhField.SetWidth(255);
		//seventhField.SetPrecision(10);
		//OGRFieldDefn eighthField("60", OFTReal);//编号，为浮点类型
		//										//eighthField.SetWidth(255);
		//eighthField.SetPrecision(10);

		//OGRFieldDefn ninthField("70", OFTReal);//编号，为浮点类型
		//									   //ninthField.SetWidth(255);
		//ninthField.SetPrecision(10);

		//OGRFieldDefn tenthField("80", OFTReal);//编号，为浮点类型
		//									   //tenthField.SetWidth(255);
		//tenthField.SetPrecision(10);

		//OGRFieldDefn eleventhField("90", OFTReal);//编号，为浮点类型
		//										  //eleventhField.SetWidth(255);
		//eleventhField.SetPrecision(10);

		//OGRFieldDefn twelfthFileld("100", OFTReal);//编号，为浮点类型
		//										   //twelfthFileld.SetWidth(255);
		//twelfthFileld.SetPrecision(10);

		////创建shp文件
		poDS1 = poDriver->Create(shpoutpath.data(), 0, 0, 0, GDT_Unknown, NULL);

		////创建图层文件，一般为1个图层，图层中添加空间参考与几何类型
		poLayer1 = poDS1->CreateLayer("DGGRID_T", &DstSPF, wkbPolygon, NULL);
		poLayer1->CreateField(&oFieldId);
		poLayer1->CreateField(&firstField);
		poLayer1->CreateField(&secondField);
		poLayer1->CreateField(&thirdField);/*
		poLayer1->CreateField(&forthField);
		poLayer1->CreateField(&fifthField);
		poLayer1->CreateField(&sixthField);
		poLayer1->CreateField(&seventhField);
		poLayer1->CreateField(&eighthField);*/
		/*poLayer->CreateField(&ninthField);
		poLayer->CreateField(&tenthField);
		poLayer->CreateField(&eleventhField);
		poLayer->CreateField(&twelfthFileld);*/
	}
	poFeature1 = OGRFeature::CreateFeature(poLayer1->GetLayerDefn());


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
			ofstream fout("grid_10.txt");
			ofstream fout_("Edge_Prox_10.txt");
			ofstream fout__("Edge_Prox_Cell_10.txt");
			while (1)
			{
				dp.nCellsAccepted++;
				dp.nCellsTested++;
				outputStatus(dp);

				DgPolygon verts(dgg);
				dgg.setVertices(*addLoc, verts, dp.nDensify);
				dgg.convert(addLoc);
				DgLocation ADDLOC = *addLoc;
				///////////遍历当前addLoc，查找其邻近单元，并输出//////////
				//vector<DgQ2DICoord>QIJ = Edge_Proximity(ADDLOC, verts, dgg.maxI(), dgg.maxJ());
				vector<DgQ2DICoord>QIJ;
				const DgGeoSphRF* geoRF = dynamic_cast<const DgGeoSphRF*>(&verts.rf());
				GeoCoord center, vertx[3];
				fout <<"3"<<endl;
				for (int i = 0; i < 3; i++)
				{
					const DgGeoCoord& v = *geoRF->getAddress(verts[i]);
					double x = v.lat();
					double y = v.lon();
					vertx[i].lat = x;
					vertx[i].lon = y;
				}
				center = sphTricenpoint(vertx);
				fout << setprecision(15) << center.lat<< "	" << center.lon<< endl;
				for (int i = 0; i < 3; i++)
				{
					const DgGeoCoord& v = *geoRF->getAddress(verts[i]);
					fout << setprecision(15) << v.lat() << "	" << v.lon() << endl;
				}
				fout << endl;
				fout_ << "3" << endl;
				fout__ << "3" << endl;
				for (int k = 0; k < QIJ.size(); k++)
				{
					DgQ2DICoord coor = QIJ[k];
					DgLocation tmp = *addLoc;
					((DgAddress<DgQ2DICoord>*)((&tmp)->address()))->setAddress(coor);
					dgg.setVertices(tmp, verts, dp.nDensify);
					const DgGeoSphRF* geoRF_ = dynamic_cast<const DgGeoSphRF*>(&verts.rf());
					for (int i = 0; i < 3; i++)
					{
						const DgGeoCoord& v = *geoRF->getAddress(verts[i]);
						double x = v.lat();
						double y = v.lon();
						vertx[i].lat = x;
						vertx[i].lon = y;
					}
					center = sphTricenpoint(vertx);
					fout_ << setprecision(15) << center.lat << "	" << center.lon << endl;
					for (int i = 0; i < 3; i++)
					{
						const DgGeoCoord& v = *geoRF->getAddress(verts[i]);
						fout__ << setprecision(15) << v.lat() << "	" << v.lon() << endl;
					}
					fout__ << endl;
				}
				fout_ << endl;
				//////////////////////////////////////////////////////////
				dgg.bndRF().incrementLocation(*addLoc);
				if (!dgg.bndRF().validLocation(*addLoc)) break;
			}
			delete addLoc;
			fout.close();
			fout_.close();
			fout__.close();
		}
		else // dp.isSuperfund
		{
			for (int q = 0; q < 12; q++)
			{
				DgHexSF baseTile(0, 0, 0, 0, true, q);
				//=baseTile.setType('P');
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
					//getValByPolygonIntesect(gridVerCor, imgData, geoTrans, Width, Height, &gridval);
					////////
					//GLC[gridval] = GLC[gridval] + 1;//统计地类个数

					//////
					//vector<Vec2D>().swap(gridVerCor);
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
					OGRLinearRing ring;
					for (int i = 0; i < gridVerCor.size(); i++)
					{
						ring.addPoint(gridVerCor[i].x, gridVerCor[i].y);
					}
						ring.closeRings();
						//环加入到polygon中
						OGRPolygon polygon;
						polygon.addRing(&ring);
						//

						poFeature1->SetGeometry(&polygon);
						poFeature1->SetField(1, 1);
						poFeature1->SetField(2, coord.i());
						poFeature1->SetField(3, coord.j());
						poLayer1->CreateFeature(poFeature1);
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
	//*GLC30 = GLC;
} // void genGrid_

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
	//remainder = "WHOLE_EARTH";
	plist->setParam(token, remainder);
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
void DGGRID_GLC30(DgGridPList &plist, char *SHP_path, int **GLC30, string polyID, char * txt_path, char * outshp_path, char * outtxt_Path)
{
	MainParam* pdp = new GridGenParam(plist);
	orientGrid(static_cast<GridGenParam&>(*pdp), plist);
	CPLSetConfigOption("SHAPE_ENCODING", "");  //解决中文乱码问题
											   //读取shp文件
	//EAF
	BDTracing(static_cast<GridGenParam&>(*pdp), SHP_path, polyID, outshp_path, outtxt_Path);
	delete pdp;
}

/////////确定三角形朝向 up or down
Vec2D FindSamePoint(vector<Vec2D> P1, vector<Vec2D> P2)
{
	for (int k = 0; k < 3; k++)
	{
		double x = P1[k].x;
		double y = P1[k].y;
		cout << setprecision(10) << x << "	" << y << endl;
		for (int j = 0; j < 3; j++)
		{
			cout << setprecision(10) << P2[j].x << "	" << P2[j].y << endl;
			if (abs(x - P2[j].x) <= 10e-10 && abs(y - P2[j].y) <= 10e-10)
			{
				cout << setprecision(10) << P2[j].x << "	" << P2[j].y << endl;
				cout<< endl;
				return P2[j];
			}
		}
	}
}
bool PolygonOrti(vector<Vec2D> gridVerCor)
{
	//三角形朝向判断 2020 05 05 修改
	//思路：预定以在第0层时，同一quam的两个三角形的公共边。
	//		按照公共边旋转到x轴正向确定旋转矩阵参数（因为本实验中国取在同一quam下，所以矩阵参数是相同的，定义在外面）
	//		每一个小grid按照参数旋转
	//		找到旋转后找到y不等于0的点，若y>0则正向；反之为负向
	int count = gridVerCor.size();
	int tmpcount1 = 0, tmpcount2 = 0;
	double RY[6] = { 0.0 };
	for (int i = 0; i < count; i++)
	{
		double cx = gridVerCor[i].x;
		double cy = gridVerCor[i].y;
		cx = cx - _AX;
		cy = cy - _AY;
		double Rx = cx * COSt + cy * SINt;
		double Ry = -cx * SINt + cy * COSt;
		RY[i] = Ry;
		//cout << Rx << " " << Ry << endl;
	}
	//cout << endl;
	RY[3] = RY[0];
	RY[4] = RY[1];
	RY[5] = RY[2];
	double RY0 = RY[0], RY2 = 0;
	int thirdP = 0;
	for (int i = 0; i < count; i++)
	{
		double deltaRY01 = abs(RY[i] - RY[i + 1]);
		double deltaRY12 = abs(RY[i + 1] - RY[i + 2]);
		double deltaRY20 = abs(RY[i + 2] - RY[i]);

		if (deltaRY01 > deltaRY12 && deltaRY20 > deltaRY12)
		{
			thirdP = i;
			break;
		}
		/*
		if (abs(RY0 - RY[i]) < 1)
		{
			thirdP = count - i;
			break;
		}*/
	}
	//cout << RY[thirdP] << " " << RY[thirdP + 1] << " " << RY[thirdP + 2] << endl;
	if (RY[thirdP] > RY[thirdP + 1] && RY[thirdP] > RY[thirdP + 2])
		return 1;
	else
		return 0;

	//double RY3 = RY[thirdP];
	//if (thirdP == 0)
	//{
	//	RY2 = RY[count - thirdP - 1];
	//}
	//else  RY2 = RY[count - thirdP];
	//if (abs(RY3) > abs(RY2))
	//{
	//	if (RY3 > 0) return 1;
	//	else return 0;
	//}
	//if (abs(RY3) < abs(RY2))
	//{
	//	if (RY3 > 0) return 0;
	//	else return 1;
	//} //0 表示 down tri 1 表示 up tri
}

bool PolygonOrti_(DgQ2DICoord coor, DgLocation ADDLOC, DgIDGG dgg, int nDensify)
{
	DgPolygon verts(dgg);
	vector<Vec2D>A1, B1, P;
	int q = coor.quadNum();
	int i = coor.coord().i();
	int j = coor.coord().j();
	DgQ2DICoord coorA1 = DgQ2DICoord(q ,DgIVec2D(i,j - 2));
	DgQ2DICoord coorB1 = DgQ2DICoord(q, DgIVec2D(i, j + 2));
	//P
	((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coor);
	dgg.setVertices(ADDLOC, verts, nDensify);
	const DgGeoSphRF* geoRF = dynamic_cast<const DgGeoSphRF*>(&verts.rf());
	int count = verts.size();
	cout << "P" << endl;

	for (int i = 0; i < count; i++)
	{
		const DgGeoCoord& v = *geoRF->getAddress(verts[i]);
		Vec2D tmp;
		tmp.x = v.lon();
		tmp.y = v.lat();
		P.push_back(tmp);
		cout << std::setprecision(10) <<v.lon() << "	" << v.lat() << endl;
	}
	//A1
	((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coorA1);
	dgg.setVertices(ADDLOC, verts, nDensify);
	geoRF = dynamic_cast<const DgGeoSphRF*>(&verts.rf());
	cout << "A1" << endl;
	for (int i = 0; i < count; i++)
	{
		const DgGeoCoord& v = *geoRF->getAddress(verts[i]);
		Vec2D tmp;
		tmp.x = v.lon();
		tmp.y = v.lat();
		A1.push_back(tmp);
		cout << std::setprecision(10)<< v.lon() << "	" << v.lat() << endl;
	}
	//B1
	cout << "B1" << endl;

	((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coorB1);
	dgg.setVertices(ADDLOC, verts, nDensify);
	geoRF = dynamic_cast<const DgGeoSphRF*>(&verts.rf());
	for (int i = 0; i < count; i++)
	{
		const DgGeoCoord& v = *geoRF->getAddress(verts[i]);
		Vec2D tmp;
		tmp.x = v.lon();
		tmp.y = v.lat();
		B1.push_back(tmp);
		cout << std::setprecision(10)<< v.lon() << "	" << v.lat() << endl;
	}
	
	Vec2D A = FindSamePoint(A1, P);
	Vec2D B = FindSamePoint(B1, P);
	//将AB旋转至x轴正向
	double Ax = A.x, Ay = A.y;
	double Bx = B.x - Ax, By = B.y - Ay;
	//计算AB的平移旋转矩阵
	double AB = sqrtf(Bx * Bx + By * By);
	double sint = By / AB;
	double cost = Bx / AB;

	//计算(i,j)的中心点坐标
	((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coor);
	dgg.setVertices(ADDLOC, verts, nDensify);
	geoRF = dynamic_cast<const DgGeoSphRF*>(&verts.rf());
	double mean_lon = 0, mean_lat = 0;
	for (int i = 0; i < count; i++)
	{
		const DgGeoCoord& v = *geoRF->getAddress(verts[i]);
		mean_lon += v.lon();
		mean_lat += v.lat();
		double xx = (v.lon()-Ax) * cost + (v.lat()-Ay) * sint;
		double yy = -(v.lon() - Ax) * sint + (v.lat() - Ay) * cost;
		cout << xx << "	" << yy << endl;
	}
	mean_lat = mean_lat / count;
	mean_lon = mean_lon / count;
	double lon = (mean_lon - Ax) * cost + (mean_lat - Ay) * sint;
	double lat = -(mean_lon - Ax) * sint + (mean_lat - Ay)* cost;
	if (lat > 0) return 1;//1 表示 up tri
	return 0;// 0 表示 down tri
}
/////////查找三角形的边邻近单元//////
vector<DgQ2DICoord> Edge_Proximity(DgLocation& Loc, int Quam, int i, int j, DgPolygon verts, int maxI, int maxJ)
{
	vector<DgQ2DICoord> COOR;
	bool UP = 1;
	//1st判断verts代表的三角形的朝向，uptri or downtri
	//bool UP = PolygonOrti(verts);
	/*const DgAddress<DgQ2DICoord> *add = (DgAddress<DgQ2DICoord>*)((Loc).address());*/
	/*int Quam = add->address().quadNum();
	int i = add->address().coord().i();
	int j = add->address().coord().j();*/
	DgQ2DICoord coor;
	DgLocation* tmp = &Loc;
	//通过ij判断DgPolygon是边界单元or非边界单元 i/j!=0/!=Imax/Jmax
	if (!(i == 0 & j == 0 || (i == 0 && j % 2 == 1) || (j == 0) || (i == maxI && j % 2 == 0) || (j == maxJ)))
	{
		//非边界单元
		if (UP)//上三角形边邻近顺序：左-右-下
		{
			DgQ2DICoord coor = DgQ2DICoord(Quam, DgIVec2D(i , j + 1));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			coor = DgQ2DICoord(Quam, DgIVec2D(i - 1, j - 1));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			coor = DgQ2DICoord(Quam, DgIVec2D(i, j - 1));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			return COOR;
		}
		else //下三角形边邻近顺序：左-右-上
		{
			coor = DgQ2DICoord(Quam, DgIVec2D(i + 1, j + 1));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			coor = DgQ2DICoord(Quam, DgIVec2D(i, j - 1));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			coor = DgQ2DICoord(Quam, DgIVec2D(i, j + 1));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			return COOR;
		}
	}
	else if(Quam <= 5)
	{
		int Quam1 = -1;
		//边界单元
		if (i == 0 && j != 0)
		{
			//右
			coor = DgQ2DICoord(Quam, DgIVec2D(i, j + 1));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			//下
			coor = DgQ2DICoord(Quam, DgIVec2D(i, j - 1));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			//左
			Quam1 = Quam;
			if (Quam1 > 5) Quam1 = Quam1 - 5;

			else if (Quam1 < 6)
			{
				Quam1 = Quam1 - 1;
				if (Quam1 == 0) Quam1 = 5;
			}
			coor = DgQ2DICoord(Quam1, DgIVec2D((maxJ-j)/2, maxJ));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			return COOR;
		}
		else if (j == 0)
		{
			//左
			coor = DgQ2DICoord(Quam, DgIVec2D(i + 1, j + 1));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			//上
			coor = DgQ2DICoord(Quam, DgIVec2D(i, j + 1));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			//右
			Quam1 = Quam;
			if (Quam1 > 5) Quam1 = Quam1 - 5;
			Quam1 = Quam1 + 4;
			if (Quam1 < 6) Quam1 = 10;
			coor = DgQ2DICoord(Quam1, DgIVec2D(i, maxJ));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			return COOR;
		}
		else if (i == maxI && j != 0 && j != maxJ)
		{
			//右
			coor = DgQ2DICoord(Quam, DgIVec2D(i, j - 1));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			//下
			coor = DgQ2DICoord(Quam, DgIVec2D(i, j + 1));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			//左
			Quam1 = Quam;
			if (Quam1 <= 5) Quam1 = Quam1 + 5;
			else if (Quam1 > 5)
			{
				Quam1 = Quam1 + 1;
				if (Quam1 > 10) Quam1 = 6;
			}
			coor = DgQ2DICoord(Quam1, DgIVec2D(0, j + 1));//可能有问题
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			return COOR;
		}
		else 
		if (j == maxJ && i != 0)
		{
			//右
			coor = DgQ2DICoord(Quam, DgIVec2D(i - 1, j - 1));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			//上
			coor = DgQ2DICoord(Quam, DgIVec2D(i, j - 1));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			//右
			Quam1 = Quam;
			if (Quam1 > 5) Quam1 = Quam1 - 5;
			Quam1 = Quam1 + 1;
			if (Quam1 > 5) Quam1 = 1;
			coor = DgQ2DICoord(Quam1, DgIVec2D(0, 2 * (maxI - i) + 1));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			return COOR;
		}
	}
	else if(Quam >= 6)
	{
		int Quam1 = -1;
		if (i == 0 && j != 0)
		{
			//左
			coor = DgQ2DICoord(Quam, DgIVec2D(i, j + 1));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			//下
			coor = DgQ2DICoord(Quam, DgIVec2D(i, j - 1));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			//右
			Quam1 = Quam;
			if (Quam1 > 5) Quam1 = Quam1 - 5;
			else if (Quam1 < 6)
			{
				Quam1 = Quam1 - 1;
				if (Quam1 == 0) Quam1 = 5;
			}
			coor = DgQ2DICoord(Quam1, DgIVec2D(maxI, j-1));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			return COOR;
		}
		else if (j == 0)
		{
			//左
			coor = DgQ2DICoord(Quam, DgIVec2D(i + 1, j + 1));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			//上
			coor = DgQ2DICoord(Quam, DgIVec2D(i, j + 1));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			//右
			Quam1 = Quam;
			if (Quam1 > 5) Quam1 = Quam1 - 5;
			Quam1 = Quam1 + 4;
			if (Quam1 < 6) Quam1 = 10;
			coor = DgQ2DICoord(Quam1, DgIVec2D(maxI, 2*(maxI - i)));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			return COOR;
		}
		else if (i == maxI && j != 0 && j != maxJ)
		{
			//右
			coor = DgQ2DICoord(Quam, DgIVec2D(i, j - 1));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			//下
			coor = DgQ2DICoord(Quam, DgIVec2D(i, j + 1));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			//左
			Quam1 = Quam;
			if (Quam1 <= 5) Quam1 = Quam1 + 5;
			else if (Quam1 > 5)
			{
				Quam1 = Quam1 + 1;
				if (Quam1 > 10) Quam1 = 6;
			}
			coor = DgQ2DICoord(Quam1, DgIVec2D((maxJ-1-j)/2, 0));//可能有问题
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			return COOR;
		}
		else if (j == maxJ && i != 0)
		{
			//右
			coor = DgQ2DICoord(Quam, DgIVec2D(i - 1, j - 1));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			//下
			coor = DgQ2DICoord(Quam, DgIVec2D(i, j - 1));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			//左
			Quam1 = Quam;
			if (Quam1 > 5) Quam1 = Quam1 - 5;
			Quam1 = Quam1 + 1;
			if (Quam1 > 5) Quam1 = 1;
			coor = DgQ2DICoord(Quam1, DgIVec2D(i, 0));
			((DgAddress<DgQ2DICoord>*)((tmp)->address()))->setAddress(coor);
			COOR.push_back(coor);
			return COOR;
		}
	}
	//DgQ2DICoord coor = DgQ2DICoord(Quam + 1, DgIVec2D(i, j));
	//DgLocation* addLoc1 = addLoc;
	//const DgAddress<DgQ2DICoord> *add1 = (DgAddress<DgQ2DICoord>*)((addLoc1)->address());
	//((DgAddress<DgQ2DICoord>*)((addLoc1)->address()))->setAddress(coor);
	////outputCellAdd2D(dp, dgg, *addLoc, verts, deg);
	//const DgGeoSphRF* geoRF = dynamic_cast<const DgGeoSphRF*>(&verts.rf());
	//out << "Quam=" << Quam << endl;
}

void Grid_located(GridGenParam& dp, unsigned short * imgData, int Width, int Height, double *geoTrans,
	OGRCoordinateTransformation *coordTrans, OGRCoordinateTransformation *coordTransInv, int **GLC30)
{
	//定义基础变量
	int total_grid = 0;
	int total_SHP = 0;//目前已生成的shp总数
	int compilte_cell = 1000 * total_SHP;//目前生成的gird总数=101*total_SHP
	int total_SHP_Max = 0;//允许生成shp的最大值 加上total_SHP的

						  //定义GDAL环境
	GDALDriver *poDriver;
	poDriver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
	OGRFeature *poFeature;
	OGRLayer *poLayer;
	////设置Dst的空间参考系
	OGRSpatialReference DstSPF;
	DstSPF.SetProjCS("UTM 45(WGS84) in northern hemisphere.");
	DstSPF.SetWellKnownGeogCS("WGS84");
	DstSPF.SetUTM(12, TRUE);
	// 创建属性字段
	OGRFieldDefn oFieldId("GLC", OFTInteger);//编号，为整数类型
	oFieldId.SetWidth(5);
	OGRFieldDefn firstField("x", OFTReal);//编号，为浮点类型
										  //firstField.SetWidth(255);
	firstField.SetPrecision(10);
	OGRFieldDefn secondField("y", OFTReal);//编号，为浮点类型
										   //secondField.SetWidth(255);
	secondField.SetPrecision(10);
	OGRFieldDefn thirdField("10", OFTReal);//编号，为浮点类型
										   //thirdField.SetWidth(255);
	thirdField.SetPrecision(10);
	OGRFieldDefn forthField("20", OFTReal);//编号，为浮点类型
										   //forthField.SetWidth(255);
	forthField.SetPrecision(10);
	OGRFieldDefn fifthField("30", OFTReal);//编号，为浮点类型
										   //fifthField.SetWidth(255);
	fifthField.SetPrecision(10);
	OGRFieldDefn sixthField("40", OFTReal);//编号，为浮点类型
										   //sixthField.SetWidth(255);
	sixthField.SetPrecision(10);
	OGRFieldDefn seventhField("50", OFTReal);//编号，为浮点类型
											 //seventhField.SetWidth(255);
	seventhField.SetPrecision(10);
	OGRFieldDefn eighthField("60", OFTReal);//编号，为浮点类型
											//eighthField.SetWidth(255);
	eighthField.SetPrecision(10);

	OGRFieldDefn ninthField("70", OFTReal);//编号，为浮点类型
										   //ninthField.SetWidth(255);
	ninthField.SetPrecision(10);

	OGRFieldDefn tenthField("80", OFTReal);//编号，为浮点类型
										   //tenthField.SetWidth(255);
	tenthField.SetPrecision(10);

	OGRFieldDefn eleventhField("90", OFTReal);//编号，为浮点类型
											  //eleventhField.SetWidth(255);
	eleventhField.SetPrecision(10);

	OGRFieldDefn twelfthFileld("100", OFTReal);//编号，为浮点类型
											   //twelfthFileld.SetWidth(255);
	twelfthFileld.SetPrecision(10);

	string shpoutpath = "E:\\program\\data\\SHP\\USA_3city\\YTZ_with_n12_40_water_local_3_line_point_18_" + to_string(total_SHP) + ".shp";

	////创建shp文件
	GDALDataset *poDS = poDriver->Create(shpoutpath.data(), 0, 0, 0, GDT_Unknown, NULL);

	////创建图层文件，一般为1个图层，图层中添加空间参考与几何类型
	poLayer = poDS->CreateLayer("DGGRID_T", &DstSPF, wkbPolygon, NULL);
	poLayer->CreateField(&oFieldId);
	poLayer->CreateField(&firstField);
	poLayer->CreateField(&secondField);
	poLayer->CreateField(&thirdField);
	poLayer->CreateField(&forthField);
	poLayer->CreateField(&fifthField);
	poLayer->CreateField(&sixthField);
	poLayer->CreateField(&seventhField);
	poLayer->CreateField(&eighthField);
	poLayer->CreateField(&ninthField);
	poLayer->CreateField(&tenthField);
	poLayer->CreateField(&eleventhField);
	poLayer->CreateField(&twelfthFileld);

	poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
	//
	double g0 = geoTrans[0],
		g1 = geoTrans[1],
		g2 = geoTrans[2],//0
		g3 = geoTrans[3],
		g4 = geoTrans[4],//0
		g5 = geoTrans[5];

	////// create the reference frames ////////

	DgRFNetwork net0;
	DgGeoSphRF geoRF(net0, dp.datum, dp.earthRadius);
	DgIDGG dgg(geoRF, dp.vert0, dp.azimuthDegs, dp.aperture, dp.actualRes,
		"DDG", dp.gridTopo, dp.projType, dp.isMixed43, dp.numAp4,
		dp.isSuperfund, dp.sfRes, dp.precision);

	// set-up to convert to degrees
	DgGeoSphDegRF deg(geoRF, geoRF.name() + "Deg");

	////open RandomPt file
	string Ptpath = "E:\\program\\data\\TXT\\USA_3city\\YTZ_with_n12_40_water_local_3_line_point.txt";
	//ofstream outfile(outpath.data());
	double A, B, C, D, E;
	double lon = 0.0, lat = 0.0;
	const int maxLine = 100;
	char line[maxLine];
	DgInputStream inFile(Ptpath, "", DgBase::Fatal);
	inFile.getline(line, maxLine);
	while (inFile.getline(line, maxLine))
	{
		//sscanf(line, "%lf, %lf, %lf, %lf, %d, %lf, %lf",&A, &B, &C, &D, &E, &lon, &lat);
		sscanf(line, ", %lf, %lf, %lf, %lf, %lf", &A, &B, &C, &lon, &lat);
		// 确定格网属性
		DgLocation* tloc = geoRF.makeLocation(DgGeoCoord(lon, lat, false));
		dgg.convert(tloc);
		DgPolygon verts(dgg);
		dgg.setVertices(*tloc, verts, 0);
		//////
		//根据tif信息，将verts投影至平面
		const DgGeoSphRF* geoRF = dynamic_cast<const DgGeoSphRF*>(&verts.rf());
		vector<Vec2D> gridVerCor = GridVercoord(verts, coordTransInv);
		//将grid投影后的平面坐标写入shp中
		OGRLinearRing ring;
		for (int i = 0; i < gridVerCor.size(); i++)
		{
			ring.addPoint(gridVerCor[i].x, gridVerCor[i].y);
		}
		ring.closeRings();
		//环加入到polygon中
		OGRPolygon polygon;
		polygon.addRing(&ring);

		poFeature->SetGeometry(&polygon);
		poLayer->CreateFeature(poFeature);
	}
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
}

bool is_In_On_Polygon(vector<Vec2D> gridVerCor, Vec2D A)
{
	int count = 0;
	int total_pt = gridVerCor.size();
	gridVerCor.push_back(gridVerCor[0]);
	Vec2D p1, p2;
	for (int i = 0; i < total_pt; i++)
	{
		p1 = gridVerCor[i];
		/*if (i == total_pt - 1)
		{
			p2 = gridVerCor[0];
		}
		else*/
		{
			p2 = gridVerCor[i + 1];
		}
		if ((A.y >= p1.y &&  A.y <= p2.y) ||
			(A.y >= p2.y &&  A.y <= p1.y))
		{
			double t = (A.y - p1.y) / (p2.y - p1.y);
			double xt = p1.x + t * (p2.x - p1.x);
			if (A.x == xt) return 1;//点在直线上
			if (A.x < xt) ++count;
		}
	}
	if (count % 2) return 1;//点在多边形内部
	else
	{
		return 0;//点在多边形外部
	}
}

void Grid_Voronoi(GridGenParam& dp, unsigned short * imgData, int Width, int Height, double *geoTrans,
	OGRCoordinateTransformation *coordTrans, OGRCoordinateTransformation *coordTransInv, int **GLC30)
{
	//定义基础变量
	int total_grid = 0;
	int total_SHP = 0;//目前已生成的shp总数
	int compilte_cell = 1000 * total_SHP;//目前生成的gird总数=101*total_SHP
	int total_SHP_Max = 0;//允许生成shp的最大值 加上total_SHP的
	double Diff_AREA = 0.0;
	double Diff_AREA1 = 0.0;
	double Diff_AREA2 = 0.0;

						  //定义GDAL环境
	GDALDriver *poDriver;
	poDriver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
	OGRFeature *poFeature;
	OGRLayer *poLayer;
	////设置Dst的空间参考系
	OGRSpatialReference DstSPF;
	DstSPF.SetProjCS("UTM 45(WGS84) in northern hemisphere.");
	DstSPF.SetWellKnownGeogCS("WGS84");
	DstSPF.SetUTM(12, TRUE);
	string shpoutpath = "D:\\data\\TestDATA\\SHP\\USA_3_city\\YTZ_with_n12_40_class_50_local_center_Vor_" + to_string(total_SHP) + ".shp";

	// 创建属性字段
	{
		OGRFieldDefn oFieldId("GLC", OFTInteger);//编号，为整数类型
		oFieldId.SetWidth(5);
		OGRFieldDefn firstField("x", OFTReal);//编号，为浮点类型
											  //firstField.SetWidth(255);
		firstField.SetPrecision(10);
		OGRFieldDefn secondField("y", OFTReal);//编号，为浮点类型
											   //secondField.SetWidth(255);
		secondField.SetPrecision(10);
		OGRFieldDefn thirdField("10", OFTReal);//编号，为浮点类型
											   //thirdField.SetWidth(255);
		thirdField.SetPrecision(10);
		OGRFieldDefn forthField("20", OFTReal);//编号，为浮点类型
											   //forthField.SetWidth(255);
		forthField.SetPrecision(10);
		OGRFieldDefn fifthField("30", OFTReal);//编号，为浮点类型
											   //fifthField.SetWidth(255);
		fifthField.SetPrecision(10);
		OGRFieldDefn sixthField("40", OFTReal);//编号，为浮点类型
											   //sixthField.SetWidth(255);
		sixthField.SetPrecision(10);
		OGRFieldDefn seventhField("50", OFTReal);//编号，为浮点类型
												 //seventhField.SetWidth(255);
		seventhField.SetPrecision(10);
		OGRFieldDefn eighthField("60", OFTReal);//编号，为浮点类型
												//eighthField.SetWidth(255);
		eighthField.SetPrecision(10);

		OGRFieldDefn ninthField("70", OFTReal);//编号，为浮点类型
											   //ninthField.SetWidth(255);
		ninthField.SetPrecision(10);

		OGRFieldDefn tenthField("80", OFTReal);//编号，为浮点类型
											   //tenthField.SetWidth(255);
		tenthField.SetPrecision(10);

		OGRFieldDefn eleventhField("90", OFTReal);//编号，为浮点类型
												  //eleventhField.SetWidth(255);
		eleventhField.SetPrecision(10);

		OGRFieldDefn twelfthFileld("100", OFTReal);//编号，为浮点类型
												   //twelfthFileld.SetWidth(255);
		twelfthFileld.SetPrecision(10);

		////创建shp文件
		GDALDataset *poDS = poDriver->Create(shpoutpath.data(), 0, 0, 0, GDT_Unknown, NULL);

		////创建图层文件，一般为1个图层，图层中添加空间参考与几何类型
		poLayer = poDS->CreateLayer("DGGRID_T", &DstSPF, wkbPolygon, NULL);
		poLayer->CreateField(&oFieldId);
		poLayer->CreateField(&firstField);
		poLayer->CreateField(&secondField);
		poLayer->CreateField(&thirdField);
		poLayer->CreateField(&forthField);
		poLayer->CreateField(&fifthField);
		poLayer->CreateField(&sixthField);
		poLayer->CreateField(&seventhField);
		poLayer->CreateField(&eighthField);
		poLayer->CreateField(&ninthField);
		poLayer->CreateField(&tenthField);
		poLayer->CreateField(&eleventhField);
		poLayer->CreateField(&twelfthFileld);
	}
	poFeature = OGRFeature::CreateFeature(poLayer->GetLayerDefn());
	//
	double g0 = geoTrans[0],
		g1 = geoTrans[1],
		g2 = geoTrans[2],//0
		g3 = geoTrans[3],
		g4 = geoTrans[4],//0
		g5 = geoTrans[5];

	////// create the reference frames ////////

	DgRFNetwork net0;
	DgGeoSphRF geoRF(net0, dp.datum, dp.earthRadius);
	DgIDGG dgg(geoRF, dp.vert0, dp.azimuthDegs, dp.aperture, dp.actualRes,
		"DDG", dp.gridTopo, dp.projType, dp.isMixed43, dp.numAp4,
		dp.isSuperfund, dp.sfRes, dp.precision);

	// set-up to convert to degrees
	DgGeoSphDegRF deg(geoRF, geoRF.name() + "Deg");

	////open RandomPt file
	string Ptpath = "D:\\data\\TestDATA\\TXT\\USA_3_city\\YTZ_with_n12_40_class_50_local_center.TXT";
	string Boundarypath = "D:\\data\\TestDATA\\TXT\\USA_3_city\\YTZ_with_n12_40_class_50_TEST_PT.txt";
	//ofstream outfile(outpath.data()); 
	double A, B, C, D, E;
	double lon = 0.0, lat = 0.0;
	const int maxLine = 100;
	char line[maxLine];
	char line_BD[maxLine];
	DgInputStream inFile(Ptpath, "", DgBase::Fatal);
	
	inFile.getline(line, maxLine);
	while (inFile.getline(line, maxLine))
	{
		Diff_AREA = 0.0;
		int total_cell = 0;
		//读取中心点对应的shp边界，存储于BD_point中
		DgInputStream inFile_BD(Boundarypath, "", DgBase::Fatal);
		inFile_BD.getline(line_BD, maxLine);
		//存储BD的起始点
		vector<Vec2D> BD_point;
		double lon_start, lat_start;
		Vec2D B_;
		while (inFile_BD.getline(line_BD, maxLine))
		{
			int total_cell = 0;
			//读取BD的坐标A,B
			//sscanf(line, "%lf, %lf, %lf, %lf, %d, %lf, %lf",&A, &B, &C, &D, &E, &lon, &lat);
			sscanf(line_BD, ", %lf, %lf, %d, %lf, %lf", &A, &B, &C, &lon, &lat);//
			B_.x = lon;
			B_.y = lat;
			BD_point.push_back(B_);
		}
		//遍历生长元
		//sscanf(line, "%lf, %lf, %lf, %lf, %d, %lf, %lf",&A, &B, &C, &D, &E, &lon, &lat);
		sscanf(line, ", %lf, %lf, %lf, %lf, %lf", &A, &B, &C, &lon, &lat);//
		//确定点所在格网属性
		//(lon,lat)==>(Q, I, J)
		int Quam = 0, i = 0, j = 0;
		DgQ2DICoord coor;
		que * Voronoi = InitQueue();
		vector<DgQ2DICoord> tmp_Vor;
		DgLocation* addLoc = geoRF.makeLocation(DgGeoCoord(lon, lat, false));
		DgPolygon verts(dgg);
		dgg.setVertices(*addLoc, verts, dp.nDensify);
		dgg.convert(addLoc);
		DgLocation ADDLOC = *addLoc;
		const DgAddress<DgQ2DICoord> *add = (DgAddress<DgQ2DICoord>*)((ADDLOC).address());
		Quam = add->address().quadNum();
		i = add->address().coord().i();
		j = add->address().coord().j();
		coor = DgQ2DICoord(Quam, DgIVec2D(i, j));
		//将DgQ2DICoord坐标coor放入vector中
		tmp_Vor.push_back(coor);
		//元素入队
		InsertQueue(Voronoi, coor);
		while (!EmptyQueue(Voronoi))
		{
			total_cell++;
			/*if (total_cell > 8000) 
				break;*/
			DeleteQueue(Voronoi, coor);
			DgLocation tmp = *addLoc;
			((DgAddress<DgQ2DICoord>*)((&tmp)->address()))->setAddress(coor);
			dgg.setVertices(tmp, verts, dp.nDensify);
			vector<DgQ2DICoord> Vor_3edge_proxi =
				Edge_Proximity(*addLoc, coor.quadNum(), coor.coord().i(), coor.coord().j(),
					verts, dgg.maxI(), dgg.maxJ());
			//将Vor_3edge_proxi中已经计算的Vor单元删除
			bool is_ok = false;
			for (auto iter1 = Vor_3edge_proxi.cbegin(); iter1 != Vor_3edge_proxi.cend();)
			{
				is_ok = false;
				//1st 是否与已有格网重复
				if (tmp_Vor.size() == 1) break;
				auto iter2 = tmp_Vor.begin();
				while (iter2 != tmp_Vor.end())
				{
					if (*iter1 == *iter2)
					{
						is_ok = true;
						iter1 = Vor_3edge_proxi.erase(iter1);
						break;
					}
					else
					{
						++iter2;
					}
				}
				//cout << tmp_Vor.size() << endl;
				//2nd 是否在多边形内
				//将QIJ==> xy
				if (is_ok) continue;
				((DgAddress<DgQ2DICoord>*)((&tmp)->address()))->setAddress(*iter1);
				DgPolygon verts_Vor(dgg);
				dgg.setVertices(tmp, verts_Vor, dp.nDensify);
				const DgGeoSphRF* geoRF = dynamic_cast<const DgGeoSphRF*>(&verts_Vor.rf());
				vector<Vec2D> gridVerCor = GridVercoord(verts_Vor, coordTransInv);
				//利用射线法与BD判断
				int tmp_is_in_count = 0, is_in_count = 0;
				int total_gridVerCor = gridVerCor.size();
				for (; tmp_is_in_count < total_gridVerCor; tmp_is_in_count++)
				{
					if(is_In_On_Polygon(BD_point, gridVerCor[tmp_is_in_count]))
					is_in_count++;
				}
				if (!is_in_count)//若is_in_count==0，即射线法判断格网在多边形外部，则用clip的方法复查
				{
					//if (!is_in)//说明gridVerCor顶点均在BD外
					{
						iter1 = Vor_3edge_proxi.erase(iter1);
						continue;
					}
				}
				if (is_in_count != total_gridVerCor)
				{
					iter1 = Vor_3edge_proxi.erase(iter1);
					continue;
				}
				
				//if (is_in_count != total_gridVerCor)
				//{
				//	double inter_AREA = clipper_intersection(gridVerCor, BD_point);
				//	Vec2D Tri_A = gridVerCor[0];
				//	Vec2D Tri_B = gridVerCor[1];
				//	Vec2D Tri_C = gridVerCor[2];
				//	cout << is_in_count << "	" << inter_AREA << endl;
				//	//Diff_AREA += AREA - inter_AREA;
				//	if (inter_AREA >= AREA_Half)
				//	{
				//		Diff_AREA += AREA - inter_AREA;//若inter_AREA>=1/2 * AREA，则总面积会多出AREA - inter_AREA
				//		Diff_AREA1 += AREA - inter_AREA;//若inter_AREA>=1/2 * AREA，则总面积会多出AREA - inter_AREA

				//	}
				//	else
				//	{
				//		Diff_AREA -= inter_AREA;//若inter_AREA<1/2 * AREA，则总面积会减少inter_AREA
				//		Diff_AREA2 -= inter_AREA;//若inter_AREA<1/2 * AREA，则总面积会减少inter_AREA
				//		iter1 = Vor_3edge_proxi.erase(iter1);
				//		continue;
				//	}
				//}
				iter1++;
			}
			if (Vor_3edge_proxi.size() != 0)
			{
				for (auto iter1 = Vor_3edge_proxi.cbegin(); iter1 != Vor_3edge_proxi.cend(); )
				{
					//元素入队
					DgQ2DICoord tc = *iter1;

					int ci = tc.coord().i();
					int cj = tc.coord().j();
					cout << ci << "	" << cj << endl;
					//if(ci==198240 && cj==501230)
					//{
					//	DgQ2DICoord tc = *iter1;

					//	int ci = tc.coord().i();
					//	int cj = tc.coord().j();
					//	cout << ci << "	" << cj << endl;

					//	//
					//	string tmpSHP = "D:\\data\\TestDATA\\SHP\\USA_3_city\\tmp.shp";
					//	////创建shp文件
					//	GDALDataset *poDS1 = poDriver->Create(tmpSHP.data(), 0, 0, 0, GDT_Unknown, NULL);

					//	////创建图层文件，一般为1个图层，图层中添加空间参考与几何类型
					//	OGRLayer *poLayer1 = poDS1->CreateLayer("DGGRID_T", &DstSPF, wkbPolygon, NULL);
					//	OGRFeature *poFeature1 = OGRFeature::CreateFeature(poLayer1->GetLayerDefn());
					//	((DgAddress<DgQ2DICoord>*)((&tmp)->address()))->setAddress(*iter1);
					//	DgPolygon verts_Vor(dgg);
					//	dgg.setVertices(tmp, verts_Vor, dp.nDensify);
					//	const DgGeoSphRF* geoRF = dynamic_cast<const DgGeoSphRF*>(&verts_Vor.rf());
					//	vector<Vec2D> gridVerCor = GridVercoord(verts_Vor, coordTransInv);
					//	//将grid投影后的平面坐标写入shp中
					//	OGRLinearRing ring;
					//	for (int i = 0; i < gridVerCor.size(); i++)
					//	{
					//		ring.addPoint(gridVerCor[i].x, gridVerCor[i].y);
					//	}
					//	ring.closeRings();
					//	//环加入到polygon中
					//	OGRPolygon polygon;
					//	polygon.addRing(&ring);

					//	poFeature1->SetGeometry(&polygon);
					//	poLayer1->CreateFeature(poFeature1);
					//	//
					//}

					InsertQueue(Voronoi, *iter1);
					tmp_Vor.push_back(*iter1);
					iter1++;
				}
			}
			/*else
			{
				break;
			}*/

			//在删除之前将其写入SHP中
			//根据tif信息，将verts投影至平面
			const DgGeoSphRF* geoRF = dynamic_cast<const DgGeoSphRF*>(&verts.rf());
			vector<Vec2D> gridVerCor = GridVercoord(verts, coordTransInv);
			//将grid投影后的平面坐标写入shp中
			OGRLinearRing ring;
			for (int i = 0; i < gridVerCor.size(); i++)
			{
				ring.addPoint(gridVerCor[i].x, gridVerCor[i].y);
			}
			ring.closeRings();
			//环加入到polygon中
			OGRPolygon polygon;
			polygon.addRing(&ring);

			poFeature->SetGeometry(&polygon);
			poLayer->CreateFeature(poFeature);

		}
	}
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
	cout << Diff_AREA << endl;
}

void LL2Coord(DgQ2DICoord coor, double lon, double lat,
	DgGeoSphRF geoRF, DgIDGG dgg, int nDensify)
{
	DgLocation* addLoc = geoRF.makeLocation(DgGeoCoord(lon, lat, false));
	DgLocation ADDLOC = *addLoc;
	const DgAddress<DgQ2DICoord> *add = (DgAddress<DgQ2DICoord>*)((ADDLOC).address());
	int Quam = add->address().quadNum();
	int I = add->address().coord().i();
	int J = add->address().coord().j();
	coor = DgQ2DICoord(Quam, DgIVec2D(I, J));
}

void QIJ2LL(DgQ2DICoord coor, DgLocation& Loc, DgPolygon verts,
	DgGeoSphRF geoRF, DgIDGG dgg, int nDensify)
{
	//DgPolygon verts(dgg);
	//DgLocation tmp = Loc;
	((DgAddress<DgQ2DICoord>*)((&Loc)->address()))->setAddress(coor);
	dgg.setVertices(Loc, verts, nDensify);
}

bool is_Same_Side(vector<Vec2D> gridVerCor, Vec2D A, Vec2D B)
{
	Vec2D AB, AP1, AP2, AP3;
	AB.x = B.x - A.x;
	AB.y = B.y - A.y;
	AP1.x = gridVerCor[0].x - A.x;
	AP2.x = gridVerCor[1].x - A.x;
	AP3.x = gridVerCor[2].x - A.x;
	AP1.y = gridVerCor[0].y - A.y;
	AP2.y = gridVerCor[1].y - A.y;
	AP3.y = gridVerCor[2].y - A.y;
	//AB与AP1叉乘
	double ABP1 = AB.x * AP1.y - AP1.x * AB.y;
	//AB与AP2叉乘
	double ABP2 = AB.x * AP2.y - AP2.x * AB.y;
	//AB与AP3叉乘
	double ABP3= AB.x * AP3.y - AP3.x * AB.y;
	//
	if ((ABP1 >= 0 && ABP2 >= 0 && ABP3 >= 0) ||
		(ABP1 <= 0 && ABP2 <= 0 && ABP3 <= 0))
		return true;
	else
	{
		return false;
	}
}

double clipper_intersection(vector<Vec2D> grid, vector<Vec2D> BD)
{
	int SCALE_ = 10000000;
	Paths subject, clip, solution;
	ClipType ct = ctIntersection;
	FillRule fr = frNonZero;
	Clipper clipper;
	int grid_size = grid.size();
	int BD_size = BD.size();
	int i = 0;
	/************************************************************************/
	ct = ctIntersection;//ctUnion;//
	fr = frEvenOdd;//frNonZero;//
				   /************************************************************************/
	subject.resize(1);
	subject[0].resize(BD_size + 1);
	for (i = 0; i < BD_size; i++)
	{
		subject[0][i].x = (int64_t)(BD[i].x * SCALE_);
		subject[0][i].y = (int64_t)(BD[i].y * SCALE_);
	}
	subject[0][BD_size].x = (int64_t)(BD[0].x * SCALE_);
	subject[0][BD_size].y = (int64_t)(BD[0].y * SCALE_);
	clip.resize(1);
	clip[0].resize(grid_size + 1);
	for (i = 0; i < grid_size; i++)
	{
		clip[0][i].x = (int64_t)(grid[i].x * SCALE_);
		clip[0][i].y= (int64_t)(grid[i].y * SCALE_);
	}
	clip[0][grid_size].x = (int64_t)(grid[0].x * SCALE_);
	clip[0][grid_size].y = (int64_t)(grid[0].y * SCALE_);

	clipper.Clear();
	clipper.AddPaths(subject, ptSubject);
	clipper.AddPaths(clip, ptClip);

	clipper.Execute(ct, solution, fr);//solution没有值 定义可能冲突了

	/*std::cout << "Result Polygon Size:" << solution.size() << std::endl;
	for (int i = 0; i < solution.size(); i++)
	{
		std::cout << "Polygon: " << i + 1 << std::endl;
		std::cout << "Point Number: " << solution[i].size() << std::endl;
		std::cout << "x\ty" << std::endl;
		for (int j = 0; j < solution[i].size(); j++)
		{
			std::cout << solution[i][j].x << "\t" << solution[i][j].y << std::endl;
		}
	}*/
	if (solution.size())
	{
		//cout << "area==" << Area(solution[0])/* / SCALE_ / SCALE_*/ << endl;
		return  (Area(solution[0]) / SCALE_ / SCALE_);//此处不对，
		//return 1;
		//相交有很多种，重合或叠加也是相交的一种
	}
	else
	{
		return 0;//表示三角形和pixel不相交
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
void SplitString_OURSGRID(string line, int **STR2I)
{
	int *str2i = new int[6];
	int n = line.size();
	for (int i = 0; i < n; ++i) {
		if (line[i] == ',') {
			line[i] = ' ';
		}
	}
	istringstream out(line);
	string str;
	int pos = 0;
	int cols = 18;
	int key = 0;
	while (pos < cols) {
		out >> str;
		istringstream tmpstr(str);
		if (pos == 6)//Q
		{
			str2i[0] = stoi(str);
			pos++;
			continue;
		}
		if (pos == 7)//I
		{
			str2i[1] = stoi(str);
			pos++;
			continue;
		}
		if (pos == 8)//J
		{
			str2i[2] = stoi(str);
			pos++;
			continue;
		}
		if (pos == 14)//wieght
		{
			str2i[3] = stoi(str);
			pos++;
			continue;
		}
		if (pos == 17)//PolyID
		{
			str2i[4] = stoi(str);
			pos++;
			break;
		}
		if (pos == 12)//type
		{
			str2i[5] = stoi(str);
			pos++;
			continue;
		}
		pos++;
	}
	*STR2I = str2i;
}
void SplitString(string line, int **STR2I)
{
	int *str2i = new int[6];
	int n = line.size();
	for (int i = 0; i < n; ++i) {
		if (line[i] == ',') {
			line[i] = ' ';
		}
	}
	istringstream out(line);
	string str;
	int pos = 0;
	int cols = 18;
	int key = 0;
	while (pos < cols) {
		out >> str;
		istringstream tmpstr(str);
		if (pos == 6)
		{
			str2i[0] = stoi(str);
			pos++;
			continue;
		}
		if (pos == 7)
		{
			str2i[1] = stoi(str);
			pos++;
			continue;
		}
		if (pos == 8)
		{
			str2i[2] = stoi(str);
			pos++;
			continue;
		}
		if (pos == 14)
		{
			str2i[3] = stoi(str);
			pos++;
			continue;
		}
		if (pos == 15)
		{
			str2i[4] = stoi(str);
			pos++;
			break;
		}
		if (pos == 12)
		{
			str2i[5] = stoi(str);
			pos++;
			continue;
		}
		pos++;
	}
	*STR2I = str2i;
}
void getEdgeGridSHPwithType(GridGenParam& dp, char *SHP_path, char *txt_path, char *outshp_path)
{
	GDALDataset   *poDS, *poDS1;
	OGRLayer  *poLayer, *poLayer1;
	OGRFeature *poFeature, *poFeature1;
	Poly_BDRasters * PBDRs, *PBDRs_Line1, *PBDRs_Line2;
	OGRGeometry *poGeometry;
	OGRwkbGeometryType pGeoType;
	OGRPolygon *BDPolygon;
	OGRLinearRing  grid_ring;
	OGRSpatialReference DstSPF;
	vector<Vec2D> gridVerCor;
	vector<DgQ2DICoord> Vctmp;

	//定义GDAL环境
	GDALDriver *poDriver;
	poDriver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
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
	DgQ2DICoord coor;
	const DgAddress<DgQ2DICoord> *add;
	DgLocation ADDLOC;
	//读取shp文件
	poDS = (GDALDataset*)GDALOpenEx(SHP_path, GDAL_OF_VECTOR, NULL, NULL, NULL);
	//创建shp文件
	poDS1 = poDriver->Create(outshp_path, 0, 0, 0, GDT_Unknown, NULL);
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
	OGRFieldDefn oFieldId("Type", OFTInteger);//多边形属性
	oFieldId.SetWidth(5);
	OGRFieldDefn firstField("Quam", OFTInteger);//编号，为浮点类型
	firstField.SetPrecision(2);
	OGRFieldDefn secondField("I", OFTInteger);//编号，为浮点类型
	secondField.SetPrecision(255);
	OGRFieldDefn thirdField("J", OFTInteger);//编号，为浮点类型
	thirdField.SetPrecision(255);
	OGRFieldDefn forthField("PolyID", OFTInteger);//编号，为浮点类型
	forthField.SetPrecision(255);
	OGRFieldDefn fifthField("Weight", OFTInteger);//编号，为浮点类型
	fifthField.SetPrecision(2);
	OGRFieldDefn sixthField("Class", OFTInteger);//编号，为浮点类型
	sixthField.SetPrecision(5);

	////创建图层文件，一般为1个图层，图层中添加空间参考与几何类型
	poLayer1 = poDS1->CreateLayer("DGGRID_T", &DstSPF, wkbPolygon, NULL);
	poLayer1->CreateField(&oFieldId);
	poLayer1->CreateField(&firstField);
	poLayer1->CreateField(&secondField);
	poLayer1->CreateField(&thirdField);
	poLayer1->CreateField(&forthField);
	poLayer1->CreateField(&fifthField);
	poLayer1->CreateField(&sixthField);

	poFeature1 = OGRFeature::CreateFeature(poLayer1->GetLayerDefn());
	addLoc = geoRF.makeLocation(DgGeoCoord(0, 0, false));
	dgg.setVertices(*addLoc, verts, dp.nDensify);
	dgg.convert(addLoc);
	ADDLOC = *addLoc;
	add = (DgAddress<DgQ2DICoord>*)((ADDLOC).address());

	//遍历TXT
	int *str2i = new int[6];
	std::ifstream in(txt_path);//
	string line;
	getline(in, line);
	//cout << line;
	while (getline(in, line))
	{
		//处理txt中的位置与qij等的对应关系。
		//FID 0 X0 Y0 X1 Y1 Q I J 0 polyid polyid+1 type area weight polyid
		//SplitString_OURSGRID(line, &str2i);
		SplitString(line, &str2i);
		//cout << str2i[0] << "," << str2i[1] << "," << str2i[2] << "," << str2i[3] << "," << str2i[4] << "," << str2i[5] << endl;

		if (str2i[5] != _TYPE_) continue;

		int Q = str2i[0], I = str2i[1], J = str2i[2];
		//int Q = 0, I = 0, J = 0;
		coor = DgQ2DICoord(Q, DgIVec2D(I, J));
		((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coor);
		dgg.setVertices(ADDLOC, verts, dp.nDensify);
		gridVerCor = GridVercoord(verts, coordTransInv);
		//将gridVerCor投影后的平面坐标写入shp中
		grid_ring.empty();
		for (int k = 0; k < 3; k++)
		{
			grid_ring.addPoint(gridVerCor[k].x, gridVerCor[k].y);
		}
		grid_ring.closeRings();
		//环加入到paddRingolygon中
		OGRPolygon polygon;
		polygon.addRing(&grid_ring);
		poFeature1->SetGeometry(&polygon);

		poFeature1->SetField(1, Q);
		poFeature1->SetField(2, I);
		poFeature1->SetField(3, J);
		poFeature1->SetField(4, str2i[3]);//polyID
		poFeature1->SetField(5, str2i[4]);//weight 
		poFeature1->SetField(6, _TYPE_);
		poLayer1->CreateFeature(poFeature1);
	}
	GDALClose(poDS);
	//GDALClose(poDS1);
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
}

bool PointInPolygon(double x, double y, Paths subject, int pathID)
{
	int count = 0;
	int total_pt = subject[pathID].size();
	Vec2D p1, p2;
	for (int i = 0; i < total_pt - 1; i++)
	{
		p1.x = subject[pathID][i].x;
		p1.y = subject[pathID][i].y;
		p2.x = subject[pathID][i+1].x;
		p2.y = subject[pathID][i+1].y;
		if ((y >= p1.y &&  y <= p2.y) ||
			(y >= p2.y &&  y <= p1.y))
		{
			double t = (y - p1.y) / (p2.y - p1.y);
			double xt = p1.x + t * (p2.x - p1.x);
			if (x == xt) return 1;//点在直线上
			if (x < xt) ++count;
		}
	}
	if (count % 2) return 1;//点在多边形内部
	else
	{
		return 0;//点在多边形外部
	}
}

void BDTracing(GridGenParam& dp, char *SHP_path, string polyID, string outshp_path, string outtxt_path)
{
	//基础变量
	int POLYCOUNT = 1;
	bool Orti = true;
	int i, i_1, j, k, tmp;
	int polycount = 0;
	int BDpointcount = 0;
	int stopcount = 0;
	double lon, lat;
	double gx, gy, gix, giy, gi_1x, gi_1y;
	double sum_x = 0, sum_y = 0;
	double angle;
	double CX[2] = { 0 }, CY[2] = { 0 };
	double dist = 0.0, min_dist = 10000000;
	double ai = 0.0;
	double sum_ai = 0;
	int edgeCount = 0;
	int SHPID = 0;

	GDALDataset   *poDS, *poDS1, *poDS2, *poDS3, *poDS4;
	OGRLayer  *poLayer, *poLayer1, *poLayer2, *poLayer3, *poLayer4;
	OGRFeature *poFeature, *poFeature1, *poFeature2, *poFeature3, *poFeature4;
	Poly_BDRasters * PBDRs, *PBDRs_Line1, *PBDRs_Line2;
	OGRGeometry *poGeometry;
	OGRwkbGeometryType pGeoType;
	OGRPolygon *BDPolygon;
	OGRLinearRing  grid_ring;
	OGRSpatialReference DstSPF;
	vector<Vec2D> gridVerCor, gridVerCor_tmp;
	vector<DgQ2DICoord> Vctmp;

	//定义GDAL环境
	GDALDriver *poDriver;
	poDriver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");

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
	DgQ2DICoord coor_i, coor_j, coor_k, coor_tmp, coor_;
	const DgAddress<DgQ2DICoord> *add;
	DgLocation ADDLOC;
	vector<DgQ2DICoord> coor;
	//读取shp文件
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
	// 创建属性字段
	OGRFieldDefn oFieldId("Type", OFTInteger);//多边形属性
	oFieldId.SetWidth(5);
	OGRFieldDefn firstField("PolyID", OFTInteger);//编号，为浮点类型
	firstField.SetWidth(255);
	OGRFieldDefn secondField("Quam", OFTInteger);//编号，为浮点类型
	secondField.SetWidth(2);
	//fifthField.SetPrecision(2);
	OGRFieldDefn thirdField("I", OFTInteger);//编号，为浮点类型
	thirdField.SetWidth(255);
	//sixthField.SetPrecision(10);
	OGRFieldDefn forthField("J", OFTInteger);//编号，为浮点类型
	forthField.SetWidth(255);
	////创建shp文件
	string SHPIDPATH = outshp_path + to_string(SHPID) + "_18.shp";
	poDS1 = poDriver->Create(SHPIDPATH.data(), 0, 0, 0, GDT_Unknown, NULL);
	////创建图层文件，一般为1个图层，图层中添加空间参考与几何类型
	poLayer1 = poDS1->CreateLayer("DGGRID_T", &DstSPF, wkbPolygon, NULL);
	poLayer1->CreateField(&oFieldId);
	poLayer1->CreateField(&firstField);
	poLayer1->CreateField(&secondField);
	poLayer1->CreateField(&thirdField);
	poLayer1->CreateField(&forthField);
	poFeature1 = OGRFeature::CreateFeature(poLayer1->GetLayerDefn());
	
	while ((poFeature = poLayer->GetNextFeature()) != NULL)
	{
		if (polycount >= 50) break;//POLYCOUNT
		cout << polycount << endl;
		polycount++;
		//get feature的type==GRIDCODE
		int n = poFeature->GetFieldCount();//get the number of field
		int GRIDCODE = poFeature->GetFieldAsDouble(1);
		if (GRIDCODE != _TYPE_) continue;
		int polyID = poFeature->GetFieldAsDouble(0);
		string TxtPath = outtxt_path + to_string(polyID) + "_" + to_string(NLEVEL) + ".txt";
		ofstream outfile;
		outfile.open(TxtPath);
		poGeometry = ConvertPolygonToPolyline(poFeature->GetGeometryRef());
		pGeoType = poGeometry->getGeometryType();
		
		if (pGeoType == wkbMultiLineString)
		{
			OGRMultiLineString *poMultiLineString = (OGRMultiLineString*)poGeometry;
			int nGeoCount = poMultiLineString->getNumGeometries();
			OGRGeometry *poLineGeometry;
			//for (int iLine = 0; iLine < 1; iLine++)
			for (int iLine = 0; iLine < nGeoCount; iLine++)
			{
				poLineGeometry = poMultiLineString->getGeometryRef(iLine);
				OGRLineString *boundary = (OGRLineString*)poLineGeometry;
				//Paths subject = PreDefineSubject((OGRLinearRing*)poLineGeometry);

				BDpointcount = boundary->getNumPoints();
				bool is_turn_P = false;
				for (int i = 0; i < BDpointcount - 1; i++)
				{
					//init the structure of arc
					vector<DgQ2DICoord> tmpV;
					//cout << i << endl;
					//if (i == 1) break;
					//1. judge the distance between the ith point and the (i+1)th one (this step is not performed)
					//2. make sure the starting grid and the end grid
					//3. trace grid based on the starting grid
					//4. generate the shape file
					//1. 
					i_1 = i + 1;
					gix = boundary->getX(i);
					giy = boundary->getY(i);
					gi_1x = boundary->getX(i_1);
					gi_1y = boundary->getY(i_1);
					//double dist = sqrtf((gix - gi_1x) * (gix - gi_1x) + (giy - gi_1y) * (giy - gi_1y));
					//if (dist < DIST) continue;
					//2.
					{
						gx = gix;
						gy = giy;
						coordTrans->Transform(1, &gx, &gy);
						//make sure the coor_i according to the latitude and lontitude of (gx, gy)
						lon = gx;
						lat = gy;
						addLoc = geoRF.makeLocation(DgGeoCoord(lon, lat, false));
						dgg.setVertices(*addLoc, verts, dp.nDensify);
						dgg.convert(addLoc);
						ADDLOC = *addLoc;
						add = (DgAddress<DgQ2DICoord>*)((ADDLOC).address());
						coor_i = DgQ2DICoord(add->address().quadNum(),
							DgIVec2D(add->address().coord().i(),
								add->address().coord().j()));
						dgg.setVertices(ADDLOC, verts, dp.nDensify);
						gridVerCor = GridVercoord(verts, coordTransInv);
						//make sure the i+1 grid
						gx = gi_1x;
						gy = gi_1y;
						coordTrans->Transform(1, &gx, &gy);
						lon = gx;
						lat = gy;
						addLoc = geoRF.makeLocation(DgGeoCoord(lon, lat, false));
						dgg.setVertices(*addLoc, verts, dp.nDensify);
						dgg.convert(addLoc);
						ADDLOC = *addLoc;
						add = (DgAddress<DgQ2DICoord>*)((ADDLOC).address());
						coor_tmp = DgQ2DICoord(add->address().quadNum(),
							DgIVec2D(add->address().coord().i(),
								add->address().coord().j()));
						//is coor_i and coor_tmp the same.
						tmpV.push_back(coor_i);
						//calcuate the intersection area
						((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coor_i);
						dgg.setVertices(ADDLOC, verts, dp.nDensify);
						/*gridVerCor = GridVercoord(verts, coordTransInv);
						ai = clipper_intersection_area2(gridVerCor[0].x, gridVerCor[0].y,
							gridVerCor[1].x, gridVerCor[1].y, gridVerCor[2].x, gridVerCor[2].y, subject);*/

						if (coor_i != coor_tmp)
						{
							tmpV.push_back(coor_tmp);
							//calcuate the intersection area
							((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coor_i);
							dgg.setVertices(ADDLOC, verts, dp.nDensify);
							/*gridVerCor = GridVercoord(verts, coordTransInv);
							ai = clipper_intersection_area2(gridVerCor[0].x, gridVerCor[0].y,
								gridVerCor[1].x, gridVerCor[1].y, gridVerCor[2].x, gridVerCor[2].y, subject);*/

						}
						else continue;
					}
					//calculate the translation and rotation matrix for AB
					double Bx = gi_1x - gix, By = gi_1y - giy;
					double AB = sqrtf(Bx * Bx + By * By);
					double sint = By / AB;
					double cost = Bx / AB;
					//3.
					while (1)
					{
						/*if (i == 3 && PBDRs_Line1[0].BDRaster_count == 5458)
						cout << endl;*/
						//break;
						((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coor_i);
						dgg.setVertices(ADDLOC, verts, dp.nDensify);
						//get the grid orientation of coor_i
						//Orti = PolygonOrti_(coor_i, ADDLOC, dgg, dp.nDensify);
						Orti = PolygonOrti(gridVerCor);
						coor = Edge_Proximity(coor_i, Orti, dgg.maxI(), dgg.maxJ());
						//delete the same grid
						Delete_same(&coor, tmpV);
						//make sure the second grid
						if (tmpV.size() == 2)
						{
							//1.does the grrid intersect with AB
							auto iter = coor.begin();
							while (iter != coor.end())
							{
								coor_ = *iter;
								((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coor_);
								dgg.setVertices(ADDLOC, verts, dp.nDensify);
								gridVerCor = GridVercoord(verts, coordTransInv);
								/*cout << coor_ << endl;
								cout << std::setprecision(10) << gridVerCor[0].x << "	" << gridVerCor[0].y << "	" << endl;
								cout << std::setprecision(10) << gridVerCor[1].x << "	" << gridVerCor[1].y << "	" << endl;
								cout << std::setprecision(10) << gridVerCor[2].x << "	" << gridVerCor[2].y << "	" << endl;
								*/
								if (!Line_Grid_intersection(gridVerCor, gix, giy, gi_1x, gi_1y, cost, sint, AB))
								{
									iter = coor.erase(iter);
									continue;
								}
								else ++iter;
							}
							for (k = 0; k < coor.size(); k++)
							{
								coor_ = coor[k];
								((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coor_);
								dgg.setVertices(ADDLOC, verts, dp.nDensify);
								gridVerCor_tmp = GridVercoord(verts, coordTransInv);
								//calculate the central coordinate of coor
								Vec2D gc = Grid_center(gridVerCor_tmp);
								dist = PointsDist(gc.x, gc.y, gi_1x, gi_1y);
								if (min_dist > dist)
								{
									min_dist = dist;
									gridVerCor = gridVerCor_tmp;
									tmp = k;
								}
							}
							if (coor.size() == 0) break;
							coor_j = coor[tmp];
							//cout << coor_j << endl;
						}
						else
						{
							min_dist = 1000000000000;
							//1.does the grrid intersect with AB
							auto iter = coor.begin();
							while (iter != coor.end())
							{
								coor_ = *iter;
								((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coor_);
								dgg.setVertices(ADDLOC, verts, dp.nDensify);
								gridVerCor = GridVercoord(verts, coordTransInv);
								/*cout << coor_ << endl;
								cout << std::setprecision(10) << gridVerCor[0].x << "	" << gridVerCor[0].y << "	" << endl;
								cout << std::setprecision(10) << gridVerCor[1].x << "	" << gridVerCor[1].y << "	" << endl;
								cout << std::setprecision(10) << gridVerCor[2].x << "	" << gridVerCor[2].y << "	" << endl;*/
								if (!Line_Grid_intersection(gridVerCor, gix, giy, gi_1x, gi_1y, cost, sint, AB))

								{
									iter = coor.erase(iter);
									continue;
								}
								else ++iter;
							}
							
							if (coor.size() == 0) break;
							for (k = 0; k < coor.size(); k++)
							{
								coor_ = coor[k];
								((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coor_);
								dgg.setVertices(ADDLOC, verts, dp.nDensify);
								gridVerCor_tmp = GridVercoord(verts, coordTransInv);
								//calculate the central coordinate of coor
								Vec2D gc = Grid_center(gridVerCor_tmp);
								dist = PointsDist(gc.x, gc.y, gi_1x, gi_1y);
								if (min_dist > dist)
								{
									min_dist = dist;
									gridVerCor = gridVerCor_tmp;
									tmp = k;
								}
							}
							coor_j = coor[tmp];
						}
						coor_i = coor_j;
						tmpV.push_back(coor_i);
						//calculate the intersection area
						((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coor_i);
						dgg.setVertices(ADDLOC, verts, dp.nDensify);
						gridVerCor = GridVercoord(verts, coordTransInv);
						/*ai = clipper_intersection_area2(gridVerCor[0].x, gridVerCor[0].y,
							gridVerCor[1].x, gridVerCor[1].y, gridVerCor[2].x, gridVerCor[2].y, subject);*/

						if (coor_tmp.coord().i() == coor_i.coord().i() &&
							coor_tmp.coord().j() == coor_i.coord().j())
							break;
					}
					for (j = 0; j < tmpV.size(); j++)//
					{
						coor_tmp = tmpV[j];
						/*ai = PBDRs_Line1[0].Ratio[j];
						sum_ai += ai;*/
						((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coor_tmp);
						dgg.setVertices(ADDLOC, verts, dp.nDensify);
						vector<Vec2D> gridVerCor = GridVercoord(verts, coordTransInv);
						//cout << coor_tmp << endl;
						//cout << verts << endl;
						//output result
						grid_ring.empty();
						for (k = 0; k < 3; k++)
						{
							grid_ring.addPoint(gridVerCor[k].x, gridVerCor[k].y);
						}
						edgeCount++;
						grid_ring.closeRings();
						OGRPolygon polygon;
						polygon.addRing(&grid_ring);
						if (edgeCount % 1000000 == 0)
						{
							//destroy feature
							OGRFeature::DestroyFeature(poFeature1);
							SHPID++;
							string SHPIDPATH = outshp_path + to_string(SHPID) + "_18.shp";
							//crate shp文件
							poDS1 = poDriver->Create(SHPIDPATH.data(), 0, 0, 0, GDT_Unknown, NULL);
							poLayer1 = poDS1->CreateLayer("DGGRID_T", &DstSPF, wkbPolygon, NULL);
							poLayer1->CreateField(&oFieldId);
							poLayer1->CreateField(&firstField);
							poLayer1->CreateField(&secondField);
							poLayer1->CreateField(&thirdField);
							poLayer1->CreateField(&forthField);
							poFeature1 = OGRFeature::CreateFeature(poLayer1->GetLayerDefn());
						}
						poFeature1->SetGeometry(&polygon);
						poFeature1->SetField(1, polyID);
						poFeature1->SetField(2, coor_tmp.quadNum());
						poFeature1->SetField(3, coor_tmp.coord().i());
						poFeature1->SetField(4, coor_tmp.coord().j());
						poLayer1->CreateFeature(poFeature1);
						outfile << polyID << "," << coor_tmp.quadNum() << "," << coor_tmp.coord().i() << "," << coor_tmp.coord().j() << endl;
					}
					vector<DgQ2DICoord>().swap(tmpV);
				}
			}
		}
		else
		{
			OGRLinearRing *boundary = (OGRLinearRing*)poGeometry;
			//Paths subject = PreDefineSubject((OGRLinearRing*)poGeometry);
			BDpointcount = boundary->getNumPoints();
			bool is_turn_P = false;
			for (int i = 0; i < BDpointcount - 1; i++)
			{
				//init the structure of arc
				vector<DgQ2DICoord> tmpV;
 				//cout << i << endl;
				//if (i == 1) break;
				//1. judge the distance between the ith point and the (i+1)th one (this step is not performed)
				//2. make sure the starting grid and the end grid
				//3. trace grid based on the starting grid
				//4. generate the shape file
				//1. 
				i_1 = i + 1;
				gix = boundary->getX(i/* + 1144*/);
				giy = boundary->getY(i /*+ 1144*/);
				gi_1x = boundary->getX(i_1 /*+ 1144*/);
				gi_1y = boundary->getY(i_1 /*+ 1144*/);
				//double dist = sqrtf((gix - gi_1x) * (gix - gi_1x) + (giy - gi_1y) * (giy - gi_1y));
				//if (dist < DIST) continue;
				//2.
				{
					gx = gix;
					gy = giy;
					coordTrans->Transform(1, &gx, &gy);
					//make sure the grid-coor_i
					lon = gx;
					lat = gy;
					addLoc = geoRF.makeLocation(DgGeoCoord(lon, lat, false));
					dgg.setVertices(*addLoc, verts, dp.nDensify);
					dgg.convert(addLoc);
					ADDLOC = *addLoc;
					add = (DgAddress<DgQ2DICoord>*)((ADDLOC).address());
					coor_i = DgQ2DICoord(add->address().quadNum(),
						DgIVec2D(add->address().coord().i(),
							add->address().coord().j()));
					dgg.setVertices(ADDLOC, verts, dp.nDensify);
					gridVerCor = GridVercoord(verts, coordTransInv);
					//make sure the grid of the i+1 point
					gx = gi_1x;
					gy = gi_1y;
					coordTrans->Transform(1, &gx, &gy);
					lon = gx;
					lat = gy;
					addLoc = geoRF.makeLocation(DgGeoCoord(lon, lat, false));
					dgg.setVertices(*addLoc, verts, dp.nDensify);
					dgg.convert(addLoc);
					ADDLOC = *addLoc;
					add = (DgAddress<DgQ2DICoord>*)((ADDLOC).address());
					coor_tmp = DgQ2DICoord(add->address().quadNum(),
						DgIVec2D(add->address().coord().i(),
							add->address().coord().j()));
					//is the coor_i and coor_tmp same
					//PBDRs_Line1[0].polyID = 0;
					tmpV.push_back(coor_i);
					//calculate the intersection area
					((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coor_i);
					dgg.setVertices(ADDLOC, verts, dp.nDensify);
					/*gridVerCor = GridVercoord(verts, coordTransInv);
					ai = clipper_intersection_area2(gridVerCor[0].x, gridVerCor[0].y,
					gridVerCor[1].x, gridVerCor[1].y, gridVerCor[2].x, gridVerCor[2].y, subject);*/
					if (coor_i != coor_tmp)
					{
						tmpV.push_back(coor_tmp);
						//calculate the intersection area
						((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coor_i);
						dgg.setVertices(ADDLOC, verts, dp.nDensify);
						/*gridVerCor = GridVercoord(verts, coordTransInv);
						ai = clipper_intersection_area2(gridVerCor[0].x, gridVerCor[0].y,
						gridVerCor[1].x, gridVerCor[1].y, gridVerCor[2].x, gridVerCor[2].y, subject);*/
					}
					else continue;
				}
				//calculate the translation and rotation matrix
				double Bx = gi_1x - gix, By = gi_1y - giy;
				double AB = sqrtf(Bx * Bx + By * By);
				double sint = By / AB;
				double cost = Bx / AB;
				//3.
				while (1)
				{
					/*if (i == 3 && PBDRs_Line1[0].BDRaster_count == 5458)
					cout << endl;*/
					//break;
					((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coor_i);
					dgg.setVertices(ADDLOC, verts, dp.nDensify);
					//claculate the orientation of coor_i
					//Orti = PolygonOrti_(coor_i, ADDLOC, dgg, dp.nDensify);
					Orti = PolygonOrti(gridVerCor);
					coor = Edge_Proximity(coor_i, Orti, dgg.maxI(), dgg.maxJ());
					//delete the same grid
					Delete_same(&coor, tmpV);
					//make sure the 2nd grid
					if (tmpV.size() == 2)
					{
						//1.does the grid intersect with the line AB
						auto iter = coor.begin();
						while (iter != coor.end())
						{
							coor_ = *iter;
							((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coor_);
							dgg.setVertices(ADDLOC, verts, dp.nDensify);
							gridVerCor = GridVercoord(verts, coordTransInv);
							//cout << coor_ << endl;
							/*cout << std::setprecision(10) << gridVerCor[0].x << "	" << gridVerCor[0].y << "	" << endl;
							cout << std::setprecision(10) << gridVerCor[1].x << "	" << gridVerCor[1].y << "	" << endl;
							cout << std::setprecision(10) << gridVerCor[2].x << "	" << gridVerCor[2].y << "	" << endl;
							*/
							if (!Line_Grid_intersection(gridVerCor, gix, giy, gi_1x, gi_1y, cost, sint, AB))
							{
								iter = coor.erase(iter);
								continue;
							}
							else ++iter;
						}
						for (k = 0; k < coor.size(); k++)
						{
							coor_ = coor[k];
							((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coor_);
							dgg.setVertices(ADDLOC, verts, dp.nDensify);
							gridVerCor_tmp = GridVercoord(verts, coordTransInv);
							//calculte the central coordinate 
							Vec2D gc = Grid_center(gridVerCor_tmp);
							dist = PointsDist(gc.x, gc.y, gi_1x, gi_1y);
							if (min_dist > dist)
							{
								min_dist = dist;
								gridVerCor = gridVerCor_tmp;
								tmp = k;
							}
						}
						if (coor.size() == 0) break;
						coor_j = coor[tmp];
						//cout << coor_j << endl;
					}
					else
					{
						min_dist = 1000000000000;
						//1.does the grid intersect with the line AB
						auto iter = coor.begin();
						while (iter != coor.end())
						{
							coor_ = *iter;
							((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coor_);
							dgg.setVertices(ADDLOC, verts, dp.nDensify);
							gridVerCor = GridVercoord(verts, coordTransInv);
							//cout << coor_ << endl;
							/*cout << std::setprecision(10) << gridVerCor[0].x << "	" << gridVerCor[0].y << "	" << endl;
							cout << std::setprecision(10) << gridVerCor[1].x << "	" << gridVerCor[1].y << "	" << endl;
							cout << std::setprecision(10) << gridVerCor[2].x << "	" << gridVerCor[2].y << "	" << endl;*/
							if (!Line_Grid_intersection(gridVerCor, gix, giy, gi_1x, gi_1y, cost, sint, AB))

							{
								iter = coor.erase(iter);
								continue;
							}
							else ++iter;
						}
						if (coor.size() == 0) break;
						for (k = 0; k < coor.size(); k++)
						{
							coor_ = coor[k];
							((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coor_);
							dgg.setVertices(ADDLOC, verts, dp.nDensify);
							gridVerCor_tmp = GridVercoord(verts, coordTransInv);
							Vec2D gc = Grid_center(gridVerCor_tmp);
							dist = PointsDist(gc.x, gc.y, gi_1x, gi_1y);
							if (min_dist > dist)
							{
								min_dist = dist;
								gridVerCor = gridVerCor_tmp;
								tmp = k;
							}
						}
						coor_j = coor[tmp];
					}
					coor_i = coor_j;
					tmpV.push_back(coor_i);
					((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coor_i);
					dgg.setVertices(ADDLOC, verts, dp.nDensify);
					gridVerCor = GridVercoord(verts, coordTransInv);
					/*ai = clipper_intersection_area2(gridVerCor[0].x, gridVerCor[0].y,
					gridVerCor[1].x, gridVerCor[1].y, gridVerCor[2].x, gridVerCor[2].y, subject);*/
					if (coor_tmp.coord().i() == coor_i.coord().i() &&
						coor_tmp.coord().j() == coor_i.coord().j())
						break;
				}
				for (j = 0; j < tmpV.size(); j++)//
				{
					edgeCount++;
					coor_tmp = tmpV[j];
					/*ai = PBDRs_Line1[0].Ratio[j];
					sum_ai += ai;*/
					((DgAddress<DgQ2DICoord>*)((&ADDLOC)->address()))->setAddress(coor_tmp);
					dgg.setVertices(ADDLOC, verts, dp.nDensify);
					vector<Vec2D> gridVerCor = GridVercoord(verts, coordTransInv);
					//cout << coor_tmp << endl;
					//cout << verts << endl;
					//output the result into a shape file
					grid_ring.empty();
					for (k = 0; k < 3; k++)
					{
						grid_ring.addPoint(gridVerCor[k].x, gridVerCor[k].y);
					}
					grid_ring.closeRings();
					OGRPolygon polygon;
					polygon.addRing(&grid_ring);
					if (edgeCount % 1000000 == 0)
					{

						//Destroyfeature
						OGRFeature::DestroyFeature(poFeature1);
						GDALClose(poDS1);
						SHPID++;
						string SHPIDPATH = outshp_path + to_string(SHPID) + "_18.shp";
						////Create shape file
						poDS1 = poDriver->Create(SHPIDPATH.data(), 0, 0, 0, GDT_Unknown, NULL);
						poLayer1 = poDS1->CreateLayer("DGGRID_T", &DstSPF, wkbPolygon, NULL);
						poLayer1->CreateField(&oFieldId);
						poLayer1->CreateField(&firstField);
						poLayer1->CreateField(&secondField);
						poLayer1->CreateField(&thirdField);
						poLayer1->CreateField(&forthField);
						poFeature1 = OGRFeature::CreateFeature(poLayer1->GetLayerDefn());
					}
					poFeature1->SetGeometry(&polygon);
					poFeature1->SetField(1, polyID);
					poFeature1->SetField(2, coor_tmp.quadNum());
					poFeature1->SetField(3, coor_tmp.coord().i());
					poFeature1->SetField(4, coor_tmp.coord().j());
					poLayer1->CreateFeature(poFeature1);

					outfile << polyID <<","<< coor_tmp.quadNum()<<","<< coor_tmp.coord().i()<<","<< coor_tmp.coord().j()<<endl;
				}
				vector<DgQ2DICoord>().swap(tmpV);
			}
		}
		outfile.close();
		//DestroyFeature
		OGRFeature::DestroyFeature(poFeature);
	}
	GDALClose(poDS);
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
}

double Angle(double pix, double piy, double pjx, double pjy)//pi-->pj
{
	//以pt为基准,计算pj--pi--pt的夹角
	double ptx = pix;
	double pty = piy;
	float theta = atan2(pix - ptx, piy - pty) - atan2(pjx - ptx, pjy - pty);
	if (theta > M_PI)
		theta -= 2 * M_PI;
	if (theta < -M_PI)
		theta += 2 * M_PI;

	theta = theta * 180.0 / M_PI;
	return theta;

}
vector<DgQ2DICoord> BDT_EdgeADJ(DgQ2DICoord cood_i, bool Orti, double Angle)
{//先不考虑格网跨边界的情况
	vector<DgQ2DICoord> coor;
	DgQ2DICoord coor_j, coor_k;
	int Quam = cood_i.quadNum();
	int i = cood_i.coord().i();
	int j = cood_i.coord().j();
	if (Orti == 1)//上三角
	{
		if(30 <= Angle && Angle <= 90)
			coor_j = DgQ2DICoord(Quam, DgIVec2D(i - 1, j - 1));//右
		else if(150 <= Angle && Angle <= 210)
			coor_j = DgQ2DICoord(Quam, DgIVec2D(i, j - 1));//下
		else if(270 <= Angle && Angle <=330)
			coor_j = DgQ2DICoord(Quam, DgIVec2D(i, j + 1));//左
		else if (Angle > 330)//左、右
		{
			coor_j = DgQ2DICoord(Quam, DgIVec2D(i, j + 1));//左
			coor_k = DgQ2DICoord(Quam, DgIVec2D(i - 1, j - 1));//右
		}
		else if (90 < Angle && Angle < 150)
		{
			coor_j = DgQ2DICoord(Quam, DgIVec2D(i, j - 1));//下
			coor_k = DgQ2DICoord(Quam, DgIVec2D(i - 1, j - 1));//右
		}
		else if (210 < Angle && Angle < 270)
		{
			coor_j = DgQ2DICoord(Quam, DgIVec2D(i, j + 1));//左
			coor_k = DgQ2DICoord(Quam, DgIVec2D(i, j - 1));//下
		}
	}
	else//下三角
	{
		if (Angle >= 330)
			coor_j = DgQ2DICoord(Quam, DgIVec2D(i, j + 1));//上
		else if(90 <= Angle && Angle <= 150)
			coor_j = DgQ2DICoord(Quam, DgIVec2D(i, j - 1));//右
		else if (210 <= Angle && Angle <= 270)
			coor_j = DgQ2DICoord(Quam, DgIVec2D(i + 1, j + 1));//左
		else if (30 < Angle && Angle < 90)
		{
			coor_j = DgQ2DICoord(Quam, DgIVec2D(i, j + 1));//上
			coor_k = DgQ2DICoord(Quam, DgIVec2D(i, j - 1));//右
		}
		else if(150 <= Angle && Angle <= 210)
		{
			coor_j = DgQ2DICoord(Quam, DgIVec2D(i, j - 1));//右
			coor_k = DgQ2DICoord(Quam, DgIVec2D(i + 1, j + 1));//左
		}
		else if (270 < Angle && Angle < 330)
		{
			coor_j = DgQ2DICoord(Quam, DgIVec2D(i, j + 1));//上
			coor_k = DgQ2DICoord(Quam, DgIVec2D(i + 1, j + 1));//左
		}
	}
	coor.push_back(coor_j);
	if (coor_k.quadNum() != -1)
	{
		coor.push_back(coor_k);
	}
	return coor;
}
double Point2LineDist(double x1, double y1, double p1x, double p1y, double p2x, double p2y)
{
	//直线两点式Ax+By+C=0
	if (p2y == p1y) return abs(y1 - p1y);
	if (p2x == p1x) return abs(x1 - p1x); 
	double A = 1.0 / (p2x - p1x);
	double B = -1.0 / (p2y - p1y);
	double C = (p1y / (p2y - p1y)) - (p1x / (p2x - p1x));
	return abs((A * x1 + B * y1 + C) / sqrt(A*A + B*B));
}
double PointsDist(double p1x, double p1y, double p2x, double p2y)
{
	return sqrt((p1x - p2x)*(p1x - p2x) + (p1y - p2y)*(p1y - p2y));
}
vector<DgQ2DICoord> BDT_AllADJ(DgQ2DICoord cood_i, bool Orti)
{
	//先不考虑格网跨边界的情况
	vector<DgQ2DICoord> coor;
	DgQ2DICoord coor_j, coor_k;
	int Quam = cood_i.quadNum();
	int i = cood_i.coord().i();
	int j = cood_i.coord().j();
	if (Orti == 1)//上三角
	{
		coor.push_back(DgQ2DICoord(Quam, DgIVec2D(i - 1, j - 1)));//右
		coor.push_back(DgQ2DICoord(Quam, DgIVec2D(i, j - 1)));//下
		coor.push_back(DgQ2DICoord(Quam, DgIVec2D(i, j + 1)));//左
		
		/*coor.push_back(DgQ2DICoord(Quam, DgIVec2D(i, j - 2)));
		coor.push_back(DgQ2DICoord(Quam, DgIVec2D(i, j + 2)));
		coor.push_back(DgQ2DICoord(Quam, DgIVec2D(i + 1, j)));
		coor.push_back(DgQ2DICoord(Quam, DgIVec2D(i + 1, j + 1)));
		coor.push_back(DgQ2DICoord(Quam, DgIVec2D(i + 1, j + 2)));
		coor.push_back(DgQ2DICoord(Quam, DgIVec2D(i - 1, j - 3)));
		coor.push_back(DgQ2DICoord(Quam, DgIVec2D(i - 1, j - 2)));
		coor.push_back(DgQ2DICoord(Quam, DgIVec2D(i - 1, j)));
		coor.push_back(DgQ2DICoord(Quam, DgIVec2D(i - 1, j + 1)));*/
	}
	else//下三角
	{
		coor.push_back(DgQ2DICoord(Quam, DgIVec2D(i, j + 1)));//上
		coor.push_back(DgQ2DICoord(Quam, DgIVec2D(i, j - 1)));//右
		coor.push_back(DgQ2DICoord(Quam, DgIVec2D(i + 1, j + 1)));//左

		/*coor.push_back(DgQ2DICoord(Quam, DgIVec2D(i, j - 2)));
		coor.push_back(DgQ2DICoord(Quam, DgIVec2D(i, j - 1)));
		coor.push_back(DgQ2DICoord(Quam, DgIVec2D(i + 1, j)));
		coor.push_back(DgQ2DICoord(Quam, DgIVec2D(i + 1, j + 1)));
		coor.push_back(DgQ2DICoord(Quam, DgIVec2D(i + 1, j + 2)));
		coor.push_back(DgQ2DICoord(Quam, DgIVec2D(i + 1, j + 3)));
		coor.push_back(DgQ2DICoord(Quam, DgIVec2D(i + 1, j - 1)));
		coor.push_back(DgQ2DICoord(Quam, DgIVec2D(i - 1, j - 1)));
		coor.push_back(DgQ2DICoord(Quam, DgIVec2D(i - 1, j - 2)));*/
	}
	return coor;
}
/////////查找三角形的边邻近单元//////
vector<DgQ2DICoord> Edge_Proximity(DgQ2DICoord cood_i, bool UP, int maxI, int maxJ)
{
	vector<DgQ2DICoord> COOR;
	int QuamTmp = 0;
	int Quam = cood_i.quadNum();
	int i = cood_i.coord().i();
	int j = cood_i.coord().j();
	DgQ2DICoord coor;
	if (!(i == 0 & j == 0 || (i == 0 && j % 2 == 1) || (j == 0) || (i == maxI && j % 2 == 0) || (j == maxJ)))
	{
		coor = DgQ2DICoord(Quam, DgIVec2D(i, j + 1));//右
		COOR.push_back(coor);
		coor = DgQ2DICoord(Quam, DgIVec2D(i, j - 1));//下
		COOR.push_back(coor);
		if (UP)
		{
			DgQ2DICoord coor = DgQ2DICoord(Quam, DgIVec2D(i - 1, j - 1));//左
			COOR.push_back(coor);
		}
		else
		{
			coor = DgQ2DICoord(Quam, DgIVec2D(i + 1, j + 1));//右
			COOR.push_back(coor);
		}
		return COOR;

		//非边界单元
	//	if (UP)//上三角形边邻近顺序：
	//	{
	//		DgQ2DICoord coor = DgQ2DICoord(Quam, DgIVec2D(i - 1, j - 1));//左
	//		COOR.push_back(coor);
	//		coor = DgQ2DICoord(Quam, DgIVec2D(i, j + 1));//右
	//		COOR.push_back(coor);
	//		coor = DgQ2DICoord(Quam, DgIVec2D(i, j - 1));//下
	//		COOR.push_back(coor);
	///*		for (int nn = 0; nn < COOR.size(); nn++)
	//			cout << COOR[nn] << endl;*/
	//		return COOR;
	//	}
	//	else //下三角形边邻近顺序：左-右-上
	//	{
	//		coor = DgQ2DICoord(Quam, DgIVec2D(i, j - 1));//左
	//		COOR.push_back(coor);
	//		coor = DgQ2DICoord(Quam, DgIVec2D(i + 1, j + 1));//右
	//		COOR.push_back(coor);
	//		coor = DgQ2DICoord(Quam, DgIVec2D(i, j + 1));//上
	//		COOR.push_back(coor);
	//		/*for (int nn = 0; nn < COOR.size(); nn++)
	//			cout << COOR[nn] << endl;*/
	//		return COOR;
	//	}
	}
	if (j == maxJ && i != 0)
	{
		//上三角
		DgQ2DICoord coor = DgQ2DICoord(Quam, DgIVec2D(i - 1, j - 1));//左
		COOR.push_back(coor);
		coor = DgQ2DICoord(Quam, DgIVec2D(i, j - 1));//下
		COOR.push_back(coor);
		if (Quam <= 5)//上菱形
		{
			if (Quam == 5) QuamTmp = 1;
			else QuamTmp = Quam + 1;
			coor = DgQ2DICoord(QuamTmp, DgIVec2D(0, 2*(maxI - i) + 1));//右
			//以下是特殊情况处理
			COOR.push_back(coor);
			coor = DgQ2DICoord(QuamTmp, DgIVec2D(0, 2 * (maxI - i) + 1+2));//右
			COOR.push_back(coor);
			return COOR;
			//////////////////////////////////////////////

		}
		else//下菱形
		{
			if (Quam == 10) QuamTmp = 1;
			else QuamTmp = Quam - 4;
			coor = DgQ2DICoord(QuamTmp, DgIVec2D(i, 0));//右
		}
		COOR.push_back(coor);
		return COOR;
	}
	if (i == maxI && j != 0 && j % 2 == 0)
	{
		//下三角
		coor = DgQ2DICoord(Quam, DgIVec2D(i, j - 1));//左
		COOR.push_back(coor);
		coor = DgQ2DICoord(Quam, DgIVec2D(i, j + 1));//上
		COOR.push_back(coor);
		if (Quam <= 5)
		{
			QuamTmp = Quam + 5;
			coor = DgQ2DICoord(QuamTmp, DgIVec2D(0, j + 1));//右
			//以下是特殊情况处理
			COOR.push_back(coor);
			coor = DgQ2DICoord(QuamTmp, DgIVec2D(0, j + 1 - 2));//右
			COOR.push_back(coor);
			return COOR;
		}
		else
		{
			if (Quam != 10) QuamTmp = 6;
			else QuamTmp = Quam + 1;
			coor = DgQ2DICoord(QuamTmp, DgIVec2D((maxJ - j - 1) / 2, 0));//右
		}
		COOR.push_back(coor);
		return COOR;
	}
	if (j == 0)
	{
		coor = DgQ2DICoord(Quam, DgIVec2D(i + 1, j + 1));//右
		COOR.push_back(coor);
		coor = DgQ2DICoord(Quam, DgIVec2D(i, j + 1));//上
		COOR.push_back(coor);
		if (Quam <= 5)
		{
			if(Quam != 1) QuamTmp = Quam + 4;
			else QuamTmp = 10;
			coor = DgQ2DICoord(QuamTmp, DgIVec2D(i, maxJ));//左
		}
		else
		{
			if (Quam != 6) QuamTmp = Quam - 1;
			else QuamTmp = 10;
			coor = DgQ2DICoord(QuamTmp, DgIVec2D(maxI, 2 * (maxI - i)));//左
		}
		COOR.push_back(coor);
		return COOR;
	}
	if (i == 0 && j % 2 == 1)
	{
		//上三角
		coor = DgQ2DICoord(Quam, DgIVec2D(i, j + 1));//右
		COOR.push_back(coor);
		coor = DgQ2DICoord(Quam, DgIVec2D(i, j - 1));//下
		COOR.push_back(coor);
		if (Quam <= 5)
		{
			if (Quam != 1) QuamTmp = Quam - 1;
			else QuamTmp = 5;
			coor = DgQ2DICoord(QuamTmp, DgIVec2D((maxJ - j) / 2, maxJ));//左
			//以下是特殊情况处理
			COOR.push_back(coor);

			coor = DgQ2DICoord(QuamTmp, DgIVec2D(((maxJ - j) / 2)-1, maxJ));//左
			COOR.push_back(coor);

			return COOR;
			/////////////////////////////////////////////////////////////////////
		}
		else
		{
			QuamTmp = Quam - 5;
			coor = DgQ2DICoord(QuamTmp, DgIVec2D(maxI, j - 1));//左
			//以下是特殊情况处理
			COOR.push_back(coor);
			coor = DgQ2DICoord(QuamTmp, DgIVec2D(maxI, j - 1 + 2));//右
			COOR.push_back(coor);
			return COOR;
		}
		COOR.push_back(coor);
		return COOR;
	}
}
void Delete_same(vector<DgQ2DICoord> *cood, vector<DgQ2DICoord> COOD)
{
	vector<DgQ2DICoord>coor = *cood;
	int i, j;
	bool is_same = false;
	DgQ2DICoord coor_tmp, coor_j;
	for (auto iter = coor.cbegin(); iter != coor.cend();)
	{
		is_same = false;
		if (COOD.size() == 1) return;
		auto ITER = COOD.cbegin();
		while (ITER != COOD.cend())
		{
			if (*iter == *ITER)
			{
				is_same = true;
				iter = coor.erase(iter);
				break;
			}
			else
			{
				++ITER;
			}
		}
		if(!is_same) ++iter;
	}
	*cood = coor;
}
bool Point_internal(double x, double y, double p1x, double p1y, double p2x, double p2y)
{
	if (p2y == p1y)
	{
		if ((p1x <= x && x <= p2x) || (p2x <= x && x <= p1x))
			return true;
	}
	
	if (p2x == p1x)
	{
		if ((p1y <= y && y <= p2y) || (p2y <= y && y <= p1y))
			return true;
	}

	if ((p1x <= x && x <= p2x) || (p2x <= x && x <= p1x)
		&& (p1y <= y && y <= p2y) || (p2y <= y && y <= p1y))
		return true;
	return false;
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
	subject[0].resize(count+1);
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
	if (solution.size())
	{
		return  (Area(solution[0]) / SCALE_ / SCALE_);//此处不对，
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


bool Line_Line_intersection(double Ax, double Ay, double Bx, double By, 
	double Px, double Py, double Qx, double Qy,double AB, double cost, double sint)
{
	//1 按A点坐标，整体平移至原点
	Px = Px - Ax;
	Qx = Qx - Ax;
	Py = Py - Ay;
	Qy = Qy - Ay;

	//2 按AB与x轴正向夹角，整体旋转
	double BxR = Bx * cost + By * sint;
	double ByR = -Bx * sint + By * cost;
	double PxR = Px * cost + Py * sint;
	double PyR = -Px * sint + Py * cost;
	double QxR = Qx * cost + Qy * sint;
	double QyR = -Qx * sint + Qy * cost;
	//3 判断PQ与AB的位置
	if (PyR * QyR > 0) return 0;//PQ位于同侧
	//
	double k = (PxR - QxR) / (PyR - QyR);
	double InterSectX = k * (QxR / k - QyR);
	if (InterSectX < 0) return 0;
	if (InterSectX > AB) return 0;
	return 1;
}
bool Line_Grid_intersection(vector<Vec2D> gridVerCor, double x1, double y1, double x2, double y2,
	double cost, double sint, double AB)
{
	//首先判断格网与MBR的关系
	double A, B, C;
	double D1, D2, D3;
	double xmin = min(x1, x2);
	double ymin = min(y1, y2);
	double xmax = max(x1, x2);
	double ymax = max(y1, y2);
	if ((gridVerCor[0].x <= xmin && gridVerCor[1].x <= xmin && gridVerCor[2].x <= xmin)
		&& (gridVerCor[0].y <= ymin && gridVerCor[1].y <= ymin && gridVerCor[2].y <= ymin)
		&& (gridVerCor[0].x >= xmax && gridVerCor[1].x >= xmax && gridVerCor[2].x >= xmax)
		&& (gridVerCor[0].y >= ymax && gridVerCor[1].y >= ymax && gridVerCor[2].y >= ymax))
		return 0;//不相交
	//其次判断线段与线段的关系
	bool LL1 = Line_Line_intersection(x1, y1, x2-x1, y2-y1, gridVerCor[0].x, gridVerCor[0].y, 
		gridVerCor[1].x, gridVerCor[1].y, AB, cost, sint);
	bool LL2 = Line_Line_intersection(x1, y1, x2 - x1, y2 - y1, gridVerCor[0].x, gridVerCor[0].y,
		gridVerCor[2].x, gridVerCor[2].y, AB, cost, sint);
	bool LL3 = Line_Line_intersection(x1, y1, x2 - x1, y2 - y1, gridVerCor[1].x, gridVerCor[1].y, 
		gridVerCor[2].x, gridVerCor[2].y, AB, cost, sint);
	if (LL1 == 1 || LL2 == 1 || LL3 == 1) return 1;
	else return 0;
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
			for (int i = 0; i<pointCount; i++)
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
			for (int i = 0; i<pointCount; i++)
			{
				x = pOGRLinearRing->getX(i);
				y = pOGRLinearRing->getY(i);
				ogrLineString.addPoint(x, y);
			}
			pOGRMultiLineString->addGeometry(&ogrLineString);

			for (int i = 0; i<innerCount; i++)
			{
				// 添加内环
				OGRLineString ogrLineString0;
				OGRLinearRing* pOGRLinearRing0 = pOGRPolygon->getInteriorRing(i);
				int pointCount = pOGRLinearRing0->getNumPoints();
				double x = 0; double y = 0;
				for (int i = 0; i<pointCount; i++)
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

		for (int i = 0; i<geometryCount; i++)
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
		for (int i = 0; i<pointCount; i++)
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

		for (int i = 0; i<geometryCount; i++)
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
int BDOri(double xi, double xj, double yi, double yj)
{
	if (yi == yj)
	{
		if (xi > xj) return 4;
		if (xi < xj) return 3;
	}
	else
	{
		if (yi > yj) return 2;
		if (yi < yj) return 1;
	}
}

void ShpSnyderFwd(GridGenParam& dp, char *SHP_path, char *txt_path)
{
	//基础变量
	int POLYCOUNT = 1;
	int i;
	int polycount = 0;
	int BDpointcount = 0;

	GDALDataset   *poDS;
	OGRLayer  *poLayer;
	OGRFeature *poFeature;
	OGRGeometry *poGeometry = NULL;
	OGRwkbGeometryType pGeoType;
	OGRPolygon *Polygon;
	OGRLinearRing  grid_ring;
	OGRLinearRing* boundary;
	OGRSpatialReference DstSPF;
	vector<Vec2D> gridVerCor, gridVerCor_tmp;
	vector<DgQ2DICoord> Vctmp;

	//定义GDAL环境
	GDALDriver *poDriver;
	poDriver = GetGDALDriverManager()->GetDriverByName("ESRI Shapefile");
	/*DstSPF.SetProjCS("UTM 12(WGS84) in northern hemisphere.");
	DstSPF.SetWellKnownGeogCS("WGS84");
	DstSPF.SetUTM(12, TRUE);*/


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
	DgQ2DICoord coor_i, coor_j, coor_k, coor_tmp, coor_;
	const DgAddress<DgQ2DICoord> *add;
	DgLocation ADDLOC;
	vector<DgQ2DICoord> coor;
	//读取shp文件
	poDS = (GDALDataset*)GDALOpenEx(SHP_path, GDAL_OF_VECTOR, NULL, NULL, NULL);

	if (poDS == NULL)
	{
		printf("Open failed.\n%s");
		exit;
	}
	poLayer = poDS->GetLayer(0); //读取层
	poLayer->ResetReading();
	ofstream fout(txt_path);   //读取文件
	if (!fout)              // 如果读取失败，打印fail
	{
		cerr << "fail" << endl;
	}
	//poFeature = poLayer->GetNextFeature();
	//生成随机数，对于图幅分块数过多的情况，利用随机数确定poFeature是否要计算
	int poFeatureID[200] = {0};
	int FeatureTotal = poLayer->GetFeatureCount(0);
	srand((int)time(0));
	/*for (int k = 0; k < 200; k++)
		poFeatureID[k] = random(FeatureTotal);
*/
	int Count = 0;
	//while ((poFeature = poLayer->GetNextFeature()) != NULL)
	//{
	//	cout << 1 << endl;
	//	poGeometry = poFeature->GetGeometryRef();
	//	pGeoType = poGeometry->getGeometryType();
	//	//OGRLinearRing *boundary = (OGRLinearRing*)poGeometry;
	//	Polygon = (OGRPolygon *)poGeometry;//获取该要素的几何形状 ;
	//	boundary = Polygon->getExteriorRing();
	//	BDpointcount = boundary->getNumPoints();
	//	OGRFeature::DestroyFeature(poFeature);

	//}
	while ((poFeature = poLayer->GetNextFeature()) != NULL)
	{
		polycount++;
		//if (Count >= 200) break;//POLYCOUNT
		bool Is_Calculat = false;
		/*for (int k = 0; k < 200; k++)
		{
			if (polycount == poFeatureID[k])
				Is_Calculat = true;
		}*/
			
		/*if (!Is_Calculat)
		{
			OGRFeature::DestroyFeature(poFeature);
			continue;
		}*/
		Count++;
		
		//cout << polycount<<"	"<<Count << endl;
		//获得第polycount个多边形顶点

		//poGeometry = ConvertPolygonToPolyline(poFeature->GetGeometryRef());	
		poGeometry = poFeature->GetGeometryRef();
		pGeoType = poGeometry->getGeometryType();
		//OGRLinearRing *boundary = (OGRLinearRing*)poGeometry;
		Polygon = (OGRPolygon *)poGeometry;//获取该要素的几何形状 ;
		boundary = Polygon->getExteriorRing();
		BDpointcount = boundary->getNumPoints();
		fout << BDpointcount <<"	" ;
		for (int i = 0; i < BDpointcount - 1; i++)
		{
			GeoCoord LL;
			LL.lon = M_PI_180 * boundary->getX(i/* + 1144*/);
			LL.lat = M_PI_180 * boundary->getY(i /*+ 1144*/);
			//调用snyder投影
			//cout << "   ll0.lon, ll0.lat: " << LL0.lon << ", " << LL0.lat << endl;
			IcosaGridPt gridpt = snyderFwd(LL, dgg.projTriRF().sphIcosa());
			/*cout << "    gridpt.triangle .x .y: " << gridpt.triangle << ", " <<
				gridpt.pt.x << ", " << gridpt.pt.y << endl;

			cout << "DgProjTriCoord: " << DgProjTriCoord(gridpt.triangle,
				DgDVec2D(gridpt.pt.x, gridpt.pt.y)) << endl;*/
			fout << gridpt.pt.x << ", " << gridpt.pt.y << ", ";
		}
		fout << endl;
		//销毁feature
		OGRFeature::DestroyFeature(poFeature);
	}
	cout << polycount<<"	"<<Count << endl;
	fout.close();
	GDALClose(poDS);
	//GDALClose(poDS1);
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

}

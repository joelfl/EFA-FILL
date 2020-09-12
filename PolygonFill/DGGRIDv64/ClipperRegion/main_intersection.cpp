#include "clipper_region.h"
#include <iostream>
//
//void SaveToSVG(const string filename, int max_width, int max_height,
//	Paths &subj, Paths &clip, Paths &solution,
//	FillRule fill_rule, bool show_coords = false)
//{
//	SVGBuilder svg;
//	svg.fill_rule = fill_rule;
//	svg.SetCoordsStyle("Verdana", 0xFF0000AA, 9);
//	svg.AddPaths(subj, false, 0x1200009C, 0xCCD3D3DA, 0.8, show_coords);
//	svg.AddPaths(clip, false, 0x129C0000, 0xCCFFA07A, 0.8, show_coords);
//	svg.AddPaths(solution, false, 0x6080ff9C, 0xFF003300, 0.8, show_coords);
//	svg.SaveToFile(filename, max_width, max_height, 80);
//}


double clipper_intersection_area(double v1_x, double v1_y, double v2_x, double v2_y,
	double v3_x, double v3_y, double upperleft_x, double upperleft_y,
	double downright_x, double downright_y)
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
	subject[0].resize(5);
	subject[0][0].x = (int64_t)(upperleft_x * SCALE_);    subject[0][0].y = (int64_t)(upperleft_y * SCALE_);
	subject[0][1].x = (int64_t)(downright_x * SCALE_);	  subject[0][1].y = (int64_t)(upperleft_y * SCALE_);
	subject[0][2].x = (int64_t)(downright_x * SCALE_);	  subject[0][2].y = (int64_t)(downright_y * SCALE_);
	subject[0][3].x = (int64_t)(upperleft_x * SCALE_);	  subject[0][3].y = (int64_t)(downright_y * SCALE_);
	subject[0][4].x = (int64_t)(upperleft_x * SCALE_);	  subject[0][4].y = (int64_t)(upperleft_y * SCALE_);
	
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

double clipper_intersection_area_D(double v1_x, double v1_y, double v2_x, double v2_y,
	double v3_x, double v3_y, double v4_x, double v4_y,
	double upperleft_x, double upperleft_y,
	double downright_x, double downright_y)
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
	subject[0].resize(5);
	subject[0][0].x = (int64_t)(upperleft_x * SCALE_);    subject[0][0].y = (int64_t)(upperleft_y * SCALE_);
	subject[0][1].x = (int64_t)(downright_x * SCALE_);	  subject[0][1].y = (int64_t)(upperleft_y * SCALE_);
	subject[0][2].x = (int64_t)(downright_x * SCALE_);	  subject[0][2].y = (int64_t)(downright_y * SCALE_);
	subject[0][3].x = (int64_t)(upperleft_x * SCALE_);	  subject[0][3].y = (int64_t)(downright_y * SCALE_);
	subject[0][4].x = (int64_t)(upperleft_x * SCALE_);	  subject[0][4].y = (int64_t)(upperleft_y * SCALE_);

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
	if (solution.size())
	{
		return  (Area(solution[0]) / SCALE_ / SCALE_);//此处不对，
	}
	else
	{
		return 0;//表示三角形和pixel不相交
	}
}

double clipper_intersection_area_H6(double v1_x, double v1_y, double v2_x, double v2_y,
	double v3_x, double v3_y, double v4_x, double v4_y,
	double v5_x, double v5_y, double v6_x, double v6_y,
	double upperleft_x, double upperleft_y,
	double downright_x, double downright_y)
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
	subject[0].resize(5);
	subject[0][0].x = (int64_t)(upperleft_x * SCALE_);    subject[0][0].y = (int64_t)(upperleft_y * SCALE_);
	subject[0][1].x = (int64_t)(downright_x * SCALE_);	  subject[0][1].y = (int64_t)(upperleft_y * SCALE_);
	subject[0][2].x = (int64_t)(downright_x * SCALE_);	  subject[0][2].y = (int64_t)(downright_y * SCALE_);
	subject[0][3].x = (int64_t)(upperleft_x * SCALE_);	  subject[0][3].y = (int64_t)(downright_y * SCALE_);
	subject[0][4].x = (int64_t)(upperleft_x * SCALE_);	  subject[0][4].y = (int64_t)(upperleft_y * SCALE_);

	clip.resize(1);
	clip[0].resize(7);
	clip[0][0].x = (int64_t)(v1_x * SCALE_);		clip[0][0].y = (int64_t)(v1_y * SCALE_);
	clip[0][1].x = (int64_t)(v2_x * SCALE_);		clip[0][1].y = (int64_t)(v2_y * SCALE_);
	clip[0][2].x = (int64_t)(v3_x * SCALE_);		clip[0][2].y = (int64_t)(v3_y * SCALE_);
	clip[0][3].x = (int64_t)(v4_x * SCALE_);		clip[0][3].y = (int64_t)(v4_y * SCALE_);
	clip[0][4].x = (int64_t)(v5_x * SCALE_);		clip[0][4].y = (int64_t)(v5_y * SCALE_);
	clip[0][5].x = (int64_t)(v6_x * SCALE_);		clip[0][5].y = (int64_t)(v6_y * SCALE_);
	clip[0][6].x = (int64_t)(v1_x * SCALE_);		clip[0][6].y = (int64_t)(v1_y * SCALE_);

	clipper.Clear();
	clipper.AddPaths(subject, ptSubject);
	clipper.AddPaths(clip, ptClip);

	clipper.Execute(ct, solution, fr);//solution没有值 定义可能冲突了
	if (solution.size())
	{
		return  (Area(solution[0]) / SCALE_ / SCALE_);//此处不对，
	}
	else
	{
		return 0;//表示三角形和pixel不相交
	}
}

double clipper_intersection_area_H5(double v1_x, double v1_y, double v2_x, double v2_y,
	double v3_x, double v3_y, double v4_x, double v4_y,
	double v5_x, double v5_y,
	double upperleft_x, double upperleft_y,
	double downright_x, double downright_y)
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
	subject[0].resize(5);
	subject[0][0].x = (int64_t)(upperleft_x * SCALE_);    subject[0][0].y = (int64_t)(upperleft_y * SCALE_);
	subject[0][1].x = (int64_t)(downright_x * SCALE_);	  subject[0][1].y = (int64_t)(upperleft_y * SCALE_);
	subject[0][2].x = (int64_t)(downright_x * SCALE_);	  subject[0][2].y = (int64_t)(downright_y * SCALE_);
	subject[0][3].x = (int64_t)(upperleft_x * SCALE_);	  subject[0][3].y = (int64_t)(downright_y * SCALE_);
	subject[0][4].x = (int64_t)(upperleft_x * SCALE_);	  subject[0][4].y = (int64_t)(upperleft_y * SCALE_);

	clip.resize(1);
	clip[0].resize(6);
	clip[0][0].x = (int64_t)(v1_x * SCALE_);		clip[0][0].y = (int64_t)(v1_y * SCALE_);
	clip[0][1].x = (int64_t)(v2_x * SCALE_);		clip[0][1].y = (int64_t)(v2_y * SCALE_);
	clip[0][2].x = (int64_t)(v3_x * SCALE_);		clip[0][2].y = (int64_t)(v3_y * SCALE_);
	clip[0][3].x = (int64_t)(v4_x * SCALE_);		clip[0][3].y = (int64_t)(v4_y * SCALE_);
	clip[0][4].x = (int64_t)(v5_x * SCALE_);		clip[0][4].y = (int64_t)(v5_y * SCALE_);
	clip[0][5].x = (int64_t)(v1_x * SCALE_);		clip[0][5].y = (int64_t)(v1_y * SCALE_);

	clipper.Clear();
	clipper.AddPaths(subject, ptSubject);
	clipper.AddPaths(clip, ptClip);

	clipper.Execute(ct, solution, fr);//solution没有值 定义可能冲突了
	if (solution.size())
	{
		return  (Area(solution[0]) / SCALE_ / SCALE_);//此处不对，
	}
	else
	{
		return 0;//表示三角形和pixel不相交
	}
}
//---------------------------------------------------------------------------
//void main()
//{
//	double area = clipper_intersection_area(1, 1, 4, 6,
//		6, 1, 2, 3,
//		-1, -6);
//}

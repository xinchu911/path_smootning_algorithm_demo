#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "OsqpEigen/OsqpEigen.h"
#include <math.h>

//(x,y)
using Point2D = std::pair<double, double>;
using namespace std;

struct spline_curve
{
    double a0;
    double a1;
    double a2;
    double a3;
    double a4;
    double a5;
    double b0;
    double b1;
    double b2;
    double b3;
    double b4;
    double b5;
    double x_start;
    double y_start;
    double x_end;
    double y_end;
};

bool readCSVInput(int N, vector<double> &input_points_x, vector<double> &input_points_y, vector<double> &anchor_points_x, vector<double> &anchor_points_y)
{
    ifstream readCSV("../src/input_data.csv", ios::in);
    int lineNum = 0;
    string lineString;
    vector<vector<string>> stringArray;
    while (getline(readCSV, lineString))
    {
        cout << lineString << endl;
        if (lineNum > 0)
        {
            stringstream ss(lineString);
            string s;
            vector<string> oneLineStrings;
            while (getline(ss, s, ','))
                oneLineStrings.push_back(s);
            if (oneLineStrings.size() != 5)
            {
                cout << "ERROR:oneLineStrings.size() != 5" << endl;
                return false;
            }
            else
            {
                input_points_x.push_back(std::stod(oneLineStrings[1]));
                input_points_y.push_back(std::stod(oneLineStrings[2]));
                anchor_points_x.push_back(std::stod(oneLineStrings[3]));
                anchor_points_y.push_back(std::stod(oneLineStrings[4]));
            }
        }
        lineNum++;
    }
    if (N == input_points_x.size())
    {
        return true;
    }
    else
    {
        input_points_x.clear();
        input_points_y.clear();
        anchor_points_x.clear();
        anchor_points_y.clear();
        cout << "ERROR:N == input_points_x.size()" << endl;
        return false;
    }
}

bool writeCSVOutput(Eigen::VectorXd const &QPSolution)
{
    ofstream writeCSV;
    writeCSV.open("../src/QPSolution.csv", ios::out);
    writeCSV << "QPSolution" << endl;
    for (int i = 1; i < QPSolution.size(); i++)
    {
        writeCSV << QPSolution[i - 1] << endl;
    }
    writeCSV << QPSolution[QPSolution.size() - 1];
    writeCSV.close();
    cout << "INFO: csv output." << endl;
    return true;
}

int main()
{
    Eigen::VectorXd c(4);
    c << 0, -0.4, 0.02, -0.004;
    int N = 6; //input point number
    int sampling_num = N - 1;
    double interval = 5;
    Point2D start_point(0, 0);
    vector<double> input_points_x;
    vector<double> input_points_y;
    vector<double> anchor_points_x;
    vector<double> anchor_points_y;
    if (!readCSVInput(N, input_points_x, input_points_y, anchor_points_x, anchor_points_y))
    {
        cout << "ERROR:readCSVInput return" << endl;
        return 1;
    }
    vector<Point2D> input_points;
    vector<Point2D> anchor_points;
    double x_ref = 0.0;
    double y_ref = 0.0;
    double ref_tan_start = 0 * c[0] + 1 * c[1] + 2 * c[2] * start_point.second + 3 * c[3] * start_point.second * start_point.second;
    // cout << "ref_tan_start: " << ref_tan_start << endl;
    double ref_tan_end = 0 * c[0] + 1 * c[1] + 2 * c[2] * (start_point.second + (N - 1) * interval) + 3 * c[3] * pow(start_point.second + (N - 1) * interval, 2);
    // cout << "ref_tan_end: " << ref_tan_end << endl;

    double bounding_box_constraint = 1e-7;
    Point2D temp_point;
    for (int i = 0; i <= N; i++)
    {
        temp_point.first = input_points_x[i];
        temp_point.second = input_points_y[i];
        input_points.push_back(temp_point);
    }
    for (int i = 0; i < sampling_num; i++)
    {
        temp_point.first = anchor_points_x[i];
        temp_point.second = anchor_points_y[i];
        anchor_points.push_back(temp_point);
    }
    //至此得到输入离散点列与位置约束的原始线上参考点(anchor points)

    Eigen::SparseMatrix<double> hessian(72, 72);
    for (int i = 0; i < hessian.rows() / N; i++)
    {
        hessian.insert(6 * i + 3, 6 * i + 3) = 36;
        hessian.insert(6 * i + 3, 6 * i + 4) = 72;
        hessian.insert(6 * i + 3, 6 * i + 5) = 120;
        hessian.insert(6 * i + 4, 6 * i + 3) = 72;
        hessian.insert(6 * i + 4, 6 * i + 4) = 192;
        hessian.insert(6 * i + 4, 6 * i + 5) = 360;
        hessian.insert(6 * i + 5, 6 * i + 3) = 120;
        hessian.insert(6 * i + 5, 6 * i + 4) = 360;
        hessian.insert(6 * i + 5, 6 * i + 5) = 720;
    }
    Eigen::SparseMatrix<double> R(72, 72);
    for (int i = 0; i < R.rows(); i++)
    {
        R.insert(i, i) = 1.0e-5;
    }
    Eigen::SparseMatrix<double> H(72, 72);
    H = hessian + R;
    int smooth_constraint_num = 4 * (N - 1) * 2;
    int position_constraint_num = (N - 1) * 2;
    int start_end_constraint_num = 2;
    int total_constraint_num = smooth_constraint_num + position_constraint_num + start_end_constraint_num;
    Eigen::VectorXd lowerBound(total_constraint_num);
    Eigen::VectorXd upperBound(total_constraint_num);
    Eigen::SparseMatrix<double> A(total_constraint_num, 12 * N); //所有约束的矩阵
    Eigen::VectorXd gradient = Eigen::VectorXd::Zero(2 * 6 * N);
    //光滑性约束共4*(N-1)*2个,由等式化为2个不等式约束
    Eigen::VectorXd T0(12);    //0阶导数向量
    Eigen::VectorXd T0_05(12); //0阶导数向量参数0.5时
    Eigen::VectorXd T1(12);    //1阶导数向量
    Eigen::VectorXd T2(12);    //2阶导数向量
    Eigen::VectorXd T3(12);    //3阶导数向量
    T0 << 1, 1, 1, 1, 1, 1, -1, 0, 0, 0, 0, 0;
    T1 << 0, 1, 2, 3, 4, 5, 0, -1, 0, 0, 0, 0;
    T2 << 0, 0, 2, 6, 12, 20, 0, 0, -2, 0, 0, 0;
    T3 << 0, 0, 0, 6, 24, 60, 0, 0, 0, -6, 0, 0;
    T0_05 << 1, 0.5, pow(0.5, 2), pow(0.5, 3), pow(0.5, 4), pow(0.5, 5), 0, 0, 0, 0, 0, 0;
    //T0光滑约束-等式约束
    for (int i = 0; i < N - 1; ++i)
    {
        //向量ai
        A.insert(i, 6 * i + 0) = T0[0];
        A.insert(i, 6 * i + 1) = T0[1];
        A.insert(i, 6 * i + 2) = T0[2];
        A.insert(i, 6 * i + 3) = T0[3];
        A.insert(i, 6 * i + 4) = T0[4];
        A.insert(i, 6 * i + 5) = T0[5];
        A.insert(i, 6 * i + 6) = T0[6];
        //向量bi
        A.insert(i + total_constraint_num / 2, 6 * i + 6 * N + 0) = T0[0];
        A.insert(i + total_constraint_num / 2, 6 * i + 6 * N + 1) = T0[1];
        A.insert(i + total_constraint_num / 2, 6 * i + 6 * N + 2) = T0[2];
        A.insert(i + total_constraint_num / 2, 6 * i + 6 * N + 3) = T0[3];
        A.insert(i + total_constraint_num / 2, 6 * i + 6 * N + 4) = T0[4];
        A.insert(i + total_constraint_num / 2, 6 * i + 6 * N + 5) = T0[5];
        A.insert(i + total_constraint_num / 2, 6 * i + 6 * N + 6) = T0[6];
        lowerBound(i) = 0;
        upperBound(i) = 0;
        lowerBound(i + total_constraint_num / 2) = 0;
        upperBound(i + total_constraint_num / 2) = 0;
    }

    //T1光滑约束-等式约束
    for (int i = 0; i < N - 1; ++i)
    {
        //向量ai
        A.insert(i + N - 1, 6 * i + 0) = T1[0];
        A.insert(i + N - 1, 6 * i + 1) = T1[1];
        A.insert(i + N - 1, 6 * i + 2) = T1[2];
        A.insert(i + N - 1, 6 * i + 3) = T1[3];
        A.insert(i + N - 1, 6 * i + 4) = T1[4];
        A.insert(i + N - 1, 6 * i + 5) = T1[5];
        A.insert(i + N - 1, 6 * i + 7) = T1[7];
        //向量bi
        A.insert(i + total_constraint_num / 2 + N - 1, 6 * i + 6 * N + 0) = T1[0];
        A.insert(i + total_constraint_num / 2 + N - 1, 6 * i + 6 * N + 1) = T1[1];
        A.insert(i + total_constraint_num / 2 + N - 1, 6 * i + 6 * N + 2) = T1[2];
        A.insert(i + total_constraint_num / 2 + N - 1, 6 * i + 6 * N + 3) = T1[3];
        A.insert(i + total_constraint_num / 2 + N - 1, 6 * i + 6 * N + 4) = T1[4];
        A.insert(i + total_constraint_num / 2 + N - 1, 6 * i + 6 * N + 5) = T1[5];
        A.insert(i + total_constraint_num / 2 + N - 1, 6 * i + 6 * N + 7) = T1[7];
        lowerBound(i + N - 1) = 0;
        upperBound(i + N - 1) = 0;
        lowerBound(i + total_constraint_num / 2 + N - 1) = 0;
        upperBound(i + total_constraint_num / 2 + N - 1) = 0;
    }

    //T2光滑约束-等式约束
    for (int i = 0; i < N - 1; ++i)
    {
        //向量ai
        A.insert(i + (N - 1) * 2, 6 * i + 0) = T2[0];
        A.insert(i + (N - 1) * 2, 6 * i + 1) = T2[1];
        A.insert(i + (N - 1) * 2, 6 * i + 2) = T2[2];
        A.insert(i + (N - 1) * 2, 6 * i + 3) = T2[3];
        A.insert(i + (N - 1) * 2, 6 * i + 4) = T2[4];
        A.insert(i + (N - 1) * 2, 6 * i + 5) = T2[5];
        A.insert(i + (N - 1) * 2, 6 * i + 8) = T2[8];
        //向量bi
        A.insert(i + total_constraint_num / 2 + (N - 1) * 2, 6 * i + 6 * N + 0) = T2[0];
        A.insert(i + total_constraint_num / 2 + (N - 1) * 2, 6 * i + 6 * N + 1) = T2[1];
        A.insert(i + total_constraint_num / 2 + (N - 1) * 2, 6 * i + 6 * N + 2) = T2[2];
        A.insert(i + total_constraint_num / 2 + (N - 1) * 2, 6 * i + 6 * N + 3) = T2[3];
        A.insert(i + total_constraint_num / 2 + (N - 1) * 2, 6 * i + 6 * N + 4) = T2[4];
        A.insert(i + total_constraint_num / 2 + (N - 1) * 2, 6 * i + 6 * N + 5) = T2[5];
        A.insert(i + total_constraint_num / 2 + (N - 1) * 2, 6 * i + 6 * N + 8) = T2[8];

        lowerBound(i + (N - 1) * 2) = 0;
        upperBound(i + (N - 1) * 2) = 0;
        lowerBound(i + total_constraint_num / 2 + (N - 1) * 2) = 0;
        upperBound(i + total_constraint_num / 2 + (N - 1) * 2) = 0;
    }
    //T3光滑约束-等式约束
    for (int i = 0; i < N - 1; ++i)
    {
        //向量ai
        A.insert(i + (N - 1) * 3, 6 * i + 0) = T3[0];
        A.insert(i + (N - 1) * 3, 6 * i + 1) = T3[1];
        A.insert(i + (N - 1) * 3, 6 * i + 2) = T3[2];
        A.insert(i + (N - 1) * 3, 6 * i + 3) = T3[3];
        A.insert(i + (N - 1) * 3, 6 * i + 4) = T3[4];
        A.insert(i + (N - 1) * 3, 6 * i + 5) = T3[5];
        A.insert(i + (N - 1) * 3, 6 * i + 9) = T3[9];
        //向量bi
        A.insert(i + total_constraint_num / 2 + (N - 1) * 3, 6 * i + 6 * N + 0) = T3[0];
        A.insert(i + total_constraint_num / 2 + (N - 1) * 3, 6 * i + 6 * N + 1) = T3[1];
        A.insert(i + total_constraint_num / 2 + (N - 1) * 3, 6 * i + 6 * N + 2) = T3[2];
        A.insert(i + total_constraint_num / 2 + (N - 1) * 3, 6 * i + 6 * N + 3) = T3[3];
        A.insert(i + total_constraint_num / 2 + (N - 1) * 3, 6 * i + 6 * N + 4) = T3[4];
        A.insert(i + total_constraint_num / 2 + (N - 1) * 3, 6 * i + 6 * N + 5) = T3[5];
        A.insert(i + total_constraint_num / 2 + (N - 1) * 3, 6 * i + 6 * N + 9) = T3[9];

        lowerBound(i + (N - 1) * 3) = 0;
        upperBound(i + (N - 1) * 3) = 0;
        lowerBound(i + total_constraint_num / 2 + (N - 1) * 3) = 0;
        upperBound(i + total_constraint_num / 2 + (N - 1) * 3) = 0;
    }

    //起止点约束-等式约束
    A.insert(smooth_constraint_num / 2, 1) = 1;
    A.insert(smooth_constraint_num / 2, 6 * N + 1) = -ref_tan_start;

    A.insert((smooth_constraint_num + total_constraint_num) / 2, 6 * (N - 1) + 1) = 1;
    A.insert((smooth_constraint_num + total_constraint_num) / 2, 6 * (N - 1) + 2) = 2;
    A.insert((smooth_constraint_num + total_constraint_num) / 2, 6 * (N - 1) + 3) = 3;
    A.insert((smooth_constraint_num + total_constraint_num) / 2, 6 * (N - 1) + 4) = 4;
    A.insert((smooth_constraint_num + total_constraint_num) / 2, 6 * (N - 1) + 5) = 5;
    A.insert((smooth_constraint_num + total_constraint_num) / 2, 6 * (2 * N - 1) + 1) = -ref_tan_end * 1;
    A.insert((smooth_constraint_num + total_constraint_num) / 2, 6 * (2 * N - 1) + 2) = -ref_tan_end * 2;
    A.insert((smooth_constraint_num + total_constraint_num) / 2, 6 * (2 * N - 1) + 3) = -ref_tan_end * 3;
    A.insert((smooth_constraint_num + total_constraint_num) / 2, 6 * (2 * N - 1) + 4) = -ref_tan_end * 4;
    A.insert((smooth_constraint_num + total_constraint_num) / 2, 6 * (2 * N - 1) + 5) = -ref_tan_end * 5;

    lowerBound(smooth_constraint_num / 2) = 0;
    upperBound(smooth_constraint_num / 2) = 0;
    lowerBound((smooth_constraint_num + total_constraint_num) / 2) = 0;
    upperBound((smooth_constraint_num + total_constraint_num) / 2) = 0;

    //位置bounding box约束-不等式约束
    for (int i = 1; i <= sampling_num; i++)
    {
        // x坐标约束
        A.insert(smooth_constraint_num / 2 + i, N * (i - 1) + 0) = T0_05[0];
        A.insert(smooth_constraint_num / 2 + i, N * (i - 1) + 1) = T0_05[1];
        A.insert(smooth_constraint_num / 2 + i, N * (i - 1) + 2) = T0_05[2];
        A.insert(smooth_constraint_num / 2 + i, N * (i - 1) + 3) = T0_05[3];
        A.insert(smooth_constraint_num / 2 + i, N * (i - 1) + 4) = T0_05[4];
        A.insert(smooth_constraint_num / 2 + i, N * (i - 1) + 5) = T0_05[5];
        // y坐标约束
        A.insert((smooth_constraint_num + total_constraint_num) / 2 + i, N * (i - 1) + 0 + 36) = T0_05[0];
        A.insert((smooth_constraint_num + total_constraint_num) / 2 + i, N * (i - 1) + 1 + 36) = T0_05[1];
        A.insert((smooth_constraint_num + total_constraint_num) / 2 + i, N * (i - 1) + 2 + 36) = T0_05[2];
        A.insert((smooth_constraint_num + total_constraint_num) / 2 + i, N * (i - 1) + 3 + 36) = T0_05[3];
        A.insert((smooth_constraint_num + total_constraint_num) / 2 + i, N * (i - 1) + 4 + 36) = T0_05[4];
        A.insert((smooth_constraint_num + total_constraint_num) / 2 + i, N * (i - 1) + 5 + 36) = T0_05[5];

        lowerBound(smooth_constraint_num / 2 + i) = anchor_points[i - 1].first - bounding_box_constraint;
        upperBound(smooth_constraint_num / 2 + i) = anchor_points[i - 1].first + bounding_box_constraint;
        lowerBound((smooth_constraint_num + total_constraint_num) / 2 + i) = anchor_points[i - 1].second - bounding_box_constraint;
        upperBound((smooth_constraint_num + total_constraint_num) / 2 + i) = anchor_points[i - 1].second + bounding_box_constraint;
        cout << i << "th x lowerBound: " << lowerBound(smooth_constraint_num / 2 + i) << "::" << anchor_points[i - 1].first - bounding_box_constraint << endl;
        cout << i << "th x upperBound: " << upperBound(smooth_constraint_num / 2 + i) << "::" << anchor_points[i - 1].first + bounding_box_constraint << endl;
        cout << i << "th y lowerBound: " << lowerBound((smooth_constraint_num + total_constraint_num) / 2 + i) << "::" << anchor_points[i - 1].second - bounding_box_constraint << endl;
        cout << i << "th y upperBound: " << upperBound((smooth_constraint_num + total_constraint_num) / 2 + i) << "::" << anchor_points[i - 1].second + bounding_box_constraint << endl;
    }
    // cout << A << endl;
    // cout << lowerBound << endl;
    // cout << "---------------------" << endl;
    // cout << upperBound << endl;

    OsqpEigen::Solver solver;

    solver.settings()->setVerbosity(false);
    solver.settings()->setWarmStart(true);

    solver.data()->setNumberOfVariables(2 * 6 * N);
    solver.data()->setNumberOfConstraints(total_constraint_num);

    if (!solver.data()->setHessianMatrix(H))
        return 1;
    if (!solver.data()->setGradient(gradient))
        return 1;
    if (!solver.data()->setLinearConstraintsMatrix(A))
        return 1;
    if (!solver.data()->setLowerBound(lowerBound))
        return 1;
    if (!solver.data()->setUpperBound(upperBound))
        return 1;

    if (!solver.initSolver())
        return 1;

    Eigen::VectorXd QPSolution;

    if (!solver.solve())
    {
        return 1;
    }
    QPSolution = solver.getSolution();
    cout << "QPSolution" << endl
         << QPSolution << endl;
    if (!writeCSVOutput(QPSolution))
    {
        cout << "ERROR" << endl;
        return 1;
    }
    spline_curve temp_curve;
    vector<spline_curve> curve_seg;
    cout << "QPSolution.rows():" << QPSolution.rows() << endl;
    Eigen::VectorXd px;
    Eigen::VectorXd py;
    for (int i = 0; i < QPSolution.rows() / N / 2; i++)
    {
        temp_curve.a0 = QPSolution(6 * i + 0);
        temp_curve.a1 = QPSolution(6 * i + 1);
        temp_curve.a2 = QPSolution(6 * i + 2);
        temp_curve.a3 = QPSolution(6 * i + 3);
        temp_curve.a4 = QPSolution(6 * i + 4);
        temp_curve.a5 = QPSolution(6 * i + 5);

        temp_curve.b0 = QPSolution(6 * (i + N) + 0);
        temp_curve.b1 = QPSolution(6 * (i + N) + 1);
        temp_curve.b2 = QPSolution(6 * (i + N) + 2);
        temp_curve.b3 = QPSolution(6 * (i + N) + 3);
        temp_curve.b4 = QPSolution(6 * (i + N) + 4);
        temp_curve.b5 = QPSolution(6 * (i + N) + 5);

        temp_curve.x_start = temp_curve.a0;
        temp_curve.y_start = temp_curve.b0;
        temp_curve.x_end = temp_curve.a0 + temp_curve.a1 + temp_curve.a2 + temp_curve.a3 + temp_curve.a4 + temp_curve.a5;
        temp_curve.y_end = temp_curve.b0 + temp_curve.b1 + temp_curve.b2 + temp_curve.b3 + temp_curve.b4 + temp_curve.b5;

        curve_seg.push_back(temp_curve);
    }

    return 0;
}
#include <iostream>
#include <cmath>
#include <ctime>
#include<vector>                   

const int population = 20;
float mutation_prob = 0.01;
float crossoover_prob = 0.8;
using namespace std;

class lw3
{
private:
    double x1[population], y1[population], x2[population], y2[population], x3[population], y3[population];
    double result1[population], result2[population], result3[population];
    short a, index = 0, index1 = 0;
    int rnd1, rnd2, rnd3, rnd11, rnd22, rnd33, iteration = 0;
    double mut_elem[1][2], cross_elem[1][2];
    int binnary_x1[9], binnary_y1[9], binnary_x2[9], binnary_y2[9];
    int binnary_x21[6], binnary_y21[6], binnary_x22[6], binnary_y22[6];
    int binnary_x31[14], binnary_y31[14], binnary_x32[14], binnary_y32[14];
    bool flag = true;
    double min1, min2, min3, new_min1, new_min2, new_min3, q = 400;
    double func1_variant[132] = { -4.5, -4.43089, -4.36178, -4.29267, -4.22356, -4.15444, -4.08533, -4.01622, -3.94711, -3.878, -3.80889, -3.73978, -3.67067, -3.60156, -3.53245, -3.46333, -3.39422, -3.32511, -3.256, -3.18689, -3.11778, -3.04867, -2.97956, -2.91045, -2.84134, -2.77222, -2.70311, -2.634, -2.56489, -2.49578, -2.42667, -2.35756, -2.28845, -2.21934, -2.15023, -2.08111, -2.012, -1.94289, -1.87378, -1.80467, -1.73556, -1.66645, -1.59734, -1.52823, -1.45912, -1.39001, -1.32089, -1.25178, -1.18267, -1.11356, -1.04445, -0.975339, -0.906228, -0.837117, -0.768006, -0.698895, -0.629784, -0.560673, -0.491562, -0.422451, -0.35334, -0.284229, -0.215118, -0.146007, -0.076896, -0.007785, 0.061326, 0.130437, 0.199548, 0.268659, 0.33777, 0.406881, 0.475992, 0.545103, 0.614214, 0.683325, 0.752436, 0.821547, 0.890658, 0.959769, 1.02888, 1.09799, 1.1671, 1.23621, 1.30532, 1.37443, 1.44355, 1.51266, 1.58177, 1.65088, 1.71999, 1.7891, 1.85821, 1.92732, 1.99643, 2.06554, 2.13466, 2.20377, 2.27288, 2.34199, 2.4111, 2.48021, 2.54932, 2.61843, 2.68754, 2.75665, 2.82577, 2.89488, 2.96399, 3.0331, 3.10221, 3.17132, 3.24043, 3.30954, 3.37865, 3.44776, 3.51688, 3.58599, 3.6551, 3.72421, 3.79332, 3.86243, 3.93154, 4.00065, 4.06976, 4.13887, 4.20799, 4.2771, 4.34621, 4.41532, 4.48443, 4.55354 };
    double func2_variant[41] = { -2, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2};
    vector<double>func3_variant;//[16384];
    double func1_cile[132] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131 };
    double func2_cile[40] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39 };
    vector<double>func3_cile;// [16384] ;
    int used_index_x1[population], used_index_y1[population], used_index_x2[population], used_index_y2[population];
    vector<double> global_min;
    vector<int>used_index_x3;
    vector<int>used_index_y3;
public:
    double func1(double x, double y)
    {
        return pow((1.5 - x + x * y), 2) + pow((2.25 - x + x * pow(y, 2)), 2) + pow((2.625 - x + x * pow(y, 3)), 2);
    }

    int func2(double x, double y)
    {
        return  (1 + pow((x + y + 1), 2) * (19 - 14 * x + 3 * x * x - 14 * y + 6 * x * y + 3 * y * y)) * (30 + pow((2 * x - 3 * y), 2) * (18 - 32 * x + 12 * x * x + 48 * y - 36 * x * y + 27 * y * y));
    }

    double func3(double x, double y)
    {
        return  -1 * (y + 47) * sin(sqrt(fabs((x / 2) + y + 47))) - x * sin(sqrt(fabs(x - y - 47)));
    }

    void rand_inp()
    {
        for (int i = 0; i < population; i++)
        {
            used_index_y3.push_back(0);
            used_index_x3.push_back(0);
        }
        int rnd, rnd1;
        for (int i = 0; i < population; i++)
        {
            if (a == 1)
            {
                rnd = 0 + rand() % (131 - 0 + 1);
                rnd1 = 0 + rand() % (131 - 0 + 1);
                x1[i] = func1_variant[rnd];
                y1[i] = func1_variant[rnd1];
                used_index_x1[i] = rnd;
                used_index_y1[i] = rnd1;
            }
            if (a == 2)
            {
                rnd = 0 + rand() % (39 - 0 + 1);
                rnd1 = 0 + rand() % (39 - 0 + 1);
                x2[i] = func2_variant[rnd];
                y2[i] = func2_variant[rnd1];
                used_index_x2[i] = rnd;
                used_index_y2[i] = rnd1;
            }
            if (a == 3)
            {
                rnd = 0 + rand() % (16383 - 0 + 1);
                rnd1 = 0 + rand() % (16383 - 0 + 1);
                x3[i] = func3_variant[rnd];
                y3[i] = func3_variant[rnd1];
                used_index_x3[i] = rnd;
                used_index_y3[i] = rnd1;
            }
        }
    }
    void set_a(short b)
    {
        a = b;
    }

    lw3 main_func()
    {
        while (flag)
        {
            for (int i = 0; i < population; i++)
            {
                for (int i = 0; i < population; i++)
                {
                    if (a == 1)
                    {
                        result1[i] = func1(x1[i], y1[i]);
                    }
                    if (a == 2)
                    {
                        result2[i] = func2(x2[i], y2[i]);
                    }
                    if (a == 3)
                    {
                        result3[i] = func3(x3[i], y3[i]);
                    }
                }

                min1 = result1[0]; min2 = result2[0]; min3 = result3[0];
                for (int i = 0; i < population; i++)
                {
                    if (a == 1)
                    {
                        if (result1[i] < min1)
                        {
                            min1 = result1[i];
                            index = i;
                        }
                    }
                    if (a == 2)
                    {
                        if (result2[i] < min2)
                        {
                            min2 = result2[i];
                            index = i;
                        }
                    }
                    if (a == 3)
                    {
                        if (result3[i] < min3)
                        {
                            min3 = result3[i];
                            index = i;
                        }
                    }
                }
                if (a == 1)
                {
                    global_min.push_back(min1);
                }
                if (a == 2)
                {
                    global_min.push_back(min2);
                }
                if (a == 3)
                {
                    global_min.push_back(min3);
                }

                if (a == 1)
                {
                    rnd1 = 0 + rand() % (19 - 0 + 1);
                    rnd11 = 0 + rand() % (19 - 0 + 1);
                }
                if (a == 2)
                {
                    rnd2 = 0 + rand() % (19 - 0 + 1);
                    rnd22 = 0 + rand() % (19 - 0 + 1);
                }
                if (a == 3)
                {
                    rnd3 = 0 + rand() % (19 - 0 + 1);
                    rnd33 = 0 + rand() % (19 - 0 + 1);
                }

                double cross = 0 + double(rand()) / RAND_MAX * (1 - 0);
                if (cross < crossoover_prob)
                {
                    if (a == 1)
                    {
                        crossover(&x1[rnd1], &y1[rnd1], &x1[rnd11], &y1[rnd11], &rnd1, &rnd11);
                    }
                    if (a == 2)
                    {
                        crossover(&x2[rnd2], &y2[rnd2], &x2[rnd22], &y2[rnd22], &rnd2, &rnd22);
                    }
                    if (a == 3)
                    {
                        crossover(&x3[rnd3], &y3[rnd3], &x3[rnd33], &y3[rnd33], &rnd3, &rnd33);
                    }
                }
                
                double getA = 0 + double(rand()) / RAND_MAX * (1 - 0);
                int ind;
                if (getA < 0.5)
                {
                    if (a == 1)
                    {
                        mut_elem[0][0] = x1[rnd1];
                        mut_elem[0][1] = y1[rnd1];
                        ind = rnd1;
                    }
                    if (a == 2)
                    {
                        mut_elem[0][0] = x2[rnd2];
                        mut_elem[0][1] = y2[rnd2];
                        ind = rnd2;
                    }
                    if (a == 3)
                    {
                        mut_elem[0][0] = x3[rnd3];
                        mut_elem[0][1] = y3[rnd3];
                        ind = rnd3;
                    }
                }
                else
                {
                    if (a == 1)
                    {
                        mut_elem[0][0] = x1[rnd11];
                        mut_elem[0][1] = y1[rnd11];
                        ind = rnd11;
                    }
                    if (a == 2)
                    {
                        mut_elem[0][0] = x2[rnd22];
                        mut_elem[0][1] = y2[rnd22];
                        ind = rnd22;
                    }
                    if (a == 3)
                    {
                        mut_elem[0][0] = x3[rnd33];
                        mut_elem[0][1] = y3[rnd33];
                        ind = rnd33;
                    }
                }
                
                double mut = 0 + rand() % (1 - 0 + 1);
                if (mut < mutation_prob)
                {
                    mutation(&mut_elem[0][0], &mut_elem[0][1], &ind);
                }
            }
            iteration++;
            int r;
            if (a == 1)
            {
                r = 10;
            }
            if (a == 2)
            {
                r = 5;
            }
            if (a == 3)
            {
                r = 4;
            }
            if (iteration % r == 0)
            {
                if (a == 1)
                {
                    new_min1 = func1(x1[0], y1[0]);
                    double minn;
                    for (int i = 0; i < population; i++)
                    {
                        minn = func1(x1[i], y1[i]);
                        if (minn < new_min1)
                        {
                            new_min1 = func1(x1[i], y1[i]);
                            index1 = i;
                        }
                    }
                    minn = global_min[0];
                    for (int i = 0; i < global_min.size(); i++)
                    {
                        if (global_min[i] < minn)
                        {
                            minn = global_min[i];
                        }
                    }
                    if (new_min1 <= minn)
                    {
                        flag = false;
                        index = index1;
                    }
                }
                if (a == 2)
                {
                    new_min2 = func2(x2[0], y2[0]);
                    double minn;
                    for (int i = 0; i < population; i++)
                    {
                        minn = func2(x2[i], y2[i]);
                        if (minn < new_min2)
                        {
                            new_min2 = func2(x2[i], y2[i]);
                            index1 = i;
                        }
                    }
                    minn = global_min[0];
                    for (int i = 0; i < global_min.size(); i++)
                    {
                        if (global_min[i] < minn)
                        {
                            minn = global_min[i];
                        }
                    }
                    if (new_min2 <= minn)
                    {
                        flag = false;
                        index = index1;
                    }
                    else if (iteration > 500)
                    {
                        flag = false;
                    }
                }
                if (a == 3)
                {
                    new_min3 = func3(x3[0], y3[0]);
                    double minn;
                    for (int i = 0; i < population; i++)
                    {
                        minn = func3(x3[i], y3[i]);
                        if (minn < new_min3)
                        {
                            new_min3 = func3(x3[i], y3[i]);
                            index1 = i;
                        }
                    }
                    minn = global_min[0];
                    for (int i = 0; i < global_min.size(); i++)
                    {
                        if (global_min[i] < minn)
                        {
                            minn = global_min[i];
                        }
                    }
                    if (new_min3 <= minn)
                    {
                        flag = false;
                        index = index1;
                    }
                    else if (iteration > 800)
                    {
                        flag = false;
                    }
                }
            }
            if (a == 1)
            {
                cout << "iteration = " << iteration << endl;
                cout << "f = " << func1(x1[index], y1[index]) << " x = " << x1[index] << " y = " << y1[index] << endl;
            }
            if (a == 2)
            {
                cout << "iteration = " << iteration << endl;
                cout << "f = " << func2(x2[index], y2[index]) << " x = " << x2[index] << " y = " << y2[index] << endl;
            }
            if (a == 3)
            {
                cout << "iteration = " << iteration << endl;
                cout << "f = " << func3(x3[index], y3[index]) << " x = " << x3[index] << " y = " << y3[index] << endl;
            }
        }
        flag = true;
        iteration = 0;
        return *this;
    }
    void mutation(double* x_1, double* y_1, int *index1)
    {
        if (a == 1)
        {
            this->get_binary(x_1, y_1, index1);
            double rnd = 0 + double(rand()) / RAND_MAX * (1 - 0);
            if (rnd <= 0.1)
            {
                if (binnary_x1[0] == 1)
                {
                    binnary_x1[0] = 0;
                }
                else
                {
                    binnary_x1[0] = 1;
                }

                if (binnary_y1[0] == 1)
                {
                    binnary_y1[0] = 0;
                }
                else
                {
                    binnary_y1[0] = 1;
                }
            }
            if (rnd > 0.1 && rnd <= 0.5)
            {
                if (binnary_x1[1] == 0)
                {
                    binnary_x1[1] = 1;
                }
                else
                {
                    binnary_x1[1] = 0;
                }

                if (binnary_y1[0] == 0)
                {
                    binnary_y1[0] = 1;
                }
                else
                {
                    binnary_y1[0] = 0;
                }
            }
            if (rnd > 0.5)
            {
                if (binnary_x1[2] == 1)
                {
                    binnary_x1[2] = 0;
                }
                else
                {
                    binnary_x1[2] = 1;
                }

                if (binnary_y1[2] == 1)
                {
                    binnary_y1[2] = 0;
                }
                else
                {
                    binnary_y1[2] = 1;
                }
            }
            double x11 = 0, y11 = 0;
            int j = 8;
            for (int i = 0; i < 9; i++)
            {
                if (binnary_x1[i] == 1)
                {
                    x11 += pow(2, j);
                }

                if (binnary_y1[i] == 1)
                {
                    y11 += pow(2, j);
                }
                j--;
            }
            for (int i = 0; i < 132; i++)
            {
                if (func1_cile[i] == x11)
                {
                    used_index_x1[*index1] = i;
                    x1[*index1] = func1_variant[i];
                    break;
                }
            }
            for (int i = 0; i < 132; i++)
            {
                if (func1_cile[i] == y11)
                {
                    used_index_y1[*index1] = i;
                    y1[*index1] = func1_variant[i];
                    break;
                }
            }
        }
        if (a == 2)
        {
            this->get_binary(x_1, y_1, index1);
            double rnd = 0 + double(rand()) / RAND_MAX * (1 - 0);
            if (rnd <= 0.1)
            {
                if (binnary_x21[0] == 1)
                {
                    binnary_x21[0] = 0;
                }
                else
                {
                    binnary_x21[0] = 1;
                }

                if (binnary_y21[0] == 1)
                {
                    binnary_y21[0] = 0;
                }
                else
                {
                    binnary_y21[0] = 1;
                }
            }
            if (rnd > 0.1 && rnd <= 0.5)
            {
                if (binnary_x21[1] == 0)
                {
                    binnary_x21[1] = 1;
                }
                else
                {
                    binnary_x21[1] = 0;
                }

                if (binnary_y21[0] == 0)
                {
                    binnary_y21[0] = 1;
                }
                else
                {
                    binnary_y21[0] = 0;
                }
            }
            if (rnd > 0.5)
            {
                if (binnary_x21[2] == 1)
                {
                    binnary_x21[2] = 0;
                }
                else
                {
                    binnary_x21[2] = 1;
                }

                if (binnary_y21[2] == 1)
                {
                    binnary_y21[2] = 0;
                }
                else
                {
                    binnary_y21[2] = 1;
                }
            }
            double x11 = 0, y11 = 0;
            int j = 5;
            for (int i = 0; i < 6; i++)
            {
                if (binnary_x21[i] == 1)
                {
                    x11 += pow(2, j);
                }

                if (binnary_y21[i] == 1)
                {
                    y11 += pow(2, j);
                }
                j--;
            }
            for (int i = 0; i < 40; i++)
            {
                if (func2_cile[i] == x11)
                {
                    used_index_x2[*index1] = i;
                    x2[*index1] = func2_variant[i];
                    break;
                }
            }
            for (int i = 0; i < 40; i++)
            {
                if (func2_cile[i] == y11)
                {
                    used_index_y2[*index1] = i;
                    y2[*index1] = func2_variant[i];
                    break;
                }
            }
        }
        if (a == 3)
        {
            this->get_binary(x_1, y_1, index1);
            double rnd = 0 + double(rand()) / RAND_MAX * (1 - 0);
            if (rnd <= 0.1)
            {
                if (binnary_x31[0] == 1)
                {
                    binnary_x31[0] = 0;
                }
                else
                {
                    binnary_x31[0] = 1;
                }

                if (binnary_y31[0] == 1)
                {
                    binnary_y31[0] = 0;
                }
                else
                {
                    binnary_y31[0] = 1;
                }
            }
            if (rnd > 0.1 && rnd <= 0.5)
            {
                if (binnary_x31[1] == 0)
                {
                    binnary_x31[1] = 1;
                }
                else
                {
                    binnary_x31[1] = 0;
                }

                if (binnary_y31[0] == 0)
                {
                    binnary_y31[0] = 1;
                }
                else
                {
                    binnary_y31[0] = 0;
                }
            }
            if (rnd > 0.5)
            {
                if (binnary_x31[2] == 1)
                {
                    binnary_x31[2] = 0;
                }
                else
                {
                    binnary_x31[2] = 1;
                }

                if (binnary_y31[2] == 1)
                {
                    binnary_y31[2] = 0;
                }
                else
                {
                    binnary_y31[2] = 1;
                }
            }
            double x11 = 0, y11 = 0;
            int j = 13;
            for (int i = 0; i < 14; i++)
            {
                if (binnary_x31[i] == 1)
                {
                    x11 += pow(2, j);
                }

                if (binnary_y31[i] == 1)
                {
                    y11 += pow(2, j);
                }
                j--;
            }
            for (int i = 0; i < 16384; i++)
            {
                if (func3_cile[i] == x11)
                {
                    used_index_x3[*index1] = i;
                    x3[*index1] = func3_variant[i];
                    break;
                }
            }
            for (int i = 0; i < 16384; i++)
            {
                if (func3_cile[i] == y11)
                {
                    used_index_y3[*index1] = i;
                    y3[*index1] = func3_variant[i];
                    break;
                }
            }
        }
    }
    void crossover(double* x_1, double* y_1, double* x_2, double* y_2, int *index1, int *index2)
    {
        if (a == 1)
        {
            this->get_binary(x_1, y_1, x_2, y_2, index1, index2);
            int temp1 = binnary_x1[0], temp2 = binnary_x1[1], temp3 = binnary_x1[2], temp4 = binnary_x1[3], temp5 = binnary_x1[4];
            binnary_x1[0] = binnary_x2[0];
            binnary_x1[1] = binnary_x2[1];
            binnary_x1[2] = binnary_x2[2];
            binnary_x1[3] = binnary_x2[3];
            binnary_x1[4] = binnary_x2[4];

            binnary_x2[0] = temp1;
            binnary_x2[1] = temp2;
            binnary_x2[2] = temp3;
            binnary_x2[3] = temp4;
            binnary_x2[4] = temp5;

            temp1 = binnary_y1[0], temp2 = binnary_y1[1], temp3 = binnary_y1[2], temp4 = binnary_y1[3], temp5 = binnary_y1[4];
            binnary_y1[0] = binnary_y2[0];
            binnary_y1[1] = binnary_y2[1];
            binnary_y1[2] = binnary_y2[2];
            binnary_y1[3] = binnary_y2[3];
            binnary_y1[4] = binnary_y2[4];

            binnary_y2[0] = temp1;
            binnary_y2[1] = temp2;
            binnary_y2[2] = temp3;
            binnary_y2[3] = temp4;
            binnary_y2[4] = temp5;

            double x11 = 0, x22 = 0, y11 = 0, y22 = 0;
            int j = 8;
            for (int i = 0; i < 9; i++)
            {
                if (binnary_x1[i] == 1)
                {
                    x11 += pow(2, j);
                }
                if (binnary_x2[i] == 1)
                {
                    x22 += pow(2, j);
                }

                if (binnary_y1[i] == 1)
                {
                    y11 += pow(2, j);
                }
                if (binnary_y2[i] == 1)
                {
                    y22 += pow(2, j);
                }
                j--;
            }

            for (int i = 0; i < 132; i++)
            {
                if (func1_cile[i] == x11)
                {
                    used_index_x1[*index1] = i;
                    x1[*index1] = func1_variant[i];
                    break;
                }
            }
            for (int i = 0; i < 132; i++)
            {
                if (func1_cile[i] == y11)
                {
                    used_index_y1[*index1] = i;
                    y1[*index1] = func1_variant[i];
                }
            }

            for (int i = 0; i < 132; i++)
            {
                if (func1_cile[i] == x22)
                {
                    used_index_x1[*index2] = i;
                    x1[*index2] = func1_variant[i];
                    break;
                }
            }
            for (int i = 0; i < 132; i++)
            {
                if (func1_cile[i] == y22)
                {
                    used_index_y1[*index2] = i;
                    y1[*index2] = func1_variant[i];
                }
            }
        }
        if (a == 2)
        {
            this->get_binary(x_1, y_1, x_2, y_2, index1, index2);
            int temp1 = binnary_x21[0], temp2 = binnary_x21[1], temp3 = binnary_x22[2];
            binnary_x21[0] = binnary_x22[0];
            binnary_x21[1] = binnary_x22[1];
            binnary_x21[2] = binnary_x22[2];

            binnary_x22[0] = temp1;
            binnary_x22[1] = temp2;
            binnary_x22[2] = temp3;

            temp1 = binnary_y21[0], temp2 = binnary_y21[1], temp3 = binnary_y21[2];
            binnary_y21[0] = binnary_y22[0];
            binnary_y21[1] = binnary_y22[1];
            binnary_y21[2] = binnary_y22[2];


            binnary_y22[0] = temp1;
            binnary_y22[1] = temp2;
            binnary_y22[2] = temp3;

            double x11 = 0, x22 = 0, y11 = 0, y22 = 0;
            int j = 5;
            for (int i = 0; i < 6; i++)
            {
                if (binnary_x21[i] == 1)
                {
                    x11 += pow(2, j);
                }
                if (binnary_x22[i] == 1)
                {
                    x22 += pow(2, j);
                }

                if (binnary_y21[i] == 1)
                {
                    y11 += pow(2, j);
                }
                if (binnary_y22[i] == 1)
                {
                    y22 += pow(2, j);
                }
                j--;
            }

            for (int i = 0; i < 40; i++)
            {
                if (func2_cile[i] == x11)
                {
                    used_index_x2[*index1] = i;
                    x2[*index1] = func2_variant[i];
                    break;
                }
            }
            for (int i = 0; i < 40; i++)
            {
                if (func2_cile[i] == y11)
                {
                    used_index_y2[*index1] = i;
                    y2[*index1] = func2_variant[i];
                }
            }

            for (int i = 0; i < 40; i++)
            {
                if (func2_cile[i] == x22)
                {
                    used_index_x2[*index2] = i;
                    x2[*index2] = func2_variant[i];
                    break;
                }
            }
            for (int i = 0; i < 40; i++)
            {
                if (func2_cile[i] == y22)
                {
                    used_index_y2[*index2] = i;
                    y2[*index2] = func2_variant[i];
                }
            }
        }
        if (a == 3)
        {
            this->get_binary(x_1, y_1, x_2, y_2, index1, index2);
            int temp1 = binnary_x31[0], temp2 = binnary_x31[1], temp3 = binnary_x31[2], temp4 = binnary_x31[3], temp5 = binnary_x31[4], temp6 = binnary_x31[5];
            binnary_x31[0] = binnary_x32[0];
            binnary_x31[1] = binnary_x32[1];
            binnary_x31[2] = binnary_x32[2];
            binnary_x31[3] = binnary_x32[3];
            binnary_x31[4] = binnary_x32[4];
            binnary_x31[5] = binnary_x31[5];

            binnary_x32[0] = temp1;
            binnary_x32[1] = temp2;
            binnary_x32[2] = temp3;
            binnary_x32[3] = temp4;
            binnary_x32[4] = temp5;
            binnary_x32[5] = temp6;

            temp1 = binnary_y31[0], temp2 = binnary_y31[1], temp3 = binnary_y31[2], temp4 = binnary_y31[3], temp5 = binnary_y31[4], temp6 = binnary_y31[5];
            binnary_y31[0] = binnary_y32[0];
            binnary_y31[1] = binnary_y32[1];
            binnary_y31[2] = binnary_y32[2];
            binnary_y31[3] = binnary_y32[3];
            binnary_y31[4] = binnary_y32[4];
            binnary_y31[5] = binnary_y32[5];

            binnary_y32[0] = temp1;
            binnary_y32[1] = temp2;
            binnary_y32[2] = temp3;
            binnary_y32[3] = temp4;
            binnary_y32[4] = temp5;
            binnary_y32[5] = temp6;

            double x11 = 0, x22 = 0, y11 = 0, y22 = 0;
            int j = 13;
            for (int i = 0; i < 14; i++)
            {
                if (binnary_x31[i] == 1)
                {
                    x11 += pow(2, j);
                }
                if (binnary_x32[i] == 1)
                {
                    x22 += pow(2, j);
                }

                if (binnary_y31[i] == 1)
                {
                    y11 += pow(2, j);
                }
                if (binnary_y32[i] == 1)
                {
                    y22 += pow(2, j);
                }
                j--;
            }

            for (int i = 0; i < 16384; i++)
            {
                if (func3_cile[i] == x11)
                {
                    used_index_x3[*index1] = i;
                    x3[*index1] = func3_variant[i];
                    break;
                }
            }
            for (int i = 0; i < 16384; i++)
            {
                if (func3_cile[i] == y11)
                {
                    used_index_y3[*index1] = i;
                    y3[*index1] = func3_variant[i];
                }
            }

            for (int i = 0; i < 16384; i++)
            {
                if (func3_cile[i] == x22)
                {
                    used_index_x3[*index2] = i;
                    x3[*index2] = func3_variant[i];
                    break;
                }
            }
            for (int i = 0; i < 16384; i++)
            {
                if (func3_cile[i] == y22)
                {
                    used_index_y3[*index2] = i;
                    y3[*index2] = func3_variant[i];
                }
            }
        }
    }
    lw3 get_binary(double* x_1, double* y_1, double* x_2, double* y_2, int* index1, int* index2)
    {
        if (a == 1)
        {
            int binnaryx1 = func1_cile[used_index_x1[*index1]], binnaryx2 = func1_cile[used_index_x1[*index2]], binnaryy1 = func1_cile[used_index_y1[*index1]], binnaryy2 = func1_cile[used_index_y1[*index2]];
            for (int i = 0; i < 9; i++)
            {
                binnary_x1[i] = binnaryx1 % 2;
                binnaryx1 = binnaryx1 / 2;
                binnary_x2[i] = binnaryx2 % 2;
                binnaryx2 = binnaryx2 / 2;
                binnary_y1[i] = binnaryy1 % 2;
                binnaryy1 = binnaryy1 / 2;
                binnary_y2[i] = binnaryy1 % 2;
                binnaryy2 = binnaryy2 / 2;
            }
        }
        if (a == 2)
        {
            int binnaryx1 = func2_cile[used_index_x2[*index1]], binnaryx2 = func2_cile[used_index_x2[*index2]], binnaryy1 = func2_cile[used_index_y2[*index1]], binnaryy2 = func2_cile[used_index_y2[*index2]];
            for (int i = 0; i < 6; i++)
            {
                binnary_x21[i] = binnaryx1 % 2;
                binnaryx1 = binnaryx1 / 2;
                binnary_x22[i] = binnaryx2 % 2;
                binnaryx2 = binnaryx2 / 2;
                binnary_y21[i] = binnaryy1 % 2;
                binnaryy1 = binnaryy1 / 2;
                binnary_y22[i] = binnaryy1 % 2;
                binnaryy2 = binnaryy2 / 2;
            }
        }
        if (a == 3)
        {
            int binnaryx1 = func3_cile[used_index_x3[*index1]], binnaryx2 = func3_cile[used_index_x3[*index2]], binnaryy1 = func3_cile[used_index_y3[*index1]], binnaryy2 = func3_cile[used_index_y3[*index2]];
            for (int i = 0; i < 14; i++)
            {
                binnary_x31[i] = binnaryx1 % 2;
                binnaryx1 = binnaryx1 / 2;
                binnary_x32[i] = binnaryx2 % 2;
                binnaryx2 = binnaryx2 / 2;
                binnary_y31[i] = binnaryy1 % 2;
                binnaryy1 = binnaryy1 / 2;
                binnary_y32[i] = binnaryy1 % 2;
                binnaryy2 = binnaryy2 / 2;
            }
        }
        return *this;
    }
    lw3 get_binary(double* x, double* y, int* index1)
    {
        if (a == 1)
        {
            int binnaryx1 = func1_cile[used_index_x1[*index1]], binnaryy1 = func1_cile[used_index_y1[*index1]];
            for (int i = 0; i < 9; i++)
            {
                binnary_x1[i] = binnaryx1 % 2;
                binnaryx1 = binnaryx1 / 2;
                binnary_y1[i] = binnaryy1 % 2;
                binnaryy1 = binnaryy1 / 2;
            }
        }
        if (a == 2)
        {
            int binnaryx1 = func2_cile[used_index_x2[*index1]], binnaryy1 = func2_cile[used_index_y2[*index1]];
            for (int i = 0; i < 6; i++)
            {
                binnary_x21[i] = binnaryx1 % 2;
                binnaryx1 = binnaryx1 / 2;
                binnary_y21[i] = binnaryy1 % 2;
                binnaryy1 = binnaryy1 / 2;
            }
        }
        if (a == 3)
        {
            int binnaryx1 = func3_cile[used_index_x3[*index1]], binnaryy1 = func3_cile[used_index_y3[*index1]];
            for (int i = 0; i < 14; i++)
            {
                binnary_x31[i] = binnaryx1 % 2;
                binnaryx1 = binnaryx1 / 2;
                binnary_y31[i] = binnaryy1 % 2;
                binnaryy1 = binnaryy1 / 2;
            }
        }
        return *this;
    }
    void interface()
    {
        short b;
        double s = q;
        bool fl = true;
        while (fl)
        {
            for (int i = 0; i < 16384; i++)
            {
                func3_cile.push_back(0);
                func3_variant.push_back(0);
            }
            for (int i = 0; i < 16384; i++)
            {
                func3_cile[i] = i;
                func3_variant[i] = s;
                s += 0.006836;
            }
            cout << "Input number of function" << endl;
            cout << "1 - Bill" << endl << "2 - Goldman-Price" << endl << "3 - Eggholder" << endl << "4 - Exit" << endl;
            cin >> b;
            this->set_a(b);
            this->rand_inp();
            switch (b)
            {
            case 1:
                this->main_func();
                break;
            case 2:
                this->main_func();
                break;
            case 3:
                this->main_func();
                break;
            case 4:
                fl = false;
                break;
            }
            global_min.clear();
            func3_cile.clear();
            func3_variant.clear();
            s = q;
        }
    }
};

int main()
{
    lw3 t;
    t.interface();
}

#include <iostream>
#include <fstream>
#include <unistd.h>
#include <sys/stat.h>
#include <stdio.h>
#include <vector>
#include <cmath>
#include <random>
#include <string.h>      //文字列の代入に使う
#include <bits/stdc++.h> //piの利用で必要(M_PI)
#include <Eigen/Eigenvalues>
//using Eigen::MatrixXd;

class ParticleInfo
{
public:
    uint32_t id;
    uint32_t type;
    /*position*/
    double posx;
    double posy;
    double posz;
    /*velocity*/
    double velx;
    double vely;
    double velz;
    /*結合*/
    uint32_t bond_pair[2];
    uint32_t bond_type[2];
    uint32_t nbond = 0;
    /*アングル*/
    uint32_t angle_pair[2][3];
    uint32_t angle_type[2];
    uint32_t nangle = 0;
    //ベシクル内の水かどうかwater1,2
    bool vesicle_water1 = false;
    bool vesicle_water2 = false;

    //sortを利用するために定義
    bool operator<(const ParticleInfo &another) const
    {
        //メンバ変数であるnum1で比較した結果を
        //この構造体の比較とする
        return id < another.id;
    }
};

class CenterOfGravity
{
public:
    /*position*/
    double x;
    double y;
    double z;
    /*数えた脂質粒子の個数*/
    int num;
    CenterOfGravity()
    {
        x = y = z = 0;
        num = 0;
    }
};

int main(int argc, char *argv[])
{
    std::vector<ParticleInfo> pinfo;
    ParticleInfo temp_info;
    /*
    
    
    
    
    座標の読み込みを行う．*/
    std::ifstream ifs0(argv[1]);
    if (!ifs0)
    {
        std::cerr << "error0" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    //いらないはじめの3行を捨てる．
    std::string delete_str[3];
    for (uint32_t i = 0; i < 3; i++)
    {
        std::getline(ifs0, delete_str[i]);
    }
    //ファイルの読み込み　粒子番号,粒子種は0から開始する．
    while (ifs0 >> temp_info.id >> temp_info.type >> temp_info.posx >> temp_info.posy >> temp_info.posz)
    {
        temp_info.id--;
        temp_info.type--;
        pinfo.push_back(temp_info);
    }
    ifs0.close();
    //はじめの文字列を読み込む
    double box_sx, box_sy, box_sz, box_ex, box_ey, box_ez, box_wt;
    sscanf(delete_str[0].c_str(), "'box_sx=%lf box_sy=%lf box_sz=%lf box_ex=%lf box_ey=%lf box_ez=%lf box_wt=%lf",
           &box_sx, &box_sy, &box_sz, &box_ex, &box_ey, &box_ez, &box_wt);
    //    std::cout <<std::setprecision(10)<< box_sx << " " << box_sy << " " << box_sz << " " << box_ex << " " << box_ey << " " << box_ez << " " << box_wt << std::endl;
    std::sort(pinfo.begin(), pinfo.end()); //classでオペレータを定義して利用している．
    /*




    ベシクル水の重心を計算する．*/
    CenterOfGravity center_water1, center_water2;
    for (int i = 0; i < pinfo.size(); i++)
    {
        if (pinfo.at(i).type == 5) //水粒子
        {
            center_water1.x += pinfo.at(i).posx;
            center_water1.y += pinfo.at(i).posy;
            center_water1.z += pinfo.at(i).posz;
            center_water1.num++;
        }
    }
    center_water1.x /= (double)center_water1.num;
    center_water1.y /= (double)center_water1.num;
    center_water1.z /= (double)center_water1.num;
    std::cout << "center of water5  " << center_water1.z << std::endl;
    for (int i = 0; i < pinfo.size(); i++)
    {
        if (pinfo.at(i).type == 6) //水粒子
        {
            center_water2.x += pinfo.at(i).posx;
            center_water2.y += pinfo.at(i).posy;
            center_water2.z += pinfo.at(i).posz;
            center_water2.num++;
        }
    }
    center_water2.x /= (double)center_water2.num;
    center_water2.y /= (double)center_water2.num;
    center_water2.z /= (double)center_water2.num;
    std::cout << "center of water6  " << center_water2.z << std::endl;
    /*



    
    ベシクル脂質の重心を計算する．*/
    CenterOfGravity center_lipid1, center_lipid2;
    for (int i = 0; i < pinfo.size(); i++)
    {
        if (pinfo.at(i).type == 1) //水粒子
        {
            center_lipid1.x += pinfo.at(i).posx;
            center_lipid1.y += pinfo.at(i).posy;
            center_lipid1.z += pinfo.at(i).posz;
            center_lipid1.num++;
        }
    }
    center_lipid1.x /= (double)center_lipid1.num;
    center_lipid1.y /= (double)center_lipid1.num;
    center_lipid1.z /= (double)center_lipid1.num;
    std::cout << "center of lipid1  " << center_lipid1.z << std::endl;
    for (int i = 0; i < pinfo.size(); i++)
    {
        if (pinfo.at(i).type == 4) //水粒子
        {
            center_lipid2.x += pinfo.at(i).posx;
            center_lipid2.y += pinfo.at(i).posy;
            center_lipid2.z += pinfo.at(i).posz;
            center_lipid2.num++;
        }
    }
    center_lipid2.x /= (double)center_lipid2.num;
    center_lipid2.y /= (double)center_lipid2.num;
    center_lipid2.z /= (double)center_lipid2.num;
    std::cout << "center of lipid4  " << center_lipid2.z << std::endl;

    /*



    
    ベシクル水脂質の重心を計算する．*/
    CenterOfGravity center_waterlipid1, center_waterlipid2;
    for (int i = 0; i < pinfo.size(); i++)
    {
        if (pinfo.at(i).type == 5 || pinfo.at(i).type == 1) //水粒子
        {
            center_waterlipid1.x += pinfo.at(i).posx;
            center_waterlipid1.y += pinfo.at(i).posy;
            center_waterlipid1.z += pinfo.at(i).posz;
            center_waterlipid1.num++;
        }
    }
    center_waterlipid1.x /= (double)center_waterlipid1.num;
    center_waterlipid1.y /= (double)center_waterlipid1.num;
    center_waterlipid1.z /= (double)center_waterlipid1.num;
    std::cout << "center of waterlipid51  " << center_waterlipid1.z << std::endl;
    for (int i = 0; i < pinfo.size(); i++)
    {
        if (pinfo.at(i).type == 6 || pinfo.at(i).type == 4) //水粒子
        {
            center_waterlipid2.x += pinfo.at(i).posx;
            center_waterlipid2.y += pinfo.at(i).posy;
            center_waterlipid2.z += pinfo.at(i).posz;
            center_waterlipid2.num++;
        }
    }
    center_waterlipid2.x /= (double)center_waterlipid2.num;
    center_waterlipid2.y /= (double)center_waterlipid2.num;
    center_waterlipid2.z /= (double)center_waterlipid2.num;
    std::cout << "center of waterlipid64  " << center_waterlipid2.z << std::endl;
    /*




    
    同系方向の粒子種を数える*/
    class Gr
    {
    public:
        double distance = 0;
        uint32_t number_of_5 = 0;
        uint32_t number_of_6 = 0;
    };
    /*自分で決めるパラメータ*/
    double dr = 0.5; //同系方向の円形ボックスの厚み．この中にいくつ粒子が入っているか数える．
    /*確定しているもの*/
    double search_r = 0.0;
    double r_now;
    double maxz = 0, minz = 100;
    for (int i = 0; i < pinfo.size(); i++)
    {
        maxz = maxz > pinfo.at(i).posz ? maxz : pinfo.at(i).posz;
        minz = minz < pinfo.at(i).posz ? minz : pinfo.at(i).posz;
    }
    //    std::cout << "粒子の最大座標と最小座標を示す " << maxz << " " << minz << std::endl;
    //    uint32_t num_of_gr = (uint32_t)((box_ez - box_sz) / dr) + 1;
    uint32_t num_of_gr = (uint32_t)((maxz - minz) / dr) + 1;

    std::vector<Gr> gr_info(num_of_gr);
    //    Gr temp_gr_info;
    for (int i = 0; i < pinfo.size(); i++)
    {
        //        std::cout << i <<" "<<pinfo.at(i).posz<<" "<<(uint32_t)(r_now / dr)<< std::endl;///////////////////////////////////////
        r_now = pinfo.at(i).posz + (-minz); //たしか，配列のインデックスがマイナスになるのを防ぐため．
        if (pinfo.at(i).type == 5)
            gr_info.at((uint32_t)(r_now / dr)).number_of_5++;
        else if (pinfo.at(i).type == 6)
            gr_info.at((uint32_t)(r_now / dr)).number_of_6++;
    }
    /*





    交点の導出*/
    uint32_t min_gr_num = 5000;
    uint32_t gr_temp = 0;
    uint32_t min_gr_id = 0;
    for (int i = (uint32_t)(40 / dr); i < (uint32_t)(55 / dr); i++)
    {
        gr_temp = std::abs(int(gr_info.at(i).number_of_5 - gr_info.at(i).number_of_6));
        //        std::cout<<"gr_temp "<<gr_temp<<" "<<i<<std::endl;
        if (min_gr_num > gr_temp)
        {
            min_gr_id = i;
            min_gr_num = gr_temp;
        }
    }
    std::cout << "min_gr_id*dr " << min_gr_id * dr << std::endl;
    /*





    導出した交点から粒子数を数える．*/
    uint32_t num_w5 = 0, num_w5out = 0, num_w6 = 0, num_w6out = 0;
    //water5を数える
    for (int i = 0; i < min_gr_id; i++)
    {
        num_w5 += gr_info.at(i).number_of_5;
        num_w6out += gr_info.at(i).number_of_6;
    }
    //water6を数える
    for (int i = num_of_gr - 1; i > min_gr_id; i--)
    {
        num_w5out += gr_info.at(i).number_of_5;
        num_w6 += gr_info.at(i).number_of_6;
    }
    std::cout << "num_w5 " << num_w5 << " num_w5out " << num_w5out << " num_w6 " << num_w6 << " num_w6out " << num_w6out << std::endl;
    /*





    慣性テンソルの導出*/
    //    double I11 = 0, I12 = 0, I13 = 0;
    //    double I21 = 0, I22 = 0, I23 = 0;
    //    double I31 = 0, I32 = 0, I33 = 0;
    Eigen::MatrixXcf A0 = Eigen::MatrixXcf::Zero(3, 3);
    Eigen::MatrixXcf A3 = Eigen::MatrixXcf::Zero(3, 3);
    Eigen::MatrixXcf D0(3, 3); //対角化後
    Eigen::MatrixXcf D3(3, 3);
    for (int i = 0; i < pinfo.size(); i++)
    {
        if (pinfo.at(i).type == 0)
        {
            A0(0, 0) += pow(pinfo.at(i).posy, 2.0) + pow(pinfo.at(i).posz, 2.0);
            A0(0, 1) -= pinfo.at(i).posx * pinfo.at(i).posy;
            A0(0, 2) -= pinfo.at(i).posx * pinfo.at(i).posz;

            A0(1, 0) -= pinfo.at(i).posy * pinfo.at(i).posx;
            A0(1, 1) += pow(pinfo.at(i).posx, 2.0) + pow(pinfo.at(i).posz, 2.0);
            A0(1, 2) -= pinfo.at(i).posy * pinfo.at(i).posz;

            A0(2, 0) -= pinfo.at(i).posz * pinfo.at(i).posx;
            A0(2, 1) -= pinfo.at(i).posz * pinfo.at(i).posy;
            A0(2, 2) += pow(pinfo.at(i).posx, 2.0) + pow(pinfo.at(i).posy, 2.0);
        }
        if (pinfo.at(i).type == 3)
        {
            A3(0, 0) += pow(pinfo.at(i).posy, 2.0) + pow(pinfo.at(i).posz, 2.0);
            A3(0, 1) -= pinfo.at(i).posx * pinfo.at(i).posy;
            A3(0, 2) -= pinfo.at(i).posx * pinfo.at(i).posz;

            A3(1, 0) -= pinfo.at(i).posy * pinfo.at(i).posx;
            A3(1, 1) += pow(pinfo.at(i).posx, 2.0) + pow(pinfo.at(i).posz, 2.0);
            A3(1, 2) -= pinfo.at(i).posy * pinfo.at(i).posz;

            A3(2, 0) -= pinfo.at(i).posz * pinfo.at(i).posx;
            A3(2, 1) -= pinfo.at(i).posz * pinfo.at(i).posy;
            A3(2, 2) += pow(pinfo.at(i).posx, 2.0) + pow(pinfo.at(i).posy, 2.0);
        }
    }
    //ソルバーの生成(http://modeling-lab.pinoko.jp/software_buster/environment/eigen.html)
    Eigen::ComplexEigenSolver<Eigen::MatrixXcf> ces;
    ces.compute(A0);
    D0 = ces.eigenvalues().asDiagonal();
    ces.compute(A3);
    D3 = ces.eigenvalues().asDiagonal();
    std::cout << "A0固有値をならべて対角行列にしたもの" << std::endl; //複素数は (実部, 虚部) と出力される．
    std::cout << D0 << std::endl
              << std::endl;
    std::cout << "A3固有値をならべて対角行列にしたもの" << std::endl; //複素数は (実部, 虚部) と出力される．
    std::cout << D3 << std::endl
              << std::endl;
    /*





    ファイルの出力*/

    //vel_file
    FILE *fpo1;
    fpo1 = fopen(argv[2], "w");
    if (fpo1 == NULL)
    {
        printf("ERROR_file_output\n");
        return -1;
    }
    for (int i = 0; i < gr_info.size(); i++)
    {
        fprintf(fpo1, "%lf   %d   %d  \n",
                dr * i,
                gr_info.at(i).number_of_5,
                gr_info.at(i).number_of_6);
    }
    fclose(fpo1);
    return 0;
}
// Copyright 2018 Zeyu Zhong
// Lincese(MIT)
// Author: Zeyu Zhong
// Date: 2018.9.2

#include "../src/Clustering.hpp"
#include <fstream>
#include <sstream>

int main() {
    std::ifstream fin("../../point.txt");
    std::string ptline;
    double x, y, z;
    std::vector<Eigen::Vector3d> points;
    while (getline(fin, ptline)) {
        std::stringstream ss(ptline);
        ss >> x >> y >> z;
        points.push_back(Eigen::Vector3d(x, y, z));
    }

    std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> lines;
    int numPoints = points.size();

    for (int i = 0; i < numPoints; i = i+2) {
        std::pair<Eigen::Vector3d, Eigen::Vector3d> line(points[i], points[i+1]);
        lines.push_back(line);
    }

    CCluster cluster;
    // Step1: local cluster, all lines into different local cluster that represent the same frag
    std::vector<std::list<int> > clusters = cluster.localCluster(lines, 4);

    std::ofstream fout("../../cluster.txt");
    int i = 0;
    for (std::vector<std::list<int> >::iterator iter1 = clusters.begin(); iter1 != clusters.end(); iter1++, i++) {
        std::list<int> linesID = *iter1;

        std::list<int>::iterator iter2;
        for (iter2 = linesID.begin(); iter2!= linesID.end(); iter2++) {
            Eigen::Vector3d a = lines[*iter2].first;
            Eigen::Vector3d b = lines[*iter2].second;
            fout << a[0] << " " << a[1] << " " << a[2] << " " << i << std::endl;
            fout << b[0] << " " << b[1] << " " << b[2] << " " << i << std::endl;
        }
    }
    return 0;
}

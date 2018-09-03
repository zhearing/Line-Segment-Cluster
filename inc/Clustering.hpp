// Copyright 2018 Zeyu Zhong
// Lincese(MIT)
// Author: Zeyu Zhong
// Date: 2018.8.5
// Update: 2018.9.2

#ifndef INC_CLUSTERING_HPP_
#define INC_CLUSTERING_HPP_

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include <utility>
#include "eigen3/Eigen/Eigen"

struct VERTEX {
    int rank_;
    int clusterID_;
    int size_;
};

class cluGraph {
 public:
    explicit cluGraph(int vertexNum): vertexNum_(vertexNum) {
        vertexs_.resize(vertexNum);

        for (int i = 0; i < vertexNum; ++i) {
            vertexs_[i].rank_ = 0;
            vertexs_[i].size_ = 1;
            vertexs_[i].clusterID_ = i;
        }
    }

    // find clusterID for vertex node
    int findCluID(int vertexID) {
        int y = vertexID;
        while (vertexID != vertexs_[vertexID].clusterID_)
            vertexID = vertexs_[vertexID].clusterID_;

        vertexs_[y].clusterID_ = vertexID;
        return vertexID;
    }

    // joins two vertexs into one cluster
    void join(int x, int y) {
        if (vertexs_[x].rank_ > vertexs_[y].rank_) {
            vertexs_[y].clusterID_ = x;
            vertexs_[x].size_ += vertexs_[y].size_;
        } else {
            vertexs_[x].clusterID_ = y;
            vertexs_[y].size_ += vertexs_[x].size_;
            if (vertexs_[x].rank_ == vertexs_[y].rank_)
                vertexs_[y].rank_++;
        }
        vertexNum_--;
    }

    int size(int x) const { return vertexs_[x].size_; }

 private:
    int vertexNum_;
    std::vector<VERTEX> vertexs_;
};

class CCluster {
 public:
    // compute the min distance from a point to a line segment
    double pt2frag(std::pair<Eigen::Vector3d, Eigen::Vector3d> frag,
                   Eigen::Vector3d c) {
        Eigen::Vector3d a = frag.first;
        Eigen::Vector3d b = frag.second;
        Eigen::Vector3d ac = c - a;
        Eigen::Vector3d ab = b - a;
        Eigen::Vector3d Iab = ab / ab.norm();
        Eigen::Vector3d ad = ac.dot(Iab)*Iab;
        Eigen::Vector3d d = a + ad;
        Eigen::Vector3d da = -ad;
        Eigen::Vector3d db = b - d;
        if (da.dot(db) < 0) {
           return ac.cross(Iab).norm();
        }

        Eigen::Vector3d bc = b - c;
        ac[1] = ac[1] * 2.0;
        bc[1] = bc[1] * 2.0;
        double lenac = ac.norm();
        double lenbc = bc.norm();

        return lenac < lenbc ? lenac : lenbc;
    }

    bool ptOnFrag(std::pair<Eigen::Vector3d, Eigen::Vector3d> frag,
                  Eigen::Vector3d c) {
        Eigen::Vector3d a = frag.first;
        Eigen::Vector3d b = frag.second;
        Eigen::Vector3d ab = b - a;
        Eigen::Vector3d ac = c - a;
        Eigen::Vector3d Iab = ab / ab.norm();    // identity vector of ab
        Eigen::Vector3d ad = ac.dot(Iab) * Iab;
        Eigen::Vector3d d = a + ad;
        Eigen::Vector3d da = -ad;
        Eigen::Vector3d db = b - d;

        if (da.dot(db) < 0)
          return true;

        return false;
    }

    // judege the similarity between two line segments
    bool line2line(const std::pair<Eigen::Vector3d, Eigen::Vector3d> lines1,
                   const std::pair<Eigen::Vector3d, Eigen::Vector3d> lines2,
                   const double distTrd) {
        Eigen::Vector2d a = Eigen::Vector2d(lines1.first[0], lines1.first[1]);
        Eigen::Vector2d b = Eigen::Vector2d(lines1.second[0], lines1.second[1]);
        Eigen::Vector2d c = Eigen::Vector2d(lines2.first[0], lines2.first[1]);
        Eigen::Vector2d d = Eigen::Vector2d(lines2.second[0], lines2.second[1]);

        double x1 = a[0]; double y1 = a[1];
        double x2 = b[0]; double y2 = b[1];
        double x3 = c[0]; double y3 = c[1];
        double x4 = d[0]; double y4 = d[1];

        double D = (y1-y2)*(x4-x3) - (y3-y4)*(x2-x1);
        double D1 = (x2*y1 - x1*y2)*(x4-x3) - (x4*y3 - x3*y4)*(x2-x1);
        double D2 = (y1-y2)*(x4*y3-x3*y4) - (y3-y4)*(x2*y1-x1*y2);
        double crossX = D1 / D;
        double crossZ = D2 / D;
        Eigen::Vector2d p = Eigen::Vector2d(crossX, crossZ);

        Eigen::Vector2d pa = a - p;
        Eigen::Vector2d pb = b - p;
        Eigen::Vector2d pc = c - p;
        Eigen::Vector2d pd = d - p;

        bool onAB = true;
        bool onCD = true;

        if (pa.dot(pb) > 0)
            onAB = false;

        if (pc.dot(pd) > 0)
            onCD = false;

        if (D == 0) {
            onAB = false;
            onCD = false;
        }

        bool isCross = false;
        double dist = 100;

        if (onAB && onCD) {    // the two line segments do cross
            isCross = true;
            return true;
            dist = 0;
        } else {                 // the two line segments does not cross
           // find the distance of a to cd, b to cd, c to ab and d to ab,
           // return the smallest one
            double distA = pt2frag(lines2, lines1.first);
            double distB = pt2frag(lines2, lines1.second);
            double distC = pt2frag(lines1, lines2.first);
            double distD = pt2frag(lines1, lines2.second);
            double minDist1 = distA < distB ? distA : distB;
            // std::cout << minDist1 << std::endl;
            double minDist2 = distC < distD ? distC : distD;
            // std::cout << minDist2 << std::endl;
            isCross = false;
            dist = minDist1 < minDist2 ? minDist1 : minDist2;
        }

        if (dist < distTrd)
            return true;

        // double angle = 0;
        // Eigen::Vector2d ab = b-a;
        // Eigen::Vector2d cd = d-c;
        // angle = acos(fabs(ab.dot(cd))/(ab.norm()*cd.norm()));
        //// std::cout << angle << std::endl;

        // if (dist < distTrd && angle < angleTrd)
        //    return true;

        return false;
    }

    // judge is the two line segments ab and cd have common region,
    // regardless oftheir horizontal offset
    bool frag2frag(const std::pair<Eigen::Vector3d, Eigen::Vector3d> ab,
                   const std::pair<Eigen::Vector3d, Eigen::Vector3d> cd) {
        Eigen::Vector3d a = ab.first;
        Eigen::Vector3d b = ab.second;
        Eigen::Vector3d c = cd.first;
        Eigen::Vector3d d = cd.second;
        bool isOnA = ptOnFrag(cd, a);
        bool isOnB = ptOnFrag(cd, b);
        bool isOnC = ptOnFrag(ab, c);
        bool isOnD = ptOnFrag(ab, d);

        if (isOnA || isOnB || isOnC || isOnD)
            return true;

        return false;
    }

    // compute the distance of pt's projection on line who go through P and direction is dir
    double pt2G(const Eigen::Vector3d P, const Eigen::Vector3d dir,
                const Eigen::Vector3d a) {
        Eigen::Vector3d pa = a - P;
        return pa.dot(dir);
    }

    // find the foot point of point pt to line who go through P and direction is dir
    Eigen::Vector3d ptProj2Line(const Eigen::Vector3d P, const Eigen::Vector3d dir,
                                const Eigen::Vector3d a) {
        Eigen::Vector3d pa = a - P;
        Eigen::Vector3d pb = pa.dot(dir)*dir;
        Eigen::Vector3d b = P + pb;
        return b;
    }

    // use a frag to represent all frags that represent the same lien segment
    // get from local cluster
    std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> mergeFrag(
        const std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> lines,
        const std::vector<std::list<int> > clusters) {
        std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> clusterRef;
        clusterRef.reserve(clusters.size());

        for (std::vector<std::list<int> >::const_iterator iter = clusters.begin();
            iter != clusters.end(); iter++) {
                std::list<int> linesIDs = *iter;
                Eigen::Vector3d P(0, 0, 0);
                Eigen::Vector3d direction(0, 0, 0);
                int n = 2*linesIDs.size();
                Eigen::MatrixXd points(3, n);
                std::vector<Eigen::Vector3d> pts;
                pts.reserve(n);

                int i = 0;
                for (std::list<int>::iterator iter2 = linesIDs.begin();
                    iter2!= linesIDs.end(); iter2++, i+=2) {
                    Eigen::Vector3d a = lines[*iter2].first;
                    Eigen::Vector3d b = lines[*iter2].second;
                    P += a;
                    P += b;
                    direction += b - a;
                    points(0, i) = a[0];
                    points(1, i) = a[1];
                    points(2, i) = a[2];
                    points(0, i+1) = b[0];
                    points(1, i+1) = b[1];
                    points(2, i+1) = b[2];
                    pts.push_back(a);
                    pts.push_back(b);
                }

                // gravity point
                P /= static_cast<double>(n);

                // use SVD to find direction
                Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n) -
                    1.0 / static_cast<double>(n)*Eigen::MatrixXd::Constant(n, n, 1.0);
                Eigen::MatrixXd Scat = points*I*points.transpose();
                Eigen::JacobiSVD<Eigen::MatrixXd> svd(Scat, Eigen::ComputeThinU);

                Eigen::MatrixXd U;
                Eigen::VectorXd S;
                U = svd.matrixU();
                S = svd.singularValues();
                int maxSValuePos;
                S.maxCoeff(&maxSValuePos);

                // Eigen::Vector3d dir = Eigen::Vector3d(U(0,maxSValuePos),
                //                                      U(1,maxSValuePos), U(2,maxSValuePos));
                Eigen::Vector3d dir = direction;
                dir.normalize();

                // for all two end point of a line segment, project it to dir and find two end of all of them
                std::vector<double> dists;
                dists.reserve(n);
                for (std::vector<Eigen::Vector3d>::iterator iter2 = pts.begin();
                    iter2 != pts.end(); iter2++) {
                    double dis = pt2G(P, dir, *iter2);
                    dists.push_back(dis);
                }

                std::vector<double>::iterator iterMin = std::min_element(dists.begin(),
                    dists.end());
                std::vector<double>::iterator iterMax = std::max_element(dists.begin(),
                    dists.end());
                int minIdx = std::distance(dists.begin(), iterMin);
                int maxIdx = std::distance(dists.begin(), iterMax);

                Eigen::Vector3d a = ptProj2Line(P, dir, pts[minIdx]);
                Eigen::Vector3d b = ptProj2Line(P, dir, pts[maxIdx]);

                clusterRef.push_back(std::pair<Eigen::Vector3d, Eigen::Vector3d>(a, b));
        }

        return clusterRef;
    }

    std::vector<std::list<int> > segCluster(std::vector<std::list<int> > clusters,
                                            const std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> clusterRef) {
        if (clusterRef.size() != clusters.size()) {
            // std::cout << "merge frag is wrong, number of clusters dismatch!";
            abort();
        }

        int numCluster = clusterRef.size();
        std::vector<std::pair<int, int> > A;

        for (int i = 0; i < numCluster; i++) {
            std::pair<Eigen::Vector3d, Eigen::Vector3d> ab = clusterRef[i];

            for (int j = i+1; j < numCluster; ++j) {
                std::pair<Eigen::Vector3d, Eigen::Vector3d> cd = clusterRef[j];
                bool isCommon = frag2frag(ab, cd);

                if (isCommon) {
                    A.push_back(std::pair<int, int>(i, j));
                    A.push_back(std::pair<int, int>(j, i));
                }
            }
        }

        // return a conjoint graph
        cluGraph graph(numCluster);

        std::vector<std::pair<int, int>>::iterator iter = A.begin();
        for (; iter != A.end(); iter++) {
            std::pair<int, int> edge = *iter;
            int a = graph.findCluID(edge.first);
            int b = graph.findCluID(edge.second);

            if (a != b)
                graph.join(a, b);
        }

        // outout the conjoint graph
        std::map<int, std::list<int> > mRegionCluster;

        for (int i = 0; i < numCluster; ++i) {
            int regionID = graph.findCluID(i);
            mRegionCluster[regionID].push_back(i);
        }

        std::vector<std::list<int> > regionClusters;
        regionClusters.resize(mRegionCluster.size());

        int i = 0;
        for (std::map<int, std::list<int> >::iterator iter = mRegionCluster.begin();
            iter != mRegionCluster.end(); iter++, i++) {
                std::list<int> clusterIDs = iter->second;

                for (std::list<int>::iterator iter2 = clusterIDs.begin();
                    iter2 != clusterIDs.end(); iter2++) {
                        // regionClusters[i].insert(regionClusters[i].end(), clusters[*iter2]);
                        regionClusters[i].insert(regionClusters[i].end(),
                            clusters[*iter2].begin(), clusters[*iter2].end());
                }
        }

        return regionClusters;
    }

    // local cluster frags to groups that all represent the same line segment
    std::vector<std::list<int> > localCluster(const std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> lines, double distTrd) {
        // get relation edge matrix A
        int num_lines = lines.size();
        std::vector<std::pair<int, int>> A;

        for (int i = 0; i < num_lines; ++i) {
            std::pair<Eigen::Vector3d, Eigen::Vector3d> ab = lines[i];

            for (int j = i+1; j < num_lines; ++j) {
                std::pair<Eigen::Vector3d, Eigen::Vector3d> cd = lines[j];
                bool isConnect = line2line(ab, cd, distTrd);

                if (isConnect) {
                    A.push_back(std::pair<int, int>(i, j));
                    A.push_back(std::pair<int, int>(j, i));
                }
            }
        }

        // return a conjoint graph
        cluGraph graph(num_lines);

        std::vector<std::pair<int, int>>::iterator iter = A.begin();
        for (; iter != A.end(); iter++) {
            std::pair<int, int> edge = *iter;
            int a = graph.findCluID(edge.first);
            int b = graph.findCluID(edge.second);

            if (a != b)
                graph.join(a, b);
        }

        // outout the conjoint graph
        std::map<int, std::list<int> > cluster2lines;

        for (int i = 0; i < num_lines; ++i) {
            int clusterID = graph.findCluID(i);
            cluster2lines[clusterID].push_back(i);
        }

        std::vector<std::list<int> > clusters;
        clusters.reserve(cluster2lines.size());

        for (std::map<int, std::list<int> >::iterator iter = cluster2lines.begin();
            iter != cluster2lines.end(); iter++) {
            clusters.push_back((*iter).second);
        }

        return clusters;
    }
};

#endif  // INC_CLUSTERING_HPP_


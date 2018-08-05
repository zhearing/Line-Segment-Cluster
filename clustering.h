//
// Created by zyzhong on 18-8-5.
//

#ifndef LINE_SEGMENT_CLUSTER_CLUSTERING_H
#define LINE_SEGMENT_CLUSTER_CLUSTERING_H

#include <vector>

struct VERTEX
{
    int rank_;
    int clusterID_;
    int size_;
};

class cluGraph
{
public:
    cluGraph(int vertexNum): vertexNum_(vertexNum)
    {
        vertexs_.resize(vertexNum);

        for(int i = 0; i < vertexNum; ++i)
        {
            vertexs_[i].rank_ = 0;
            vertexs_[i].size_ = 1;
            vertexs_[i].clusterID_ = i;
        }

    }

    // find clusterID for vertex node
    int findCluID(int vertexID)
    {
        int y = vertexID;
        while(vertexID != vertexs_[vertexID].clusterID_)
            vertexID = vertexs_[vertexID].clusterID_;

        vertexs_[y].clusterID_ = vertexID;
        return vertexID;
    }

    // joins two vertexs into one cluster
    void join(int x, int y)
    {
        if(vertexs_[x].rank_ > vertexs_[y].rank_)
        {
            vertexs_[y].clusterID_ = x;
            vertexs_[x].size_ += vertexs_[y].size_;
        }
        else
        {
            vertexs_[x].clusterID_ = y;
            vertexs_[y].size_ += vertexs_[x].size_;
            if(vertexs_[x].rank_ == vertexs_[y].rank_)
                vertexs_[y].rank_++;
        }
        vertexNum_--;
    }

    int size(int x) const { return vertexs_[x].size_; }

private:
    int vertexNum_;
    std::vector<VERTEX> vertexs_;
};

#endif //LINE_SEGMENT_CLUSTER_CLUSTERING_H

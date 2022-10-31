#include <algorithm>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>
#include <tuple>
#include <vector>
#include "utils/types.h"
#include "utils/utils.h"
#include "graph/graph.h"

Graph::Graph()
: edge_count_(0)
, vlabel_count_(0)
, elabel_count_(0)
, neighbors_{}
, elabels_{}
, updates_{}
, vlabels_{}
{}

void Graph::AddVertex(uint id, uint label)
{
    if (id >= vlabels_.size())
    {
        vlabels_.resize(id + 1, NOT_EXIST);  //如果插入的点超出了vector原定的尺寸，则扩大vector的尺寸，并用NOT_EXIST这个数据来填充那些扩大的空间
        vlabels_[id] = label;
        neighbors_.resize(id + 1);
        elabels_.resize(id + 1);
    }
    else if (vlabels_[id] == NOT_EXIST)
    {
        vlabels_[id] = label;
    }
    
    vlabel_count_ = std::max(vlabel_count_, label + 1);  //TODO: ？？？
    // print graph
    /*std::cout << "labels: ";
    for (uint i = 0; i < vlabels_.size(); i++)
    {
        std::cout << i << ":" << vlabels_[i] << " (";
        for (uint j = 0; j < neighbors_[i].size(); j++)
        {
            std::cout << neighbors_[i][j] << ":" << elabels_[i][j] << " ";
        }
        std::cout << ")" << std::endl;
    }*/
}

void Graph::RemoveVertex(uint id)
{
    vlabels_[id] = NOT_EXIST;
    neighbors_[id].clear();
    elabels_[id].clear();
}



void Graph::AddEdge(uint v1, uint v2, uint label)
{
    auto lower = std::lower_bound(neighbors_[v1].begin(), neighbors_[v1].end(), v2);   //通过lower_bound函数，使v1的邻居全是有序排列的
    if (lower != neighbors_[v1].end() && *lower == v2) return;   //如果在v1邻居中找到了大于等于v2的位置，并且v2此时已经在v1的nei中了，所以return
    
    size_t dis = std::distance(neighbors_[v1].begin(), lower);  //计算从v1的nei的起点~lower这个范围内的元素的个数
    neighbors_[v1].insert(lower, v2);  //在lower的地方插入v2 
    elabels_[v1].insert(elabels_[v1].begin() + dis, label);  //因为dis其实也是一个迭代器，所以可以和begin相加
    
    lower = std::lower_bound(neighbors_[v2].begin(), neighbors_[v2].end(), v1);  //无向图，所以更新完v1后更新v2
    dis = std::distance(neighbors_[v2].begin(), lower);
    neighbors_[v2].insert(lower, v1);
    elabels_[v2].insert(elabels_[v2].begin() + dis, label);  

    edge_count_++;     //更新两个节点的邻居，G图中的边则会多加一条
    elabel_count_ = std::max(elabel_count_, label + 1);
    // print graph
    /*std::cout << "labels: ";
    for (uint i = 0; i < vlabels_.size(); i++)
    {
        std::cout << i << ":" << vlabels_[i] << " (";
        for (uint j = 0; j < neighbors_[i].size(); j++)
        {
            std::cout << neighbors_[i][j] << ":" << elabels_[i][j] << " ";
        }
        std::cout << ")" << std::endl;
    }*/
}

void Graph::RemoveEdge(uint v1, uint v2)
{
    auto lower = std::lower_bound(neighbors_[v1].begin(), neighbors_[v1].end(), v2);
    if (lower == neighbors_[v1].end() || *lower != v2)   //如果没有在v1的邻居中找到v2，则无法删除边
    { 
        std::cout << "deletion error" << std::endl;
        exit(-1);
    }
    neighbors_[v1].erase(lower);  //从v1的邻居中删除v2
    elabels_[v1].erase(elabels_[v1].begin() + std::distance(neighbors_[v1].begin(), lower));  //从边标签中删除v2的标签
    
    lower = std::lower_bound(neighbors_[v2].begin(), neighbors_[v2].end(), v1);  //反向同理
    if (lower == neighbors_[v2].end() || *lower != v1)
    {
        std::cout << "deletion error" << std::endl;
        exit(-1);
    }
    neighbors_[v2].erase(lower);
    elabels_[v2].erase(elabels_[v2].begin() + std::distance(neighbors_[v2].begin(), lower));

    edge_count_--;    //执行成功后，边的总数减一
}

uint Graph::GetVertexLabel(uint u) const
{
    return vlabels_[u];
}

const std::vector<uint>& Graph::GetNeighbors(uint v) const   //返回的是一个vector的引用
{
    return neighbors_[v];
}

const std::vector<uint>& Graph::GetNeighborLabels(uint v) const
{
    return elabels_[v];
}

std::tuple<uint, uint, uint> Graph::GetEdgeLabel(uint v1, uint v2) const
{
    uint v1_label, v2_label, e_label;
    v1_label = GetVertexLabel(v1);
    v2_label = GetVertexLabel(v2);

    //从这里往下，都是再找与v1 v2相连的那条边的label
    const std::vector<uint> *nbrs;
    const std::vector<uint> *elabel;
    uint other;
    if (GetDegree(v1) < GetDegree(v2))   //GetDegree就是返回 neighbors_[v].size();
    {
        nbrs = &GetNeighbors(v1);
        elabel = &elabels_[v1];
        other = v2;
    }
    else
    {
        nbrs = &GetNeighbors(v2);
        elabel = &elabels_[v2];
        other = v1;
    }
    
    long start = 0, end = nbrs->size() - 1, mid;
    while (start <= end)    //因为nei在存储的时候，是通过顺序存储的，所以在这里用二分法去查找，速度更快
    {
        mid = (start + end) / 2;
        if (nbrs->at(mid) < other)
        {
            start = mid + 1;
        }
        else if (nbrs->at(mid) > other)
        {
            end = mid - 1;
        }
        else
        {
            e_label = elabel->at(mid);
            return {v1_label, v2_label, e_label};  //如果成功找到，与v1 v2相连的边的label返回
        }
    }
    return {v1_label, v2_label, -1};   //如果失败，则只返回v1 v2的label
}




uint Graph::GetDegree(uint v) const
{
    return neighbors_[v].size();
}




//返回直径：通过BFS获取图G的最长的路径
uint Graph::GetDiameter() const
{
    uint diameter = 0;
    for (uint i = 0u; i < NumVertices(); i++)       //NumVertices() ：return vlabels_.size()
    if (GetVertexLabel(i) != NOT_EXIST)   //GetVertexLabel(i) ：return vlabels_[u];
    {
        std::queue<uint> bfs_queue;      //因为BFS是一圈一圈往外遍历的，所以增加一个queue来实现一圈一圈的遍历效果
        std::vector<bool> visited(NumVertices(), false);
        uint level = UINT_MAX;
        bfs_queue.push(i);
        visited[i] = true;
        while (!bfs_queue.empty())
        {
            level++;   //因为是BFS，所有是一层一层的往外扩张的，所以每遍历一层，level+1，最后如果找到了比diameter还长的值，则说明最大半径不止diameter
            uint size = bfs_queue.size();
            for (uint j = 0u; j < size; j++)
            {
                uint front = bfs_queue.front();
                bfs_queue.pop();

                const auto& nbrs = GetNeighbors(front);
                for (const uint nbr: nbrs)
                {
                    if (!visited[nbr])
                    {
                        bfs_queue.push(nbr);
                        visited[nbr] = true;
                    }
                }
            }
        }
        if (level > diameter) diameter = level;
    }
    return diameter;
}


void Graph::LoadFromFile(const std::string &path)
{
    if (!io::file_exists(path.c_str()))
    {
        std::cout << "Failed to open: " << path << std::endl;
        exit(-1);
    }
    std::ifstream ifs(path);

    char type;
    while (ifs >> type)
    {
        if (type == 't')
        {
            char temp1;
            uint temp2;
            ifs >> temp1 >> temp2;
        }
        else if (type == 'v')
        {
            uint vertex_id, label;
            ifs >> vertex_id >> label;
            AddVertex(vertex_id, label);
        }
        else
        {
            uint from_id, to_id, label;
            ifs >> from_id >> to_id >> label;
            AddEdge(from_id, to_id, label);
        }
    }
    ifs.close();
}

void Graph::LoadUpdateStream(const std::string &path)
{
    if (!io::file_exists(path.c_str()))
    {
        std::cout << "Failed to open: " << path << std::endl;
        exit(-1);
    }
    std::ifstream ifs(path);

    std::string type;
    while (ifs >> type)
    {
        if (type == "v" || type == "-v")
        {
            uint vertex_id, label;
            ifs >> vertex_id >> label;
            updates_.emplace('v', type == "v", vertex_id, 0u, label);
        }
        else
        {
            uint from_id, to_id, label;
            ifs >> from_id >> to_id >> label;
            updates_.emplace('e', type == "e", from_id, to_id, label);
        }
    }
    ifs.close();
}

void Graph::PrintMetaData() const
{
    std::cout << "# vertices = " << NumVertices() <<
        "\n# edges = " << NumEdges() << std::endl;
}
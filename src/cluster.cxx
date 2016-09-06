#include "cluster.hh"

//std/stl
#include <string>
#include <iostream>
#include <complex>

using namespace std;
using namespace vmm;

//--------------------------------------------------------------//
// Constructor
Cluster::Cluster() :
    m_chamber_no(0),
    n_duplicate_strips(0),
    // pdo
    m_pdo(0),
    m_pdo_calibrated(0.0),
    //position
    m_position_centroid(-1),
    m_position_centroid_calibrated(-1)
{
    m_hits.clear();
}
struct clusHitStripSmaller {
    bool operator()(vmm::Hit a, vmm::Hit b) { return a.strip() < b.strip(); }
} byClusHitStrip;
//--------------------------------------------------------------//
void Cluster::addHit(Hit hit)
{
    m_hits.push_back(hit);
    m_pdo += hit.pdo();
}
//--------------------------------------------------------------//
//bool Cluster::operator==(Cluster& rhs)
//{
//    bool is_same = false;
//    if(m_clus_pdo == rhs.pdo()) is_same = true;
//    return is_same;
//}
//--------------------------------------------------------------//
void Cluster::sortClusterHitsByChannel()
{
    //sort the hits in order of increasing strip number
    sort(m_hits.begin(), m_hits.end(), byClusHitStrip); 
}
//--------------------------------------------------------------//
void Cluster::removeDuplicateStrips()
{
    n_duplicate_strips = 0;
    if(m_hits.size()>1) {
        for(int i = 0; i < (int)m_hits.size(); i++) {
            for(int j = i+1; j < (int)m_hits.size(); j++) {
                if( ((fabs(m_hits.at(i).strip() - m_hits.at(j).strip())>1) &&
                    (m_hits.at(i).pdo() == m_hits.at(j).pdo()) &&
                    (m_hits.at(i).tdo() == m_hits.at(j).tdo()) &&
                    (m_hits.at(i).bcid() == m_hits.at(j).bcid())) ||
                    (m_hits.at(i).strip() == m_hits.at(j).strip())) {

                        m_hits.erase(m_hits.begin() + j);
                        j--;
                        n_duplicate_strips++;

                    } // duplicate strip
            } // j
        } // i
    }
}
//--------------------------------------------------------------//
void Cluster::print()
{
    cout << " ---------------- Cluster::print ------------------ " << endl;
    cout << "  cluster on chamber: " << chamber() << ", size: " << size() << endl;
    for(int i = 0; i < (int)m_hits.size(); i++)
        m_hits.at(i).print();
}


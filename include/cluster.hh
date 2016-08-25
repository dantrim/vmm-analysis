#ifndef CLUSTER_HH
#define CLUSTER_HH

//vmm
#include "hit.hh"

//std/stl
#include <vector>

namespace vmm {

class Cluster
{

    public :
        Cluster();
        virtual ~Cluster(){};

        void addHit(Hit h);
        int size() { return (int)m_hits.size(); }

        std::vector<Hit>& hits() { return m_hits; }

        void sortClusterHitsByChannel();
        void removeDuplicateStrips();
        int duplicateStrips() { return n_duplicate_strips; }
        
        int pdo() { return m_clus_pdo; }

        void print();

    private :
        std::vector<Hit> m_hits;
        int m_clus_pdo;
        int n_duplicate_strips;

}; // class

} // namespace vmm

#endif
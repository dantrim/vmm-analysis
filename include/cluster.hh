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

        void setChamberNo(int number) { m_chamber_no = number; }

        void addHit(Hit h);
        int size() { return (int)m_hits.size(); }

        std::vector<Hit>& hits() { return m_hits; }

        void sortClusterHitsByChannel();
        void removeDuplicateStrips();
        int duplicateStrips() { return n_duplicate_strips; }
        
        int pdo() { return m_clus_pdo; }
        int chamber() { return m_chamber_no; }

        void print();
        std::vector<Hit> m_hits;

    private :
        int m_clus_pdo;
        int m_chamber_no;
        int n_duplicate_strips;

}; // class

} // namespace vmm

#endif

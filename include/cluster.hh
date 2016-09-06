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
        
        int chamber() { return m_chamber_no; }

        void print();

        ////////////////////////////////////////
        // pdo
        ////////////////////////////////////////

        int pdo() { return m_pdo; } // set automatically as we add hits

        void setPdoCalibrated(int pdo) { m_pdo_calibrated = pdo; }
        double pdoCalibrated() { return m_pdo_calibrated; }

        ////////////////////////////////////////
        // position
        ////////////////////////////////////////
        void setPositionCentroid(double pos) { m_position_centroid = pos; }
        double positionCentroid() { return m_position_centroid; }

        void setPositionCentroidCalibrated(double pos) { m_position_centroid_calibrated = pos; }
        double positionCentroidCalibrated() { return m_position_centroid_calibrated; }

    private :
        std::vector<Hit> m_hits;
        int m_chamber_no;
        int n_duplicate_strips;

        // pdo
        int m_pdo;
        double m_pdo_calibrated;

        // position
        double m_position_centroid;
        double m_position_centroid_calibrated;

}; // class

} // namespace vmm

#endif

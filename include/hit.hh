#ifndef HIT_H
#define HIT_H

namespace vmm {

class Hit
{
    public :
        Hit();
        virtual ~Hit(){};

        void fill(int chip_no, int bcid, int gray, int channel, int strip,
                    int threshold, int pdo, int tdo);
        void setChamberNo(int chamberno) { m_chamber_no = chamberno; }

        bool operator==(Hit&); 
        bool isCrossFire(Hit& rhs);

        void print();

        int chamber() { return m_chamber_no; }
        int chip() { return m_chip_no; }
        int bcid() { return m_bcid; }
        int gray() { return m_gray; }
        int channel() { return m_channel; }
        int strip() { return m_strip; }
        int threshold() { return m_threshold; }
        int pdo() { return m_pdo; }
        int tdo() { return m_tdo; }

        void setStrip(int strip) { m_strip = strip; }

    private :
        int m_chamber_no;
        int m_chip_no;
        int m_bcid;
        int m_gray;
        int m_channel;
        int m_strip;
        int m_threshold;
        int m_pdo;
        int m_tdo;


}; // class Hit

} // namespace vmm




#endif

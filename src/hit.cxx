#include "hit.hh"

#include <string>
#include <iostream>
#include <complex>

using namespace std;

using namespace vmm;

//--------------------------------------------------------------//
// Constructor
Hit::Hit() :
    m_chip_no(-1),
    m_bcid(-1),
    m_gray(-1),
    m_channel(-1),
    m_strip(-1),
    m_threshold(-1),
    m_pdo(-1),
    m_tdo(-1)
{
}
//--------------------------------------------------------------//
// fill the hit attributes
void Hit::fill(int chip, int bcid, int gray, int channel, int strip,
    int threshold, int pdo, int tdo)
{
    m_chip_no = chip;
    m_bcid = bcid;
    m_gray = gray;
    m_channel = channel;
    m_strip = strip;
    m_threshold = threshold;
    m_pdo = pdo;
    m_tdo = tdo;
}
//--------------------------------------------------------------//
// equivalence operator
bool Hit::operator==(Hit& rhs)
{
    bool is_same = false;
    if(this != &rhs) {
        if( m_strip==rhs.strip() &&
            m_pdo==rhs.pdo() &&
            m_tdo==rhs.tdo() &&
            m_bcid==rhs.bcid()) is_same = true;
    }
    return is_same;
}
//--------------------------------------------------------------//
// is cross fire with RHS Hit
bool Hit::isCrossFire(Hit& rhs)
{
    bool is_cross_fire = false;
    if(this != &rhs) {
        if(m_chip_no==rhs.chip() &&
            (fabs(m_channel-rhs.channel())==1) &&
            m_pdo==rhs.pdo() &&
            m_tdo==rhs.tdo() &&
            m_bcid==rhs.bcid()) is_cross_fire = true;
    }
    return is_cross_fire;
}

//--------------------------------------------------------------//
// print
void Hit::print()
{
    cout << " ---------------------- Hit::Print ---------------------- " << endl;
    cout << "  chip: " << m_chip_no << "  channel: " << m_channel << "  strip: " << m_strip << endl;
    cout << "  bcid: " << m_bcid << "  gray   : " << m_gray <<  "   threshold: " << m_threshold << endl;
    cout << "  pdo : " << m_pdo  << "  tdo    : " << m_tdo << endl;
    
}


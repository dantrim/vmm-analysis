#!/usr/bin/env python

class Cluster :
    def __init__(self) :
        self.strip_hits = []
        self.chamber_number = -1

        # event data for this cluster
        self.cluster_pdo = -1
        self.cluster_pdo_calibrated = -1

        self.position_centroid = -1
        self.position_centroid_calibrated = -1

    def setChamber(self, chamber_no) :
        self.chamber_number = chamber_no

    def chamber(self) :
        return self.chamber_number

    def addHit(self, in_hit) :
        self.strip_hits.append(in_hit)

    def hits(self) :
        return self.strip_hits

    def size(self) :
        return len(self.strip_hits)

    ##############################
    # cluster charge data
    def setPdo(self, pdo_) :
        self.cluster_pdo += pdo_

    def pdo(self) :
        return self.cluster_pdo

    def setPdoCalibrated(self, pdo_) :
        self.cluster_pdo_calibrated += pdo_

    def pdo_calibrated(self) :
        return self.cluster_pdo_calibrated

    ##############################
    # cluster position data
    def setPositionCentroid(self, centroid) :
        self.position_centroid = centroid

    def positionCentroid(self) :
        return self.position_centroid

    def setPositionCentroidCalibrated(self, centroid) :
        self.position_centroid_calibrated

    def positionCentroidCalibrated(self) :
        return self.position_centroid_calibrated

    


    ##############################
    # sort methods
    def sortHitsByChamberStrip(self) :
        """
        Sort the individual chamber hits in order
        of increasing strip number
        """
        self.strip_hits.sort(key=lambda x: x.strip_number, reverse=False)


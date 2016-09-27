#!/usr/bin/env python

class Hit :
    def __init__(self) :
        self.chamber_number = -1 # chamber to which this hit belongs
        self.chip_number = -1 # chip that is registering this hit

        # VMM data
        self.bcid = -1
        self.gray = -1
        self.channel = -1
        self.strip = -1
        self.threshold = 0
        self.pdo = 0
        self.tdo = 0
        self.filled_ok = False

    def Print(self) :
        out = "Hit::Print    "
        out += "  chip: %d  channel: %d  strip: %d | bcid: %d  gray: %d  threshold: %d"%(self.chip_number, self.channel, self.strip, self.bcid, self.gray, self.threshold)
        out += "  pdo: %d  tdo: %d"%(self.pdo, self.tdo)
        print out

    def filled(self) :
        return (self.filled_ok) and (self.chamber_number>=0)

    def setChamber(self, chamber_no_) :
        self.chamber_number = chamber_no_

    def fill(self, chip_, bcid_, gray_, channel_, strip_, threshold_, pdo_, tdo_) :
        self.chip_number = chip_
        self.bcid = bcid_
        self.gray = gray_
        self.channel = channel_
        self.strip = strip_
        self.threshold = threshold_
        self.pdo = pdo_
        self.tdo = tdo_
        self.filled_ok = True

    def __eq__(self, rhs) :
        strip_eq = (self.strip==rhs.strip)
        pdo_eq = (self.pdo==rhs.pdo)
        tdo_eq = (self.tdo==rhs.tdo)
        bcid_eq = (self.bcid==rhs.bcid)
        return (strip_eq and pdo_eq and tdo_eq and bcid_eq)

    def crossfire(self, rhs) :
        chip_eq = (self.chip_number==rhs.chip_number)
        chan_offset = (abs(self.channel-rhs.channel)==1)
        pdo_eq = (self.pdo==rhs.pdo)
        tdo_eq = (self.tdo==rhs.tdo)
        bcid_eq = (self.bcid==rhs.bcid)
        return (chip_eq and chan_offset and pdo_eq and tdo_eq and bcid_eq)
        
        

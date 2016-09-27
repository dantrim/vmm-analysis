#!/usr/bin/env python

import ROOT as r
r.PyConfig.IgnoreCommandLineOptions = True
r.gROOT.SetBatch(True)



class VMMCanvas :
    def __init__(self, name_="", w_=800, h_=600) :
        self.name = name_
        self._canvas = r.TCanvas("c_%s"%name_, "", w_, h_)

    def canvas(self) :
        return self._canvas

   # def divide(self, rows_, columns_) :
   #     self.canvas.Divide(rows_, columns_)

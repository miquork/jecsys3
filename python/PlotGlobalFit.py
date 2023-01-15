import os, ROOT
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
import ctypes
from collections import OrderedDict
from tdrstyle_JERC import *
import tdrstyle_JERC as TDR
TDR.extraText  = 'Preliminary'

class PlotGlobalFit():
    def __init__(self, run='Run2Test', year='Run2', algo='AK4 CHS'):
        self.run = run
        self.year = year
        self.algo = algo
        self.inputPath = '../rootfiles/'
        self.outputPath = './pdfs/'
        os.system('mkdir -p '+self.outputPath)
        self.filename = 'output'+run+'.root'
        TDR.cms_lumi = TDR.commonScheme['legend'][self.year]+', '+TDR.commonScheme['lumi'][self.year]+' fb^{-1}'
        TDR.cms_energy = TDR.commonScheme['energy'][self.year]
        TDR.extraText3 = []
        TDR.extraText3.append(self.algo)
        TDR.extraText3.append('|#eta| < 1.3')

    def LoadInputs(self):
        self.inputfile = ROOT.TFile(os.path.join(self.inputPath,self.filename))
        self.infos = OrderedDict()

        self.infos['herr']    = {'objName':'herr',               'leg':1, 'legName': 'herr',     'legStyle': 'f',   'plotinfo': {'opt':'e3', 'fcolor':ROOT.kCyan+1, 'lcolor':ROOT.kCyan+1} }
        self.infos['hjesref'] = {'objName':'herr_ref',           'leg':1, 'legName': 'herr_ref', 'legStyle': 'f',   'plotinfo': {'opt':'e3', 'fcolor':ROOT.kYellow+1, 'lcolor':ROOT.kYellow+1} }
        self.infos['jesFit']  = {'objName':'jesFit_func_Resp',   'leg':1, 'legName': None,       'legStyle': '',    'plotinfo': { 'lcolor' : ROOT.kBlack, 'lstyle' : ROOT.kSolid, 'lwidth' : 1,} }
        self.infos['zjet']    = {'objName':'Resp_zjet_mpf',      'leg':1, 'legName': 'Z+jet',    'legStyle': 'lep', 'plotinfo': {'opt':'p', 'mcolor':ROOT.kRed+1} }

        pf_colors = OrderedDict([('Resp',ROOT.kBlack), ('chf',ROOT.kRed+1), ('nhf',ROOT.kGreen+2), ('nef',ROOT.kAzure+2)])
        self.infos_PF = OrderedDict()
        for mode, color in pf_colors.items():
            self.infos_PF[mode+'_prefit'] = {'objName':mode+'_zjet_mpf_prefit', 'leg':1, 'legName': mode, 'legStyle': 'lep', 'plotinfo': {'opt':'Pz', 'mcolor':color, 'marker':ROOT.kFullCircle} }
            self.infos_PF[mode+'_postfit'] = {'objName':mode+'_zjet_mpf_postfit', 'leg':1, 'legName': mode, 'legStyle': 'lep', 'plotinfo': {'opt':'Pz', 'mcolor':color, 'marker':ROOT.kOpenCircle} }
            self.infos_PF[mode+'_variation'] = {'objName':mode+'_zjet_mpf_variation', 'leg':1, 'legName': mode, 'legStyle': 'lep', 'plotinfo': {'opt':'e3', 'fcolor':color, 'lcolor':color} }

        for name, info in self.infos.items():
            if 'Resp_' in info['objName']:
                for mode in ['prefit', 'postfit']:
                    info['obj_'+mode] = self.inputfile.Get(info['objName']+'_'+mode)
            else:
                info['obj'] = self.inputfile.Get(info['objName'])

        for name, info in self.infos_PF.items():
            info['obj'] = self.inputfile.Get(info['objName'])

    def FixXAsisPartition(self, shift=None):
        GettdrCanvasHist(self.canv).GetXaxis().SetNoExponent(True)
        latex = ROOT.TLatex()
        latex.SetTextFont(42)
        latex.SetTextSize(0.05)
        latex.SetTextAlign(23)
        if shift==None:
            YMin, YMax = (GettdrCanvasHist(self.canv).GetYaxis().GetXmin(), GettdrCanvasHist(self.canv).GetYaxis().GetXmax())
            shift = YMin-0.018*(YMax-YMin)
        for xbin in [30,100,300,1000, 3000]:
            latex.DrawLatex(xbin,shift,str(xbin))

    def CreateCanvasGlobalFit(self, canvName):
        XMin, XMax = (15, 4500)
        YMin, YMax = (0.95,1.07)
        yName = 'jet response ratio'
        if canvName =='raw': yName = 'raw '+yName
        if canvName =='postfit': yName = 'post-fit '+yName
        yName = yName.capitalize()
        self.canv = tdrCanvas('PlotGlobalFit'+self.year+canvName, XMin, XMax, YMin, YMax, 'p_{T} [GeV]', yName, square=kSquare, isExtraSpace=True)
        self.canv.SetLogx(True)
        self.FixXAsisPartition()
        l,r,t,b = (self.canv.GetLeftMargin(),self.canv.GetRightMargin(),self.canv.GetTopMargin(),self.canv.GetBottomMargin())
        self.legs = {}
        xpos = l+0.04*(1-l-r)
        ypos = b+0.04*(1-t-b)
        self.legs[0] = tdrLeg(xpos,ypos,xpos+0.2,ypos+0.04*4, textSize=0.04)
        xpos = 1-r-0.04*(1-l-r)-0.1
        ypos = 1-t-0.04*(1-t-b)
        # self.legs[1] = tdrLeg(xpos-0.2,ypos-0.04*2,xpos,ypos, textSize=0.04)
        self.legs[1] = tdrLeg(xpos-0.2,ypos-0.04*4,xpos,ypos, textSize=0.04)
        self.line = ROOT.TLine(XMin, 1, XMax, 1)
        tdrDrawLine(self.line, lcolor=ROOT.kBlack, lstyle=ROOT.kDashed, lwidth=1)

    def PlotResponse(self, canvName):
        self.CreateCanvasGlobalFit(canvName)
        for info in self.infos.values():
            obj = info['obj_'+canvName if 'obj_'+canvName in info else 'obj']
            if not 'jesFit_func' in info['objName']:
                tdrDraw(obj, **info['plotinfo'])
            elif canvName=='postfit':
                tdrDrawLine(obj, **info['plotinfo'])
            if info['legName']:
                self.legs[info['leg']].AddEntry(obj, info['legName'], info['legStyle'])

        self.line.Draw('same')
        fixOverlay()
        self.canv.SaveAs(os.path.join(self.outputPath, self.filename.replace('.root', '_'+canvName+'.pdf')))
        self.canv.Close()

    def PlotPFVariations(self):
        XMin, XMax = (15, 4500)
        YMin, YMax = (-2+0.001,2.5-0.001)
        canvName = 'PFVariations'+self.year
        yName = 'PF composition change (10^{-2})'
        self.canv = tdrCanvas(canvName, XMin, XMax, YMin, YMax, 'p_{T} [GeV]', yName, square=kSquare, isExtraSpace=True)
        self.canv.SetLogx(True)
        self.FixXAsisPartition()
        line = ROOT.TLine(XMin,0, XMax, 0)
        tdrDrawLine(line, lcolor=ROOT.kBlack, lstyle=ROOT.kDashed, lwidth=1)

        for info in self.infos_PF.values():
            obj = info['obj']
            tdrDraw(obj, **info['plotinfo'])
            if info['legName']:
                self.legs[info['leg']].AddEntry(obj, info['legName'], info['legStyle'])

        fixOverlay()
        self.canv.SaveAs(os.path.join(self.outputPath,canvName+'.pdf'))
        self.canv.Close()

    def PlotAll(self):
        self.LoadInputs()
        self.PlotResponse('prefit')
        self.PlotResponse('postfit')
        self.PlotPFVariations()
        self.inputfile.Close()

def main():
    PlotGlobalFit().PlotAll()

if __name__ == '__main__':
    main()

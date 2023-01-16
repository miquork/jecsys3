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

        self.functionforms = OrderedDict([
            ('track',   {'leg': 'tracks (p0)',    'color': rt.kGreen+2,  'npar':4, 'func': 'ftd'}),
            ('pho',     {'leg': 'photons (p1)',   'color': rt.kCyan+1,   'npar':1, 'func': 'fp'}),
            ('had',     {'leg': 'hadrons (p2)',   'color': rt.kRed+1,    'npar':4, 'func': 'fhx'}),
            ('hcal',    {'leg': 'had. hcal (p3)', 'color': rt.kAzure-2,  'npar':4, 'func': 'fhh'}),
            ('ecal',    {'leg': 'had. ecal (p4)', 'color': rt.kViolet+2, 'npar':4, 'func': 'feh'}),
            ('shower',  {'leg': 'shower (p5)',    'color': rt.kOrange,   'npar':6, 'func': 'fhw'}),
            ('offset',  {'leg': 'offset (p6)',    'color': rt.kTeal+2,   'npar':3, 'func': 'fl1-1'}),
            ('mctrack', {'leg': 'MC tracks (p7)', 'color': rt.kGreen-1,  'npar':4, 'func': 'ftd-ftm'}),
            ('flav',    {'leg': 'flavor (p8)',    'color': rt.kOrange+1, 'npar':3, 'func': 'f1q3-1'}),
            ])

        self.infos = OrderedDict()

        self.infos['herr']    = {'objName':'herr',               'legName': 'herr',            'legStyle': 'f',   'plotinfo': {'opt':'e3', 'fcolor':ROOT.kCyan+1, 'lcolor':ROOT.kCyan+1, 'msize':0} }
        self.infos['hjesref'] = {'objName':'herr_ref',           'legName': 'herr_ref',        'legStyle': 'f',   'plotinfo': {'opt':'e3', 'fcolor':ROOT.kYellow+1, 'lcolor':ROOT.kYellow+1, 'msize':0} }
        self.infos['jesFit']  = {'objName':'jesFit_func_Resp',   'legName': None,              'legStyle': '',    'plotinfo': { 'lcolor' : ROOT.kBlack, 'lstyle' : ROOT.kSolid, 'lwidth' : 1,} }

        self.infos['zjet']    = {'objName':'Resp_zjet_mpf',      'legName': 'Z+jet',           'legStyle': 'lep', 'plotinfo': {'opt':'p', 'msize':1.2, 'marker': ROOT.kFullStar,   'mcolor':ROOT.kRed+1} }
        self.infos['gjet']    = {'objName':'Resp_gamjet_mpf',    'legName': '#gamma+jet',      'legStyle': 'lep', 'plotinfo': {'opt':'p', 'msize':0.8, 'marker': ROOT.kFullSquare, 'mcolor':ROOT.kBlue+1} }
        self.infos['hadw']    = {'objName':'Resp_hadw_mpf',      'legName': 'W#rightarrow qq', 'legStyle': 'lep', 'plotinfo': {'opt':'p', 'msize':0.9, 'marker': ROOT.kFullCircle, 'mcolor':ROOT.kGreen+2} }

        pf_colors = OrderedDict([('Resp',ROOT.kBlack), ('chf',ROOT.kRed+1), ('nhf',ROOT.kGreen+2), ('nef',ROOT.kAzure+2)])
        self.infos_PF = OrderedDict()
        for mode, color in pf_colors.items():
            self.infos_PF[mode+'_prefit'] = {'objName':mode+'_zjet_mpf_prefit', 'legName': mode, 'legStyle': 'lep', 'plotinfo': {'opt':'Pz', 'mcolor':color, 'marker':ROOT.kFullCircle} }
            self.infos_PF[mode+'_postfit'] = {'objName':mode+'_zjet_mpf_postfit', 'legName': mode, 'legStyle': 'lep', 'plotinfo': {'opt':'Pz', 'mcolor':color, 'marker':ROOT.kOpenCircle} }
            # self.infos_PF[mode+'_variation'] = {'objName':mode+'_zjet_mpf_variation', 'legName': mode, 'legStyle': 'lep', 'plotinfo': {'opt':'e3', 'fcolor':color, 'lcolor':color} }

        for name, info in self.infos.items():
            if 'Resp_' in info['objName']:
                for mode in ['prefit', 'postfit', 'raw']:
                    info['obj_'+mode] = self.inputfile.Get(info['objName']+'_'+mode)
            else:
                info['obj'] = self.inputfile.Get(info['objName'])

        for name, info in self.infos_PF.items():
            info['obj'] = self.inputfile.Get(info['objName'])

    def FixXAsisPartition(self, shift=None):
        self.canv.SetLogx(True)
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
        self.FixXAsisPartition()
        l,r,t,b = (self.canv.GetLeftMargin(),self.canv.GetRightMargin(),self.canv.GetTopMargin(),self.canv.GetBottomMargin())
        # self.leg = {}
        xpos = l+0.04*(1-l-r)
        ypos = b+0.04*(1-t-b)
        # self.leg = tdrLeg(xpos,ypos,xpos+0.2,ypos+0.04*4, textSize=0.04)
        xpos = 1-r-0.04*(1-l-r)-0.1
        ypos = 1-t-0.04*(1-t-b)
        self.leg = tdrLeg(xpos-0.2,ypos-0.04*2,xpos,ypos, textSize=0.04)
        self.leg = tdrLeg(xpos-0.2,ypos-0.04*4,xpos,ypos, textSize=0.04)
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
                self.leg.AddEntry(obj, info['legName'], info['legStyle'])

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
        self.FixXAsisPartition()
        line = ROOT.TLine(XMin,0, XMax, 0)
        tdrDrawLine(line, lcolor=ROOT.kBlack, lstyle=ROOT.kDashed, lwidth=1)

        for info in self.infos_PF.values():
            obj = info['obj']
            tdrDraw(obj, **info['plotinfo'])
            if info['legName']:
                self.leg.AddEntry(obj, info['legName'], info['legStyle'])

        fixOverlay()
        self.canv.SaveAs(os.path.join(self.outputPath,canvName+'.pdf'))
        self.canv.Close()

    def PlotShapes(self, mode):
        XMin, XMax = (15, 4500)
        # YMin, YMax = (-2+0.001,2.5-0.001)
        YMin, YMax = (-3,3)
        canvName = 'GlobalFitShapes_'+mode+self.year
        self.canv = tdrCanvas(canvName, XMin, XMax, YMin, YMax, 'p_{T} [GeV]', 'JES change (%)', square=kSquare, isExtraSpace=True)
        self.FixXAsisPartition()
        leg = tdrLeg(0.4, 0.9- 0.035*(len(self.functionforms)+(3 if mode=='postfit' else 1))/2, 0.92, 0.9, 0.035)
        leg.SetNColumns(2)
        self.shapes = OrderedDict()
        sum='0'
        for shape, info in self.functionforms.items():
            self.shapes[shape] = self.inputfile.Get('shape_'+mode+'_'+info['func'])
            sum += '+('+self.shapes[shape].GetTitle()+')'
            if mode !='input':
                for par in range(self.shapes[shape].GetNpar()):
                    sum = sum.replace('['+str(par)+']',str(self.shapes[shape].GetParameter(par)))
                    # print(self.shapes[shape].GetNpar(), par, self.shapes[shape].GetParameter(par))
            tdrDrawLine(self.shapes[shape], lcolor=info['color'], lstyle=rt.kSolid, lwidth=2)
            leg.AddEntry(self.shapes[shape], info['leg'], 'l')
        if mode=='prefit': self.sum = sum
        self.shapes['sum'] = rt.TF1('sum', sum,15,4500)
        tdrDrawLine(self.shapes['sum'], lcolor=rt.kBlack, lstyle=rt.kSolid, lwidth=2)
        leg.AddEntry(self.shapes['sum'], 'sum', 'l')
        if mode=='postfit':
            var_sum = sum+'-('+self.sum+')'
            self.shapes['var_sum'] = rt.TF1('var_sum', var_sum,15,4500)
            tdrDrawLine(self.shapes['var_sum'], lcolor=rt.kBlack, lstyle=rt.kDashed, lwidth=2)
            leg.AddEntry(self.shapes['var_sum'], 'post-pre', 'l')
        fixOverlay()
        self.canv.SaveAs(os.path.join(self.outputPath, canvName+'.pdf'))
        self.canv.Close()

    def PlotAll(self):
        self.LoadInputs()
        self.PlotResponse('raw')
        self.PlotResponse('prefit')
        self.PlotResponse('postfit')
        # self.PlotPFVariations()
        self.PlotShapes('input')
        self.PlotShapes('prefit')
        self.PlotShapes('postfit')
        self.inputfile.Close()

def main():
    PlotGlobalFit().PlotAll()

if __name__ == '__main__':
    main()

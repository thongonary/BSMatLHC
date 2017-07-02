import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *

generator = cms.EDFilter("Pythia8GeneratorFilter",
    maxEventsToPrint = cms.untracked.int32(1),
    pythiaPylistVerbosity = cms.untracked.int32(1),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(13000.),
    PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CUEP8M1SettingsBlock,
        processParameters = cms.vstring(
            'HiddenValley:ffbar2Zv = on',
                
                'HiddenValley:spinFv = 0 ! zero flavor in hadronization',
                'HiddenValley:FSR = on ! turn on FSR',
                'HiddenValley:fragment = on ! turn on fragmentation',
                
                'HiddenValley:ffbar2Zv = on',
                '4900023:m0 = 3000',
                '4900023:mWidth = 300',
                '4900023:onMode = 0',
                '4900023:onIfAny = 4900101 ',
                
                '4900111:m0 = 20  ! spin0 meson mass (need to fix decay)',
                '4900113:m0 = 20  ! spin1 meson mass (need to fix decay)',
                '4900101:m0 = 10  ! hidden valley quark mass ',
                '4900111:onMode = on',
                '4900113:onMode = on',

                '4900111:0:bRatio = 0.0000012',
                '4900111:1:bRatio = 0.0000005',
                '4900111:2:bRatio = 0.0004920',
                '4900111:3:bRatio = 0.0793840',
                '4900111:4:bRatio = 0.8679713',
                '4900111:5:bRatio = 0.0000000',
                '####',
                '4900111:6:bRatio = 0.43000000',
                '####',
                '4900111:7:bRatio = 0.0001840',
                '4900111:8:bRatio = 0.0000000',
                '4900111:9:bRatio = 0.0519670',
                '4900111:10:bRatio = 0.0000000',
                
                '4900113:0:bRatio = 0.0500000',
                '4900113:1:bRatio = 0.2000000',
                '4900113:2:bRatio = 0.0500000',
                '4900113:3:bRatio = 0.2000000',
                '4900113:4:bRatio = 0.0500000',
                '4900113:5:bRatio = 0.1500000',
                '####',
                '4900113:6:bRatio = 0.4300000',
                '####',
                '4900113:7:bRatio = 0.1500000',
                '4900113:8:bRatio = 0.0000000',
                '4900113:9:bRatio = 0.1500000',
                '4900113:10:bRatio = 0.000000',

                'HiddenValley:alphaOrder = 1 ! turn on running',
                'HiddenValley:Ngauge  = 2 ! Nc for the gauge group SU(Nc)',
                'HiddenValley:Lambda = 10. ! confinement scale.',
                '#HiddenValley:NBFlavRun = 0 ! number of boson in the running',
                'HiddenValley:NFlav = 2 ! number of flavors in the running',
                'HiddenValley:probVector = .75 ! ratio of spin1 vs spin0 meson production',
            ),
        parameterSets = cms.vstring('pythia8CommonSettings',
                                    'pythia8CUEP8M1Settings',
                                    'processParameters',
                                    )
    )
)


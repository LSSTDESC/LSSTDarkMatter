import math

from lsstDebug import getDebugFrame
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable
from lsst.meas.astrom import AstrometryTask, displayAstrometry, createMatchMetadata,\
    LoadAstrometryNetObjectsTask
from lsst.obs.base import ExposureIdInfo
import lsst.daf.base as dafBase
from lsst.afw.math import BackgroundList
from lsst.afw.table import IdFactory, SourceTable
from lsst.meas.algorithms import SourceDetectionTask
#from lsst.meas.astrom import AstrometryTask, displayAstrometry, createMatchMetadata
from lsst.meas.base import SingleFrameMeasurementTask, ApplyApCorrTask, CatalogCalculationTask
from lsst.meas.deblender import SourceDeblendTask
#from .photoCal import PhotoCalTask
import lsst.afw.image
import lsst.afw.detection
import numpy as np

__all__ = ["ProcessSimConfig", "ProcessSimTask"]

class ProcessSimConfig(pexConfig.Config):
    detection = pexConfig.ConfigurableField(
        target=SourceDetectionTask,
        doc="Detect sources"
    )
    doDeblend = pexConfig.Field(
        dtype=bool,
        default=True,
        doc="Run deblender input exposure"
    )
    deblend = pexConfig.ConfigurableField(
        target=SourceDeblendTask,
        doc="Split blended sources into their components"
    )
    measurement = pexConfig.ConfigurableField(
        target=SingleFrameMeasurementTask,
        doc="Measure sources"
    )


class ProcessSimTask(pipeBase.Task):
    _DefaultName = "processSim"
    ConfigClass = ProcessSimConfig
    
    def __init__(self, **kwds):
        pipeBase.Task.__init__(self, **kwds)
        
        # add schema here
        self.schemaMapper = None
        self.schema = afwTable.SourceTable.makeMinimalSchema()
        
        self.makeSubtask("detection", schema=self.schema)
        if self.config.doDeblend:
            self.makeSubtask("deblend", schema=self.schema)
        self.makeSubtask("measurement", schema=self.schema)
    
    
    def run(self, exposure):
        
        table = SourceTable.make(self.schema)
        background = None
        
        detRes = self.detection.run(table=table, exposure=exposure)
        sourceCat = detRes.sources
        
        if background is None:
            background = BackgroundList()
        if detRes.fpSets.background:
            background.append(detRes.fpSets.background)
        
        if self.config.doDeblend:
            self.deblend.run(exposure=exposure, sources=sourceCat)
        
        self.measurement.run(
            measCat=sourceCat,
            exposure=exposure
        )
        
        #results2 = self.deblend.run(...)
        #results3 = self.measurement.run(...)
        
        return pipeBase.Struct(
            exposure = exposure,
            background = background,
            sourceCat = sourceCat
        )

if __name__ == "__main__":
    import sys
    import pyfits
    fits = pyfits.open(sys.argv[1])
    data = fits[0].data
    exposure = lsst.afw.image.ExposureF(data.shape[1], data.shape[0])
    exposure.getMaskedImage().getImage().getArray()[:,:] = data
    exposure.getMaskedImage().getVariance().getArray()[:,:] = np.var(data)
    exposure.setPsf(lsst.afw.detection.GaussianPsf(39, 39, 0.6))
    
    config = ProcessSimTask.ConfigClass()
    task = ProcessSimTask(config=config)
    result = task.run(exposure)
    result.sourceCat.writeFits(sys.argv[2])

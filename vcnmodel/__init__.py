#!/usr/bin/env python

# Use Semantic Versioning, http://semver.org/
version_info = (0, 2, 0, 'a')
__version__ = "%d.%d.%d%s" % version_info

#print ("apparent version: ", __version__)

import vcnmodel.model_run
import vcnmodel.model_params
import vcnmodel.cell_config
import vcnmodel.generate_run
import vcnmodel.NoiseTrainingGen
import vcnmodel.IVPlots
import vcnmodel.analyze_run
import vcnmodel.cellInitialization

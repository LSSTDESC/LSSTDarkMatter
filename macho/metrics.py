from lsst.sims.maf.metrics import BaseMetric
import numpy as np


# Define our class, inheriting from BaseMetric
class massMetric(BaseMetric):
    # Add a doc string to describe the metric.
    """
    From a mass in solar masses, use the time = 1.2*sqrt(mass/10.) years scaling
    to determine 
        1) the time scale
        2) number of visits
        3) is N_visits > 30?
        4) in year long blocks starting at first epoch, are there more than 10 visits per year?
        5) if 3) is True and 4) is True, return True, else return False
    This is a model for detectability of a microlensing event of the given timescale 
    """
    # Add our "__init__" method to instantiate the class.
    # **kwargs allows additional values to be passed to the BaseMetric that you 
    #     may not have been using here and don't want to bother with. 
    def __init__(self, mass, metricName='massmetric', **kwargs):
        self.mass = mass
        self.time_scale = 1.2*np.sqrt(mass/10.)
        # Set the values we want to keep for our class.
        self.MJDcol = 'expMJD'
        self.RAcol = 'fieldRA'
        self.DECcol = 'fieldDec'
        cols=[self.MJDcol, self.RAcol, self.DECcol]
        # Now we have to call the BaseMetric's __init__ method, to get the "framework" part set up.
        # We currently do this using 'super', which just calls BaseMetric's method.
        # The call to super just basically looks like this .. you must pass the columns you need, and the kwargs.
        super(massMetric, self).__init__(col=cols, metricName=metricName, **kwargs)
    # Now write out "run" method, the part that does the metric calculation.
    def run(self, dataSlice, slicePoint=None):
        result = detect_mass(dataSlice[self.MJDcol], dataSlice[self.RAcol], dataSlice[self.DECcol], self.time_scale)
        return result

def detect_mass(mjd_obs, ra_obs, dec_obs, time_scale) :
    # Add a doc string to describe the metric.
    """
    From a mass in solar masses, use the time = 1.2*sqrt(mass/10.) years scaling
    to determine 
        1) the time scale
        2) number of visits
        3) is N_visits > 30?
        4) in year long blocks starting at first epoch, are there more than 10 visits per year?
        5) if 3) is True and 4) is True, return True, else return False
    This is a model for detectability of a microlensing event of the given timescale. 

    There is no cumulative time counting here: if the range of MJDs is 5 years while 
    the timescale is 1 year, as long as any one year passes criteria 4, then criteria 4 is passed
    """
    verbose = False
    #verbose = True
    min_visits = 30
    min_visits_year = 10
    detectable = True

    # criteria 3
    nvisits = len(mjd_obs)
    if nvisits < min_visits:
        detectable = False
    if verbose: print "nvisits", nvisits, ">", min_visits, detectable

    # criteria 4
    if verbose: print "years", mjd_obs.min()/365.25, mjd_obs.max()/365.25
    nyears = (mjd_obs.max() - mjd_obs.min())/365.25
    counter = 0
    year_long_visits = []
    for i in np.arange(0,nyears,1):
        year_long_visits.append( [i,i+1])

    for i in np.arange(0,nyears,1):
        y_visit = year_long_visits[int(i)]
        ix = (mjd_obs >= mjd_obs.min() + y_visit[0]*365.25) &  \
            (mjd_obs <= mjd_obs.min() + y_visit[1]*365.25)
        if verbose: print "year min",mjd_obs.min()/365.25+y_visit[0], "year max",mjd_obs.min()/365.25+y_visit[1],
        if verbose: print "nvisits in a year",len(mjd_obs[ix]),">", min_visits_year
        if len(mjd_obs[ix]) < min_visits_year:
            detectable = False
        else :
            counter += 1

    if verbose: print "\t", detectable
    if counter < time_scale:
        detectable = False

    if verbose: print "\t", detectable
    return detectable

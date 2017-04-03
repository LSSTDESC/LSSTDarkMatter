from __future__ import print_function
import httplib2
import urllib
try:
	import BeautifulSoup
	Parser = BeautifulSoup.RobustHTMLParser
except:
	from bs4 import BeautifulSoup as Parser
	
import collections
import copy
import read_girardi
import tempfile

# Code to download girardi isochrones

def sender(http, postget, url, headers, body=None):
	if postget=='POST':
		headers['Content-type']='application/x-www-form-urlencoded; charset=UTF-8'
	else:
		if 'Content-type' in headers:
			del headers['Content-type']
	
	try:
		urlencode = urllib.urlencode
	except:
		urlencode = urllib.parse.urlencode
	if body is None:
		response, content = http.request(url, postget, headers=headers)
	else:
		response, content = http.request(url, postget, headers=headers,
		            body=urlencode(body))
	if response['status'] not in ('200', '302'):
		print ( url,'\n',response,'\n',headers,'\n',content)
		raise Exception()
	return response,content

url = 'http://stev.oapd.inaf.it/cgi-bin/cmd_2.7'

__params  = [
('submit_form','Submit'),
('cmd_version',2.6),
('isoc_kind','parsec_CAF09_v1.2S'),
("eta_reimers",0.2),
("kind_interp",1),
("kind_tpagb",0),
("kind_pulsecycle",0),
("kind_postagb",-1),
("photsys_version","yang"),
("photsys_file","tab_mag_odfnew/tab_mag_decam.dat"),
("kind_cspecmag","aringer09"),
("dust_sourceM","nodustM"),
("dust_sourceC","nodustC"),
("kind_mag",2),
("kind_dust",0),
("extinction_av",0.0),
("imf_file","tab_imf/imf_chabrier_lognormal.dat"),
("isoc_val",0),
("isoc_age",1.0e9),
("isoc_zeta",0.019),
("isoc_zeta0",0.008),
("isoc_lage0",6.6),
("isoc_lage1",10.13),
("isoc_dlage",0.05),
("isoc_age0",12.7e9),
("isoc_z0",0.0001),
("isoc_z1",0.03),
("isoc_dz",0.0001),
("output_kind",0),
("output_evstage",1),
("lf_maginf",20),
("lf_magsup",-20),
("lf_deltamag",0.2),
(".cgifields",'output_evstage,kind_cspecmag.output_gzip,output_kind,dust_sourceC,isoc_val,dust_sourceM,photsys_version,isoc_kind')]

paramsdict = collections.OrderedDict(__params)

__filters = [ ("tab_mag_odfnew/tab_mag_2mass_spitzer_wise.dat","2MASS + Spitzer (IRAC+MIPS) + WISE"),
("tab_mag_odfnew/tab_mag_2mass.dat","2MASS JHKs"),
("tab_mag_odfnew/tab_mag_ogle_2mass_spitzer.dat","OGLE + 2MASS + Spitzer (IRAC+MIPS)"),
("tab_mag_odfnew/tab_mag_2mass_spitzer_wise_washington_ddo51.dat","2MASS+Spitzer+WISE+Washington+DDO51"),
("tab_mag_odfnew/tab_mag_ubvrijhk.dat","UBVRIJHK"), # (cf. Maiz-Apellaniz 2006 + Bessell 1990)"),
("tab_mag_odfnew/tab_mag_bessell.dat","UBVRIJHKLMN"),# (cf. Bessell 1990 + Bessell & Brett 1988)"),
("tab_mag_odfnew/tab_mag_akari.dat","AKARI"),
("tab_mag_odfnew/tab_mag_batc.dat","BATC"),
("tab_mag_odfnew/tab_mag_megacam_wircam.dat","CFHT Megacam + Wircam (all ABmags)"),
("tab_mag_odfnew/tab_mag_wircam.dat","CFHT Wircam"),
("tab_mag_odfnew/tab_mag_megacam.dat","CFHT/Megacam u*g'r'i'z'"),
("tab_mag_odfnew/tab_mag_ciber.dat","CIBER"),
("tab_mag_odfnew/tab_mag_dcmc.dat","DCMC"),
("tab_mag_odfnew/tab_mag_decam.dat","DECAM (ABmags)"),
("tab_mag_odfnew/tab_mag_decam_vista.dat","DECAM ugrizY (ABmags) + VISTA ZYJHK (Vegamags)"),
("tab_mag_odfnew/tab_mag_denis.dat","DENIS"),
("tab_mag_odfnew/tab_mag_dmc14.dat","DMC 14 filters"),
("tab_mag_odfnew/tab_mag_dmc15.dat","DMC 15 filters"),
("tab_mag_odfnew/tab_mag_eis.dat","ESO/EIS (WFI UBVRIZ + SOFI JHK)"),
("tab_mag_odfnew/tab_mag_wfi.dat","ESO/WFI"),
("tab_mag_odfnew/tab_mag_wfi_sofi.dat","ESO/WFI+SOFI"),
("tab_mag_odfnew/tab_mag_wfi2.dat","ESO/WFI2"),
("tab_mag_odfnew/tab_mag_galex.dat","GALEX FUV+NUV (Vegamag) + Johnsons UBV"),
("tab_mag_odfnew/tab_mag_galex_sloan.dat","GALEX FUV+NUV + SDSS ugriz (all ABmags) "),
("tab_mag_odfnew/tab_mag_gaia.dat","Gaias G, G_BP and G_RP (Vegamags)"),
("tab_mag_odfnew/tab_mag_UVbright.dat","HST+GALEX+Swift/UVOT UV filters"),
("tab_mag_odfnew/tab_mag_acs_hrc.dat","HST/ACS HRC"),
("tab_mag_odfnew/tab_mag_acs_wfc.dat","HST/ACS WFC"),
("tab_mag_odfnew/tab_mag_nicmosab.dat","HST/NICMOS AB"),
("tab_mag_odfnew/tab_mag_nicmosst.dat","HST/NICMOS ST"),
("tab_mag_odfnew/tab_mag_nicmosvega.dat","HST/NICMOS vega"),
("tab_mag_odfnew/tab_mag_stis.dat","HST/STIS imaging mode"),
("tab_mag_odfnew/tab_mag_wfc3ir.dat","HST/WFC3 IR channel (final throughputs)"),
("tab_mag_odfnew/tab_mag_wfc3uvis1.dat","HST/WFC3 UVIS channel, chip 1 (final throughputs)"),
("tab_mag_odfnew/tab_mag_wfc3uvis2.dat","HST/WFC3 UVIS channel, chip 2 (final throughputs)"),
("tab_mag_odfnew/tab_mag_wfc3_wideverywide.dat","HST/WFC3 all W+LP+X filters (UVIS1+IR, final throughputs)"),
("tab_mag_odfnew/tab_mag_wfc3_verywide.dat","HST/WFC3 long-pass and extremely wide filters (UVIS1, final throughputs)"),
("tab_mag_odfnew/tab_mag_wfc3_medium.dat","HST/WFC3 medium filters (UVIS1+IR, final throughputs)"),
("tab_mag_odfnew/tab_mag_wfc3_wide.dat","HST/WFC3 wide filters (UVIS1+IR, final throughputs)"),
("tab_mag_odfnew/tab_mag_wfpc2.dat","HST/WFPC2 (Vegamag, cf. Holtzman et al. 1995)"),
("tab_mag_odfnew/tab_mag_int_wfc.dat","INT/WFC (Vegamag)"),
("tab_mag_odfnew/tab_mag_iphas.dat","IPHAS"),
("tab_mag_odfnew/tab_mag_kepler.dat","Kepler + SDSS griz + DDO51 (in ABmags)"),
("tab_mag_odfnew/tab_mag_kepler_2mass.dat","Kepler + SDSS griz + DDO51 (in ABmags) + 2MASS (~Vegamag)"),
("tab_mag_odfnew/tab_mag_lbt_lbc.dat","LBT/LBC (Vegamag)"),
("tab_mag_odfnew/tab_mag_lsst.dat","LSST ugrizY, March 2012 total filter throughputs (all ABmags)"),
("tab_mag_odfnew/tab_mag_noao_ctio_mosaic2.dat","NOAO/CTIO/MOSAIC2 (Vegamag)"),
("tab_mag_odfnew/tab_mag_ogle.dat","OGLE-II"),
("tab_mag_odfnew/tab_mag_panstarrs1.dat","Pan-STARRS1"),
("tab_mag_odfnew/tab_mag_sloan.dat","SDSS ugriz"),
("tab_mag_odfnew/tab_mag_sloan_2mass.dat","SDSS ugriz + 2MASS JHKs"),
("tab_mag_odfnew/tab_mag_sloan_ukidss.dat","SDSS ugriz + UKIDSS ZYJHK"),
("tab_mag_odfnew/tab_mag_swift_uvot.dat","SWIFT/UVOT UVW2, UVM2, UVW1,u (Vegamag) "),
("tab_mag_odfnew/tab_mag_spitzer.dat","Spitzer IRAC+MIPS"),
("tab_mag_odfnew/tab_mag_stroemgren.dat","Stroemgren-Crawford"),
("tab_mag_odfnew/tab_mag_suprimecam.dat","Subaru/Suprime-Cam (ABmags)"),
("tab_mag_odfnew/tab_mag_TESS_2mass_kepler.dat","TESS + 2MASS (Vegamags) + Kepler + SDSS griz + DDO51 (in ABmags)"),
("tab_mag_odfnew/tab_mag_tycho2.dat","Tycho VTBT"),
("tab_mag_odfnew/tab_mag_ukidss.dat","UKIDSS ZYJHK (Vegamag)"),
("tab_mag_odfnew/tab_mag_visir.dat","VISIR"),
("tab_mag_odfnew/tab_mag_vista.dat","VISTA ZYJHKs (Vegamag)"),
("tab_mag_odfnew/tab_mag_vphas.dat","VPHAS+ (ABmags)"),
("tab_mag_odfnew/tab_mag_vst_omegacam.dat","VST/OMEGACAM (Vegamag)"),
("tab_mag_odfnew/tab_mag_vilnius.dat","Vilnius"),
("tab_mag_odfnew/tab_mag_washington.dat","Washington CMT1T2"),
("tab_mag_odfnew/tab_mag_washington_ddo51.dat","Washington CMT1T2 + DDO51"),
("tab_mag_odfnew/tab_mag_jwst_wide.dat","planned JWST wide filters")]

filtersdict = dict([_[::-1] for _ in __filters])

headers = collections.OrderedDict()
headers['User-Agent']='Mozilla/5.0 (X11; Linux x86_64; rv:31.0) Gecko/20100101 Firefox/31.0'



def getit(logage, feh, filters=None):
	""" get the isochrone of a given log age and feh
	if the logage is a tuple then it is interpreted as a grid parameters
	(minlogage,maxlogage,steplogage)
	filters is the keyword specifying the photometric system
	"""
	http = httplib2.Http()
	curparams = copy.copy(paramsdict)
	#print filtersdict,f
	curparams['photsys_file'] = filtersdict[filters]	
	agegrid = False
	
	try : 
		nages = len(logage)
		if nages != 3:
			raise
		lage1,lage2,lage_step = logage
		assert(lage2>lage1)
		agegrid=True
	except:
		pass

	zzero = 0.0152
	z = zzero * 10**feh

		
	if agegrid:
		curparams['isoc_lage0']=lage1
		curparams['isoc_lage1']=lage2
		curparams['isoc_dlage']=lage_step
		curparams['isoc_zeta0']=z
		curparams['isoc_val']=1		
	else:
		curparams['isoc_age'] = 10**logage
		curparams['isoc_zeta'] =  z

	response,content=sender(http, 'POST' , url, headers, body=curparams)
	doc = Parser(content)
	ref1=doc.findAll('a')[0]['href']
	ref0=doc.findAll('base')[0]['href']
	ref0=ref0[:ref0.rfind('/')]
	isourl = ref0+'/'+ref1
	try:
		urlopen = urllib.urlopen
	except:
		urlopen = urllib.request.urlopen
	with tempfile.NamedTemporaryFile() as fd:
		fname = fd.name
		fcontents=urlopen(isourl).read()
		fd.file.write(fcontents)
		fd.file.flush()
		#print content
		#import os
		#os.system('cp %s /tmp/__iso.dat'%fname)
		isodat= read_girardi.read_girardi(fname)
	return isodat
		
		
def filters():
	return sorted(filtersdict.keys())
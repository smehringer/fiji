import re
import os.path

from ij import ImagePlus, IJ, Prefs, WindowManager
from ij.io import OpenDialog, FileSaver, Opener, ImportDialog, DirectoryChooser
from ij.plugin import ZProjector
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable
from ij.plugin.filter import ParticleAnalyzer
from ij.plugin.filter import Binary, MaximumFinder
from java.awt import Rectangle
from ij.process import ImageProcessor
from ij.gui import GenericDialog,WaitForUserDialog

# ==============================================================================
# Classes
# ==============================================================================

# ------------------------------------------------------------------------------
# Class ImageHolder
# ------------------------------------------------------------------------------
class ImageHolder(object):
	""" Hold images that are important for further analysis."""
	def __init__(self):
		self.niba = None
		self.bf   = None
		self.cfp  = None
		self.wu   = None
		self.mrna = None

	def checkStatus(self):
		for img in [self.niba, self.bf, self.wu, self.cfp, self.mrna]:
			if img == None:
				return False
		return True

# ------------------------------------------------------------------------------
# Class My_table
# ------------------------------------------------------------------------------
class My_table(object):
	""" Provide (result-)table that is easy to fill, access and alter"""
	def __init__(self,column_names):
		self.rows  = []
		self.keys  = column_names
		self.n     = len(column_names)
		self.count = 0

	def addRow(self,r): # input is a list of attributes
		if len(r) != len(self.keys):
    		# exception method
			return(False)
		row = {}
		for k in range(self.n):
			row[self.keys[k]] = r[k]
		self.rows.append(row)
		self.count = self.count + 1

	def updateRow(self,i,r):
		if len(r) != len(self.keys):
    		# exception method
			return(False)
		for k in range(self.n):
			self.rows[i][self.keys[k]] = r[k]

	def getEntry(self,i,name):
		return(self.rows[i][name])

	def setEntry(self,i,name,entry):
		self.rows[i][name] = entry
		return(True)

	def delRow(self,i):
		del self.rows[i]
		self.count = self.count - 1

	def getColumn(self,name):
		c = []
		for r in range(self.count):
			c = c + [self.rows[r][name]]
		return(c)

	def getIndexByEntry(self,name,entry):
		c = []
		for r in range(self.count):
			if self.rows[r][name]==entry:
				c = c + [r]
		return(c)

	def getRow(self, i):
		return self.rows[i]

	def __repr__(self):
		string = "row_n\t"
		for k in self.keys:
			string = string + k + "\t"
		string = string + "\n0\t"
		for r in range(self.count):
			for k in self.keys:
				string = string + str(self.rows[r][k]) + "\t"
			string = string + "\n" + str(r+1) + "\t"
		return string

# ==============================================================================
# Functions
# ==============================================================================

# ------------------------------------------------------------------------------
# maxZProjection()
# ------------------------------------------------------------------------------
def maxZprojection(stackimp,start,stop):
	"""
	Return a maximum projection of param:'stackimp' (stacked image) 
	from slice param:'start' to param:'stop'.
	[Implements Fiji>Image>Stacks>Z Project]
	"""
	zp = ZProjector(stackimp)
	zp.setMethod(ZProjector.MAX_METHOD)
	zp.setStartSlice(start)
	zp.setStopSlice(stop)
	zp.doProjection()
	zpimp = zp.getProjection()
	return zpimp

# ------------------------------------------------------------------------------
# calculateThresholdValue()
# ------------------------------------------------------------------------------
def calculateThresholdValue(imp,percentage):
	"""
	Return grey value (0-255) (threshold value) where param:'percentage' % 
	of the param:'img' is chosen (black area).
	[Mirrors Histogramm in Fiji>Image>Adjust>Threshold]
	"""
	# calculates a histogram of pixel values
	ma = imp.getProcessor().getMax()
	mi = imp.getProcessor().getMin()
	hist = [0 for i in range(256)]
	steps = round((ma-mi)/256 ,3 )
	for x in range(imp.width):
		for y in range(imp.height):
			index = int( (imp.getPixel(x,y)[0]-mi)/steps )
			if index >= 256: index = 255
			hist[index] = hist[index] + 1
	i = 1
	val = 0 # percent
	while val < percentage: 
		val = float(sum(hist[-i:])*100.0)/float(imp.width*imp.height)
		i = i + 1
	return(int((256-i+2)*steps) + mi)

# ------------------------------------------------------------------------------
# getInitialROIs()
# ------------------------------------------------------------------------------
def getInitialROIs(niba, rm, pa):
	"""
	Process image param:'niba' for automatic segmentation using the Particle Analyser 
	(param:'pa'). Resulting ROIs are stored in the Roi-Manager (param:'rm').
	Return processed image.
	[User interaction needed]
	"""
	# image 'niba' must be open when given to this function
	if niba.getNSlices() != 1:
		# todo:: identify "good" slices?
		Zniba = maxZprojection(niba, 7, niba.getNSlices()-5)
		print "hello"
	else: Zniba = niba.duplicate()

	IJ.run(Zniba, "Enhance Contrast","saturated Pixel=0.5")
	IJ.run(Zniba, "Subtract Background...","radius = 100")
	IJ.run(Zniba, "Bandpass Filter...","filter_large=150 ; filter_small=10")
	IJ.run(Zniba, "Smooth","")
	IJ.run(Zniba, "Gaussian Blur...","radius=10")

	Zniba.show()
	# manual thresholding through user
	IJ.run(Zniba, "Threshold...","")
	WaitForUserDialog("Segmentation Threshold", 
					  "Please alter threshold before you proceed with OK").show()
	
	Zniba.show()
	IJ.run(Zniba, "Make Binary","")
	IJ.run(Zniba, "Dilate","") # slightly enlarge cells
	IJ.run(Zniba, "Watershed","")

	if pa.analyze(Zniba):
		print "All ok"
	else:
		print "There was a problem in analyzing"

	rm.show()
	return(Zniba)

# ------------------------------------------------------------------------------
# countNuclei()
# ------------------------------------------------------------------------------
def countNuclei(rm, wu, pa, roi_table, niba):
	"""
	Return table of analyzed nuclei. Update param:'roi_table' to assign each 
	nuclei to its cell.
	param:'rm' Roi Manager
	param:'wu' image that is is processed for nuclei segmentation.
	param:'pa' Particle Analyser used for segmentation.
	param:'roi_table' table with info for each cell(roi) 
	param:'niba' image to measure whi5 itensities for each nuclei
	"""
	if wu.getNSlices() != 1: 
		Zwu = maxZprojection(wu,10,15)
	else: 
		Zwu = wu.duplicate()
		Zwu.killRoi()# if any roi appears on image ignore it

	# analyze wu image
	IJ.run(Zwu, "Enhance Contrast","saturated Pixel=0.5")
	IJ.run(Zwu, "Smooth","")
	IJ.run(Zwu, "Gaussian Blur...","radius=12")
	Zwu.show()
	lower_threshold	= calculateThresholdValue(Zwu,3.3)
	Zwu.getProcessor().setThreshold(lower_threshold, Zwu.getProcessor().getMax(), 
									ImageProcessor.NO_LUT_UPDATE)
	IJ.run(Zwu, "Threshold","dark background")
	old_roi_count = rm.getCount()
	rois = rm.getRoisAsArray()
	pa.analyze(Zwu)

	# prepare for measurements
	nuclei_table = My_table(["area","whi5","X","Y","width","height"])
	Zniba = maxZprojection(niba, 7, niba.getNSlices()-5);
	rt = ResultsTable.getResultsTable()
	rt.reset()
	IJ.run("Set Measurements...",
		   "area ; centroid ;  integrated density ; bounding rectangle")
	rm.deselect()

	# store measurements in nuclei_table
	for i in range(old_roi_count,rm.getCount()):
		rm.select(Zniba,i)
		IJ.run(Zniba,"Measure","")
		i = i - old_roi_count # result table starts with 0
		nuclei_table.addRow([rt.getValue('Area',i),rt.getValue('IntDen',i),
						     rt.getValue('X',i),rt.getValue('Y',i),
						     rt.getValue('Width',i),rt.getValue('Height',i)])

	# assign each roi to a cell 
	# by storing its id (row number in nuclei_table) in roi_table>n_id
	for i in range(nuclei_table.count):
		x = int(nuclei_table.getEntry(i,'X'))
		y = int(nuclei_table.getEntry(i,'Y'))
		for roi in rois:
			if roi.contains(x,y): 
				roi_table.setEntry(rm.getRoiIndex(roi),"nuclei", 
					roi_table.getEntry(rm.getRoiIndex(roi),"nuclei") + 1)
				roi_table.setEntry(rm.getRoiIndex(roi),"n_id", 
					roi_table.getEntry(rm.getRoiIndex(roi),"n_id") + [i])

	rt.reset()		
	rm.setSelectedIndexes(range(old_roi_count,rm.getCount()))
	rm.runCommand("delete") # delete nuclei rois
	return(nuclei_table)

# ------------------------------------------------------------------------------
# countAndMeasureSPB()
# ------------------------------------------------------------------------------
def countAndMeasureSPB(rm, cfp, pa, roi_table, mean_grey_value):
	"""
	Return table of analyzed spindle pole bodies (spb). Update param:'roi_table' to 
	assign each spb to its cell.
	param:'rm' Roi Manager
	param:'cfp' image that is is processed for spb segmentation.
	param:'pa' Particle Analyser used for segmentation.
	param:'roi_table' table with info for each cell(roi) 
	param:'mean:grey_value' mean grey value of the background in 'niba' image
	"""
	old_roi_count = rm.getCount()
	rois = rm.getRoisAsArray()

	# analyze cfp image
	Zcfp = maxZprojection(cfp,1,cfp.getNSlices())
	lower_threshold	= calculateThresholdValue(Zcfp,0.19)
	Zcfp.getProcessor().setThreshold(lower_threshold, Zcfp.getProcessor().getMax(), 
									 ImageProcessor.NO_LUT_UPDATE)
	IJ.run(Zcfp, "Threshold","dark background")
	IJ.run(Zcfp, "Despeckle","")
	pa.analyze(Zcfp)

	# prepare for measurements
	spb_table = My_table(["area","spb_intensity","high_intensity","X","Y"])
	rt = ResultsTable.getResultsTable()
	rt.reset()
	IJ.run("Set Measurements...",
		   "area ; min & max grey value; mean grey value; centroid ;  integrated density")
	rm.deselect()

	# store measurements in spb_table
	for i in range(old_roi_count,rm.getCount()):
		rm.select(Zcfp,i)
		IJ.run(Zcfp,"Measure","")
		i = i - old_roi_count # result table starts with 0
		spb_table.addRow([rt.getValue('Area',i),rt.getValue('IntDen',i),"no",
						  rt.getValue('X',i),rt.getValue('Y',i)])
	# evaluate each spb, wether it might be noise
	to_delete = []
	for i in range(spb_table.count):
		rm.deselect()
		rt.reset()
		if spb_table.getEntry(i,"area") < 9:
			rm.select(Zcfp, old_roi_count + i)
			IJ.run(Zcfp, "Enlarge...", "enlarge=3")
			rm.runCommand("Update")
			touch = False
			# if spb is now touching another one, it is probably noise
			rs = range(old_roi_count,rm.getCount())
			rs.remove(i+old_roi_count)
			for r in rs:
				if touchingRoi(rm.getRoi(i+old_roi_count),rm.getRoi(r)):
					touch = True
			if touch: to_delete = [i] + to_delete
			else:
				IJ.run(Zcfp,"Measure","")
				spb_table.setEntry(i,"spb_intensity",rt.getValue('IntDen',0))
	for i in to_delete:
		spb_table.delRow(i)

	# assign each roi to a cell 
	# by storing its id (row number in spb_table) in roi_table>spb_id
	for i in range(spb_table.count):
		x = int(spb_table.getEntry(i,'X'))
		y = int(spb_table.getEntry(i,'Y'))
		for roi in rois:
			# because cast to int causes round off, check for all possible points
			if (roi.contains(x,y) or roi.contains(x+1,y) or 
			   roi.contains(x,y+1) or roi.contains(x+1,y+1)):  
				roi_table.setEntry(rm.getRoiIndex(roi),"spb", 
					roi_table.getEntry(rm.getRoiIndex(roi),"spb") + 1)
				roi_table.setEntry(rm.getRoiIndex(roi),"spb_id",
					roi_table.getEntry(rm.getRoiIndex(roi),"spb_id") + [i])

	rt.reset()		
	rm.setSelectedIndexes(range(old_roi_count,rm.getCount()))
	rm.runCommand("delete") # delete spindle pole body rois

	# ensure that every cell has no more than two spb's
	for i in range(roi_table.count):
		if roi_table.getEntry(i,"spb") > 2:
			spbs = sorted([(spb_table.getEntry(j,'area'),j) for j in roi_table.getEntry(i,"spb_id")])
			spbs = spbs[-2:] # only take the biggest two spindle pole bodys for the cell
			roi_table.setEntry(i,"spb",2)
			roi_table.setEntry(i,"spb_id",[spbs[j][1] for j in [0,1]])

	# identify spb's with a twice as high intensity
	spbs = sum(roi_table.getColumn("spb_id"),[])
	intensities = sorted( [ spb_table.getEntry(j,"spb_intensity") for j in spbs] )
	if len(intensities) %2 == 1:
		av_intensity = intensities[((len(intensities)+1)/2)-1]
	else:
		av_intensity = float(sum(intensities[(len(intensities)/2)-1:(len(intensities)/2)+1]))/2.0
	for i in spbs:
		if spb_table.getEntry(i,"spb_intensity") >= 2.1*av_intensity:
			spb_table.setEntry(i,"high_intensity","yes")

	return(spb_table)

# ------------------------------------------------------------------------------
# getRoiMeasurements()
# ------------------------------------------------------------------------------
def getRoiMeasurements(rm,roi_table,niba):
	"""
	Return the mean grey value of the background of param:'niba'.
	Store roi measurements in param:'roi_table'.
	"""
	# prepare for measurements
	Zniba = maxZprojection(niba, 7, niba.getNSlices()-5);Zniba.show()
	rt = ResultsTable.getResultsTable()
	rt.reset()
	IJ.run("Set Measurements...",
		   "area ; mean grey value ; centroid ;  integrated density ; bounding rectangle")
	
	# get mean grey value
	rm.deselect()
	rm.setSelectedIndexes(range(rm.getCount()))
	IJ.run(Zniba,"Make Inverse","")
	IJ.run(Zniba,"Measure","")
	mean_grey_value = rt.getValue('Mean',0)
	rt.reset()
	rm.deselect()

	# get measurements and store them in roi_table
	for i in range(rm.getCount()):
		rm.select(Zniba,i)
		IJ.run(Zniba,"Measure","")
		roi_table.addRow([rm.getName(i)[:4],rt.getValue('Area',i),0,[],0,[],
						 int(rt.getValue('X',i)),int(rt.getValue('Y',i)),
						 "no",rt.getValue('IntDen',i),rt.getValue('Width',i),
						 rt.getValue('Height',i)])

	rt.reset()
	rm.deselect() 
	Zniba.hide()
	return(mean_grey_value)

# ------------------------------------------------------------------------------
# touchingRoi()
# ------------------------------------------------------------------------------
def touchingRoi(roi,roi2):
	""" Return True if param:'roi' and param:'roi2' are touching. """
	r = roi.getBounds()
	for x in range(r.width):
		for y in range(r.height):
			if roi.contains(r.x+x , r.y+y):
				if roi2.contains(r.x+x , r.y+y):
					return(True)	
	return(False)

# ------------------------------------------------------------------------------
# shrinkRoi()
# ------------------------------------------------------------------------------
def shrinkRoi(roi,touching_rois):
	""" 
	Shrink param:'roi' until it only touches the closest roi(s) in 
	param:'touching_rois' with as few pixels as possible.
	"""
	n = len(touching_rois)
	index = rm.getRoiIndex(roi)
	which = [1 for i in range(n)]
	s = sum(which)
	change = []
	while (s >= 1):
		for i in change: which[i] = 0
		IJ.run("Enlarge...","enlarge=-1") # shrink roi
		rm.runCommand("Update")
		roi = rm.getRoi(index) # get shrinked roi
		for r in range(n):
			if which[r] == 1:
				if touchingRoi(roi,rm.getRoi(touching_rois[r])) == False :
					s = s - 1
					change = change + [r]

	return([touching_rois[i] for i in [x for x,y in enumerate(which) if y == 1]])

# ------------------------------------------------------------------------------
# combineTwoRois()
# ------------------------------------------------------------------------------
def combineTwoRois(index,index2,roi_table,rm):
	""" 
	Cobine Measurements of rois at param:'index' and param:'index2'. Delete
	roi that has the bigger index.
	"""
	rm.deselect()
	rm.setSelectedIndexes([index,index2])
	# "Or" and "Update" will combine both rois 
	# and save the resulting new roi at smaller index
	rm.runCommand("OR")
	rm.runCommand("Update") 
	rm.deselect()

	if index < index2: 
		# roi at index2 will be deleted -> store information in roi at index
		rm.select(index2)
		rm.runCommand("delete")
		roi_table.setEntry(index,"nuclei",roi_table.getEntry(index2,"nuclei")+roi_table.getEntry(index,"nuclei")) 
		roi_table.setEntry(index,"n_id",roi_table.getEntry(index2,"n_id")+roi_table.getEntry(index,"n_id"))
		roi_table.setEntry(index,"spb",roi_table.getEntry(index2,"spb")+roi_table.getEntry(index,"spb"))
		roi_table.setEntry(index,"spb_id",roi_table.getEntry(index2,"spb_id")+roi_table.getEntry(index,"spb_id"))
		roi_table.setEntry(index,"area",roi_table.getEntry(index2,"area")+roi_table.getEntry(index,"area"))
		roi_table.setEntry(index,"whi5",roi_table.getEntry(index2,"whi5")+roi_table.getEntry(index,"whi5"))
		roi_table.setEntry(index,"name",roi_table.getEntry(index2,"name"))
		roi_table.delRow(index2)
	else: 
		# roi at index will be deleted -> store information in roi at index2
		rm.select(index)
		rm.runCommand("delete")
		roi_table.setEntry(index2,"nuclei",roi_table.getEntry(index2,"nuclei")+roi_table.getEntry(index,"nuclei"))
		roi_table.setEntry(index2,"n_id",roi_table.getEntry(index2,"n_id")+roi_table.getEntry(index,"n_id"))
		roi_table.setEntry(index2,"spb",roi_table.getEntry(index,"spb")+roi_table.getEntry(index2,"spb"))
		roi_table.setEntry(index2,"spb_id",roi_table.getEntry(index,"spb_id")+roi_table.getEntry(index2,"spb_id"))
		roi_table.setEntry(index2,"area",roi_table.getEntry(index2,"area")+roi_table.getEntry(index,"area"))
		roi_table.setEntry(index2,"whi5",roi_table.getEntry(index2,"whi5")+roi_table.getEntry(index,"whi5"))
		roi_table.delRow(index)
	print "Combined",index,"and",index2
		
# ------------------------------------------------------------------------------
# evaluate_noiseOrBud()
# ------------------------------------------------------------------------------
def evaluate_noiseOrBud(index, rm, roi_table):
	"""
	Alter roi_table entries according to the evaluation of param:'index'
	Evaluation:  
	* Enlarge roi at param:'index'
	* If roi does not touch any other roi its considered NOISE is deleted
	* If roi touches any roi with a nuclei (identified cell) its considered a BUD
	  ** Count touching rois
	  ** If there are more than one potential mother cell: evaluate_Bud()
	  ** If its only one mother cell:
	     *** If mother cell has 2 spb: evaluate_G2_vs_PM()
	     *** If mother cell has 1 spb: ANAPHASE
	     *** If mother cell has 0 spb: LATE S
	  ** combine mother cell and bud
	"""
	rm.deselect() 
	rm.select(index)
	
	# enlarge particle to half the size of an average cell
	av = sum(roi_table.getColumn('area'))/rm.getCount() # average area
	roi_area = roi_table.getEntry(index,"area")
	pix = int((av/2 - roi_area)/100) # enlarge by 1 pixel enlarges circular area by ~100
	if pix <= 1: pix = 2
	enlarge_by_pix	= "enlarge=" + str(pix)
	rm.select(index)
	IJ.run("Enlarge...",enlarge_by_pix)
	rm.runCommand("Update")
	roi = rm.getRoi(index) # get updated enlarged roi

	# check if other rois are even near the investigated roi
	rois_to_evaluate = roi_table.getIndexByEntry("eval","no")
	rois_in_range = []
	for check_roi in rois_to_evaluate:
		x1 =  roi_table.getEntry(index,'X')
		x2 =  roi_table.getEntry(check_roi,'X')
		y1 =  roi_table.getEntry(index,'Y')
		y2 =  roi_table.getEntry(check_roi,'Y')
		w1 =  roi_table.getEntry(index,'width')
		w2 =  roi_table.getEntry(check_roi,'width')
		h1 =  roi_table.getEntry(index,'height')
		h2 =  roi_table.getEntry(check_roi,'height')
		if ( (abs(x1-x2) < w1 + w2) and (abs(y1-y2) < h1 + h2) ):
			rois_in_range = rois_in_range + [check_roi]
	
	# check if remaining other rois touch investigated roi
	touching_rois = []
	for check_roi in rois_in_range:
		if touchingRoi(roi,rm.getRoi(check_roi)):
			if roi_table.getEntry(check_roi,"nuclei") > 0:
				touching_rois = touching_rois + [check_roi]
				print "Found a touching roi ",check_roi,"for roi " + str(index)

	# evlaute investigated roi
	if len(touching_rois) == 0: # roi is noise
		rm.runCommand("Delete")
		roi_table.delRow(index)
		return(False,index) 
	else: # roi is a bud
		touching_rois = shrinkRoi(roi,touching_rois)
		IJ.run("Enlarge...","enlarge=1")
		rm.runCommand("Update")
		roi = rm.getRoi(index) # get updated enlarged roi
		
		if len(touching_rois) > 1:
			# evaluate touching rois returning one 'mother cell'
			touching_rois = evaluate_Bud(index,roi,touching_rois)
			if touching_rois == None: return False
		else:
			# evaluate cell cycle phase
			if roi_table.getEntry(touching_rois[0],"spb")==2:
				evaluate_G2_vs_PM(touching_rois[0],roi_table,rm,nuclei_table,spb_table)
			elif roi_table.getEntry(index,"spb")==1:
				# probably nuclei analysis was wrong because roi (index) has a
				# spb but no nucleus (very unlikely)
				roi_table.setEntry(touching_rois[0],"name","ANA")
			else:
				roi_table.setEntry(touching_rois[0],"name","Late S")
		
		roi_table.setEntry(index,"eval","yes")
		print "~~~ Evaluated ROI",index,roi_table.getEntry(index,"name")
		roi_table.setEntry(touching_rois[0],"eval","yes")
		print "~~~ Evaluated ROI",touching_rois[0],roi_table.getEntry(touching_rois[0],"name")
		
		combineTwoRois(index,touching_rois[0],roi_table,rm)
	return True

# ------------------------------------------------------------------------------
# evaluate_watershed()
# ------------------------------------------------------------------------------
def evaluate_watershed(index,rm,roi_table):
	"""
	Return list of rois that were seperated from param:'index' by watershed.
	Evaluation:
	* rois that are only 2 pixels apart are considered as "watershedded"
	"""
	# rois that are only 2 pixels apart can only (with a very high probability) 
	# have been seperated by watersehd
	# -> check for touching rois in an 2 pixels wider area:
	rm.deselect() 
	rm.select(index)
	IJ.run("Enlarge...","enlarge=2")
	rm.runCommand("Update")
	roi = rm.getRoi(index) # get updated enlarged roi
	rois = roi_table.getIndexByEntry("eval","no")
	rois.remove(index)
	rois2 = rois + [] # copy instead of reference

	# check if other rois are even near the investigated roi (faster)
	for check_roi in rois:
		x1 =  roi_table.getEntry(index,'X')
		x2 =  roi_table.getEntry(check_roi,'X')
		y1 =  roi_table.getEntry(index,'Y')
		y2 =  roi_table.getEntry(check_roi,'Y')
		w1 =  roi_table.getEntry(index,'width')
		w2 =  roi_table.getEntry(check_roi,'width')
		h1 =  roi_table.getEntry(index,'height')
		h2 =  roi_table.getEntry(check_roi,'height')
		# + 5 to encounter any rounding mistakes to be sure 
		if ( (abs(x1-x2) > (w1/2 + w2/2 + 5)) or (abs(y1-y2) > (h1/2 + h2/2 + 5)) ): 
			rois2.remove(check_roi)

	# now check if remaining other rois touch investigated roi
	touching_rois = []
	for check_roi in rois2:
		if touchingRoi(roi,rm.getRoi(check_roi)):
			touching_rois = touching_rois + [check_roi]

	IJ.run("Enlarge...","enlarge=-2") # shrink roi again
	rm.runCommand("Update")
	
	return(touching_rois)

# ------------------------------------------------------------------------------
# evaluate_Bud()										   [evaluate_NoiseOrBud]
# ------------------------------------------------------------------------------
def evaluate_Bud(index,roi,touching_rois):
	"""
	Return cell from param:'touching_rois' that suits best as the 'mother cell'.
	Evaluation:
	* if roi has no spb
	  ** potential mother cells = cells with two spbs
	  ** if there is one potential cells
	  	 *** FOUND mother cell. evaluate cell cycle phase
	  ** if there are more thann one potential cells
	  	 *** take the cell that has the spb with the minimum distance to the bud
	  	 *** FOUND mother cell.
	  ** if there is no potential mother cell
	  	 *** there might have been an ERROR in analyzing the images
	* if roi has a spb
	  ** take the cell that has the spb with the minimum distance to the buds spb
	  ** FOUND mother cell.
	"""
	touching_rois_spb = [ (i,roi_table.getEntry(i,"spb")) for i in touching_rois ]

	# evaluate roi at index
	if roi_table.getEntry(index,"spb")==0:
		# potential mother cells are those with two spb's
		potential_mothercells = [i for i,x in touching_rois_spb if x == 2]

		if len(potential_mothercells)==1:
			evaluate_G2_vs_PM(potential_mothercells[0], roi_table, rm, nuclei_table, spb_table)
			return(potential_mothercells)
		elif len(potential_mothercells)>1 :
			x = roi_table.getEntry(index,"X")
			y = roi_table.getEntry(index,"Y")
			# calculate distance to bud (roi at index) for each spb
			distances = sum([[ (((x-spb_table.getEntry(j,"X"))**2+(y-spb_table.getEntry(j,"Y"))**2)**(0.5),i) for j in roi_table.getEntry(i,"spb_id")] for i in potential_mothercells],[])
			evaluate_G2_vs_PM(sorted( distances )[0][1],roi_table,rm,nuclei_table,spb_table)
			return([ sorted( distances )[0][1] ])
		else:
			print "Could not evaluate bud. Leaving it to be manually evaluated"
			return (None)
	else:
		x = spb_table.getEntry(roi_table.getEntry(index,"spb_id")[0],"X")
		y = spb_table.getEntry(roi_table.getEntry(index,"spb_id")[0],"Y")
		# take the cell that has the spb with the minimum distance to the buds spb
		spbs = sum([[ (((x-spb_table.getEntry(j,"X"))**2+(y-spb_table.getEntry(j,"Y"))**2)**(0.5),i) for j in roi_table.getEntry(i,"spb_id")] for i in touching_rois],[])
		roi_table.setEntry(sorted( spbs )[0][1],"name","P/M")
		return([ sorted( spbs )[0][1] ])

# ------------------------------------------------------------------------------
# evaluate_G2_vs_PM()
# ------------------------------------------------------------------------------
def evaluate_G2_vs_PM(index, roi_table, rm, nuclei_table, spb_table):
	"""
	Alter roi_table entries according to the cell cycle evaluation of param:'index'
	Evaluation:  
	* get distance d between spbs and diameter D of the nucleus
	* if d > D cell is in P/M
	* if d < D cell is in G2   ('the official rules for segmentation')
	"""
	
	# get spb koordinates and (euclidean-)distance d between them
	sk = [(spb_table.getEntry(i,"X"),spb_table.getEntry(i,"Y")) for i in roi_table.getEntry(index,"spb_id")]
	d  = ( (sk[0][0]-sk[1][0])**2  + (sk[0][1]-sk[1][1])**2 )**(0.5)

	# get nuclei diameter D (mean of height and width)
	D = ( nuclei_table.getEntry(roi_table.getEntry(index,"n_id")[0],"width") + 
		  nuclei_table.getEntry(roi_table.getEntry(index,"n_id")[0],"height") )/2

	# evaluate cell
	if d > D: roi_table.setEntry(index,"name","P/M")
	else :    roi_table.setEntry(index,"name","G2")

# ------------------------------------------------------------------------------
# overlay_area()
# ------------------------------------------------------------------------------
def overlay_area(r, r2, rm):
	""" Return the area (number of pixels) that roi r and r2 have in common."""
	roi  = rm.getRoi(r)
	roi2 = rm.getRoi(r2)
	r = roi.getBounds()
	count = 0
	for x in range(r.width):
		for y in range(r.height):
			if roi.contains(r.x+x , r.y+y):
				if roi2.contains(r.x+x , r.y+y):
					count = count +1
	return(count)

# ------------------------------------------------------------------------------
# high_whi5()
# ------------------------------------------------------------------------------
def high_whi5(index,nucleus_id,roi_table):
	"""Return TRUE if whi5 intensity is considered high (whi5_diff <= -10) """
	whi5_diff = ((roi_table.getEntry(index,"whi5")/
				 (roi_table.getEntry(index,"area"))) -
				 (nuclei_table.getEntry(nucleus_id,"whi5")/
				 (nuclei_table.getEntry(nucleus_id,"area")))
				)
	if whi5_diff <= -10: return True
	else: return False	

# ------------------------------------------------------------------------------
# userDialog()
# ------------------------------------------------------------------------------
def userDialog(rm):
	""" Communicate with the user if he wants to alter rois or their cell cycle"""

	WaitForUserDialog("Cellsegmentation was finished", 
					 ("Please look at your images and make any neccessary "
					  "changes with the ROI Manager. \n You can delete ROIs "
					  "or add new ones using Fiji. \n When you press OK a next "
					  "window will let you change the cell cycle phases.")).show()

	# create dialog with a text box for each cell and its cell cycle phase
	gd = GenericDialog("Cell Cycle")  
	for r in range(rm.getCount()):
		if len(re.split("\ -\ ",rm.getName(r)))==1: gd.addStringField(("roi#"+str(r+1))," ")
		else:	gd.addStringField(("roi#"+str(r+1)), re.split("\ -\ ",rm.getName(r))[1] )   
	gd.setCancelLabel("Return to Roi Manager")
	gd.setOKLabel("OK and Proceed")
	gd.showDialog()  

	# process user input
	if gd.wasCanceled():
		# read out options anyway
		for r in range(rm.getCount()):
			name = gd.getNextString()
			rm.deselect()
			rm.select(r)
			rm.runCommand("Rename",(str(r+1)+ " - "+ name))
		return False	
	else:
		# Read out the options  
		for r in range(rm.getCount()):
			name = gd.getNextString()
			correct = name in ["G1", "g1", "G2", "g2",
							   "Ana", "ANA", "ana",
							   "PM", "pm", "p/m", "P/M",
							   "early s", "Early S", "EARLY S", "early_s", "Early_S",
							   "late s", "Late S", "LATE S", "late_s", "Late_S",
							   "tc","TC","t/c","T/C"]
			if name == "" or correct == False:
				WaitForUserDialog("Cellcylcle names were incomplete or false", 
								 ("You have not filled in all cell cycle names "
								  "or no correct ones. \n Please check again and "
								  "correct your input")).show()
				return False
		return True

# ------------------------------------------------------------------------------
# getImages()
# ------------------------------------------------------------------------------
def getImages():
	"""
	Return images saved by openImage().
	Call openImage() as long as the dialog wasn't canceled and until all images 
	have bin selected.
	"""
	# close all active windows that might interfere with plugin
	wins = WindowManager.getIDList()
	if wins != None:
		for w in wins:
			WindowManager.getImage(w).close()

	# get images
	print "i have been here"
	images = ImageHolder()
	o = False
	
	while (o == False):
		WaitForUserDialog("Choose images", 
						 ("Please choose the following images:"
						  "\nNIBA, CFP, WU, NG/CYS, BF")).show()
		if not openImages(images):
			return None
		if images.checkStatus():
			o = True

	return images

# ------------------------------------------------------------------------------
# openImages()													   [getImages()]
# ------------------------------------------------------------------------------
def openImages(images):
	""" Open user dialog and save images in an ImageHolder object."""
	o = Opener()
	o.openMultiple()
	wins = WindowManager.getIDList()
	if wins == None:
		return False
	for w in wins:
		img = WindowManager.getImage(w)
		name = img.getTitle()
		if not ("TIF" in name or "tif" in name):
			print "ERROR: No .tif image selected" 
		if "NIBA" in name:
			images.niba = img
		elif "WU" in name:
			images.wu = img
		elif "CFP" in name:
			images.cfp = img
		elif "BF" in name:
			images.bf = img
		elif "NG" in name: # or maybe other names ?
			images.mrna = img
	return True

# ------------------------------------------------------------------------------
# evaluate_no_nuclei()
# ------------------------------------------------------------------------------
def evaluate_no_nuclei(roi_table, rm):
	""" Call 'evaluate_noiseOrBud()' for each cell without an nucleus."""
	# reversed list because then if a roi is deleted, no shift in roi indices occurs
	cells_without_nuclei = roi_table.getIndexByEntry("nuclei", 0)[::-1]
	for index in cells_without_nuclei:
		print "~~~ Evaluate roi", index, roi_table.getEntry(index, "name")
		roi_table.setEntry(index, "eval", "yes")
		evaluate_noiseOrBud(index, rm, roi_table)
	return True

# ------------------------------------------------------------------------------
# evaluate_two_spbs()
# ------------------------------------------------------------------------------
def evaluate_two_spbs(roi_table, rm, nuclei_table, spb_table):
	"""
	Alter roi_table for each cell that has two spb's.
	Evaluation:
	* if cell has one nucleus: evaluate_G2_vs_PM()
	* if cell has more than one nucleus: ANAPHASE
	(* cell can't have no nucleus because it would have been already evaluated!)
	"""
	rois_to_evaluate = roi_table.getIndexByEntry("eval", "no")
	cells_with_two_spbs = [c for c in roi_table.getIndexByEntry("spb", 2)[::-1] if c in rois_to_evaluate]
	for index in cells_with_two_spbs:
		print "~~~ Evaluate roi", index, roi_table.getEntry(index, "name")
		if roi_table.getEntry(index, "nuclei")==1:
			evaluate_G2_vs_PM(index, roi_table, rm, nuclei_table, spb_table)
			roi_table.setEntry(index, "eval", "yes")
		else:
			roi_table.setEntry(index, "name", "ANA")
			roi_table.setEntry(index, "eval", "yes")
	return True

# ------------------------------------------------------------------------------
# evaluate_highIntensity()
# ------------------------------------------------------------------------------
def evaluate_highIntensity_spb(roi_table, spb_table):
	""" Alter roi_table for each cell that has a spb with high intensity."""
	rois_to_evaluate = roi_table.getIndexByEntry("eval", "no")
	cells_with_high_intensity_spb = [c for c in roi_table.getIndexByEntry("spb",1)[::-1] if roi_table.getEntry(c,"spb_id")[0] in spb_table.getIndexByEntry("high_intensity","yes") and c in roi_table.getIndexByEntry("eval","no")]
	for index in cells_with_high_intensity_spb:
		print "~~~ Evaluate roi", index, roi_table.getEntry(index, "name")
		roi_table.setEntry(index,"name","Late S")
		roi_table.setEntry(index,"eval","yes")
	return True

# ------------------------------------------------------------------------------
# evaluate_many_neighbours()
# ------------------------------------------------------------------------------
def evaluate_many_neighbours(roi_table, rm):
	"""
	Alter roi_table for each cell that has more than one neighbour.
	Assumption: If cells still stick together and are not evauated yet
				they are considered to be one cell.
	Evaluation:
	* choose another cell as a partner that has the most overlaying area
	* combine those cells
	* evaluate cell cycle phase based on whi5 intensity
	  ** high whi5 intensity: T/C
	  ** low whi5 intensity: ANAPHASE
	"""
	# if a cell has very near neighbours those cells have probably been seperated by watershed
	remaining_cells = [(c,evaluate_watershed(c,rm,roi_table)) for c in roi_table.getIndexByEntry("eval","no") ]
	cells_with_too_many_neighbours = [(c,w) for c,w in remaining_cells if len(w)>=2][::-1]
	for c,w in cells_with_too_many_neighbours:
		print "~~~ Evaluate roi", c, roi_table.getEntry(c, "name")
		rm.select(c)
		IJ.run("Enlarge...","enlarge=1")
		rm.runCommand("Update")
		chosen_cell = sorted([ (overlay_area(c,c1,rm),c1) for c1 in w ])[-1:][0][1]
		combineTwoRois(c,chosen_cell,roi_table,rm)
		
		#evaluate cell
		min_index = min(c,chosen_cell)
		nucleus_id = roi_table.getEntry(min_index,"n_id")[0]
		if high_whi5(min_index,nucleus_id,roi_table):
			roi_table.setEntry(min_index,"name","T/C")
		else:
			roi_table.setEntry(min_index,"name","ANA")
		roi_table.setEntry(min_index,"eval","yes")

# ------------------------------------------------------------------------------
# evaluate_one_neighbour()
# ------------------------------------------------------------------------------
def evaluate_one_neighbour(roi_table, rm):
	""" 
	Combine each cell with its one neighbour.
	Assumption: If cells still stick together and are not evauated yet
				they are considered to be one cell.
	"""
	remaining_cells = [(c, evaluate_watershed(c,rm,roi_table)) for c in roi_table.getIndexByEntry("eval", "no") ]
	cells_with_one_neighbour = [(c, w) for c,w in remaining_cells if len(w)==1][::-1]

	# remove duplicates ( if (a,b) are neighbouring cells, than also (b,a) -> remove one tuple )
	for i in range(len(cells_with_one_neighbour)/2) :
		cells_with_one_neighbour.remove((cells_with_one_neighbour[i][1][0], [cells_with_one_neighbour[i][0]]))

	for i in range(len(cells_with_one_neighbour)):
		print "~~~ Evaluate roi", cells_with_one_neighbour[i][0], roi_table.getEntry(cells_with_one_neighbour[i][0], "name")
		rm.select(cells_with_one_neighbour[i][0])
		IJ.run("Enlarge...", "enlarge=1")
		rm.runCommand("Update")
		combineTwoRois(cells_with_one_neighbour[i][0], cells_with_one_neighbour[i][1][0], roi_table, rm)

# ------------------------------------------------------------------------------
# evaluate_remaining()
# ------------------------------------------------------------------------------
def evaluate_remaining(roi_table, rm):
	""" 
	Alter roi_table for each remaining single cell.
	Evaluation:
	* if cell has one nucleus
	  ** if nucleus has a high whi5 intensity: G1
	  ** if nucleus has a low  whi5 intensity: EARLY S
	* if cell has more than one nucleus
	  ** if nucleus has a high whi5 intensity: T/C
	  ** if nucleus has a low  whi5 intensity: ANAPHASE
	"""
	rois_to_evaluate = roi_table.getIndexByEntry("eval", "no")
	for index in rois_to_evaluate:
		print "~~~ Evaluate roi", index, roi_table.getEntry(index, "name")
		nucleus_id = roi_table.getEntry(index, "n_id")[0]
		if roi_table.getEntry(index, "nuclei") == 1:
			if high_whi5(index, nucleus_id, roi_table):
				roi_table.setEntry(index, "name", "G1")
			else:
				roi_table.setEntry(index, "name", "Early S")
		else:	
			if high_whi5(index, nucleus_id, roi_table):
				roi_table.setEntry(index, "name", "T/C")
			else:
				roi_table.setEntry(index, "name", "ANA")
		roi_table.setEntry(index, "eval", "yes")

# ------------------------------------------------------------------------------
# renameRois()
# ------------------------------------------------------------------------------
def renameRois(rm, roi_table):
	""" Rename Rois according to roi_table entries."""
	rm.deselect()
	for index in range(rm.getCount()):
		rm.select(index)
		rm.runCommand("Rename",(str(index+1)+ " - "+ roi_table.getEntry(index, "name")))
		rm.deselect()

# ------------------------------------------------------------------------------
# createMask()
# ------------------------------------------------------------------------------
def createMask(mask, rm):
	""" Creates an RGB cell mask that can be processed by IDLmerger."""
	ip = mask.getProcessor()
	ip.setValue(255)
	ip.fill() # fills the whole picture black

	# now fill each roi with a different shade of grey
	rois = rm.getRoisAsArray()
	for index in range(len(rois)):
		colour = index*255/rm.getCount() # because cell numbers are low, not two should get the same colour
		ip.setValue(colour)
		ip.fill(rois[index])

# ------------------------------------------------------------------------------
# saveMask()
# ------------------------------------------------------------------------------
def saveMask(mask, rm):
	""" Ask user for a directory and store param:'mask' in there."""
	gdSave = GenericDialog("Save mask")
	gdSave.addMessage(("Press OK if the mask image should be stored in the same "
					   "directory where your input files came from."
			  		   "\nOr choose another directory."))
	gdSave.setCancelLabel("Choose directory...")
	gdSave.showDialog()
	if not gdSave.wasCanceled():
		path = IJ.getDirectory("current")
	else:
		path = DirectoryChooser("Choose directory to store mask file").getDirectory()
		if path == None: # user canceled directory dialog
			path = IJ.getDirectory("current")

	fs = FileSaver(mask)
	maskFile = path + images.niba.getShortTitle() + "_mask_cells.tif"
	fs.saveAsTiff(maskFile)
	return (path)

# ------------------------------------------------------------------------------
# proceed_with_IDLmerger()
# ------------------------------------------------------------------------------
def proceed_with_IDLmerger(path):
	""" If chosen, ask user for merger.py directory, save last preferences and 
	launch IDLmerger. """
	pathToMerger = DirectoryChooser("Choose directory to lastprefs file").getDirectory()
	# try to already set preferences with optained mask
	# try to set preferences with optained mask
	lastprefs = "last_preferences.pref"
	text = ("(dp0\n"
		   "S'locpath'\n"
		   "p1\n"
		   "S''\n"
		   "p2\n"
		   "sS'channeltokens'\n"
		   "p3\n"
		   "S''\n"
		   "p4\n"
		   "sS'mskpath'\n"
		   "p5\n"
		   "S'"+path+"'\n"
		   "p6\n"
		   "sS'outpath'\n"
		   "p7\n"
		   "S''\n"
		   "p8\n"
		   "s.")
	os.system("echo "+ '"'+text+'"'+ " > "+ lastprefs)
	#launch merger
	os.system("python "+ pathToMerger+ "merger.py &")

# ==============================================================================
# Main
# ==============================================================================

# get images through user Dialog
images = getImages()

# prepare imageJ instances (roi manager, particle analyzer, results table)
rm = RoiManager.getInstance() 
if rm != None: rm.reset()
else: rm = RoiManager()
table = ResultsTable()
options = ParticleAnalyzer.ADD_TO_MANAGER | ParticleAnalyzer.EXCLUDE_EDGE_PARTICLES | ParticleAnalyzer.INCLUDE_HOLES
pa = ParticleAnalyzer(options, ParticleAnalyzer.AREA, table, 0, 100000000)

# get initial rois
mask = getInitialROIs(images.niba, rm, pa) 

#get roi measurements such like area, koordinates, grey value, and whi5 intensity
roi_table = My_table(["name","area","nuclei","n_id","spb","spb_id","X","Y","eval","whi5","width","height"])
mean_grey_value = getRoiMeasurements(rm, roi_table, images.niba)

# get nuclei info for each cell
nuclei_table = countNuclei(rm, images.wu, pa, roi_table, images.niba) # update roi_table with nuclei counts

# get spindle pole body info for each cell
spb_table    = countAndMeasureSPB(rm, images.cfp, pa, roi_table, mean_grey_value)


# now evaluate each roi
# the order of evaluation is crucial!
print "\n======================= Start Evaluation =======================\n"
print "\n======================= Cells without nuclei are evaluated. "
evaluate_no_nuclei(roi_table, rm)
print "\n======================= Cells with two spindle pole bodies are evaluated. "
evaluate_two_spbs(roi_table, rm, nuclei_table, spb_table)
print "\n======================= Cells with a spb of high intensity are evaluated. "
evaluate_highIntensity_spb(roi_table, spb_table)
print "\n======================= Cells with too many neighbours are evaluated. "
evaluate_many_neighbours(roi_table, rm)
print "\n======================= Cells with one neighbour are evaluated. "
evaluate_one_neighbour(roi_table, rm)
print "\n======================= Cells than remain are evaluated. "
evaluate_remaining(roi_table, rm)
print "\n======================= Done Evaluation=======================\n"

renameRois(rm, roi_table)

wins = WindowManager.getIDList()
for w in wins:
	WindowManager.getImage(w).hide()
images.mrna.show()
rm.runCommand("Show All")
images.cfp.show()
rm.runCommand("Show All")
images.bf.show()
rm.runCommand("Show All")
Zniba = maxZprojection(images.niba, 7, images.niba.getNSlices()-5)
Zniba.show()
rm.runCommand("Show All")
IJ.run("Tile")
print "### Done."

u = userDialog(rm)
while(u == False): u = userDialog(rm)

rm = RoiManager.getInstance()
rm.runCommand("Show None")
wins = WindowManager.getIDList()
for w in wins:
	WindowManager.getImage(w).close()

# create RGB cell mask and save it
createMask(mask, rm)
path = saveMask(mask, rm)

# ask user if she/he wants to proceed using the mask with IDL merger
gd = GenericDialog("Choose merger.py directory")
gd.addMessage("Do you want to proceed with the optained cell mask "
			  "and launch the IDL merger?"
			  "\nBe aware that you need a loc file in addition.")
gd.setOKLabel("Proceed with IDL merger")
gd.showDialog()
if gd.wasOKed(): proceed_with_IDLmerger(path)
print "Done for real."
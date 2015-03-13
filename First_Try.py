import re
import os.path
from ij import ImagePlus, IJ, Prefs
from ij.io import OpenDialog
from ij.plugin import ZProjector
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable
from ij.plugin.filter import ParticleAnalyzer
from ij.plugin.filter import Binary, MaximumFinder
from java.awt import Rectangle
from ij.process import ImageProcessor

import sys
sys.path.insert(0,'/home/basar/Personal/Svenja/Scripts/')
from table import My_table

def maxZprojection(stackimp,start,stop):
    zp = ZProjector(stackimp)
    zp.setMethod(ZProjector.MAX_METHOD)
    zp.setStartSlice(start)
    zp.setStopSlice(stop)
    zp.doProjection()
    zpimp = zp.getProjection()
    return zpimp

def getVal(table, i):
	x = table.getValue('X', i)
	y = table.getValue('Y', i)
	return [x, y] 

def calculateThresholdValue(imp):
	##
	## implements tha manually definable percentage of pixel when using Image>>Adjust>>Threshold
	##
	# calculates a histogram of pixel values
	ma = imp.getProcessor().getMax()
	mi = imp.getProcessor().getMin()
	hist = [0 for i in range(256)]
	steps = round((ma-mi)/256,3 )
	
	for x in range(imp.width):
		for y in range(imp.height):
			index = int( (imp.getPixel(x,y)[0]-mi)/steps )
			if index >= 256: index = 255
			hist[index] = hist[index] + 1
	
	i = 1
	val = 0 # percent
	while val < 0.19: 
		val = float(sum(hist[-i:])*100.0)/float(imp.width*imp.height)
		i = i + 1
	
	return(int((256-i)*steps) + mi)


def getInitialROIs(niba,pa):
	
	Zniba = maxZprojection(niba,10,15) # identify "good" slices instead of setting the start(10) and stop(15) slice?
	
	IJ.run(Zniba, "Enhance Contrast","saturated Pixel=0.5")
	IJ.run(Zniba, "Smooth","")
	IJ.run(Zniba, "Bandpass Filter...","filter_large=200 ; filter_small=10")
	IJ.run(Zniba, "Gaussian Blur...","radius=10")
	IJ.run(Zniba, "Threshold","Triangle")
	
	#table = ResultsTable()
	rm = RoiManager.getInstance() 
	if rm != None: rm.reset()
	else: rm = RoiManager()

	#options = ParticleAnalyzer.ADD_TO_MANAGER | ParticleAnalyzer.EXCLUDE_EDGE_PARTICLES
	#pa = ParticleAnalyzer(options, ParticleAnalyzer.AREA, table, 0, 100000000)
	#pa.setHideOutputImage(True)
	if pa.analyze(Zniba):
		print "All ok"
	else:
		print "There was a problem in analyzing"
	
	rm = RoiManager.getInstance()
	Zniba.show()
	#rm.show()
	#rois = rm.getRoisAsArray()
	
	for e in range(rm.getCount()):
		rm.deselect()
		rm.select(e)
		IJ.run("Enlarge...","enlarge=6") 
		IJ.run(Zniba,"Fill","")
	
	rm.deselect()
	#pause = OpenDialog("Choose Track Data...", "")
	rm.reset()
	IJ.run(Zniba, "Watershed","")
	
	if pa.analyze(Zniba):
		print "All ok"
	else:
		print "There was a problem in analyzing"

	rm = RoiManager.getInstance()
	Zniba.show()
	rm.show()
	return([Zniba,rm,table])

def countNuclei(rm,wu,raw_wu,pa,roi_table):

	if raw_wu: Zwu = maxZprojection(wu,10,15)
	else: Zwu = wu
	
	IJ.run(Zwu, "Enhance Contrast","saturated Pixel=0.5")
	IJ.run(Zwu, "Smooth","")
	IJ.run(Zwu, "Gaussian Blur...","radius=5 ; accuracy=0.01")
	IJ.run(Zwu, "Threshold","dark background")
	IJ.run(Zwu, "Convert to Mask", "")

	old_roi_count = rm.getCount()
	rois = rm.getRoisAsArray()
	
	pa.analyze(Zwu)	
	nuclei_table = My_table(["area","nuclei_intensity","X","Y","width","height"])
	
	rt = ResultsTable.getResultsTable()
	rt.reset()
	IJ.run("Set Measurements...","area ; centroid ;  integrated density ; bounding rectangle")
	rm.deselect()
	for i in range(old_roi_count,rm.getCount()):
		rm.select(Zwu,i)
		IJ.run(Zwu,"Measure","")
		i = i - old_roi_count # result table starts with 0
		nuclei_table.addRow([rt.getValue('Area',i),rt.getValue('IntDen',i),
						  rt.getValue('X',i),rt.getValue('Y',i),rt.getValue('Width',i),rt.getValue('Height',i)])
	av = sum(nuclei_table.getColumn("area"))/nuclei_table.count
			
	for i in range(nuclei_table.count):
		# if nuclei is smaller than 10% of the average nuclei area it just might be noise
		if nuclei_table.getEntry(i,"area") >= 0.1*av: 
			x = int(nuclei_table.getEntry(i,'X'))
			y = int(nuclei_table.getEntry(i,'Y'))
			for roi in rois:
				if roi.contains(x,y): 
					roi_table.setEntry(rm.getRoiIndex(roi),"nuclei", roi_table.getEntry(rm.getRoiIndex(roi),"nuclei") + 1)
					roi_table.setEntry(rm.getRoiIndex(roi),"nuclei_id", roi_table.getEntry(rm.getRoiIndex(roi),"nuclei_id") + [i])

	rt.reset()		
	rm.setSelectedIndexes(range(old_roi_count,rm.getCount()))
	rm.runCommand("delete") # delete spindle pole body rois
	
	return(nuclei_table)


def countAndMeasureSPB(rm,cfp,pa,roi_table): # spindel-pole-bodies

	old_roi_count = rm.getCount()
	rois = rm.getRoisAsArray()
	
	Zcfp = maxZprojection(cfp,1,cfp.getNSlices()) # take all slices
	
	lower_threshold	= calculateThresholdValue(Zcfp)
	print lower_threshold
	Zcfp.getProcessor().setThreshold(lower_threshold, Zcfp.getProcessor().getMax(), ImageProcessor.NO_LUT_UPDATE)
	IJ.run(Zcfp, "Threshold","dark background")
	IJ.run(Zcfp, "Despeckle","")
	Zcfp.show()
	
	pa.analyze(Zcfp)	
	spb_table = My_table(["area","spb_intensity","high_intensity","X","Y"])
	
	rt = ResultsTable.getResultsTable()
	rt.reset()
	IJ.run("Set Measurements...","area ; min & max grey value; mean grey value; centroid ;  integrated density")
	rm.deselect()
	for i in range(old_roi_count,rm.getCount()):
		rm.select(Zcfp,i)
		IJ.run(Zcfp,"Measure","")
		i = i - old_roi_count # result table starts with 0
		spb_table.addRow([rt.getValue('Area',i),rt.getValue('IntDen',i),"no",
						  rt.getValue('X',i),rt.getValue('Y',i)])

	to_delete = []
	for i in range(spb_table.count):
		rm.deselect()
		rt.reset()
		if spb_table.getEntry(i,"area") < 3:
			rm.select(Zcfp,old_roi_count + i)
			IJ.run("Enlarge...","enlarge=3")
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
	
	for i in range(spb_table.count):
		x = int(spb_table.getEntry(i,'X'))
		y = int(spb_table.getEntry(i,'Y'))
		for roi in rois:
			if roi.contains(x,y): 
				roi_table.setEntry(rm.getRoiIndex(roi),"spb", roi_table.getEntry(rm.getRoiIndex(roi),"spb") + 1)
				roi_table.setEntry(rm.getRoiIndex(roi),"spb_id", roi_table.getEntry(rm.getRoiIndex(roi),"spb_id") + [i])

	rt.reset()		
	rm.setSelectedIndexes(range(old_roi_count,rm.getCount()))
	rm.runCommand("delete") # delete spindle pole body rois
	
	for i in range(roi_table.count):
		if roi_table.getEntry(i,"spb") > 2:
			spbs = sorted([(spb_table.getEntry(j,'area'),j) for j in roi_table.getEntry(i,"spb_id")])
			spbs = spbs[-2:] # only take the biggest two spindle pole bodys for the cell
			roi_table.setEntry(i,"spb",2)
			roi_table.setEntry(i,"spb_id",[spbs[j][1] for j in [0,1]])
	
	spbs = sum(roi_table.getColumn("spb_id"),[])
	intensities = sorted( [ spb_table.getEntry(j,"spb_intensity") for j in spbs] )
	# calculate median
	if len(intensities) %2 == 1: av_intensity = intensities[((len(intensities)+1)/2)-1]
	else: av_intensity = float(sum(intensities[(len(intensities)/2)-1:(len(intensities)/2)+1]))/2.0
	
	for i in spbs:
		if spb_table.getEntry(i,"spb_intensity") >= 2*av_intensity:
			spb_table.setEntry(i,"high_intensity","yes")
			
	return(spb_table)
	
def getRoiMeasurements(rm,roi_table):
	rois = rm.getRoisAsArray()
	n = len(rois)
	indexes = rm.getSelectedIndexes() # save cuurent selection
	rm.deselect()

	IJ.run("Set Measurements...","area ; centroid ; bounding rectangle")
	
	rt = ResultsTable.getResultsTable()
	rt.reset()
	for i in range(n):
		rm.select(i)
		IJ.run("Measure","")
		roi_table.addRow([rm.getName(i),rt.getValue('Area',i),0,[],0,[],rt.getValue('X',i),rt.getValue('Y',i),"no",
			rt.getValue('BX',i),rt.getValue('BY',i),rt.getValue('Width',i),rt.getValue('Height',i)])
		
	rt.reset()
	rm.deselect() 
	rm.setSelectedIndexes(indexes) # restore selection
	return(True)

def touchingRoi(roi,roi2):
	r = roi.getBounds()
	for x in range(r.width):
		for y in range(r.height):
			if roi.contains(r.x+x , r.y+y):
				if roi2.contains(r.x+x , r.y+y):
					return(True)	
	return(False)
	
def shrinkRoi2(roi,touching_rois):
	n = len(touching_rois)
	index = rm.getRoiIndex(roi)
	while (n>=1):
		IJ.run("Enlarge...","enlarge=-1") # shrink roi
		rm.runCommand("Update")
		roi = rm.getRoi(index) # get shrinked roi
		r = 0
		while( r < n and n>=1):
			if touchingRoi(roi,touching_rois[r]) == False :
				if n == 1 : 
					return(touching_rois)
				else:
					touching_rois = touching_rois[:r] + touching_rois[r+1:]
					n = n - 1
					r = r - 1
			r = r + 1

def shrinkRoi(roi,touching_rois):
	# shrink roi until it touches other cells only slightly
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

def combineTwoRois(index,index2,roi_table,rm):
	rm.deselect()
	rm.setSelectedIndexes([index,index2])
	rm.runCommand("OR")
	rm.runCommand("Update") # "Or" and "Update" will combine the two rois and save this new roi in the first roi (smallest index)
	rm.deselect()
	
	if index < index2: 
		# roi at index2 will be deleted -> store information in roi at index
		rm.select(index2)
		rm.runCommand("delete")
		roi_table.setEntry(index,"nuclei",roi_table.getEntry(index2,"nuclei")+roi_table.getEntry(index,"nuclei")) 
		roi_table.setEntry(index,"nuclei_id",roi_table.getEntry(index2,"nuclei_id")+roi_table.getEntry(index,"nuclei_id"))
		roi_table.setEntry(index,"spb",roi_table.getEntry(index2,"spb")+roi_table.getEntry(index,"spb"))
		roi_table.setEntry(index,"spb_id",roi_table.getEntry(index2,"spb_id")+roi_table.getEntry(index,"spb_id"))
		roi_table.setEntry(index,"area",roi_table.getEntry(index2,"area")+roi_table.getEntry(index,"area"))
		roi_table.setEntry(index,"name",roi_table.getEntry(index2,"name"))
		roi_table.delRow(index2)
		return(True)
	else: 
		# roi at index will be deleted -> store information in roi at index2
		rm.select(index)
		rm.runCommand("delete")
		roi_table.setEntry(index2,"nuclei",roi_table.getEntry(index2,"nuclei")+roi_table.getEntry(index,"nuclei"))
		roi_table.setEntry(index2,"nuclei_id",roi_table.getEntry(index2,"nuclei_id")+roi_table.getEntry(index,"nuclei_id"))
		roi_table.setEntry(index2,"spb",roi_table.getEntry(index,"spb")+roi_table.getEntry(index2,"spb"))
		roi_table.setEntry(index2,"spb_id",roi_table.getEntry(index,"spb_id")+roi_table.getEntry(index2,"spb_id"))
		roi_table.setEntry(index2,"area",roi_table.getEntry(index2,"area")+roi_table.getEntry(index,"area"))
		roi_table.delRow(index)
		
	
def evaluate_noiseOrBud(index,roi,av,roi_area,rois_to_evaluate,rm,roi_table):
	print "Evaluate: Noise or Bud" 
	rm.deselect() 
	rm.select(index)
	
	# enlarge particle to half the size of an average cell 
	# enlarge by 1 pixel enlarges area by about 100 (if cell is round) -> divide by 100
	pix = int((av/2 - roi_area)/100)
	if pix <= 1: pix = 2
	s_pix	= "enlarge=" + str(pix)
	
	rm.select(index)
	IJ.run("Enlarge...",s_pix)
	rm.runCommand("Update")
	roi = rm.getRoi(index) # get updated enlarged roi
	
	# check if rois are even near the investigated roi
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
	
	# check if rois in range touch investigated roi
	touching_rois = []
	for check_roi in rois_in_range:
		if touchingRoi(roi,rm.getRoi(check_roi)):
			if roi_table.getEntry(check_roi,"nuclei") > 0:
				touching_rois = touching_rois + [check_roi]
				print "found a touching roi for roi " + str(index)
			
	if len(touching_rois) == 0: 
		rm.runCommand("Delete")
		roi_table.delRow(index)
		return(False,index) 
	else:
		touching_rois = shrinkRoi(roi,touching_rois)
		IJ.run("Enlarge...","enlarge=1")
		rm.runCommand("Update")
		roi = rm.getRoi(index) # get updated enlarged roi
		
		if len(touching_rois) > 1: touching_rois = evaluate_Bud(index,roi,touching_rois)
		else: 
			if roi_table.getEntry(touching_rois[0],"spb")==2:
				evaluate_G2_vs_PM(touching_rois[0],roi_table,rm,nuclei_table,spb_table)
			else: roi_table.setEntry(index,"name","cc: Late S")
		
		print touching_rois
		combineTwoRois(index,touching_rois[0],roi_table,rm)
		return(True,touching_rois[0])

def evaluate_watershed(index,roi,rm,roi_table):
	# rois that are only 2 pixels apart can only (with a very high probability) 
	# been seperated by watersehd
	# -> check for touching rois in an 2 pixels wider area:
	rm.deselect() 
	rm.select(index)
	
	IJ.run("Enlarge...","enlarge=2")
	rm.runCommand("Update")
	roi = rm.getRoi(index) # get updated enlarged roi

	rois = rm.getRoisAsArray()
	rois_without_current_roi = rois[:index] + rois[index+1:]
	
	# check if rois are even near the investigated roi
	rois_in_range = []
	for check_roi in rois_without_current_roi:
		x1 =  roi_table.getEntry(index,'X')
		x2 =  roi_table.getEntry(rm.getRoiIndex(check_roi),'X')
		y1 =  roi_table.getEntry(index,'Y')
		y2 =  roi_table.getEntry(rm.getRoiIndex(check_roi),'Y')
		w1 =  roi_table.getEntry(index,'width')
		w2 =  roi_table.getEntry(rm.getRoiIndex(check_roi),'width')
		h1 =  roi_table.getEntry(index,'height')
		h2 =  roi_table.getEntry(rm.getRoiIndex(check_roi),'height')
		if ( (abs(x1-x2) < w1/2 + w2/2) and (abs(y1-y2) < h1/2 + h2/2) ):
			rois_in_range = rois_in_range + [check_roi]
	
	# check if rois in range touch investigated roi
	touching_rois = []
	for check_roi in rois_in_range:
		if touchingRoi(roi,check_roi):
			touching_rois = touching_rois + [check_roi]

	IJ.run("Enlarge...","enlarge=-2")
	rm.runCommand("Update")
	
	return(touching_rois)

def evaluate_Bud(index,roi,touching_rois):
	touching_rois_spb = [ (i,roi_table.getEntry(i,"spb")) for i in touching_rois ]
	if roi_table.getEntry(index,"spb")==0:
		# potential mother cells are those with two spb's
		potential_mothercells = [i for i,x in touching_rois_spb if x == 2]
		if len(potential_mothercells)==1:
			evaluate_G2_vs_PM(potential_mothercells[0],roi_table,rm,nuclei_table,spb_table)
			return(potential_mothercells)
		if len(potential_mothercells)>1 :
			x = roi_table.getEntry(index,"X")
			y = roi_table.getEntry(index,"Y")
			# take the cell that has the spb with the minimum distance to the bud (bud centroid)
			spbs = sum([[ (((x-spb_table.getEntry(j,"X"))**2+(y-spb_table.getEntry(j,"Y"))**2)**(0.5),i) for j in roi_table.getEntry(i,"spb_id")] for i in potential_mothercells],[])
			evaluate_G2_vs_PM(sorted( spbs )[0][1],roi_table,rm,nuclei_table,spb_table)
			return([ sorted( spbs )[0][1] ])
		else:
			print "Could not evaluate bud. Leaving it to be manually evaluated"
	else:
		x = spb_table.getEntry(roi_table.getEntry(index,"spb_id")[0],"X")
		y = spb_table.getEntry(roi_table.getEntry(index,"spb_id")[0],"Y")
		# take the cell that has the spb with the minimum distance to the buds spb
		spbs = sum([[ (((x-spb_table.getEntry(j,"X"))**2+(y-spb_table.getEntry(j,"Y"))**2)**(0.5),i) for j in roi_table.getEntry(i,"spb_id")] for i in touching_rois],[])
		roi_table.setEntry(sorted( spbs )[0][1],"name","cc: P/M")
		return([ sorted( spbs )[0][1] ])

def evaluate_G2_vs_PM(index,roi_table,rm,nuclei_table,spb_table):
	## if cell has two spindle-pole-bodys but one nucleus
	## this function decides wether the cell cycle phase is G2 or P/M
	## the decision is made by comparing the distance of the spb's to the nucleus diameter
	
	# get spb koordinates
	sk = [(spb_table.getEntry(i,"X"),spb_table.getEntry(i,"Y")) for i in roi_table.getEntry(index,"spb_id")]
	print sk,roi_table.getEntry(index,"spb_id")
	d  = ( (sk[0][0]-sk[1][0])**2  + (sk[0][1]-sk[1][1])**2 )**(0.5) # euclidean distance

	# get nuclei diameter (mean of height and width)
	D = ( nuclei_table.getEntry(roi_table.getEntry(index,"nuclei_id")[0],"width") + 
		  nuclei_table.getEntry(roi_table.getEntry(index,"nuclei_id")[0],"height") )/2
	print "spb distance:",d,"nucleus diameter:",D
	if d > D: roi_table.setEntry(index,"name","cc: P/M")
	else :    roi_table.setEntry(index,"name","cc: G2")
	

############## main #######################

# get images
op = OpenDialog("Choose Track Data...", "")
path = op.getDirectory()
imgName = op.getFileName()

# make regular expression out of imgName
whichIMG = re.split("w",imgName)
whichIMG = whichIMG[0]

cfp  = IJ.openImage(path + whichIMG + "w1CFP.TIF")  # spindle pole bodies
niba = IJ.openImage(path + whichIMG + "w2NIBA.TIF") # TF Whi5
ng   = IJ.openImage(path + whichIMG + "w4NG.TIF")   # mRNA spotting
print os.path.exists( path + "MAX_" + whichIMG + "w5WU.tif"), path + "MAX_" + whichIMG + "w5WU.tif"
if os.path.exists( path + "MAX_" + whichIMG + "w5WU.tif"):
	raw_wu = False
	wu   = IJ.openImage(path + "MAX_" + whichIMG + "w5WU.tif")
	print "Use user input max projection of WU.tif"
else:
	raw_wu = True
	wu   = IJ.openImage(path + whichIMG + "w5WU.TIF")   # DAPI staining of nuclei
bf   = IJ.openImage(path + whichIMG + "w6BF.TIF")   # BrightField

table = ResultsTable()
options = ParticleAnalyzer.ADD_TO_MANAGER | ParticleAnalyzer.EXCLUDE_EDGE_PARTICLES | ParticleAnalyzer.INCLUDE_HOLES
pa = ParticleAnalyzer(options, ParticleAnalyzer.AREA, table, 0, 100000000)


# get rois
initROI = getInitialROIs(niba,pa)
initROI[0].show() # pic
initROI[1].show() # manager
rm = initROI[1]
Zniba = initROI[0]

roi_table = My_table(["name","area","nuclei","nuclei_id","spb","spb_id","X","Y","eval","BX","BY","width","height"])

getRoiMeasurements(rm,roi_table)

	
nuclei_table = countNuclei(rm,wu,raw_wu,pa,roi_table) # update roi_table with nuclei counts
spb_table    = countAndMeasureSPB(rm,cfp,pa,roi_table)

print spb_table
print nuclei_table

# now evaluate each roi

av = sum(roi_table.getColumn('area'))/rm.getCount() # average area

cells_without_nuclei = roi_table.getIndexByEntry("nuclei",0)[::-1] 
# reversed list -> if a roi is deleted no shift in roi indices occurs

for index in cells_without_nuclei:
	roi = rm.getRoi(index)
	roi_area = roi_table.getEntry(index,"area")
	print "### Evaluate roi" ,rm.getName(index),"with no nuclei and an area of",roi_area,":"
	# if roi has no nucleus and is smaller than half an average cell (do not depend merly on nuclei analysis)
	# it is either noise or a bud
	if roi_area <= av/2:
		roi_table.setEntry(index,"eval","yes")
		rois_to_evaluate = roi_table.getIndexByEntry("eval","no")
		# was roi disconnected from (an)other cell(s) by watershed?
		watershedded = evaluate_watershed(index,roi,rm,roi_table)
		if len(watershedded) != 0: print "was watershedded"
		
		bud,mothercell_index = evaluate_noiseOrBud(index,roi,av,roi_area,rois_to_evaluate,rm,roi_table)
		if mothercell_index != index : roi_table.setEntry(mothercell_index,"eval","yes")

cells_with_two_spbs = [c for c in roi_table.getIndexByEntry("spb",2)[::-1] if c in roi_table.getIndexByEntry("eval","no")]

for index in cells_with_two_spbs:
	print "### Evaluate roi" ,rm.getName(index),"with two spindle poly bodies"
	if roi_table.getEntry(index,"nuclei")==1:
		evaluate_G2_vs_PM(index,roi_table,rm,nuclei_table,spb_table)
		roi_table.setEntry(index,"eval","yes")
	if roi_table.getEntry(index,"nuclei")>1:
		roi_table.setEntry(index,"name","cc: ANA")
		roi_table.setEntry(index,"eval","yes")
		

print roi_table
print "Still to evaluate:",roi_table.getIndexByEntry("eval","no")
print "### Done."
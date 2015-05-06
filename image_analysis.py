import re
import os.path
from ij import ImagePlus, IJ, Prefs, WindowManager
from ij.io import OpenDialog, FileSaver, Opener, ImportDialog
from ij.plugin import ZProjector
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable
from ij.plugin.filter import ParticleAnalyzer
from ij.plugin.filter import Binary, MaximumFinder
from java.awt import Rectangle
from ij.process import ImageProcessor
from ij.gui import GenericDialog,WaitForUserDialog

import sys
sys.path.insert(0,'/home/basar/Personal/Svenja/Scripts/')
from class_ImageHolder import ImageHolder
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

def calculateThresholdValue(imp,percentage):
	##
	## implements tha manually definable percentage of pixel 
	## when thresholding (using Image>>Adjust>>Threshold)
	##
	
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

def getInitialROIs(niba, rm, pa):
	# in order to work properly the image niba must be open when given to this function
	if niba.getNSlices() != 1:
		# todo: identify "good" slices instead of setting the start(7) and stop(n-5) slice?
		Zniba = maxZprojection(niba, 7, niba.getNSlices()-5)
		print "hello"
	else: Zniba = niba.duplicate()

	IJ.run(Zniba, "Enhance Contrast","saturated Pixel=0.5")
	IJ.run(Zniba, "Subtract Background...","radius = 100")
	IJ.run(Zniba, "Bandpass Filter...","filter_large=150 ; filter_small=10")
	IJ.run(Zniba, "Smooth","")
	IJ.run(Zniba, "Gaussian Blur...","radius=10")

	Zniba.show()
	IJ.run(Zniba, "Threshold...","")
	WaitForUserDialog("Segmentation Threshold", "Please alter threshold before you proceed with OK").show()
	
	Zniba.show()
	IJ.run(Zniba, "Make Binary","")
	IJ.run(Zniba, "Watershed","")

	if pa.analyze(Zniba):
		print "All ok"
	else:
		print "There was a problem in analyzing"

	#Zniba.show()
	rm.show()
	#Zniba.hide()
	#Zniba.close()
	return([Zniba,rm,table])

def countNuclei(rm, wu, pa, roi_table, niba):
	
	if wu.getNSlices() != 1: Zwu = maxZprojection(wu,10,15)
	else: 
		Zwu = wu.duplicate()
		Zwu.killRoi()# if any roi appears on image ignore it. if no roi appears this will have no influence
	
	IJ.run(Zwu, "Enhance Contrast","saturated Pixel=0.5")
	IJ.run(Zwu, "Smooth","")
	IJ.run(Zwu, "Gaussian Blur...","radius=12")
	Zwu.show()
	
	lower_threshold	= calculateThresholdValue(Zwu,3.3)
	
	Zwu.getProcessor().setThreshold(lower_threshold, Zwu.getProcessor().getMax(), ImageProcessor.NO_LUT_UPDATE)
	IJ.run(Zwu, "Threshold","dark background")
	old_roi_count = rm.getCount()
	rois = rm.getRoisAsArray()
	
	pa.analyze(Zwu)
	
	nuclei_table = My_table(["area","whi5","X","Y","width","height"])
	Zniba = maxZprojection(niba, 7, niba.getNSlices()-5);#Zniba.show()
	rt = ResultsTable.getResultsTable()
	rt.reset()
	IJ.run("Set Measurements...","area ; centroid ;  integrated density ; bounding rectangle")
	rm.deselect()
	
	for i in range(old_roi_count,rm.getCount()):
		rm.select(Zniba,i)
		IJ.run(Zniba,"Measure","")
		i = i - old_roi_count # result table starts with 0
		nuclei_table.addRow([rt.getValue('Area',i),rt.getValue('IntDen',i),
						  rt.getValue('X',i),rt.getValue('Y',i),rt.getValue('Width',i),rt.getValue('Height',i)])
	av = sum(nuclei_table.getColumn("area"))/nuclei_table.count
			
	for i in range(nuclei_table.count):
		x = int(nuclei_table.getEntry(i,'X'))
		y = int(nuclei_table.getEntry(i,'Y'))
		for roi in rois:
			if roi.contains(x,y): 
				roi_table.setEntry(rm.getRoiIndex(roi),"nuclei", roi_table.getEntry(rm.getRoiIndex(roi),"nuclei") + 1)
				roi_table.setEntry(rm.getRoiIndex(roi),"n_id", roi_table.getEntry(rm.getRoiIndex(roi),"n_id") + [i])

	rt.reset()		
	rm.setSelectedIndexes(range(old_roi_count,rm.getCount()))
	rm.runCommand("delete") # delete spindle pole body rois
	#Zniba.hide()
	return(nuclei_table)


def countAndMeasureSPB(rm,cfp,pa,roi_table,mean_grey_value): # spindel-pole-bodies

	old_roi_count = rm.getCount()
	rois = rm.getRoisAsArray()
	
	Zcfp = maxZprojection(cfp,1,cfp.getNSlices()) # take all slices
	
	lower_threshold	= calculateThresholdValue(Zcfp,0.19)
	Zcfp.getProcessor().setThreshold(lower_threshold, Zcfp.getProcessor().getMax(), ImageProcessor.NO_LUT_UPDATE)
	IJ.run(Zcfp, "Threshold","dark background")
	IJ.run(Zcfp, "Despeckle","")
	#Zcfp.show()

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
	
	for i in range(spb_table.count):
		x = int(spb_table.getEntry(i,'X'))
		y = int(spb_table.getEntry(i,'Y'))
		for roi in rois:
			# because cast to int causes round off, check for possible points
			if roi.contains(x,y) or roi.contains(x+1,y) or roi.contains(x,y+1) or roi.contains(x+1,y+1):  
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
		if spb_table.getEntry(i,"spb_intensity") >= 2.1*av_intensity:
			spb_table.setEntry(i,"high_intensity","yes")

	#Zcfp.hide()
	rt.reset()
	return(spb_table)
	
def getRoiMeasurements(rm,roi_table,niba):
	
	rois = rm.getRoisAsArray()
	n = len(rois)
	indexes = rm.getSelectedIndexes() # save current selection
	rm.deselect()
	Zniba = maxZprojection(niba, 7, niba.getNSlices()-5);Zniba.show()
	rt = ResultsTable.getResultsTable()
	rt.reset()
	
	IJ.run("Set Measurements...","area ; mean grey value ; centroid ;  integrated density ; bounding rectangle")
	
	# get mean grey value
	rm.deselect()
	rm.setSelectedIndexes(range(rm.getCount()))
	IJ.run(Zniba,"Make Inverse","")
	IJ.run(Zniba,"Measure","")
	mean_grey_value = rt.getValue('Mean',0)
	rt.reset()
	rm.deselect()

	# get measurements
	for i in range(n):
		rm.select(Zniba,i)
		IJ.run(Zniba,"Measure","")
		roi_table.addRow([rm.getName(i)[:4],rt.getValue('Area',i),0,[],0,[],int(rt.getValue('X',i)),int(rt.getValue('Y',i)),
		"no",rt.getValue('IntDen',i),rt.getValue('Width',i),rt.getValue('Height',i)])
		
	rt.reset()
	rm.deselect() 
	rm.setSelectedIndexes(indexes) # restore selection
	Zniba.hide()
	return(mean_grey_value)

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
	print "Combined",index,"and",index2
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
		roi_table.setEntry(index,"n_id",roi_table.getEntry(index2,"n_id")+roi_table.getEntry(index,"n_id"))
		roi_table.setEntry(index,"spb",roi_table.getEntry(index2,"spb")+roi_table.getEntry(index,"spb"))
		roi_table.setEntry(index,"spb_id",roi_table.getEntry(index2,"spb_id")+roi_table.getEntry(index,"spb_id"))
		roi_table.setEntry(index,"area",roi_table.getEntry(index2,"area")+roi_table.getEntry(index,"area"))
		roi_table.setEntry(index,"whi5",roi_table.getEntry(index2,"whi5")+roi_table.getEntry(index,"whi5"))
		roi_table.setEntry(index,"name",roi_table.getEntry(index2,"name"))
		roi_table.delRow(index2)
		return(True)
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
				print "Found a touching roi ",check_roi,"for roi " + str(index)
			
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
			else:
				if roi_table.getEntry(index,"spb")==1:
					roi_table.setEntry(touching_rois[0],"name","ANA")
					# probably nuclei analysis was wrong
				else:
					roi_table.setEntry(touching_rois[0],"name","Late S")
		
		roi_table.setEntry(index,"eval","yes")
		print "~~~ Evaluated ROI",index,roi_table.getEntry(index,"name")
		roi_table.setEntry(touching_rois[0],"eval","yes")
		print "~~~ Evaluated ROI",touching_rois[0],roi_table.getEntry(touching_rois[0],"name")
		
		combineTwoRois(index,touching_rois[0],roi_table,rm)
		return(True,touching_rois[0])

def evaluate_watershed(index,rm,roi_table):
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
	rois2 = rois + [] # copy
	# check if rois are even near the investigated roi
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
		if ( (abs(x1-x2) > (w1/2 + w2/2 + 5)) or (abs(y1-y2) > (h1/2 + h2/2 + 5)) ): rois2.remove(check_roi)
	# check if rois in range touch investigated roi
	touching_rois = []
	for check_roi in rois2:
		if touchingRoi(roi,rm.getRoi(check_roi)):
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
		roi_table.setEntry(sorted( spbs )[0][1],"name","P/M")
		return([ sorted( spbs )[0][1] ])

def evaluate_G2_vs_PM(index,roi_table,rm,nuclei_table,spb_table):
	##
	## if cell has two spindle-pole-bodys but one nucleus
	## this function decides wether the cell cycle phase is G2 or P/M
	## the decision is made by comparing the distance of the spb's to the nucleus diameter
	##
	
	# get spb koordinates and (euclidean-)distance d between them
	sk = [(spb_table.getEntry(i,"X"),spb_table.getEntry(i,"Y")) for i in roi_table.getEntry(index,"spb_id")]
	d  = ( (sk[0][0]-sk[1][0])**2  + (sk[0][1]-sk[1][1])**2 )**(0.5)

	# get nuclei diameter D (mean of height and width)
	D = ( nuclei_table.getEntry(roi_table.getEntry(index,"n_id")[0],"width") + 
		  nuclei_table.getEntry(roi_table.getEntry(index,"n_id")[0],"height") )/2

	# evaluate cell
	if d > D: 
		roi_table.setEntry(index,"name","P/M")
	else :    
		roi_table.setEntry(index,"name","G2")
	
def overlay_area(r,r2,rm):
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

def high_whi5(index,nucleus_id,roi_table):
	whi5_diff = (roi_table.getEntry(index,"whi5")/(roi_table.getEntry(index,"area")))-(nuclei_table.getEntry(nucleus_id,"whi5")/(nuclei_table.getEntry(nucleus_id,"area")) )
	if whi5_diff <= -10: return True
	else: return False	

def userDialog(rm):
	WaitForUserDialog("Cellsegmentation was finished", "Please look at your images and make any neccessary changes with the ROI Manager. \n You can delete ROIs or add new ones using Fiji. \n When you press OK a next window will let you change the cell cycle phases.").show()
	
	gd = GenericDialog("Cell Cycle")  
	for r in range(rm.getCount()):
		if len(re.split("\ -\ ",rm.getName(r)))==1: gd.addStringField(("roi#"+str(r+1))," ")
		else:	gd.addStringField(("roi#"+str(r+1)), re.split("\ -\ ",rm.getName(r))[1] )   
	gd.showDialog()  

	if gd.wasCanceled(): 
		print "User canceled dialog!"
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
								"tc","TC","t/c","T/C"
								]
			if name == "" or correct== False: 
				WaitForUserDialog("Cellcylcle names were incomplete or false", "You have not filled in all cellcycle names or no correct ones. \n Please check again and correct your input").show()
				return False
		return True

def getImages():
	# close all active windows that might interfere with plugin
	wins = WindowManager.getIDList()
	if wins != None:
		for w in wins:
			WindowManager.getImage(w).close()
	
	# get images
	o = Opener()
	o.openMultiple()
	images = ImageHolder()
	wins = WindowManager.getIDList()
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

	return images

############################################# main ####################################################

# get images through user Dialog
images = getImages()
if not images.checkStatus():
	print "ERROR: Didn't get all neccessary images.\nPlease Select the following Image files: NIBA, WU, CFP, BF, NG"

# prepate imagej instances (roi manager, particle analyzer, results table)
rm = RoiManager.getInstance() 
if rm != None: rm.reset()
else: rm = RoiManager()
table = ResultsTable()
options = ParticleAnalyzer.ADD_TO_MANAGER | ParticleAnalyzer.EXCLUDE_EDGE_PARTICLES | ParticleAnalyzer.INCLUDE_HOLES
pa = ParticleAnalyzer(options, ParticleAnalyzer.AREA, table, 0, 100000000)

# get initial rois
initROI = getInitialROIs(images.niba, rm, pa)

#get roi measurements such like area, koordinates, grey value, and whi5 intensity
roi_table = My_table(["name","area","nuclei","n_id","spb","spb_id","X","Y","eval","whi5","width","height"])
mean_grey_value = getRoiMeasurements(rm,roi_table,images.niba)

# get nuclei info for each cell
nuclei_table = countNuclei(rm, images.wu, pa, roi_table, images.niba) # update roi_table with nuclei counts
initROI[0].show()

# get spindle pole body info for each cell
spb_table    = countAndMeasureSPB(rm, images.cfp, pa, roi_table, mean_grey_value)


# now evaluate each roi
av = sum(roi_table.getColumn('area'))/rm.getCount() # average area

cells_without_nuclei = roi_table.getIndexByEntry("nuclei",0)[::-1] 
# reversed list because if a roi is deleted no shift in roi indices occurs
print "\n======================= Start Evaluation=======================\n"
print "\n======================= Cells to be Evaluated with no nuclei: ",cells_without_nuclei
for index in cells_without_nuclei:
	roi = rm.getRoi(index)
	roi_area = roi_table.getEntry(index,"area")
	print "\n### Evaluate roi" ,rm.getName(index),"with no nuclei and an area of",roi_area,":"
	# if roi has no nucleus and is smaller than half an average cell (do not depend merly on nuclei analysis)
	# it is either noise or a bud
	#if roi_area <= av/2:
	roi_table.setEntry(index,"eval","yes")
	print "~~~ Evaluated ROI",index,roi_table.getEntry(index,"name")
	rois_to_evaluate = roi_table.getIndexByEntry("eval","no")
				
	bud,mothercell_index = evaluate_noiseOrBud(index,roi,av,roi_area,rois_to_evaluate,rm,roi_table)
			

cells_with_two_spbs = [c for c in roi_table.getIndexByEntry("spb",2)[::-1] if c in roi_table.getIndexByEntry("eval","no")]
print "\n======================= Cells to be Evaluated with twp spbs: ",cells_with_two_spbs
for index in cells_with_two_spbs:
	print "\n### Evaluate roi" ,rm.getName(index),"with two spindle poly bodies"
	if roi_table.getEntry(index,"nuclei")==1:
		evaluate_G2_vs_PM(index,roi_table,rm,nuclei_table,spb_table)
		roi_table.setEntry(index,"eval","yes")
	if roi_table.getEntry(index,"nuclei")>1:
		roi_table.setEntry(index,"name","ANA")
		roi_table.setEntry(index,"eval","yes")
	print "~~~ Evaluated ROI",index,roi_table.getEntry(index,"name")


cells_with_high_intensity_spb = [c for c in roi_table.getIndexByEntry("spb",1)[::-1] if roi_table.getEntry(c,"spb_id")[0] in spb_table.getIndexByEntry("high_intensity","yes") and c in roi_table.getIndexByEntry("eval","no")]
print "\n======================= Cells to be Evaluated with a spb of high intensity: ",cells_with_high_intensity_spb
for index in cells_with_high_intensity_spb:
	roi_table.setEntry(index,"name","Late S")
	roi_table.setEntry(index,"eval","yes")


# if a cell had any near neighbours they would have been seperated by watershed
remaining_cells = [(c,evaluate_watershed(c,rm,roi_table)) for c in roi_table.getIndexByEntry("eval","no") ]
cells_with_too_many_neighbours = [(c,w) for c,w in remaining_cells if len(w)>=2][::-1]
print "\n======================= Cells to be Evaluated with too many neighbours: ",cells_with_too_many_neighbours
for c,w in cells_with_too_many_neighbours:
	
	rm.select(c)
	IJ.run("Enlarge...","enlarge=1")
	rm.runCommand("Update")
	#choose those cell as a partner that has the most overlaying area
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

#print roi_table
#WaitForUserDialog("Cellsegmentation was finished", "Please look at your images and make any neccessary changes with the ROI Manager. \n You can delete ROIs or add new ones using Fiji. \n When you press OK a next window will let you change the cell cycle phases.").show()

remaining_cells = [(c,evaluate_watershed(c,rm,roi_table)) for c in roi_table.getIndexByEntry("eval","no") ]
cells_with_one_neighbour = [(c,w) for c,w in remaining_cells if len(w)==1][::-1]
print "\n======================= Cells to be Evaluated with one neighbour: ",cells_with_one_neighbour

for i in range(len(cells_with_one_neighbour)/2) :
	cells_with_one_neighbour.remove((cells_with_one_neighbour[i][1][0],[cells_with_one_neighbour[i][0]]))
print "\n======================= Cells to be Evaluated with one neighbour: ",cells_with_one_neighbour
# remove doubles (if cell A and B are neighbours they should only be combined once!)
for i in range(len(cells_with_one_neighbour)):
	rm.select(cells_with_one_neighbour[i][0])
	IJ.run("Enlarge...","enlarge=1")
	rm.runCommand("Update")
	combineTwoRois(cells_with_one_neighbour[i][0],cells_with_one_neighbour[i][1][0],roi_table,rm)
	


remaining_cells = [ c for c in roi_table.getIndexByEntry("eval","no") ]
print "\n======================= Cells to be Evaluated (remaining): ",remaining_cells
for index in remaining_cells:
	nucleus_id = roi_table.getEntry(index,"n_id")[0]
	if roi_table.getEntry(index,"nuclei") == 1:
		if high_whi5(index,nucleus_id,roi_table):
			roi_table.setEntry(index,"name","G1")
		else:
			roi_table.setEntry(index,"name","Early S")
	else:	
		if high_whi5(index,nucleus_id,roi_table):
			roi_table.setEntry(index,"name","T/C")
		else:
			roi_table.setEntry(index,"name","ANA")
	roi_table.setEntry(index,"eval","yes")

print "\n======================= Done Evaluation=======================\n"

def renameRois(rm, roi_table):
	rm.deselect()
	for index in range(rm.getCount()):
		rm.select(index)
		rm.runCommand("Rename",(str(index+1)+ " - "+ roi_table.getEntry(index, "name")))
		rm.deselect()

renameRois(rm, roi_table)
print roi_table
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
print rm.getCount()
rm.select(Zniba,0)
rm.deselect()
rm.runCommand("Show None")
rm.runCommand("Save", "/home/svenja/Documents/data2/RoiSet.zip")
#fs = FileSaver(processed_Zniba) 
#fs.saveAsTiff("/home/svenja/Documents/data2/label.tiff")
#fs2 = FileSaver(Zniba) 
#fs2.saveAsTiff("/home/svenja/Documents/data2/wekaInput.tiff")

wins = WindowManager.getIDList()
for w in wins:
	WindowManager.getImage(w).close()

import os.path
import sys
from time import sleep

from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable
from ij import ImagePlus, IJ
from ij.gui import Roi
from ij.io import OpenDialog
from java.awt import Rectangle
from ij.gui import GenericDialog, WaitForUserDialog

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

def openUserDialog():
	gd = GenericDialog("Hough Segmentation")
	gd.enableYesNoCancel("OK","Choose Directory")
	gd.centerDialog(True)
	gd.addMessage(("Please Fill out all Parameter Settings\n"+
                   "according to your image!"))
	gd.addNumericField("MINIMUM Cell Radius (Pixel)",6 ,0)
	gd.addNumericField("MAXIMUM Cell Radius (Pixel)",20,0)
	gd.addNumericField("Gradient Threshold  (Pixel)",50,0)
	gd.addNumericField("Filter Radius       (Pixel)",14,0)
	gd.addNumericField("Multi Radii         (Pixel)",0.5,1)
	gd.showDialog()
	return(gd)

def getUserInput():
	gd = openUserDialog()
	if gd.wasCanceled():
		return False

	# if user did not cancel he must have chosen "Chose directory" or
	# falsly pressed OK (without choosing directory)
	condition = True
	while(condition):
		od = OpenDialog("Choose image")
		imp = IJ.openImage(od.getPath())
		gd = openUserDialog()
		condition = (gd.wasCanceled()==False and gd.wasOKed()==False)

	if gd.wasCanceled():
		return False
	if gd.wasOKed():
		minR = gd.getNextNumber()
		maxR = gd.getNextNumber()
		grdth = gd.getNextNumber()
		filterR = gd.getNextNumber()
		multiR = gd.getNextNumber()
		args = [minR, maxR, grdth, filterR, multiR]
	return([imp, od.getPath(), args])

def callMatlab(path, args):
	sys.stdout.write("Starting Hough Transform..")
	os.system(('matlab -nodisplay -nodesktop -nojvm -r '+
	          '" cd('+"'~/Projects/wolfgang/Scripts/');"+
	          'seg_oneImage('+ "'"+ path + "'," + 
	          str(args[0])+","+str(args[1])+","+ str(args[2])+
	          ","+str(args[3])+'); exit" &'
	          ))
	# wait for matlab to finish before proceeding
	while(os.path.exists("/tmp/example_rois.txt")==False):
		sys.stdout.write(".")
		sys.stdout.flush()
		sleep(1)
	print "Done"

def getRois(rm,roi_table):
	"""
	"""
	fileIn = open("/tmp/example_rois.txt")

	# prepare for measurements
	IJ.run("Set Measurements...", "area")
	rt = ResultsTable.getResultsTable()
	rt.reset()
	rm.reset()
	i = 0

	for line in fileIn.readlines():
		x,y,r = line.split("\t")
		r = float(r)
		x = float(x)
		y = float(y)
		roi = Roi(x-r, y-r, 2*r, 2*r, 200) # 200 is max -> circle no rectangle
		rm.addRoi(roi)
		#i = rm.getRoiIndex(roi) # does not work somehow. returns -1
		rm.select(i)
		rm.runCommand("Measure")
		roi_table.addRow([rm.getName(i)[:4], rt.getValue('Area',i), x, y, r])
		rm.deselect()
		i = i+1

	rt.reset()
	rm.deselect()
	return(True)

# ----------------------------------------------------------------
# Part EdgeDetection: look for Border in between two cells
# ----------------------------------------------------------------
def hasEdgeInBetween(imp_original, rm, roi1, roi2, roi_table):
	imp = imp_original.duplicate()
	IJ.run(imp, "Enhance Contrast","0.8")
	imp.show()
	# edge detectiono via fiji
	# create rectangular roi between the two rois
	roi = getPolygonRoi(roi1, roi2, roi_table)
	edge_table = []

def getPolygonRoi(roi1, roi2, roi_table):
	"""
	Computes (vector theory) a trapezoid between two circles 
	and stores it as a PolygonRoi.
	"""
	r1 = float(roi_table.getEntry(roi1,"r"))
	r2 = float(roi_table.getEntry(roi2,"r"))
	x1 = float(roi_table.getEntry(roi1,"x"))
	y1 = float(roi_table.getEntry(roi1,"y"))
	x2 = float(roi_table.getEntry(roi2,"x"))
	y2 = float(roi_table.getEntry(roi2,"y"))

	vx = x2-x1
	vy = y2-y1
	
	alpha11 =  ((r1**2/(vx**2+vy**2))**(0.5))
	alpha12 = -((r1**2/(vx**2+vy**2))**(0.5))
	alpha21 =  ((r2**2/(vx**2+vy**2))**(0.5))
	alpha22 = -((r2**2/(vx**2+vy**2))**(0.5))
	
	p11_x = x1 + alpha11*(-vx)
	p11_y = y1 + alpha11*(vy)
	p12_x = x1 + alpha12*(-vx)
	p12_y = y1 + alpha12*(vy)
	
	p21_x = x2 + alpha21*(-vx)
	p21_y = y2 + alpha21*(vy)
	p22_x = x2 + alpha22*(-vx)
	p22_y = y2 + alpha22*(vy)
	
	xPoints = [p11_x, p21_x, p22_x, p12_x]
	yPoints = [p11_y, p21_y, p22_y, p12_y]
	roi = PolygonRoi(xPoints, yPoints, Roi.POLYGON)
	return roi

# -- end Part EdgeDetection ---------------------------------------

################################ main ######################################

# build Generic Dialog
userInput = getUserInput()
if (not (userInput == False)):
	
	imp, path, args = userInput

	# process with matlab
	callMatlab(path, args)
	imp.show()
	
	# get circular rois
	rm = RoiManager.getInstance()
	rm.show()
	roi_table = My_table(["name","area","x","y","r"])
	getRois(rm, roi_table)



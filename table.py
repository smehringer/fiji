class My_table(object):
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


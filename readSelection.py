from pymol import cmd

def readSelection (fileName):
	'''
	DESCRIPTION

	reads a csv file containing a column of selection names and a column of
	selection strings and creates selection objects for each
	'''
	line = []
	#read in the selection file
	with open(fileName) as f:
		for line in f:
			sel = line.split(",")
			print sel[0] + " " + sel[1]
			#cmd.select(selection[0], selection[1])

	f.close()

cmd.extend('readSelection',readSelection)





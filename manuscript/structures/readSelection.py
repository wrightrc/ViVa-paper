from pymol import cmd, stored

def readSelection (fileName):
	'''
	DESCRIPTION

	reads a csv file containing a column of selection names and a column of
	selection strings and creates selection objects for each
	'''
	#read in the selection file
	with open(fileName) as f:
		lines = f.readlines()
		for line in lines:
			sel = line.split(",")
			cmd.select(sel[0], sel[1])

	f.close()

cmd.extend('readSelection',readSelection);





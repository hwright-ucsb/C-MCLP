
def parse(filepath):

	configfile = open(filepath, "r")
	flag = True
	dat = []
	while flag:
		l = configfile.readline()
		if (len(l)==0):
			flag = False
		elif(not(l[0]=="#" or l[0]=="\n")):
			dat.append([int(i) for i in l.split()])


	#num sides
	border=[]
	for i in range(0,dat[0[0]]):
		border.append(dat[i])

	return border


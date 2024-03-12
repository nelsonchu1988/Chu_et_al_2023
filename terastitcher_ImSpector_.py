## This python code is adopted from terastitcher and perform stitching alignment of 3D tile images https://abria.github.io/TeraStitcher/

import sys
import os
from xml.dom import minidom
from subprocess import call
import subprocess
from datetime import datetime
import re
import csv

##################################################################################
# change settings below
#*********************************************************************************

##################################################################################
# WORK_MODE = 1 => start within ImSpector and run terastitcher
# WORK_MODE = 2 => start within ImSpector and export parameter in file work_path/imspector_parameters.txt
# WORK_MODE = 3 => start without ImSpector with an absolute path to imspector_parameters.txt as argument and stitch
WORK_MODE = 1

#*********************************************************************************
# work_path: output of terastitcher configuration and stitched images
# select_path_in_dialog: select output path in gui-dialog
work_path='D:\\Stitched' + os.sep + datetime.now().strftime('%Y%m%d_%H%M%S') + os.sep
# work_path='/home/markus/myData/stitched' + os.sep + datetime.now().strftime('%Y%m%d_%H%M%S') + os.sep
select_path_in_dialog = False

#*********************************************************************************
# tool_path: absolute path of cli terastitcher
tool_path='D:\\TeraStitcher\\bin\\terastitcher.exe'
# tool_path='/home/markus/myPrograms/terastitcher_bin/bin/terastitcher'

# set overlap to False if the overlap should be defined manually
readOverlap = True
############# reference_time und reference_channel aus den ImageSeriesViewer-Slidern auslesen? #############
# reference_time and reference_channel: sutable time and channel for stitching.
# displacements will be calculated only once for this time and channel and saved
# terastitcher will use saved displacements etc. for stitching of other channels
reference_point_from_sliders = True

# will be set to True if TKinter was found
foundTKinter = False

class Parameters(object):
	def __init__(self):
		self.axesNames = {}
		self.reference_time = 0
		self.reference_channel = 0
		self.path_to_stack = str()
		self.files = []
		self.stack_cols = 0
		self.stack_rows = 0
		self.stack_slices = 0
		self.width_px = 0
		self.height_px = 0
		self.width_phy = 0.0
		self.height_phy = 0.0
		self.voxel_depth = 0.0
		self.overlap_I = 0
		self.sV = str(100)
		self.sH = str(100)
		self.sD = str(1)

#*********************************************************************************
# predefined axes names
############# die Namen der Treiber koennen unterschiedlich sein... z.B. table mit kleinem t... ##########################
x_axis_name='xyz-Table X'
y_axis_name='xyz-Table Y'
z_axis_name='xyz-Table Z'
xRes_name = "xyz-Table X Resolution"
yRes_name = "xyz-Table Y Resolution"
zRes_name = "xyz-Table Z Resolution"
ZStart_name = "xyz-Table Z Start position"
ZStop_name = "xyz-Table Z Stop position"
XYOverlap_name = "xyz-Table XY Overlap"
##################################################################################

#*********************************************************************************
# special_channel_axis_name: not splitted measurement of Ultra II
special_channel_axis_name='UltraII Filter'
##################################################################################

if WORK_MODE != 3: import lvbt
else:
	try:
		import Tkinter, tkFileDialog # for possible selecting of the output folder
		foundTKinter = True
	except ImportError:
		foundTKinter = False
		print('Tkinter module not found.')
	if len(sys.argv) != 2:
		print("provide an absolute path to imspector_parameters.txt as the only argument")
		raise SystemExit
	if not os.path.exists(sys.argv[1]):
		print('could not create imspector_parameters.txt path. path was ' + sys.argv[1])
		raise SystemExit

global parameters
parameters = Parameters()

def saveParameters(filePath):
	with open(filePath, "wb+") as f:
		writer = csv.writer(f, delimiter="|")
		writer.writerow( ('numAxes', len(parameters.axesNames)) )
		for i in range(len(parameters.axesNames)): writer.writerow( (i, parameters.axesNames[i]) )
		writer.writerow( ('reference_time', parameters.reference_time) )
		writer.writerow( ('reference_channel', parameters.reference_channel) )
		writer.writerow( ('path_to_stack', parameters.path_to_stack) )
		writer.writerow( ('stack_cols', parameters.stack_cols) )
		writer.writerow( ('stack_rows', parameters.stack_rows) )
		writer.writerow( ('stack_slices', parameters.stack_slices) )
		writer.writerow( ('width_px', parameters.width_px) )
		writer.writerow( ('height_px', parameters.height_px) )
		writer.writerow( ('width_phy', parameters.width_phy) )
		writer.writerow( ('height_phy', parameters.height_phy) )
		writer.writerow( ('voxel_depth', parameters.voxel_depth) )
		writer.writerow( ('overlap_I', parameters.overlap_I) )
		writer.writerow( ('numFiles', len(parameters.files)) )
		for i in range(len(parameters.files)): writer.writerow( (i, parameters.files[i]) )

def loadParameters(filePath):
	global parameters
	with open(filePath, "r") as f:
		reader = csv.reader(f, delimiter="|")
		parameters.axesNames = {}
		parameters.files = []
		numAxes = int(next(reader)[1])
		for i in range(numAxes): parameters.axesNames[i] = str(next(reader)[1])
		parameters.reference_time = int(next(reader)[1])
		parameters.reference_channel = int(next(reader)[1])
		parameters.path_to_stack = str(next(reader)[1])
		parameters.stack_cols = int(next(reader)[1])
		parameters.stack_rows = int(next(reader)[1])
		parameters.stack_slices = int(next(reader)[1])
		parameters.width_px = int(next(reader)[1])
		parameters.height_px = int(next(reader)[1])
		parameters.width_phy = float(next(reader)[1])
		parameters.height_phy = float(next(reader)[1])
		parameters.voxel_depth = float(next(reader)[1])
		if readOverlap:
			parameters.overlap_I = float(next(reader)[1])
		else:
			savedOverlap = float(next(reader[1]))
		numFiles = int(next(reader)[1])
		for i in range(numFiles): parameters.files.append(str(next(reader)[1]))

if select_path_in_dialog:
	if WORK_MODE != 3: work_path = lvbt.io.getPath() + os.sep + datetime.now().strftime('%Y%m%d_%H%M%S') + os.sep
	else:
		if foundTKinter:
			root = Tkinter.Tk()
			work_path = tkFileDialog.askdirectory(parent = root, initialdir="/", title = 'Please select a directory') + os.sep + datetime.now().strftime('%Y%m%d_%H%M%S') + os.sep
		print(work_path)
if not os.path.exists(work_path):
	os.makedirs(work_path)
if not os.path.exists(work_path):
	print('could not create output path. path was ' + work_path)
	raise SystemExit

if not os.path.exists(tool_path):
	print('terastitcher is not found. path is ' + tool_path)
	raise SystemExit

def readReferencePoinFromSliders():
	time_axis_names = ['Time Time', 'custom name for time axis']
	channel_axis_names = ['UltraII', 'TriMScope', 'custom name for channel axis', 'UltraII Filter']
	for i in range(int(lvbt.viewer().numSlider())):
		slider = lvbt.viewer.slider(i)
		for j in range(len(time_axis_names)):
			if time_axis_names[j] == slider.name():
				parameters.reference_time = slider.posCurrent()
		for j in range(len(channel_axis_names)):
			if channel_axis_names[j] == slider.name():
				parameters.reference_channel = slider.posCurrent()
	print('reference_time: ' + str(parameters.reference_time))
	print('reference_channel: ' + str(parameters.reference_channel))

##################################################################################
##################################################################################
# begin
##################################################################################
##################################################################################

def readParametersFromImSpector():
	global parameters
	parameters.axesNames = {}
	parameters.files = []
	if reference_point_from_sliders: readReferencePoinFromSliders()
	# get access to the measurement
	measurement = lvbt.measurement("Measurement 1")
	# copy list of all properties [internal_name:value]
	allProps = measurement.getProperties()
	# read axes names
	parameters.axesNames = readoutAxesNames(allProps)
	files = lvbt.viewer().files()
	# read path to series
	parameters.path_to_stack = lvbt.viewer().directory()
	# read all names of files and remove stack path from names
	files = lvbt.viewer().files()
	for i in range(len(files)):
		for j in range(len(files[i])):
			fname = files[i][j]
			fname = fname.replace(parameters.path_to_stack, '')
			parameters.files.append(fname)
	parameters.stack_cols = measurement.getProperty(xRes_name)
	parameters.stack_rows = measurement.getProperty(yRes_name)
	parameters.stack_slices = measurement.getProperty(zRes_name)
	stacks = lvbt.viewer().stacks()
	if len(stacks) <= 0:
		print('invalid stack size (no stacks found)')
		raise SystemExit
	stack = stacks[0]
	res = stack.res
	parameters.width_px = res[0]
	parameters.height_px = res[1]
	phy_len = stack.length
	parameters.width_phy = phy_len[0]
	parameters.height_phy = phy_len[1]
	z_start = measurement.getProperty(ZStart_name)
	z_stop = measurement.getProperty(ZStop_name)
	parameters.voxel_depth = abs(z_stop - z_start) / parameters.stack_slices;
	# overlap in % in ImSpector
	if readOverlap:
		parameters.overlap_I = measurement.getProperty(XYOverlap_name)

def main():
	print('main start')
	if WORK_MODE == 3: loadParameters(sys.argv[1])
	else:
		readParametersFromImSpector()
		saveParameters(work_path + 'imspector_parameters.txt')
		if WORK_MODE == 2:
			print('parameters are saved in ' + work_path + 'imspector_parameters.txt')
			raise SystemExit
	# create ordered list of files
	files_T_C_Y_X_Z = structure_from_properties_and_filenames()
	print('t steps: ' + str(len(files_T_C_Y_X_Z)))
	print('channels: ' + str(len(files_T_C_Y_X_Z[0])))
	rows = len(files_T_C_Y_X_Z[0][0])
	cols = len(files_T_C_Y_X_Z[0][0][0])
	print('rows: ' + str(rows))
	print('cols: ' + str(cols))
	# create subdirectories for each timestamp and channel
	for t in files_T_C_Y_X_Z.keys():
		for c in files_T_C_Y_X_Z[t].keys():
			subpath = make_subpath(t, c)
			os.makedirs(work_path + subpath)
	doc = make_import_xml(files_T_C_Y_X_Z)
	prettyxmlstr = doc.toprettyxml(indent="", newl="\n", encoding="UTF-8")
	refsubpath = make_subpath(parameters.reference_time, parameters.reference_channel)
	with open(work_path + refsubpath + os.sep + 'xml_import.xml', "wb") as f:
		f.write(prettyxmlstr)
	# stitch reference frame and create all configuration files
	stitch_initial(work_path + refsubpath)
	refdoc = minidom.parse(work_path + refsubpath + os.sep + 'xml_merging.xml')
	#print refdoc.toprettyxml(indent="", newl="\n", encoding="UTF-8")
	for t in files_T_C_Y_X_Z.keys():
		for c in files_T_C_Y_X_Z[t].keys():
			if t == parameters.reference_time and c == parameters.reference_channel: continue
			stacks = refdoc.getElementsByTagName('Stack')
			for stack in stacks:
				row = int(stack.attributes['ROW'].value)
				col = int(stack.attributes['COL'].value)
				stack.setAttribute('IMG_REGEX', files_T_C_Y_X_Z[t][c][row][col])
			prettyxmlstr = refdoc.toprettyxml(indent="  ", newl="\n", encoding="UTF-8")
			subpath = make_subpath(t, c)
			path = work_path + subpath
			with open(path + os.sep + 'xml_merging.xml', "wb") as f:
				f.write(prettyxmlstr)
			delete_bad()
			print("merging time " + str(t) + " and channel " + str(c))
			print("calling merging")
			subprocess.call([tool_path, '--merge', '--projin=' + path + 'xml_merging.xml', '--imout_depth=16', '--volout=' + path])
			print("done")
	print('main end')

# finds and reads out property by internal name (unique name)
def findPropertyByName(allProps, name):
	if allProps.has_key(name):
		return [1, allProps.get(name)]
	else:
		return [0, '']

# read out and check axes names from measurement properties
def readoutAxesNames(allProps):
	axes_names = {}
	axes_names[0] = findPropertyByName(allProps, 'FirAxis')
	axes_names[1] = findPropertyByName(allProps, 'SecAxis')
	axes_names[2] = findPropertyByName(allProps, 'ThdAxis')
	axes_names[3] = findPropertyByName(allProps, 'FthAxis')
	if axes_names[2][1] == 'None':
		axes_names[2] = findPropertyByName(allProps, 'ThrAxis')
	emptyAxisFound = False
	clean_axes_names = {}
	for i in range(0, 4):
		if axes_names[i][0] == 0 or axes_names[i][1] == 'None':
			emptyAxisFound = True
		elif emptyAxisFound == True:
			print('invalid or not supported axes configuration')
			raise SystemExit
		else: clean_axes_names[i] = axes_names[i][1]
	return clean_axes_names

def structure_from_properties_and_filenames():
	global parameters
	time_axis = findTimeAxis(parameters.axesNames)
	maybe_have_special_channel_name = findSpecialChannelAxis(parameters.axesNames)
	if maybe_have_special_channel_name:
		print("found special channel axis in properties")
	have_special_channel_name = False
	both_xy_are_available = areBothXYAvailable(parameters.axesNames)
	xy_are_swapped = False
	if both_xy_are_available:
		xy_are_swapped = areXandYswapped(parameters.axesNames)
	else:
		idxX = getIndex(parameters.axesNames, x_axis_name)
		if idxX != -1: xy_are_swapped = True
		else:
			idxY = getIndex(parameters.axesNames, y_axis_name)
			if idxY == -1:
				print('found no x and no y axis')
				raise SystemExit
	z_index = getIndex(parameters.axesNames, z_axis_name)
	files_T_C_Y_X_Z = dict()
	for i in range(len(parameters.files)):
		fname = parameters.files[i]
		if i == 0:
			print('first file: ' + fname)
			if maybe_have_special_channel_name:
				have_special_channel_name = checkSpecialChannelAxisInFileName(fname)
				if have_special_channel_name:
					print("found special channel axis in the first file name")
		t = 0
		if time_axis:
			t = get4DigitIndex(fname, time_axis)
		c = 0
		if have_special_channel_name:
			c = get4DigitIndex(fname, special_channel_axis_name)
		else:
			c = getChannelIndex(fname)
		C_pos = findCPos(fname)
		row = 0
		col = 0
		if both_xy_are_available:
			row, col = get_row_and_column(fname, C_pos, xy_are_swapped)
		else:
			row = get_row_or_column(fname, C_pos)
			if xy_are_swapped:
				col = row
				row = 0
		if not t in files_T_C_Y_X_Z: files_T_C_Y_X_Z[t] = dict()
		if not c in files_T_C_Y_X_Z[t]: files_T_C_Y_X_Z[t][c] = dict()
		if not row in files_T_C_Y_X_Z[t][c]: files_T_C_Y_X_Z[t][c][row] = dict()
		if not col in files_T_C_Y_X_Z[t][c][row]:
			fname = escape_string(fname)
			z = 0
			if z_index > -1:
				z = get4DigitIndex(fname, z_axis_name)
				_Z_pos = fname.find(z_axis_name)
				digitPos = _Z_pos + len(z_axis_name)
				fname = fname[:digitPos] + '[0-9]+' + fname[digitPos + 4:]
			files_T_C_Y_X_Z[t][c][row][col] = fname
	#endfor i
	return files_T_C_Y_X_Z

def findTimeAxis(axes_names):
	found_time = False
	time_axis = ""
	for i in range(len(axes_names)):
		name = axes_names[i]
		if name.find(' Time') == (len(name) - len(' Time')):
			if found_time == True:
				print('can not handle multiple time axes')
				raise SystemExit
			else:
				found_time = True
				time_axis = name
	return time_axis

def get4DigitIndex(fname, axis):
	pos = fname.find(axis)
	if pos == -1:
		print('axis is not found in file name.')
		raise SystemExit
	digitPos = pos + len(axis)
	numDigits = 0
	while (digitPos + numDigits < len(fname)) and fname[digitPos + numDigits].isdigit():
		numDigits = numDigits + 1
	if (digitPos + numDigits < len(fname)) and (numDigits == 4):
		return int(fname[digitPos : digitPos + 4])
	else:
		return int(0)

def findSpecialChannelAxis(axes_names):
	return findAxisByName(axes_names, special_channel_axis_name)

def checkSpecialChannelAxisInFileName(fname):
	found = False
	pos = fname.find(special_channel_axis_name)
	if pos > -1 and pos < len(fname) - 4:
		found = True
	return found

# search first occurance of _Cxx, x is a digit
def getChannelIndex(fname):
	pos = fname.find('_C')
	if pos != -1:
		numDigits = 0
		while (pos + 2 + numDigits < len(fname)) and fname[pos + 2 + numDigits].isdigit():
			numDigits = numDigits + 1
		if ((numDigits == 2) and ((pos + 2 + numDigits) < len(fname))):
			c = fname[pos + 2 + numDigits]
			if (c == '_') or (c == '.'):
				return int(fname[pos + 2:pos + 2 + numDigits])
	print('can not find channel index')
	raise SystemExit

# check if measurement have both x and y axes
def areBothXYAvailable(axes_names):
	cnt = 0
	if findAxisByName(axes_names, x_axis_name) == True: cnt = cnt + 1
	if findAxisByName(axes_names, y_axis_name) == True: cnt = cnt + 1
	return cnt == 2

def findAxisByName(axes_names, axis_name):
	found = False
	for i in range(len(axes_names)):
		name = axes_names[i]
		pos = name.find(axis_name)
		if pos != -1 and pos == (len(name) - len(axis_name)):
			if found == True:
				print('can not handle multiple axes with identical names')
				raise SystemExit
			else:
				found = True
	return found

def areXandYswapped(axes_names):
	idxX = getIndex(axes_names, x_axis_name)
	idxY = getIndex(axes_names, y_axis_name)
	if idxX == -1 or idxY == -1:
		print('predefined x- or y-axis is not found')
		raise SystemExit
	return idxX > idxY

def getIndex(axes_names, axis_name):
	idx = -1
	for i in range(len(axes_names)):
		name = axes_names[i]
		if name.find(axis_name) == (len(name) - len(axis_name)):
			idx = i
			break
	return idx

def findCPos(fname):
	C_pos = fname.find('_C')
	if(C_pos == -1):
		print('invalid file name, _C is not found')
		raise SystemExit
	return C_pos

# cut row and col from string of the form "cccccccccc[aa x bb]_Ccccccc", c - character, aa - row, bb - col
def get_row_and_column(fname, _C_pos, xy_are_swapped):
	if not xy_are_swapped:
		xPos = _C_pos - 5
		row = int(fname[xPos - 3 : xPos - 1])
		col = int(fname[xPos + 2 : xPos + 4])
	else:
		xPos = _C_pos - 5
		col = int(fname[xPos - 3 : xPos - 1])
		row = int(fname[xPos + 2 : xPos + 4])		
	return row, col

# cut index from string of the form "cccccccccc[aa]_Ccccccc", c - character, aa - row or col
def get_row_or_column(fname, _C_pos):
	row = int(fname[_C_pos - 3 : _C_pos - 1])
	return row

# escape special characters
def escape_string(str):
	str = str.replace('\\', '\\\\')
	for sc in '[].':
		str = str.replace(sc, '\\'+sc)
	return str

# tXXX\cXXX\
def make_subpath(t, c):
	return 't' + str(t).zfill(3) + os.sep + 'c' + str(c).zfill(3) + os.sep

def make_import_xml(files_T_C_Y_X_Z):
	doc = minidom.Document();
	dt = minidom.getDOMImplementation('').createDocumentType('TeraStitcher', '', "TeraStitcher.DTD")
	doc.insertBefore(dt, doc.documentElement)
	fill_2d_series(doc, files_T_C_Y_X_Z)
	return doc

# fill xml structure
def fill_2d_series(doc, files_T_C_Y_X_Z):
	global parameters
	mainNode = doc.createElement('TeraStitcher')
	mainNode.setAttribute('volume_format', 'TiledXY|2Dseries')
	doc.appendChild(mainNode)
	node = doc.createElement('stacks_dir')
	node.setAttribute('value', parameters.path_to_stack)
	mainNode.appendChild(node)
	node = doc.createElement('dimensions')
	node.setAttribute('stack_columns', str(parameters.stack_cols))
	node.setAttribute('stack_rows', str(parameters.stack_rows))
	node.setAttribute('stack_slices', str(parameters.stack_slices))
	mainNode.appendChild(node)
	voxel_width = parameters.width_phy / parameters.width_px
	voxel_height = parameters.height_phy / parameters.height_px
	node = doc.createElement('voxel_dims')
	node.setAttribute('D', str(parameters.voxel_depth))
	node.setAttribute('H', str(voxel_width))
	node.setAttribute('V', str(voxel_height))
	mainNode.appendChild(node)
	# the origin of the stage. Keep it at 0
	origin_x = 0
	origin_y = 0
	origin_z = 0
	node = doc.createElement('origin')
	node.setAttribute('D', str(origin_z))
	node.setAttribute('H', str(origin_x))
	node.setAttribute('V', str(origin_y))
	mainNode.appendChild(node)
	# factor that defines the position
	overlap_factor = 1.0 - (parameters.overlap_I / 100.0)
	# offset of between two stacks
	mech_displ_x = parameters.width_phy * overlap_factor
	mech_displ_y = parameters.height_phy * overlap_factor
	node = doc.createElement('mechanical_displacements')
	node.setAttribute('H', str(mech_displ_x))
	node.setAttribute('V', str(mech_displ_y))
	mainNode.appendChild(node)
	stacksNode = doc.createElement('STACKS')
	for row in files_T_C_Y_X_Z[parameters.reference_time][parameters.reference_channel].keys():
		for col in files_T_C_Y_X_Z[parameters.reference_time][parameters.reference_channel][row].keys():
			fname = files_T_C_Y_X_Z[parameters.reference_time][parameters.reference_channel][row][col]
			stackNode = doc.createElement('Stack')
			stackNode.setAttribute('ROW', str(row))
			stackNode.setAttribute('COL', str(col))
			abs_x = mech_displ_x * col
			abs_y = mech_displ_y * row
			abs_z = 0
			stackNode.setAttribute('ABS_D', str(abs_z))
			stackNode.setAttribute('ABS_H', str(abs_x))
			stackNode.setAttribute('ABS_V', str(abs_y))
			stackNode.setAttribute('DIR_NAME', '')
			stackNode.setAttribute('IMG_REGEX', fname)
			stackNode.setAttribute('STITCHABLE', 'no')
			stackNode.setAttribute('Z_RANGES', '[0,' + str(parameters.stack_slices) + ')')
			stackNode.appendChild(doc.createElement('NORTH_displacements'))
			stackNode.appendChild(doc.createElement('EAST_displacements'))
			stackNode.appendChild(doc.createElement('SOUTH_displacements'))
			stackNode.appendChild(doc.createElement('WEST_displacements'))
			stacksNode.appendChild(stackNode)
		# endfor col
	# endfor row
	mainNode.appendChild(stacksNode) 


# stitch reference frame
def stitch_initial(path):
	global parameters
	print("path: ", path)
	oH = str(int((parameters.overlap_I / 100.0) * parameters.width_px))
	oV = str(int((parameters.overlap_I / 100.0) * parameters.height_px))
	print('width_px: ' + str(parameters.width_px) + ', height_px: ' + str(parameters.height_px))
	print('overlap_I: ' + str(parameters.overlap_I))
	print('oH: ' + str(oH) + ', oV: ' + str(oV))
	print("processing time " + str(parameters.reference_time) + " and channel " + str(parameters.reference_channel))
	delete_bad()
	print("calling displcompute")
	print(path,oV,oH,parameters.sV,parameters.sH,parameters.sD)
	subprocess.call([tool_path, '--displcompute', '--rescan', '--projin=' + path + 'xml_import.xml', '--projout=' + path + 'xml_displcomp.xml', '--oV=' + oV, '--oH=' + oH, '--sV=' + parameters.sV, '--sH=' + parameters.sH, '--sD=' + parameters.sD])
	print("done")
	if not os.path.exists(path + 'xml_displcomp.xml'):
		print('xml_displcomp.xml' + 'does not exists!')
		raise SystemExit
	delete_bad()
	print("calling displproj")
	subprocess.call([tool_path, '--displproj', '--projin=' + path + 'xml_displcomp.xml', '--projout=' + path + 'xml_displproj.xml'])
	print("done")
	if not os.path.exists(path + 'xml_displproj.xml'):
		print('xml_displproj' + 'does not exists!')
		raise SystemExit
	delete_bad()
	print("calling displthres")
	subprocess.call([tool_path, '--displthres', '--threshold=0', '--projin=' + path + 'xml_displproj.xml', '--projout=' + path + 'xml_displthres.xml'])
	print("done")
	if not os.path.exists(path + 'xml_displthres.xml'):
		print('xml_displthres' + 'does not exists!')
		raise SystemExit
	delete_bad()
	print("calling placetiles")
	subprocess.call([tool_path, '--placetiles', '--projin=' + path + 'xml_displthres.xml', '--projout=' + path + 'xml_merging.xml'])
	print("done")
	if not os.path.exists(path + 'xml_merging.xml'):
		print('xml_merging' + 'does not exists!')
		raise SystemExit
	delete_bad()
	print("calling merging")
	subprocess.call([tool_path, '--merge', '--projin=' + path + 'xml_merging.xml', '--imout_depth=16', '--volout=' + path])
	print("done")

# delete bad files
def delete_bad():
	bad_file_path = parameters.path_to_stack + os.sep + 'mdata.bin'
	############# hier auch noch null.xml loeschen? #############
	if os.path.exists(bad_file_path):
		os.remove(bad_file_path)

if __name__=="__main__":
	main()

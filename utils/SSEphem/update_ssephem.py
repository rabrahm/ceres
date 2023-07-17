import os
import numpy as np

def LeapSecUpdate():
	#os.system('wget http://maia.usno.navy.mil/ser7/leapsec.dat')
	os.system('wget https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat -O leapsec.dat')
	if os.access('leapsec.tab', os.F_OK):
		os.system('mv leapsec.tab leapsec_old.tab')
	if True:
		f = open('leapsec.dat','r')
		lines = f.readlines()
		fo = open('leapsec.tab','w')
		for line in lines:
			if not '#' in line:
				cos = line.split()
				#if cos[1] == 'JAN':
				#	date = cos[0]+' January 1 '
				#else:
				#	date = cos[0]+' July 1    '
				#jd = float(cos[4])
				if cos[2] == '1':
					date = cos[0]+' January 1 '
				else:
					date = cos[0]+' July 1    '
				mjd = float(cos[0])
				leap = float(cos[4])
				leap = int(np.around(leap))
				#nline = '    '+str(int(np.around(jd-2400000.5)))+'       '+date+'TAI-UTC = '+str(leap)+'.0\n'
				nline = '    '+str(int(np.around(mjd))) + '       '+date+'TAI-UTC = '+str(leap)+'.0\n'
				fo.write(nline)
		f.close()
		fo.close()
		os.system('rm leapsec.dat')
	else:#except:
		print 'No luck...'
		os.system('rm leapsec.dat')
		os.system('mv leapsec_old.tab leapsec.tab')

def SSEphemDownload():
	if os.access('DEc403',os.F_OK):
		os.system('mv DEc403 DEc403_old')
	if os.access('ascp2000.403', os.F_OK):
		os.system('mv ascp2000.403 ascp2000_old.403')
	os.system('wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de403/ascp2000.403')
	if os.access('header.403', os.F_OK):
		os.system('mv header.403 header_old.403')
	os.system('wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de403/header.403')
	if os.access('ascp2000.403', os.F_OK) == False or os.access('header.403', os.F_OK) == False:
		print 'No Luck...'
		os.system('mv ascp2000_old.403 ascp2000.403')
		os.system('mv header_old.403 header.403')
	else:
		os.system('cat header.403 ascp2000.403 > asc_cat.403')
		os.system('./asc2bin asc_cat.403 2451544 2484394')
	if os.access('DEc403',os.F_OK) == False:
		print 'No Luck...'
		if os.access('DEc403_old',os.F_OK):
			os.system('mv DEc403_old DEc403')

def IersUpdate():
	if os.access('finals2000A.data',os.F_OK):
		os.system('mv finals2000A.data finals2000A_old.data')
	#os.system('wget http://maia.usno.navy.mil/ser7/finals2000A.data')
	os.system('wget --no-proxy https://datacenter.iers.org/data/10/finals2000A.data')
	if os.access('finals2000A.data',os.F_OK) == False:
		print 'one'
		print 'No Luck...'
		os.system('mv finals2000A_old.data finals2000A.data')
		pass
	else:
		if os.access('iers.tab',os.F_OK):
			os.system('mv iers.tab iers_old.tab')
		try:
			output = open('iers.tab','w')
			finaldata = open('finals2000A.data','r')
			for line in finaldata:
				mj = line[7:15]
				if len(line.split()) > 5:
					c1 = line[18:27]
					c2 = line[37:46]
					c3 = line[58:68]
					l  = ' '+mj+' '+c1+' '+c2+' '+c3+' '+'\n'
					output.write(l)
			finaldata.close()
			output.close()
		except:
			print 'two'
			print 'No Luck...'
			if os.access('iers_old.tab',os.F_OK):
				os.system('mv iers_old.tab iers.tab')
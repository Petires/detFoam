# Windows, Cantera
from cantera import *
from SDToolbox import *
thisfilename = os.path.abspath('lookupY_05.py')
write = 1
directory = '/home/piotr/'
writetofilename = 'cTable_fpT_07'
logfilename = directory + writetofilename + '.log' # full logfile name

fn = directory + writetofilename + '.csv' #full filename

h2=[0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5, 0.52, 0.54, 0.56, 0.58, 0.6]
nh2=len(h2) # H2 mole fraction 
p=[0.1*1e5, 0.9*1e5, 1.1*1e5, 2*1e5, 5*1e5] + [10*1e5, 20*1e5, 30*1e5, 40*1e5] + [60*1e5, 80*1e5, 100*1e5, 125*1e5, 150*1e5]
nP =len(p)
T=[250.0, 270.0, 290.0, 310.0, 330.0, 350.0, 370.0, 390.0, 410.0, 430.0, 450.0, 470.0, 490.0, 510.0, 530.0, 550.0, 570.0, 590.0, 610.0, 630.0, 650.0, 670.0, 690.0, 710.0, 730.0, 750.0, 770.0, 790.0, 810.0, 830.0, 850.0, 870.0, 890.0, 910.0, 930.0, 950.0, 970.0, 990.0, 1010.0, 1030.0, 1050.0, 1070.0, 1090.0, 1110.0, 1130.0, 1150.0, 1170.0, 1190.0, 1210.0, 1230.0, 1250.0, 1270.0, 1290.0, 1310.0, 1330.0, 1350.0, 1370.0, 1390.0, 1410.0, 1430.0, 1450.0, 1470.0, 1490.0, 1510.0, 1530.0, 1550.0, 1570.0, 1590.0]
nT=len(T)
o2=list(h2)
n2=list(h2)

fH = [0.0]*len(h2)

lenght = nh2*nP*nT
print(lenght)
reply = raw_input('Continue? Y/N [Y]: ')
if reply=='N':
    sys.exit('program stopped by user')

mech = 'OConaire.cti'
gas0 = Solution(mech) 
nsp =  gas0.n_species
ynames=gas0.species_names

# find fuel, nitrogen, and oxygen indices
fuel = raw_input('Fuel\n')
ih2 = gas0.species_index(fuel)
io2  = gas0.species_index('O2')
in2  = gas0.species_index('N2')
#ih2o = gas0.species_index('H2O')
#ioh  = gas0.species_index('OH')
#ih  = gas0.species_index('H')

x = [0]*nsp
Yout = np.zeros((nh2,nP,nT,nsp))*np.nan
Tout = np.zeros((nh2,nP,nT))*np.nan
# log file:
if write:
       logid = open(logfilename, 'w')
       #fid.write(logid, 'Induction times for constant-volume explosion\n');
       logid.write('burned gas composition\n')
       logid.write('fH, p, Tu')
       for i in range (0,len(ynames)):
         logid.write(', ')
         logid.write(ynames[i])
       logid.write('\n')


for i in range (0,nh2):
    for j in range (0,nP):
        for k in range (0,nT):
               print([i, j, k])
               x=[0.0]*nsp
               x[ih2]=h2[i]
               o2[i]=0.21*(1-h2[i])
	       x[io2]=o2[i]
               n2[i]=1-h2[i]-o2[i]
	       x[in2]=n2[i]
               gas0.TPX = T[k], p[j], x
               mtemp=gas0.Y
	       fH[i]=mtemp[ih2]
	       gas0.equilibrate('HP')
               mtemp=gas0.Y
	       Yout[i][j][k]=mtemp
               Tout[i][j][k]=gas0.T
	       if write:
                  logid.write('%1.5f, %1.3e, %4.1f' %(fH[i], p[j], T[k]))
                  for l in range (0,len(mtemp)):
                      logid.write(' %1.5e' %mtemp[l])
                  logid.write('\n' %mtemp[l])
if write:
       logid.close # close log file
print('writing file')
#disp(fn);

#d = date;
	
fid = open(fn, 'w')
#fid.write('Shock and Detonation Toolbox\n');
fid.write('/*--------------------------------*- C++ -*---------------------------------*\\\n')
fid.write('| =========                |                                                 |\n')
fid.write('| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n')
fid.write('|  \\    /   O peration     | Version:  1.7.1                                 |\n')
fid.write('|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |\n')
fid.write('|    \\/     M anipulation  |                                                 |\n')
fid.write('\\*--------------------------------------------------------------------------*/\n')

fid.write('\nFoamFile\n');
fid.write('{\n');
fid.write('version     2.0;\n')
fid.write('format      ascii;\n')
fid.write('class       dictionary;\n')
fid.write('location    "constant";\n')
fid.write('object      %s;\n' % writetofilename)
fid.write('}\n');
fid.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n')



fid.write('\n// created with %s' % thisfilename)
fid.write('\n// %s mechanism '% mech[0:len(mech)-4])

fid.write('\n\noptions\n');
fid.write('{\n');
fid.write('\tlogarithmic false;\n')
fid.write('}\n');

fid.write('\nfields\n')
fid.write('3\n')
fid.write('(\n')

fid.write('\n')
fid.write(' {\n')
fid.write('     name            fH;\n');
fid.write('     min             %1.5f;\n'% min(fH))
fid.write('     max             %1.5f;\n'% max(fH))
fid.write('     N               %d;\n'% int(len(fH)-1))
fid.write(' }\n')
#fid.write('\n');

fid.write('\n')
fid.write(' {\n')
fid.write('     name            p;\n');
fid.write('     min             %1.3e;\n'% min(p))
fid.write('     max             %1.3e;\n'% max(p))
fid.write('     N               %d;\n'% int(len(p)-1))
fid.write(' }\n')
#fid.write('\n');

fid.write('\n')
fid.write(' {\n')
fid.write('     name            Tu;\n');
fid.write('     min             %4.1f;\n'% min(T))
fid.write('     max             %4.1f;\n'% max(T))
fid.write('     N               %d;\n'% int(len(T)-1))
fid.write(' }\n')
#fid.write('\n');
fid.write(');\n\n');


fid.write('\noutput\n')
fid.write('%d\n'% nsp)
fid.write('(\n')
for i in range(0,nsp):
    fid.write(' {    name \t %s; }\n'% ynames[i])
    # fid.write(' {    name \t O2; }');
fid.write(')\n')
fid.write(';\n\n')

fid.write('values\n')
fid.write('%d\n'% (3+nsp))
fid.write('(\n')

fid.write('\n%d  // fH\n( '% lenght)

for i in range (0,nh2):    
    for j in range (0,nP):        
        for k in range (0,nT):
		fid.write('\n%1.5f '% fH[i])
fid.write('\n)\n');

fid.write('\n%d  // p\n( '% lenght)
for i in range (0,nh2):    
    for j in range (0,nP):        
        for k in range (0,nT):
		fid.write('\n%1.3e '% p[j])
fid.write('\n)\n');

fid.write('\n%d  // Tu\n( '% lenght)
for i in range (0,nh2):    
    for j in range (0,nP):        
        for k in range (0,nT):
		fid.write('\n%4.1f '% T[k])
fid.write('\n)\n');


for spec in range (0,nsp):
	fid.write('\n%d //%s\n( '%(lenght, ynames[spec]))
	for i in range (0,nh2):    
	    for j in range (0,nP):        
		for k in range (0,nT):
		    fid.write('\n%1.4e '% Yout[i][j][k][spec])
	fid.write('\n)\n')
fid.write('\n)\n')
fid.write(';\n')

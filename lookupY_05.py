# Windows, Cantera
from cantera import *
from SDToolbox import *
thisfilename = os.path.abspath('lookupY_05.py')
write = 1
writetofilename = 'cTable_fpT_07'
logfilename = writetofilename + '.log' # full logfile name

fn = writetofilename + '.csv' #full filename

h2=np.append(np.arange(0,31,1),np.arange(32,61,2))
nh2=len(h2) # H2 mole fraction 
p=np.append(np.append([0.1, 0.9, 1.1, 2, 5],np.arange(10,40,10)),[60, 80, 100, 125, 150])*1e5
nP =len(p)
T=np.arange(250,1601,20)
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
ih2 = gas0.species_index('H2')
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
                  logid.write('\n')
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

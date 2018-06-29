# Windows, Cantera
# tabulation of log10(tIgn)
from cantera import *
from SDToolbox import *
thisfilename = os.path.abspath('lookupIGN_06.py')
write = 1
#directory = raw_input('Path to result directory\n')
directory = '/home/piotr/'
writetofilename = 'LOG_tignTable_fpT_09'
logfilename = directory + writetofilename + '.log' # full logfile name

fn = directory + writetofilename + '.csv' #full filename

h2=[1*1e-2, 5*1e-2, 10*1e-2, 20*1e-2, 30*1e-2, 40*1e-2, 50*1e-2, 60*1e-2] 
nh2=len(h2); # H2 mole fraction 
fH = [0.0]*len(h2)
  
p=[0.1*1e5, 0.9*1e5, 1.1*1e5, 2*1e5, 5*1e5, 10*1e5, 20*1e5, 30*1e5, 50*1e5, 75*1e5, 100*1e5]
nP =len(p)


T=np.arange(800,2001,20)
nT=len(T)
o2=list(h2)
n2=list(h2)
lenght = nh2*nP*nT
print(lenght)

#mech = raw_input('Reaction mechanism name\n')
mech = 'OConaire.cti'
gas0 = Solution(mech) 
nsp =  gas0.n_species
# find fuel, nitrogen, and oxygen indices
#fuel = raw_input('Fuel\n')
fuel = 'H2'
ih2 = gas0.species_index(fuel)
io2  = gas0.species_index('O2')
in2  = gas0.species_index('N2')
#ih2o = gas0.species_index('H2O')
#ioh  = gas0.species_index('OH')
#ih  = gas0.species_index('H')

x = [0.0]*nsp
tIgn=np.zeros((nh2,nP,nT))*np.nan

fig_num = 2 #0=no plots
while fig_num == 2:
    check = raw_input('Would you like to create plots?(y/n)\n')
    if check == 'y':
        fig_num = 1
    if check == 'n':
        fig_num = 0

# log file:
if(write):
       logid = open(logfilename, 'w')
       logid.write('Induction times\n')
       logid.write('fH, p, T, tign, log10(tign)\n')


dtdT = 0 # maximum error dt_ign/dT
for i in range (0,nh2):
    for j in range (0,nP):
        for k in range (0,nT):
               print([i, j, k])
               x=[0.0]*nsp
               x[ih2]=h2[i]
               o2[i]=0.21*(1.0-h2[i])
               x[io2]=o2[i]
               n2[i]=1.0-h2[i]-o2[i]
               x[in2]=n2[i]               
               gas0.TPX = T[k], p[j], x
               mtemp=gas0.Y
               fH[i]=mtemp[ih2]
               if not fig_num == 0:
                   fig_num=i*nP*nT+j*nT+k+1
               out = explosion(gas0,fig_num) # constant volume explosion
               #[out] = explosionHP(gas0,fig_num); # constant pressure explosion
               tIgn[i][j][k] = out.ind_time
               if(tIgn[i][j][k]>0.9):
                   tIgn[i][j][k]=1e1  
               print(log10(tIgn[i][j][k]))
               if(write):
                  logid.write('%1.5f, %1.3e, %4.1f, %1.5e, %1.5f\n'% (fH[i], p[j], T[k], tIgn[i][j][k], log10(tIgn[i][j][k])));
               
               if(tIgn[i][j][k]<0.1):
                   if(k>0):
                      dtdT_temp = tIgn[i][j][k] - tIgn[i][j][k-1]/(T[k]-T[k-1])
                      if(dtdT_temp > dtdT):
                          dtdT = dtdT_temp
                          dtdT_info = [i,j,k]


if(write):
       logid.close # close log file
print('writing file')
print(fn)

	
fid = open(fn, 'w')
fid.write('/*--------------------------------*- C++ -*---------------------------------*\\\n')
fid.write('| =========                |                                                 |\n')
fid.write('| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n')
fid.write('|  \\    /   O peration     | Version:  1.7.1                                 |\n')
fid.write('|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |\n')
fid.write('|    \\/     M anipulation  |                                                 |\n')
fid.write('\\*--------------------------------------------------------------------------*/\n')

fid.write('\nFoamFile\n')
fid.write('{\n');
fid.write('version     2.0;\n')
fid.write('format      ascii;\n')
fid.write('class       dictionary;\n')
fid.write('location    "constant";\n')
fid.write('object      %s;\n'%writetofilename)
fid.write('}\n');
fid.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n')
fid.write('\n// created with %s'%thisfilename)
fid.write('\n// %s mechanism '% mech[0:len(mech)-4])


fid.write('\n\noptions\n')
fid.write('{\n')
fid.write('    logarithmic true;\n')
fid.write('}\n')

fid.write('\nfields\n')
fid.write('3\n')
fid.write('(\n')

fid.write('\n');
fid.write(' {\n')
fid.write('     name            fH;\n')
fid.write('     min             %1.5f; // xH2=%1.3f\n'%(min(fH),min(h2)))
fid.write('     max             %1.5f; // xH2=%1.3f\n'%(max(fH),max(h2)))
fid.write('     N               %d;\n'%(len(fH)-1))
fid.write(' }\n')
#fid.write('\n')

fid.write('\n')
fid.write(' {\n')
fid.write('     name            p;\n')
fid.write('     min             %1.3e;\n'%min(p))
fid.write('     max             %1.3e;\n'%max(p))
fid.write('     N               %d;\n'%(len(p)-1))
fid.write(' }\n')
#fid.write('\n')

fid.write('\n')
fid.write(' {\n')
fid.write('     name            T;\n')
fid.write('     min             %4.1f;\n'%min(T))
fid.write('     max             %4.1f;\n'%max(T))
fid.write('     N               %d;\n'%(len(T)-1))
fid.write(' }\n')
#fid.write('\n')
fid.write(');\n\n')


fid.write('\noutput\n')
fid.write('%d\n' % 1)
fid.write('(\n')
fid.write(' {    name \t tIgn; }')
fid.write('\n)\n')
fid.write(';\n\n')

fid.write('values\n')
fid.write('4\n')
fid.write('(\n')

fid.write('\n%d  // fH\n( '%lenght)

for i in range (0,nh2):
    for j in range (0,nP):
        for k in range (0,nT):
            fid.write('\n%1.5f '%fH[i])
fid.write('\n)\n');

fid.write('\n%d  // p\n( '%lenght)
for i in range (0,nh2):
    for j in range (0,nP):
        for k in range (0,nT):
            fid.write('\n%1.3e '%p[j])
fid.write('\n)\n');

fid.write('\n%d  // T\n( '%lenght)
for i in range (0,nh2):
    for j in range (0,nP):
        for k in range (0,nT):
            fid.write('\n%4.1f '%T[k])
fid.write('\n)\n');


fid.write('\n%d // log10(tIgn) \n( '%lenght)
for i in range (0,nh2):
    for j in range (0,nP):
        for k in range (0,nT):
            fid.write('\n%1.5f '%log10(tIgn[i][j][k]))
fid.write('\n)\n')

fid.write('\n)\n')
fid.write(';\n')
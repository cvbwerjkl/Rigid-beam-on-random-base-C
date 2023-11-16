## used for comparison C and Octave realisations
## octave packege need to be installed

import os

os.system('make clean')
os.system('make')
os.system('./bin/my_prog -t')
os.system('octave tests.m')
os.system('tail -n +6 solution_Octave.txt > solution_Oct_clean.txt')
os.system('rm solution_Octave.txt')

octave = []
c_gsl = []

with open('solution_Oct_clean.txt') as file:
    while line := file.readline():
        if line != '\n':
            octave.append(float(line.rstrip()))

with open('testsolutionVector.txt') as file:
    while line := file.readline():
        if line != '\n':
            c_gsl.append(float(line.rstrip()))

if len(octave) != len(c_gsl):
    print("Error, solutions are not equal\n")
    exit()

for i in range(len(octave)):
    if abs(octave[i] - c_gsl[i]) > 0.1:
        print("Error, solutions are not equal\n")
        break
else:
    print("Ok,solutions are equal\n") 

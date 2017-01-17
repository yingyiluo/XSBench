from optparse import OptionParser
import os
import re
def get_energy(filepath):
    file = open(filepath)
    lines = file.readlines()[1:]
    for line in lines:
	if re.match('RAPL_CPU1_ENERGY_J', line):
            e = float(line[19:])       
    print e
    return e

def main():
    parser = OptionParser()
    parser.add_option('-i', '--ifile', dest="inputfile", help="give path where the data files locate", metavar="FILE/FOLDER")
#    dutero_avg = get_avg("dutero-arria10-idle.txt", 2)
    opts, args = parser.parse_args()
    if not opts.inputfile:   # if filename is not given
        parser.error('Inputfile not given')
#    ifile = os.listdir(opts.inputfile)
    
    lookups = 15000000; 
    
    e = get_energy(opts.inputfile)

    print "kernel energy on average is "+ str(e)+" J"
    efficiency = (e*1e6)/(15000000)
    print "energy efficiency is "+ str(efficiency) + " uJ/lookup"
if __name__ == "__main__":
    main()

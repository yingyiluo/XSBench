from optparse import OptionParser
import os
import re
def get_avg(filepath, pos):
    file = open(filepath)
    lines = file.readlines()[1:]

    counter = 0
    sum = 0

    for line in lines:
        sum = sum + float(line.split(" ")[pos])
        counter = counter + 1

    return (sum/counter)

def get_time(filepath):
    file = open(filepath)
    lines = file.readlines()[1:]
    for line in lines:
	if re.match('OCL_KERNEL_SEC', line):
            time = float(line[15:])       
    #print time
    return time

def main():
    parser = OptionParser()
    parser.add_option('-i', '--ifile', dest="inputfile", help="give path where the data files locate", metavar="FILE/FOLDER")
#    dutero_avg = get_avg("dutero-arria10-idle.txt", 2)
    opts, args = parser.parse_args()
    if not opts.inputfile:   # if filename is not given
        parser.error('Inputfile not given')
#    ifile = os.listdir(opts.inputfile)
    
    lookups = 15000000;
    card_avg = get_avg(opts.inputfile, 2)
    dutero_avg = 128.24
    card_idle = 29
    print opts.inputfile
    card_kr_avg = card_idle - dutero_avg + card_avg
    print "average power is " + str(card_kr_avg)
    
    time = get_time('../xsperf/'+opts.inputfile)

    card_kr_engy = card_kr_avg * time
    
    print "kernel energy on average is "+ str(card_kr_engy)+" J"
    efficiency = (card_kr_engy*1e6)/(15000000)
    print "energy efficiency is "+ str(efficiency) + " uJ/lookup"
if __name__ == "__main__":
    main()

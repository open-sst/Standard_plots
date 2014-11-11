import numpy as np
import matplotlib.pyplot as plt

class ts:

    def __init__(self,years,months,data):
        self.years = years
        self.months = months
        self.data = data
        self.make_time_axis()

    def make_time_axis(self):
        self.time_axis = self.years[:]
        for i in range(0,len(self.time_axis)):
            self.time_axis[i] += (self.months[i]-1.)/12.

    def plot_ts(self):
        plt.plot(self.time_axis,self.data)
        plt.plot(self.time_axis,np.zeros(len(self.time_axis)),color="Black")
        plt.plot([2014,2014],[-5,5],color="Red")
        plt.show()

    def rebaseline(self,year1,year2):
        print

def slices(s, *args):
    position = 0
    splits = []
    for length in args:
        splits.append(s[position:position + length])
        position += length
    return splits
    
def read_soi():
    f = open('Data/soi', 'r')

    year = []
    month = []
    data = []

    for i in range(1,5):
        f.readline()

    readon = 1

    for line in f:
        if readon == 1:
            line = line.strip()
            line = slices(line,4,6,6,6,6,6,6,6,6,6,6,6,6)

            for i in range(1,13):
                if float(line[i]) != -999.9:
                    year.append(float(line[0]))
                    month.append(float(i))
                    data.append(float(line[i]))

            if line[0] == "2014":
                readon = 0
                    

    f.close()
    soi = ts(year,month,data)
    return soi

def read_nino1():
    return read_psd('Data/nina1.data')

def read_nino4():
    return read_psd('Data/nina4.data')

def read_nino34():
    return read_psd('Data/nina34.data')

def read_psd(filename):

    f = open(filename,'r')

    line = f.readline()
    line = line.strip()
    line = line.split()

    year1 = float(line[0])
    year2 = float(line[1])

    year = []
    month = []
    data = []

    readon = 1

    for line in f:
        if readon == 1:
            line = line.strip()
            line = line.split()

            for i in range(1,13):
                if float(line[i]) != -99.99:
                    year.append(float(line[0]))
                    month.append(float(i))
                    data.append(float(line[i]))

            if float(line[0]) == year2:
                readon = 0


    f.close()
    
    psd = ts(year,month,data)
    return psd

def read_nao():
    f= open('Data/nao_index.tim', 'r')

    year = []
    month = []
    data = []

    for i in range(1,10):
        f.readline()

    for line in f:
        line = line.strip()
        line = line.split()
        year.append(float(line[0]))
        month.append(float(line[1]))
        data.append(float(line[2]))

    f.close()

    nao = ts(year,month,data)
    return nao

def read_ao():
    f = open('Data/monthly.ao.index.b50.current.ascii', 'r')

    year = []
    month = []
    data = []

    for line in f:
        line = line.strip()
        line = line.split()
        year.append(float(line[0]))
        month.append(float(line[1]))
        data.append(float(line[2]))

    f.close()

    ao = ts(year,month,data)
    return ao

def read_aao():
    f = open('Data/monthly.aao.index.b79.current.ascii', 'r')

    year = []
    month = []
    data = []

    for line in f:
        line = line.strip()
        line = line.split()
        year.append(float(line[0]))
        month.append(float(line[1]))
        data.append(float(line[2]))

    f.close()

    aao = ts(year,month,data)
    return aao

nino1 = read_nino1()
nino1.plot_ts()

nino4 = read_nino4()
nino4.plot_ts()

nino34 = read_nino34()
nino34.plot_ts()

aao = read_aao()
aao.plot_ts()

soi = read_soi()
soi.plot_ts()

ao = read_ao()
ao.plot_ts()

nao = read_nao()
nao.plot_ts()



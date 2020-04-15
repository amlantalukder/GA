import sys, random, math, re

# -----------------------------------------------
def readFile(file_name):

    fl = open(file_name)
    data = fl.readlines()
    fl.close()
    return data

# -----------------------------------------------
def writeFile(file_name, data, mode='w'):

    fl = open(file_name, mode)
    fl.write(data)
    fl.close()

# -----------------------------------------------
def readDataTable(filename, delim=','):
    
    return [item.strip().split(delim) for item in readFile(filename)]

# -----------------------------------------------
def writeDataTable(data, filename, mode='w'):
    
    writeFile(filename, formatDataTable(data, ",", "\n"), mode)

# -----------------------------------------------
def formatDataTable(data, col_sep=",", row_sep="\n"):

    return row_sep.join([col_sep.join([str(item1) for item1 in item]) for item in data])

# -----------------------------------------------
def printDec(msg):
    horizontal_border = '_' * 50
    vertical_border = '|'

    l = len(horizontal_border)

    print(horizontal_border)
    print(' ' * l)
    msg_part = msg.strip()
    while len(msg_part) >= l - 4:
        print(vertical_border + ' ' + msg_part[:l - 4] + ' ' + vertical_border)
        msg_part = msg_part[l - 4:].strip()
    print(vertical_border + ' ' + msg_part + ' ' * (l - 3 - len(msg_part)) + vertical_border)
    print(horizontal_border)
    print("")

#-----------------------------------------------
def getSettings(file_name):

    data = [re.sub(r'[ ]+', ' ', item.strip()) for item in readFile(file_name)]
    
    #-----------------------------------------------
    # Get GA settings
    #-----------------------------------------------
    return dict([item.split(' :') for item in data[:18]])

# -----------------------------------------------
def showPercBar(counter, size, perc, perc_inc=10):
    num_prints = 100 // perc_inc
    progress = int(counter * 10 / size) * 10
    if progress >= perc:
        sys.stdout.write('=' * (int((progress - perc) // num_prints) + 1))
        sys.stdout.flush()
        if progress >= 100:
            print('100%')
        perc = progress + perc_inc
    return perc

# -----------------------------------------------
class Queue:
    def __init__(self, qlen):
        self.items = []
        self.max_len = qlen

    def isEmpty(self):
        return self.items == []

    def isFull(self):
        return len(self.items) >= self.max_len

    def enqueue(self, item):
        self.items.insert(0,item)

    def dequeue(self):
        return self.items.pop()

    def size(self):
        return len(self.items)
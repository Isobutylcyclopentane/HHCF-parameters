'''
AFM parameter extraction program
'''
from matplotlib import pyplot as plt
from scipy import optimize as opt

import math
import numpy as np
import statistics
import os
import glob
import ntpath


def read_file (file_name):

    result_x = []
    result_y = []
    result = []

    file = open(file_name, 'r')
    for i in range(3):
        file.readline()
    for line in file:
        tempLst1 = line.strip().split('  ')
        for i in range(len(tempLst1)):
            if 'e' in tempLst1[i]:
                tempInt = tempLst1[i].strip().split('e')
                tempInt = float(tempInt[0]) * (10 ** int(tempInt[1]))
            elif 'E' in tempLst1[i]:
                tempInt = tempLst1[i].strip().split('E')
                tempInt = float(tempInt[0]) * (10 ** int(tempInt[1]))
            else:
                tempInt = float(tempLst1[i])

            if (i == 0):
                result_x.append(tempInt)
            elif (i == 1):
                result_y.append(tempInt)

        result.append((result_x[-1],result_y[-1]))

    return result_x, result_y, result


def find_saturation(quadrum):

    temp_y = []
    temp = []
    auto_range = round(len(quadrum) / 10)
    #print("DEBUG auto_range",auto_range)
    for i in range(auto_range):
        tmp_val = quadrum[-(i+1)]
        temp_y.append(tmp_val[1])
        temp.append(tmp_val)

    for index in range((len(quadrum)-auto_range),-1,-1):
        tmp_val = quadrum[index]

        stdv = statistics.stdev(temp_y)
        avg = statistics.mean(temp_y)
        upper_bound = avg + stdv
        lower_bound = avg - stdv

        if not (lower_bound < tmp_val[1] < upper_bound):
            break
        temp.append(tmp_val)
        temp_y.append(tmp_val[1])
    avg = statistics.mean(temp_y)
    stdv = statistics.stdev(temp_y)
    result = (avg,stdv)
    return result


def get_first_derivative(data):
    result = []
    for i in range(1,len(data)):
        dx = data[i][0] - data[i-1][0]
        dy = data[i][1] - data[i-1][1]
        dydx = dy/dx
        temp = data[i-1],dydx
        result.append(temp)
    return result


def get_second_derivative(firstDerivative):
    result = []
    for i in range(1,len(firstDerivative)):
        d2y = firstDerivative[i][1] - firstDerivative[i-1][1]
        dx2 = firstDerivative[i][0][0] - firstDerivative[i-1][0][0]
        d2ydx2 = d2y / dx2
        temp = firstDerivative[i-1][0],d2ydx2
        result.append(temp)
    return result

def get_xi(firstDerivative):
    temp = []
    tempDxdy = []
    for i in range(3):
        dydx = firstDerivative[i][1]
        temp.append(firstDerivative[i])
        tempDxdy.append(dydx)

    for i in range(3,len(firstDerivative)):
        dydx = firstDerivative[i][1]

        avg = statistics.mean(tempDxdy)
        stdv = statistics.stdev(tempDxdy)
        upperBond = avg + stdv
        lowerBond = avg - stdv

        if not (lowerBond <= dydx <= upperBond):
            break

        tempDxdy.append(dydx)
        temp.append(firstDerivative[i])
    return temp[-1][0][0]


def log_scale_points(data):
    result = []
    for x,y in data:
        if x*y == 0:
            continue
        x = math.log(x)
        y = math.log(y)
        result.append((x,y))
    return result

def find_alpha(logData):
    logFirstDerivative = get_first_derivative(logData)
    tempy = []
    for i in range(3):
        tempy.append(logFirstDerivative[i][1])

    for i in range(3,len(logFirstDerivative)):
        dydx = logFirstDerivative[i][1]
        avg = statistics.mean(tempy)
        stdv = statistics.stdev(tempy)
        upperBond = avg + stdv
        lowerBond = avg + stdv

        if not(lowerBond <= dydx <= upperBond):
            break
        tempy.append(dydx)

    return avg,stdv

def log_scale(data):

    result = []
    for i in data:
        if i <= 0:
           continue
        i = math.log(i)
        result.append(i)
    return result

def HHCF(x,w,xi,alpha):
    return 2*(w**2)*(1-np.exp(-(abs(x)/xi)**(2*alpha)))

def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

def error_check(guesses,x_,y_):
    w,xi,alpha = guesses
    temp = (y_ - HHCF(x_, w, xi, alpha))**2
    temp *= (1000/(x_+1))
    #temp *= np.exp(-x)
    temp_err = statistics.mean(list(temp))
    fit_error.append(temp_err)
    return temp

def log_scale_plot(x,y):
    plt.plot(log_scale(x),log_scale(y))


def main(filename):
    global fit_error
    fit_error = list()


    dataX, dataY, data = read_file(filename)
    data_raw = data.copy()
    dataX_raw = dataX.copy()
    dataY_raw = dataY.copy()
    data = [(x[0] * 10 ** 10, x[1] * 10 ** 20) for x in data]
    dataX = [x * 10 ** 10 for x in dataX]
    dataY = [x * 10 ** 20 for x in dataY]

    x = np.array(dataX)
    y = np.array(dataY)

    x_raw = np.array(dataX_raw)
    y_raw = np.array(dataY_raw)

    saturation, sat_error = find_saturation(data)

    firstDerivative = get_first_derivative(data)
    guessedXi = get_xi(firstDerivative)

    logData = log_scale_points(data)
    guessedAlpha = find_alpha(logData)

    guessedAlpha, aUncert = guessedAlpha
    guessedW = math.sqrt(saturation / 2)
    guessedWUncert = math.sqrt(sat_error / 2)

    guesses = [guessedW, guessedXi, guessedAlpha]
    guesses_demonstrate = guesses.copy()
    guesses_demonstrate[0] *= 10 ** -10
    guesses_demonstrate[1] *= 10 ** -10

    alphaUpper = guessedAlpha + aUncert + 0.5 * guessedAlpha
    alphaLower = guessedAlpha - aUncert - 0.5 * guessedAlpha
    satUpper = guessedW + guessedWUncert + 0.5 * guessedW
    satLower = guessedW - guessedWUncert - 0.5 * guessedW

    boundaries = [(satLower, -np.inf, alphaLower), (satUpper, np.inf, alphaUpper)]

    a0 = opt.least_squares(error_check, guesses[:], bounds=boundaries, \
                           args=(x, y), xtol=1 * 10 ** -15, \
                           ftol=1 * 10 ** -15)
    popt = a0.x
    perr = a0.jac

    w, xi, alpha = list(popt)

    w_raw = w * 10 ** -10
    xi_raw = xi * 10 ** -10
    alpha_raw = alpha

    fname = path_leaf(filename).replace(".txt", "")
    bpath = "C:\\Users\\J. Cheng\\Dropbox\\AFM_image\\0923_AFM\\result"

    print(fname)
    print("approximated value: [w, xi, alpha]\n ", \
          guesses_demonstrate)
    print("RMS roughness: ", w_raw)
    print("lateral correlation length: ", xi_raw)
    print("local roughness: ", alpha_raw)
    # print("uncertainties: ",perr)
    print()
    estimateY = []
    for i in x_raw:
        tempY = HHCF(i, w_raw, xi_raw, alpha_raw)
        estimateY.append(float(tempY))

    estY = np.array(estimateY)

    fit_error = np.array(fit_error)



    plt.plot(x_raw, y_raw, label="observation")
    plt.plot(x_raw, estY, label="model")
    plt.title("Heigh-Height Correlation Function")
    plt.xlabel("r (m)")
    plt.ylabel("H(r) (m)")
    plt.grid(True)
    ppath = os.path.join(bpath, fname+"_HHCF" )
    plt.savefig(ppath)
    plt.close()

    plt.plot(x_raw, y_raw, label="observation")
    plt.plot(x_raw, estY, label="model")
    plt.title("Height-Height Correlation Function in log scale")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("r (m)")
    plt.ylabel("H(r) (m^2)")
    plt.grid(True)
    ppath = os.path.join(bpath, fname + "_log")
    plt.savefig(ppath)
    plt.close()

    # log_scale_plot(x_raw,y_raw)
    # log_scale_plot(x_raw,estY)
    # plt.show(block=True)

    fit_error_raw = np.array(fit_error) * 10 ** -20
    plt.plot(fit_error_raw)
    plt.title("Error Function")
    plt.xlabel("attemps(times)")
    plt.ylabel("errors")
    plt.grid(True)
    ppath = os.path.join(bpath, fname + "_error")
    plt.savefig(ppath)
    plt.close()


if __name__ == "__main__":
    path = "C:\\Users\\J. Cheng\\Dropbox\\AFM_image\\0923_AFM"
    for pathname in glob.glob(os.path.join(path, '*.txt')):
        main(pathname)
'''
    global fit_error
    fit_error = list()

    filename = input("File name? => ")
    dataX,dataY,data = read_file(filename)
    data_raw = data.copy()
    dataX_raw = dataX.copy()
    dataY_raw = dataY.copy()
    data = [(x[0]*10**10,x[1]*10**20) for x in data]
    dataX = [x *10 **10 for x in dataX]
    dataY = [x*10**20 for x in dataY]

    x = np.array(dataX)
    y = np.array(dataY)

    x_raw = np.array(dataX_raw)
    y_raw = np.array(dataY_raw)

    saturation,sat_error = find_saturation(data)

    firstDerivative = get_first_derivative(data)
    guessedXi= get_xi(firstDerivative)

    logData = log_scale_points(data)
    guessedAlpha = find_alpha(logData)

    guessedAlpha,aUncert = guessedAlpha
    guessedW = math.sqrt(saturation/2)
    guessedWUncert = math.sqrt(sat_error/2)

    guesses = [guessedW,guessedXi,guessedAlpha]
    guesses_demonstrate = guesses.copy()
    guesses_demonstrate[0] *= 10 ** -10
    guesses_demonstrate[1] *= 10 ** -10

    alphaUpper = guessedAlpha + aUncert + 0.5*guessedAlpha
    alphaLower = guessedAlpha - aUncert - 0.5*guessedAlpha
    satUpper = guessedW + guessedWUncert + 0.5*guessedW
    satLower = guessedW - guessedWUncert - 0.5*guessedW

    boundaries = [(satLower,-np.inf,alphaLower),(satUpper,np.inf,alphaUpper)]

    a0 = opt.least_squares(error_check,guesses[:],bounds=boundaries,\
                           args=(x,y),xtol=1*10**-15,\
                           ftol=1*10**-15)
    popt = a0.x
    perr = a0.jac

    w,xi,alpha = list(popt)

    w_raw = w * 10 **-10
    xi_raw = xi * 10 ** -10
    alpha_raw = alpha

    print("approximated value: [w, xi, alpha]\n ",\
          guesses_demonstrate)
    print("RMS roughness: ",w_raw)
    print("lateral correlation length: ",xi_raw)
    print("local roughness: ",alpha_raw)
    #print("uncertainties: ",perr)

    estimateY = []
    for i in x_raw:
        tempY = HHCF(i,w_raw,xi_raw,alpha_raw)
        estimateY.append(float(tempY))

    estY = np.array(estimateY)

    fit_error = np.array(fit_error)
    bpath = "C:\\Users\\J. Cheng\\Dropbox\\AFM_image\\0923_AFM\\result"

    plt.plot(x_raw, y_raw,label="observation")
    plt.plot(x_raw,estY,label = "model")
    plt.title("Heigh-Height Correlation Function")
    plt.xlabel("r (m)")
    plt.ylabel("H(r) (m)")
    plt.grid(True)
    
    plt.savefig()

    plt.plot(x_raw, y_raw, label="observation")
    plt.plot(x_raw, estY, label="model")
    plt.title("Height-Height Correlation Function in log scale")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("r (m)")
    plt.ylabel("H(r) (m^2)")
    plt.grid(True)
    plt.show(block=True)

    #log_scale_plot(x_raw,y_raw)
    #log_scale_plot(x_raw,estY)
    #plt.show(block=True)

    fit_error_raw = np.array(fit_error) * 10 ** -20
    plt.plot(fit_error_raw)
    plt.title("Error Function")
    plt.xlabel("attemps(times)")
    plt.ylabel("errors")
    plt.grid(True)
    plt.show(block = True)
'''

'''
    #Alternative plot method
    
    plt.plot(x_raw,y_raw)
    plt.plot(x_raw,estY_raw)
    plt.show(block=True)

    log_scale_plot(x_raw,y_raw)
    log_scale_plot(x_raw,estY_raw)
    plt.show(block=True)


    fit_error_raw = np.array(fit_error) * 10 ** -20
    print("DEBUG fit_error",(fit_error_raw[1]/len(y)))

'''
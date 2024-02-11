import math
import numpy as np


def bisect(f, a, b, it=20, tol=1.e-5):
    """
    Bisection method to estimate the root of f(x)
    Inputs:
        f => function
        a => lower bracket
        b => upper bracket
        it => iterations (default: 20)
        tol => tolerance (default: 1.e-5)
    Output: 
        c => root estimate
            OR...
        Error message if guesses do not bracket root
    """

    prev = a  # sets initial value for tolerance calculation

    if f(a)*f(b) > 0:  # checks the initial brackets
        return "Initial guesses do not bracket solution."

    for i in range(it):
        c = (a+b)/2  # calculates the midpoint
        if abs(c - prev) < tol:
            break  # checks the midpoint compared to tolerance
        if f(c)*f(a) > 0:  # if true(+), root does not lie in this interval
            a = c  # creates new lower bound => old midpoint
        else:
            b = c  # creates new upper bound => old midpoint
        prev = c  # current value of c => prev to calc new error

    return c


def newtraph(f, df, x, it=100, tol=1.e-5):
    """
    This function uses the Newton-Raphson method to estimate the root of f(x)
    Inputs:
        f => main function
        df => derivative of the function
        x => initial guess
        it => iterations (Default: 100)
        tol => tolerance of error (Default: 1.e-5)
    Output:
        xN => root estimate
        i => number of iterations
        error => relative approximate error
    """

    for i in range(it):
        xN = x - (f(x) / df(x))
        error = abs(xN - x)
        if error < tol:
            break
        x = xN

    return xN, i, error


def strlineregr(x, y):
    '''
    Given x and y values, the function returns best fit straight line
    regression parameters.

    Inputs:
        Paired x and y values

    Outputs:
        a0 => line intercept
        a1 => line slope
        Rsq => R squared value (coefficient of determination)
        SE => standard error
    '''

    # Error handling
    if len(x) != len(y):
        return "X and Y must have same length."

    # Stores number of values
    n = len(x)

    # Sum of x-values
    sumX = np.sum(x)
    # Average of sumX
    xBar = sumX / n

    # Sum of y-values
    sumY = np.sum(y)
    # Average of sumY
    yBar = sumY / n

    # Initializing sum of x^2
    sumSqX = 0
    # Initializing sum of XY
    sumXY = 0
    # Initializing error variable => empty matrix of length n to hold error vals
    e = np.zeros((n))
    # Initializing SST and SSE
    SST = 0
    SSE = 0

    # Loop for calculating sum of squared x-vals and sum of x*y-vals
    for i in range(n):
        # Taking each element in list of x-values and adding it to summed squared
        sumSqX += x[i]**2
        # Taking each element in list of x & y-values and adding it to XY sum
        sumXY += x[i]*y[i]
        # Equation from slides to give the slope of equation
        a1 = (n*sumXY - sumX*sumY) / (n*sumSqX - sumX**2)

    # Equation of the intercept
    a0 = yBar - a1*xBar

    # Loop for calculating the errors and squared errors
    for i in range(n):
        # Difference observed and predicted
        e[i] = y[i] - (a0 + a1*x[i])
        # Total squared error of the data
        SST += (y[i] - yBar)
        # Calculates the sum of squared error
        SSE += e[i]**2

    # Calculating the portiong of accounted for variability in model
    SSR = SST - SSE
    # Calculating the R^2 value (regression value)
    Rsq = SSR / SST
    # Standard Error
    SE = np.sqrt(SSE / (n-2))

    return a0, a1, Rsq, SE


def newtint(x, y, xx):
    '''
    Function that uses a Newton Interpolating Polynomial
    Uses n-1 order polynomial given an n amount of data points to
    return a value of the dependent variable => yint at a provided
    independent variable => xx

    Input:
        x => array of independent variables
        y => array of dependent variables
        xx => desired x-value to interpolate

    Output:
        yint => interpolated value at xx
    '''

    # Check if the data is the same length
    n = len(x)
    if len(y) != n:
        return 'X and Y must be the same length.'

    # Initialize holder for calculations -> used to fill with dependent variables
    b = np.zeros((n, n))
    # Assign dependent variables to the 1st column of b for all rows in that column
    b[:, 0] = np.transpose(y)  # switch from input row to input column in b

    # Looping through each value in array to calculate the divided differences
    for j in range(1, n):
        for i in range(n-j):
            # Equation for divided differences based on next and prev points (see slides)
            b[i, j] = (b[i+1, j-1] - b[i, j-1]) / (x[i+j] - x[i])

    xt = 1
    yint = b[0, 0]

    # Solves for yint based on xx and divided differences
    for j in range(n-1):
        xt *= (xx - x[j])
        yint += b[0, j+1] * xt

    return yint


def lagrange(x, y, xx):
    '''
    Lagrange Interpolating Polynomial
    Uses n-1 order lagrange interpolating polynomial based on n number of data points to 
    return a value of the dependent variable yint given the independent variable xx

    Input:
        x => array of independent variable values
        y => array of dependent variable values
        xx => desired independent variable to interpolate

    Output:
        yint => interpolated value
    '''

    # Checking to see if everything is the same length
    n = len(x)
    if len(y) != n:
        return 'X and Y must be the same length'

    # Creating a placeholder
    s = 0

    for i in range(n):
        product = y[i]

        for j in range(n):
            if i != j:
                # This is the weighting equation => L (see slides)
                product *= (xx - x[j]) / (x[i] - x[j])
        s += product

    yint = s
    return yint

# Defining the trapezoidal rule quadrature


def trapes(f, a, b, n=100):
    '''
    Composite Trapezoidal Rule Quadrature:

    Inputs:
        f => name of function to be integrated
        a, b => integration limits (the interval)
        n => number of segment intervals (Default: 100)

    Outputs:
        I => estimate of integration (the result)
    '''

    # Error checking to ensure interval bounds are correct
    if b <= a:
        return 'Upper bound must be greater than lower bound.'
    # Assigning lower bound to first value (x) which will be then calculated for in f(x)
    x = a
    # Width of each segment
    h = (b-a) / n
    # Summation of the values => from the formula I =  ...
    s = f(a)  # funciton evaluated at the lower bound
    # Looping to find inner interval values for each segment
    for i in range(n-1):
        x += h
        s += 2*f(x)
    # Adds functin evaluated at the upper bound
    s += f(b)
    # Composite Trapezoidal Rule
    I = ((b-a)/(n*2))*s

    return I

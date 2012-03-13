def average(l):
    # Calculate average
    return sum(l)/float(len(l))
    
def covariance(xs,ys):
    # Calculate covariance
    prod=[]
    if len(xs)!=len(ys):
        raise Exception('xs is not same length as ys')
    for count in range(len(xs)):
        prod.append(xs[count]*ys[count])
    covar=average(prod)-average(xs)*average(ys)
    return covar
    
def variance(l):
    # Calculate the variance
    import math
    var=[]
    avg=average(l)
    for val in l:
        var.append(math.pow(val-avg,2))
    return average(var)
    
def stddev(l):
    # Calculate the standard deviation
    import math
    return math.sqrt(variance(l))
    
def correlation(xs,ys):
    # Calculate Pearson's correlation coefficient
    if len(xs)==1:
        return 1.0
    return covariance(xs,ys)/(stddev(xs)*stddev(ys))
import sys
import numpy as np
import math
from hmmlearn import hmm
from sklearn.externals import joblib

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt

#Adapted from Steve Busan's generateReactivityProfiles.py script in the Shapemapper pipeline.
#Altered to process multiple replicates, one replicate per column.
#profile contains +modRate - -modRate for each nucleotide
#maxBackground for SHAPEMaP is 0.05, but icSHAPE reactivity are much higher; 0.3 equivocal.
#bgSignal is the background (minus control) modification rate.
#HMM functions added to handle 5' drop-off rates declining before a gap in total read coverage. (Artifact of RT-stop probing).
#INPUT must be in the format: bgCov_rep1...bgCov_repN rxCov_rep1...rxCov_repN
#A crucial assumption in this script is that there are the same number of treated and background control replicates.

global confminDepth, confmaxBackground, hmmDepth
confminDepth=180 #positions below this depth will be reported as -999 reactivity ("no data")
confmaxBackground=0.3
hmmDepth=180 #depth of total coverage below which it is considered a potential "gap" in hmm
#if len(sys.argv)>4:
#    confmaxBackground=float(sys.argv[4])

#find the local maxima that occurs before each gap. Set state=0 in the region between localMax and gap.
def backtrackGap(states,depth): #states and depth are 2d np arrays of the same length. One rep per col.
    if np.ndim(depth)==1: #in case given data has only one replicate, calling depth[i][j] would error
        depth = depth.reshape(-1,1) #reshapes [1,2,3...] to [[1],[2],[3],..]
        #states is guaranteed to be 2d because I used unflatten() on it in HMM()
    for j in range(0,len(states[0])):
        lastGap = 0
        i=0
        while i < (len(states)):
            if states[i][j]==0:
                localMax = 0
                localMaxPosition = i
                position = i-1
                while position>lastGap:
                    if depth[position][j] >= localMax:
                        localMax = depth[position][j]
                        localMaxPosition = position
                    position -= 1
                #Change the state in rxStates to 0 for the positions between localMaxPosition and i
                for k in range(localMaxPosition+1,i):
                    states[k][j] = 0
                #Fast forward i until we pass the current gap region
                gap = True
                while gap and i<len(states):
                    if states[i][j]!=0:
                        gap = False
                        i-=1
                    else:
                        i+=1
                lastGap = i
            i+=1
    #plotting states
    '''
    x=np.arange(0,len(states))
    aveCov = np.mean(depth[:,0])
    temp = states*depth
    plt.plot(x, temp[:,0],'b-', label="Coverage")
    plt.plot(x, states[:,0]*aveCov,'g-', label="HMM State")
    plt.legend(loc='best')
    plt.ylabel('Coverage or HMM')
    plt.xlabel('Mapping Position')
    plt.savefig('newHMMplot.png')
    plt.clf() 
    '''
    return states*depth #converts gap regions in depth array to 0

def unflatten(array1d,lengths): #this only works if the flattened arrays had the same length, which they should here
    array2d = np.zeros((int(len(array1d)/len(lengths)),len(lengths)))
    #print(array1d, len(array1d))
    for i in range(0,len(array1d)):
        j = int(i/lengths[0])
        array2d[i-j*lengths[j]][j] = array1d[i]
    return array2d

def HMM(bgDepth,rxDepth): #Determine positions of too-low coverage in bg and rx with a HMM
    SCRIPTPATH = __file__.rstrip("normalizeSHAPE.py")
    model = joblib.load(SCRIPTPATH+"coverageHMM.pkl") #load the trained HMM model
    rxLengths = [len(rxDepth)]*len(rxDepth[0])
    bgLengths = [len(bgDepth)]*len(bgDepth[0])
    #print(rxDepth[1:10])
    rxDepth1d = rxDepth.flatten('F').reshape(-1,1) #converting replicate columns into one long column
    bgDepth1d = bgDepth.flatten('F').reshape(-1,1) #reshape(-1,1) converts from [1,2,3...] to [[1],[2],[3],...] which hmmlearn requires
    rxDepth1d[rxDepth1d<hmmDepth] = 0
    rxDepth1d[rxDepth1d>=hmmDepth] = 1
    bgDepth1d[bgDepth1d<hmmDepth] = 0
    bgDepth1d[bgDepth1d>=hmmDepth] = 1
    bgStates = model.predict(bgDepth1d,lengths=bgLengths) #predict gap or no_gap states based on model
    rxStates = model.predict(rxDepth1d,lengths=rxLengths)
    bgStates = (np.array(bgStates)-1)*-1 #Gap regions are encoded as '1' in the states arrays,
    rxStates = (np.array(rxStates)-1)*-1 #flipping this to be encoded by '0',
    #this is so when depth matrix is multiplied by a states matrix, gap regions become 0
    #and non-gap regions keep the same depth.
    #print states and depths to check:
    '''
    states = open("states.txt", 'w')
    depths = open("depths.txt", 'w')
    for i in range(0,len(bgStates)):
        states.write(str(bgStates[i]) +"\t" + str(rxStates[i]) + "\n")
        depths.write(str(bgDepth1d[i]) +"\t" + str(rxDepth1d[i]) + "\n")
    states.close()
    depths.close()
    '''
    bgStates2d = unflatten(bgStates,bgLengths) #converting states frm 1d to 2d: [[stateRep1,stateRep2],.]
    rxStates2d = unflatten(rxStates,rxLengths)
    
    #plotting states
    '''
    x=np.arange(0,len(rxStates[0:rxLengths[0]]))
    aveCov = np.mean(rxDepth[:,0])
    #print(rxDepth[1:10])
    plt.plot(x, bgDepth[:,0],'b-', label="Coverage")
    plt.plot(x, bgStates2d[:,0]*aveCov,'g-', label="HMM State")
    plt.legend(loc='best')
    plt.ylabel('Coverage or HMM')
    plt.xlabel('Mapping Position')
    plt.savefig('newHMMplot.png')
    plt.clf()  
    '''
    newrxDepth = backtrackGap(rxStates2d,rxDepth)
    newbgDepth = backtrackGap(bgStates2d,bgDepth)

    return newbgDepth,newrxDepth

def filterProfile(profile,bgSignal,rxDepth,bgDepth,maxBackground=confmaxBackground, minDepth=confminDepth, trimEnds=True):
# filter reactivity profile -
# exclude nucs with high background mutation rates
# exclude nucs with low read coverage
# exclude nucs with high standard errors (optional)
#There are a lot of outliers with icSHAPE data. Trying filtering out profiles with -0.2>=rates>=0.2
    filteredProfile = np.zeros((len(profile),len(profile[0])))
    bgDepth, rxDepth = HMM(bgDepth, rxDepth) #filters out positions that are inacurrate due to loss of adequate coverage
    profileStart=0
    profileStop=len(profile)-1
    if trimEnds:
        profileStart = 2
        profileStop -= 20
    for i in range(0,len(profile)):
        for j in range(0,len(profile[0])):
            goodNuc = True
            #if bgSignal[i][j] > maxBackground:
            #    goodNuc = False
            #if profile[i][j] <= -0.01 or profile[i][j]>=0.01:
            #    goodNuc = False
            if bgDepth[i][j] < 1 or rxDepth[i][j] < 1:
                goodNuc = False

            if profile[i][j] > -500 and goodNuc == True and i<profileStop and i>=profileStart: #**End trimming note
                filteredProfile[i][j] = profile[i][j]

            else:
                filteredProfile[i][j] = -999.0
    return filteredProfile
#**5' coverage at the 3' end is nonexistant for technical reasons, even though total coverage can be hi
#5' coverage at the very 5' ends can still be high though.
#Thus I only need to "end trim" (zero out total depth profiles) for the last 20-30 at the 3' end.                               
#modified from Gregg Rice's boxplot normalization script
#
# 1.find the scaling factor by ranking all of the shape values
# 2. take 1.5*abs(Q1-Q3) as a cutoff
# 3. remove either the top 10% of the RNA or the positions above this cutoff, whichever is smaller
# 4. Average the next 10% from the original length of the RNA --> this is the scaling factor
def calcQuartile(x,q,qtype=7):
    #source: http://adorio-research.org/wordpress/?p=125
    # x = array, q = quartile (in % as a decimal)
    y=x
    n = len(y)
    abcd = [(0,   0, 1, 0), # inverse empirical distrib.function., R type 1
      (0.5, 0, 1, 0), # similar to type 1, averaged, R type 2
      (0.5, 0, 0, 0), # nearest order statistic,(SAS) R type 3
      (0,   0, 0, 1), # California linear interpolation, R type 4
      (0.5, 0, 0, 1), # hydrologists method, R type 5
      (0,   1, 0, 1), # mean-based estimate(Weibull method), (SPSS,Minitab), type 6 
      (1,  -1, 0, 1), # mode-based method,(S, S-Plus), R type 7
      (1.0/3, 1.0/3, 0, 1), # median-unbiased ,  R type 8
      (3/8.0, 0.25, 0, 1)   # normal-unbiased, R type 9.
     ]                               
    a, b, c, d = abcd[qtype-1]
    g, j = math.modf( a + (n+b) * q -1)
    if j<0:
        return x[0]
    elif j>=n:
        return x[n-1]
    j = int(math.floor(j))
    if g == 0:
        return x[j]
    else:
        return y[j] + (y[j+1]- y[j])* (c + d * g)
                               
def findBoxplotFactor(array,rxDepth,bgDepth,minDepth): #input could be 2D if more than one replicate
    #in shapemapper it looks like positions with poor depth are not considered in the normFactor calculation, but positions with too high neg mod rate are considered...as 0s...?? But why? Comment in the code says:
    # Following line is behavior that normalization and structure modeling were optimized with,
    # although this behavior is probably not ideal
    # x = [n if n>-500 else 0 for n in array[:,rep]]                          
    normFactor = np.zeros(len(array[0]))
    for rep in range(0,len(array[0])):
        x,o,a = [],[],0
        for n in range(0,len(array[:,rep])):
            if array[n][rep]>-500:
                x.append(array[n][rep])
            #else:
            #    if rxDepth[n][rep]>=minDepth and bgDepth[n][rep]>=minDepth:
            #        x.append(0.)
                    
        if len(x)/10 < 1:
            normFactor[rep] = 1.0
        else:
            x.sort()
            tenPct = int(len(x)/10)
            fivePct = int(len(x)/20)
            #calculate the interquartile range *1.5
            qLimit = 1.5*abs(calcQuartile(x,0.25)-calcQuartile(x,0.75))
            tenLimit = x[len(x)-1 - tenPct]
            fiveLimit = x[len(x)-1 - fivePct]
            #choose the cutoff that eliminates the fewest points
            limit = max(qLimit,tenLimit)
            if len(x)<100:
                limit = max(qLimit,fiveLimit)
            #make new list without the outliers
            for i in range(len(x)):
                if x[i]<limit:
                    o.append(x[i])
            #avg next ten percent
            try:
                for i in range(-tenPct,0):
                    a = o[i] + a
                normFactor[rep] = a/tenPct
                if normFactor[rep] == 0:
                    normFactor[rep] = 1.0
            except IndexError:
                normFactor[rep] = 1.0
    #print "normFactor:", normFactor
    return normFactor
                               
def normalizeData(array, normFactor): #normFactor should be an array of length r, the number of replicates
    newArray = np.zeros((len(array),len(array[0])))
    for i in range(0,len(array)):
        for j in range(0,len(array[0])):
            # -999 indicates high background/no data
            if array[i][j] > -500:
                newArray[i][j] = array[i][j]/normFactor[j] #different reps will have diff norm factors
            else:
                newArray[i][j] = -999.0
    return newArray

def normalizeSHAPEmain(cov,cov5,trimEnds=True):
    rates = cov5/cov #neg1Rate neg2Rate... negnRate pos1Rate treat2Rate...treatnRate
    #print(rates[35:40])    
    numSamples = int(len(rates[0])/2) #assumed to be same number of neg (bg) and treated (rx) samples
    profiles = ((rates[:,numSamples:] - rates[:,0:numSamples])) #+0.025)/20 #/4.45040 #make numbers closer to SHAPE-MaP dists. On average icSHAPE mod rates (in + sample) are 4.77 times higher than in SHAPE-MaP. raw reactivity distributions are spread between -.2 and .2, but SHAPEMaP vals are from -.01 to .01
    #print(profiles[35:40])
    bgSignals = rates[:,0:numSamples]
    rxDepths = cov[:,numSamples:]
    bgDepths = cov[:,0:numSamples]
    
    filteredProfiles = filterProfile(profiles,bgSignals,rxDepths,bgDepths,confmaxBackground,confminDepth,trimEnds)
    #print(filteredProfiles[35:40]
    normFactors = findBoxplotFactor(filteredProfiles,rxDepths,bgDepths, confminDepth)
    #print(normFactors)                      
    normProfile = normalizeData(filteredProfiles, normFactors)
    #print(normProfile[35:40])
    #print(np.mean(normProfile[normProfile>-999]))
    #print(np.median(normProfile[normProfile>-999]))
    #print(np.std(normProfile[normProfile>-999]))
    return normProfile

if __name__ == "__main__":
    #needs 5coverage and coverage files with an equal number of - and + samples,
    #neg are the first columns in .5coverage and .coverage files
    cov5File = open(sys.argv[2],'r')
    covFile = open(sys.argv[1], 'r')
    cov5Lines = cov5File.read().splitlines()
    covLines = covFile.read().splitlines()
    cov5File.close()
    covFile.close()

    cov5 = np.zeros((len(cov5Lines),len(cov5Lines[0].split())), dtype='float')
    cov = np.zeros((len(covLines),len(covLines[0].split())), dtype='int')
    for i in range(0,len(cov5Lines)):
        thisline = cov5Lines[i].split()
        numbers = []
        for j in thisline:
            numbers.append(float(j)) #need to convert coverage values to float
        cov5[i] = numbers
    for i in range(0,len(covLines)):
        thisline = covLines[i].split()
        numbers = []
        for j in thisline:
            numbers.append(int(float(j)))
        cov[i] = numbers
    #print cov5[0:10]
    #print cov[0:10]        
    normProfile = normalizeSHAPEmain(cov,cov5)
    
    outfile = open(sys.argv[3], 'w')
    for i in normProfile:
        outstr=""
        for j in i:
            outstr+=str(j)+"\t"
        outstr = outstr.rstrip('\t')
        outfile.write(outstr+"\n")

    outfile.close()    
    

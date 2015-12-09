'''
jackknife.py

Author: Amanda Buyan
Written: March2013-October2013

Runs an interative jackknife analysis for crossing angle and contacts on a set of TM helix dimer
simulations - see full description below in the main method for a more detailed explanation
'''

import os,sys,MDAnalysis,numpy,shutil,math,multiprocessing,subsets,time,MDAnalysis.analysis.distances,ctypes,gc,tables
import scipy.weave as weave
from decimal import *
from tables import *
from random import randint
from multiprocessing import Process,Manager,Value,Array

'''
factorial

Function taken from this website:

http://pythonstarter.blogspot.co.uk/2009/08/writing-factorial-function-in-python.html

This basically calculates x! and returns the value

x: integer used to calculate the factorial 
'''
def factorial(x):
    if x==0:
        return 1
    else:
        return x*factorial(x-1)

'''
frange

function from
 
http://code.activestate.com/recipes/66472/

This function is exactly like the range function, but returns float increments instead.

start: starting number - provide as a float
end: ending number - provide as a float
inc: increment (for increasing or decreasing)  
'''
def frange(start,end=None,inc=None):
    L=[]
    if end==None:
            end=start+0.0
            start=0.0
    else:
            start+=0.0
    if inc==None:
            inc=1.0
    while 1:
            next=start+len(L)*inc
            if inc>0 and next>=end:
                    break
            elif inc<0 and next<=end:
                    break
            L.append(next)
    return L

'''
mean

Calculates the mean of a dataset and returns it.

nums: dataset used to calculate the mean
'''
def mean(nums):
    total=0.0
    length=0
    for i in range(len(nums)):
        if nums[i]!=-99.0:  
            total=total+nums[i]
            length+=1 
    average=float(total)/float(length)
    return average

'''
stddev

Calculates the standard deviation of a dataset and returns it.

nums: dataset used to calculate mean and standard deviation
average: mean of the nums
'''
def stddev(nums,average):
    total=0.0
    length=0
    for i in range(len(nums)):
        if nums[i]!=-99.0:
            total=total+(nums[i]-average)**2
            length+=1
    standarddev=math.sqrt(float(total)/float(length))
    return standarddev

'''
splitandmean

Splits the xanglearray into right-handed and left-handed crossing angles, calculates the mean of 
both datasets and returns the means.

xanglearray: array to be sorted
'''
def splitandmean(xanglearray):
    lharray=[]
    rharray=[]
    for xangle in xanglearray:
        if float(xangle)>0:
            lharray.append(float(xangle))
        else:
            rharray.append(float(xangle))
    if len(lharray)>0:
        lhmean=mean(lharray)
    else:
        lhmean=-99.00
    if len(rharray)>0:
        rhmean=mean(rharray)
    else:
        rhmean=-99.00
    return lhmean,rhmean

'''
boxplotstats

This function takes a dataset and returns the following statistics: minimum, the first quartile,
the median, the third quartile, and the maximum.

dataset: the data used to get stats for the boxplots
'''
def boxplotstats(dataset):

    #check if there are any -99.0 in the dataset
    while -99.0 in dataset:
        dataset.remove(-99.0)

    #sort the dataset
    dataset.sort()

    #find the median and separate into two arrays
    if len(dataset)%2!=0:
            index_median=((len(dataset)+1)/2)-1
            median=dataset[index_median]
            list1=dataset[0:index_median]
            list2=dataset[index_median+1:]
    else:
            lower=(len(dataset))/2-1
            upper=(len(dataset))/2
            median=(dataset[lower]+dataset[upper])/2
            list1=dataset[0:lower+1]
            list2=dataset[upper:]   

    #find the first and third quartile of the dataset for the first dataset
    if len(list1)%2!=0:
            index_median=((len(list1)+1)/2)-1
            first_quartile=list1[index_median]
    else:
            lower=(len(list1))/2-1
            upper=(len(list1))/2
            first_quartile=(list1[lower]+list1[upper])/2
    
    #find the first and third quartile of the dataset for the first dataset
    if len(list2)%2!=0:
            index_median=((len(list2)+1)/2)-1
            third_quartile=list2[index_median]
    else:
            lower=(len(list2))/2-1
            upper=(len(list2))/2
            third_quartile=(list2[lower]+list2[upper])/2

    #find min and max
    minimum=list1[0]
    maximum=list2[-1]

    #return the variables
    return minimum,first_quartile,median,third_quartile,maximum

'''
rotmatxyz

rotates a specified matrix around a combination of the x,y,and z axis by angles specified in the
input, and returns it.

alpha: angle to rotate about the x-axis
beta: angle to rotate about the y-axis
gamma: angle to rotate about the z-axis
inmat: the matrix that will be rotated
'''
def rotmatxyz(alpha,beta,gamma,inmat):
    rotmatx=[[1,        0,             0        ],
             [0,math.cos(alpha),-math.sin(alpha)],
             [0,math.sin(alpha),math.cos(alpha)]]
    rotmaty=[[math.cos(beta),0,math.sin(beta)],
             [     0,         1,       0],
             [-math.sin(beta),0,math.cos(beta)]]
    rotmatz=[[math.cos(gamma),-math.sin(gamma),0],
             [math.sin(gamma),math.cos(gamma), 0],
             [     0,               0,         1]]
    temprotmat=numpy.dot(rotmatx,rotmaty)
    wholerotmat=numpy.dot(temprotmat,rotmatz)
    outmat=numpy.dot(wholerotmat,inmat)
    return outmat

'''
centroid

calculates the centroid of a set of points and returns it.

coords: the coordinates for calculating the centroid
numcoords: the number of residues in the helix
'''
def centroid(coords,numcoords):
    xcoord=0
    ycoord=0
    zcoord=0
    centroid=[]
    for i in range(0,numcoords):
        xcoord = xcoord+coords[i][0]
        ycoord = ycoord+coords[i][1]
        zcoord = zcoord+coords[i][2]
    centroid.append(xcoord/numcoords)
    centroid.append(ycoord/numcoords)
    centroid.append(zcoord/numcoords)
    return centroid

'''
COM_helix_vector

Calculates the two vectors through the helices by using the COM of the first four and the last four
residues in each helix.  Returns vector through helix 1 and vector through helix 2

h1coords: coordinates of helix 1
h2coords: coordinates of helix 2
lengthh1: length of helix 1
lengthh2: length of helix 2
'''
def COM_helix_vector(h1coords,h2coords,lengthh1,lengthh2):

    #get first four and last four coordinates for each helix
    h1first4 = h1coords[:4]
    h1last4 = h1coords[lengthh1-4:lengthh1]
    h2first4 = h2coords[:4]
    h2last4 = h2coords[lengthh2-4:lengthh2]

    #get the COM of each vector
    COMh1first4 = numpy.array(centroid(h1first4,4))
    COMh1last4 = numpy.array(centroid(h1last4,4))
    COMh2first4 = numpy.array(centroid(h2first4,4))
    COMh2last4 = numpy.array(centroid(h2last4,4))

    #get vectors through helix
    temp1 = COMh1last4-COMh1first4
    temp2 = COMh2last4-COMh2first4

    #hormalize the vectors
    v1 = temp1/math.sqrt(numpy.dot(temp1,temp1))
    v2 = temp2/math.sqrt(numpy.dot(temp2,temp2))

    #return the unit vectors
    return v1,v2

'''
crossing_angle

This crossing angle function does the following: 

1. Reads in the coordinates of each helix from the universe given
2. Does either a Singular Value Decomposition (SVD) through each helix to find the eigenvector 
   through each helix, or finds the Center of Mass (COM) of the first four Calphas and the last 
   four residues of the Calphas, then makes an eigenvector with this
3. Rotates helix 1 coordinates (using eigenvector from previous step to calculate angles) so they
   align with the POSITIVE x-axis, and rotate helix 2 by the same angles.  These new rotated 
   coordinates are used to calculate new eigenvectors through the helices (using either SVD or COM).
4. Calculate cross-product of new vectors and align it with the POSITIVE x-axis.  Rotate the eigen-
   vectors by these angles as well.
5. Calculate crossing angle using the original eigenvectors, then determine if the angle between 
   the two helices is right-handed or left-handed (more details below in the code).
6. Calculates distance as well, and puts a "1" in the matrix where the distance and the crossing
   angle corresponds
7. Returns the matrix to a Queue object, where it will be processed in the main function

NOTE ABOUT SELECTIONS: When dealing with MARTINI for the forcefield, make sure you have the right
martini option chosen (either martini-phill or martini), and make sure you're selecting the ATOMS
and NOT the residues - will not work if you go by the number of backbone particles.  With the 
other forcefields, you should be able to use the residue selection and the default option (reading
backbone atoms as CA only).  Also, should work for both AT and CG simulations.

NOTE ABOUT THE EIGENVECTOR OPTION: You can calculate the crossing angle in two ways: by defining 
the vector through the helix with the COM's of the ends of the helix, or by doing a SVD to 
calculate the "best fit line" through the helix.  Either way should work, but the first option
is the default.

To prepare gro and xtc files so they just contain protein (so program runs faster):
#create an index file to be read into gromacs 
make_ndx -f md.tpr -n index.ndx

#make a trajectory file containing just the protein
#selection for output is protein with all backbone atoms (I just chose 1)
echo '1\n' | trjconv -f md.xtc -s md.tpr -n index.ndx -pbc mol -o protein.xtc

#create a gro file to go with the trajectory file
#selection is same as above (1)
echo '1\n' | editconf -f md.gro -n index.ndx -o protein.gro
'''
def crossing_angle(universe):

    ############################################################################
    #Get helix coordinates, geometric center, and length for helix 1 and 2
    ############################################################################

    #get the CA selection, CA coordinates, and geometric center for 1st helix
    TM1CAs=universe.selectAtoms("bynum %d:%d and (name %s)" % (options["-res"].value[0], options["-res"].value[1], options["-ca"].value))
    helix_1_CA_coordinates = TM1CAs.coordinates()
    helix_1_CA_geometric_center = TM1CAs.centerOfGeometry()
    lengthh1 = len(TM1CAs)

    #get the CA selection, CA coordinates, and geometric center for 2nd helix
    TM2CAs=universe.selectAtoms("bynum %d:%d and (name %s)" % (options["-res"].value[2],options["-res"].value[3], options["-ca"].value))
    helix_2_CA_coordinates = TM2CAs.coordinates()
    helix_2_CA_geometric_center = TM2CAs.centerOfGeometry()
    lengthh2 = len(TM2CAs)

    ############################################################################
    #If option for eigenvectors is false:
    #Calculate the vectors from the COM of the ends of the helix, ELSE:
    #Do a SVD of both vectors to get eigenvectors of the original helices
    ############################################################################

    if (options["-eigen"].value):

        uu, dd, vv = numpy.linalg.svd(helix_1_CA_coordinates-helix_1_CA_geometric_center)
        helix_1_vector = vv[0]
        magh1=math.sqrt(numpy.dot(helix_1_vector,helix_1_vector))

        uu, dd, vv = numpy.linalg.svd(helix_2_CA_coordinates-helix_2_CA_geometric_center)
        helix_2_vector = vv[0]
        magh2=math.sqrt(numpy.dot(helix_2_vector,helix_2_vector))

    else:

        helix_1_vector,helix_2_vector=COM_helix_vector(helix_1_CA_coordinates,helix_2_CA_coordinates,lengthh1,lengthh2)
        magh1=math.sqrt(numpy.dot(helix_1_vector,helix_1_vector))
        magh2=math.sqrt(numpy.dot(helix_2_vector,helix_2_vector))

    ############################################################################
    #Rotate helix 1 until it is aligned with the POSITIVE x-axis - use these 
    #matrices to rotate helix 2 as well
    ############################################################################

    #calculate the angle needed to rotate about z-axis so y-component is 0 - figure out gamma
    gamma=0
    if (helix_1_vector[1]<0):
        if (helix_1_vector[0]<0):
            gamma=math.radians(180)-math.atan(helix_1_vector[1]/helix_1_vector[0])
        elif (helix_1_vector[0]==0):
            if helix_1_vector[1]>0:
                gamma=-math.radians(90)
            else:
                gamma=math.radians(90)
        else:
            gamma=-math.atan(helix_1_vector[1]/helix_1_vector[0])
    else:
        if (helix_1_vector[0]<0):
            gamma=-(math.radians(180)+math.atan(helix_1_vector[1]/helix_1_vector[0]))
        elif (helix_1_vector[0]==0):
            if helix_1_vector[1]>0:
                gamma=math.radians(90)
            else:
                gamma=-math.radians(90)
        else:
            gamma=-math.atan(helix_1_vector[1]/helix_1_vector[0])

    #Rotate helix 1 vector about z-axis
    temph1=rotmatxyz(0,0,gamma,helix_1_vector)

    #calculate the angle needed to rotate about y-axis - figure out beta
    beta=math.atan(temph1[2]/temph1[0])

    #set dummy variable for the helix
    helix_1_CA_coords_rotated=[]
    helix_2_CA_coords_rotated=[]

    #Rotate entire helices about z-axis, then about y-axis
    for i in range(0,lengthh1):
        newcoordh1 = rotmatxyz(0,beta,gamma,helix_1_CA_coordinates[i])
        helix_1_CA_coords_rotated.append(newcoordh1)
    for i in range(0,lengthh2):
        newcoordh2 = rotmatxyz(0,beta,gamma,helix_2_CA_coordinates[i])
        helix_2_CA_coords_rotated.append(newcoordh2)

    #calculate the Center of Mass (COM) of the two helices
    center1=centroid(helix_1_CA_coords_rotated,lengthh1)
    center2=centroid(helix_2_CA_coords_rotated,lengthh2)

    #if the value of eigenvector is true, then use SVD to calculate the new vectors.  
    #if not, then use the same technique as above (ends of helices) 
    if (options["-eigen"].value):

        #do a SVD of the newly rotated coordinates to get the new eigenvectors
        svdvector1=[]
        svdvector2=[]
        for i in range(0,lengthh1):
            svdvector1.append(helix_1_CA_coords_rotated[i]-center1)
        for i in range(0,lengthh2):
            svdvector2.append(helix_2_CA_coords_rotated[i]-center2)
        uu,dd,vv = numpy.linalg.svd(svdvector1)
        newh1 = vv[0]
        uu,dd,vv = numpy.linalg.svd(svdvector2)
        newh2 = vv[0]

    else:

        newh1,newh2=COM_helix_vector(helix_1_CA_coords_rotated,helix_2_CA_coords_rotated,lengthh1,lengthh2)
        magnewh1=math.sqrt(numpy.dot(newh1,newh1))
        magnewh1=math.sqrt(numpy.dot(newh2,newh2))

    ############################################################################
    #take the cross product of the two new vectors and align it with the POSITIVE 
    #x-axis
    ############################################################################

    #calculate cross-product
    cross=numpy.cross(newh1,newh2)
    magcross=math.sqrt(numpy.dot(cross,cross))

    #align with positive x-axis - rotate around z-axis
    gamma=0
    if (cross[1]<0):
        if (cross[0]<0):
            gamma=math.radians(180)-math.atan(cross[1]/cross[0])
        else:
            if (cross[0]<0.005 and cross[0]>-0.005):
                gamma=math.radians(90)
            else:
                gamma=-math.atan(cross[1]/cross[0])
    else:
        if (cross[0]<0):
            gamma=-(math.radians(180)+math.atan(cross[1]/cross[0]))
        else:
            if (cross[0]<0.005 and cross[0]>-0.005):
                   gamma=math.radians(90)
            else:
                   gamma=-math.atan(cross[1]/cross[0])

    #rotate the vector
    temp=rotmatxyz(0,0,gamma,cross)

    #align with x-axis - rotate around y-axis
    beta=math.atan(temp[2]/temp[0])
    newcross=rotmatxyz(0,beta,0,temp)

    #rotate the eigenvectors of the helices
    newh1again=rotmatxyz(0,beta,gamma,newh1)
    newh2again=rotmatxyz(0,beta,gamma,newh2)

    #rotate the actual helices by the same angles
    helix_1_CA_coords_rotagain=[]
    helix_2_CA_coords_rotagain=[]

    #Rotate entire helices about z-axis, then about y-axis
    for i in range(0,lengthh1):
        newcoordh1 = rotmatxyz(0,beta,gamma,helix_1_CA_coords_rotated[i])
        helix_1_CA_coords_rotagain.append(newcoordh1)
    for i in range(0,lengthh2):
        newcoordh2 = rotmatxyz(0,beta,gamma,helix_2_CA_coords_rotated[i])
        helix_2_CA_coords_rotagain.append(newcoordh2)

    ###########################################################################################
    #Calculate crossing angle with original vectors, then determine if the angle between the 
    #two helices is right-handed or left-handed based on the newest eigenvectors.  If 
    #helix 1's COM x-component < helix 2's COM x-component and the x-component of the cross-
    #product between the two new eigenvectors is < 0, then the angle is right-handed.  If 
    #helix 1's COM x-component > helix 2's COM x-component and the x-component of the cross-
    #product between the two new eigenvectors is > 0, then the angle is also right-handed. Any
    #other cases are left-handed angles.
    ###########################################################################################

    #calculate crossing angle
    cross=numpy.cross(helix_1_vector,helix_2_vector)
    magcross=math.sqrt(numpy.dot(cross,cross))
    if (magcross/(magh1*magh2))>1 or (magcross/(magh1*magh2))<-1:
        xangle=90.0000000
    else:
        xangle=math.degrees(math.asin(magcross/(magh1*magh2)))

    #calculate final geometric center
    helix_1_rotated_geometric_center=centroid(helix_1_CA_coords_rotagain,lengthh1)
    helix_2_rotated_geometric_center=centroid(helix_2_CA_coords_rotagain,lengthh2)

    #figure out if the angle between the two helices is right-handed or left-handed  
    if (helix_1_rotated_geometric_center[0]<helix_2_rotated_geometric_center[0]):
        cross2=numpy.cross(newh1again, newh2again)
        if (cross2[0]<0):
            xangle=-1*xangle
    else:
        cross2=numpy.cross(newh1again,newh2again)
        if (cross2[0]>0):
            xangle=-1*xangle

    #return crossing angle
    return xangle

'''
mode

finds the range(s) that the most crossing angle values reside in for both the left-handed and the 
right-handed crossing angles based on the probabilities, and returns the values of the lower and 
upper bounds of the crossing angle for both right-handed and left-handed values. 

xanglematrix: the matrix with all of the different probabilities of the xangle
'''
def mode(xanglematrix):

    #initalize variables
    tempxangle=options["-bounds"].value[0]
    lowermodeslh=[]
    uppermodeslh=[]
    lowermodesrh=[]
    uppermodesrh=[]
    templhmodes=[]
    temprhmodes=[]
    lhprobstemp=[]
    rhprobstemp=[]

    #loop over xanglematrix to figure out the modes
    for i in range(len(xanglematrix)):
        if i==len(xanglematrix)/2:
            if len(lowermodesrh)==0:
                lowermodesrh.append(-99.00)
                uppermodesrh.append(-99.00)
                rhprobstemp.append(xanglematrix[i])
        elif i>3 and i<len(xanglematrix)-3 and i<len(xanglematrix)/2:
            inflectioncheck=checkifinflection(xanglematrix[i-3],xanglematrix[i-2],xanglematrix[i-1],xanglematrix[i],xanglematrix[i+1],xanglematrix[i+2],xanglematrix[i+3])
            if (inflectioncheck):
                lowermodesrh.append(tempxangle)
                uppermodesrh.append(tempxangle+options["-binsize"].value)
                rhprobstemp.append(xanglematrix[i])
        elif i>3 and i<len(xanglematrix)-3 and i>len(xanglematrix)/2 and i!=len(xanglematrix)-1:
            inflectioncheck=checkifinflection(xanglematrix[i-3],xanglematrix[i-2],xanglematrix[i-1],xanglematrix[i],xanglematrix[i+1],xanglematrix[i+2],xanglematrix[i+3])
            if (inflectioncheck):
                lowermodeslh.append(tempxangle)
                uppermodeslh.append(tempxangle+options["-binsize"].value)
                lhprobstemp.append(xanglematrix[i])
        elif i==len(xanglematrix)-1:
            if len(lowermodeslh)==0:
                lowermodeslh.append(-99.00) 
                uppermodeslh.append(-99.00)
                lhprobstemp.append(xanglematrix[i])
        tempxangle+=options["-binsize"].value

    #average the lower and upper bounds for left-handed modes
    for i in range(len(lowermodeslh)):
        mode=(float(lowermodeslh[i])+float(uppermodeslh[i]))/2
        templhmodes.append(mode)

    #average the lower and upper bounds for right-handed modes
    for i in range(len(lowermodesrh)):
        mode=(float(lowermodesrh[i])+float(uppermodesrh[i]))/2
        temprhmodes.append(mode)

    #sort the lists
    lhprobs=len(lhprobstemp)*[None]
    rhprobs=len(rhprobstemp)*[None]
    lhmodes=len(templhmodes)*[None]
    rhmodes=len(temprhmodes)*[None]
    highprob=0.0
    for i in range(len(lhprobstemp)):
        if lhprobstemp[i]>highprob:
            for j in reversed(range(0,len(lhprobstemp)-1,1)):
                lhprobs[j+1]=lhprobs[j]
                lhmodes[j+1]=lhmodes[j]
            lhprobs[0]=lhprobstemp[i]
            lhmodes[0]=templhmodes[i]
            highprob=lhprobstemp[i]
        else:
            lhprobs[i]=lhprobstemp[i]
            lhmodes[i]=templhmodes[i]
            highprob=lhprobstemp[i]
    highprob=0.0	
    for i in range(len(rhprobstemp)):
        if rhprobstemp[i]>highprob:
            for j in reversed(range(0,len(rhprobstemp)-1,1)):
                rhprobs[j+1]=rhprobs[j]
                rhmodes[j+1]=rhmodes[j]
            rhprobs[0]=rhprobstemp[i]
            rhmodes[0]=temprhmodes[i]
            highprob=rhprobstemp[i]
        else:
            rhprobs[i]=rhprobstemp[i]
            rhmodes[i]=temprhmodes[i]
            highprob=rhprobstemp[i]

    #return the modes
    return lhmodes,rhmodes,lhprobs,rhprobs

'''
checkifinflection

This function checks if the point that you are at is an inflection point (i.e. maximum)

twobelow: the value two bins below
onebelow: the value one bin below
value: the value you want to check
oneabove: the value one bin above
twoabove: the value two bins above
'''
def checkifinflection(threebelow,twobelow,onebelow,value,oneabove,twoabove,threeabove):

    #check if inflection point - if so, return true
    if threebelow<value and twobelow<value and onebelow<value and oneabove<value and twoabove<value and threeabove<value:
        return True

    #default return
    return False

'''
put_crossing_angle

put_crossing_angle goes through the whole crossing angle array given, and bins them - if the 
crossing angle is in a certain bin, then it increases the number in that bin by one.  The 
binsize and bounds are specified in the options.

xanglearraytemp: the raw xangle data
R: number of rows in the xangle matrix
'''
def put_crossing_angle(xanglearraytemp,R):

    #convert the array to a ctype and R to an int
    xanglearray=(ctypes.c_float * len(xanglearraytemp))(*xanglearraytemp)
    R=int(R)

    #C code that will be run
    c_code="""

    #include <Python.h>

    //initialize variables
    int tempR=0;
    int SIZE=xanglearray.size();
    float xanglematrix[R];
    float returnxanglematrix[R];

    //set all values of xanglematrix to 0 for initialization
    for (int i=0; i<R; i++)
    {
        xanglematrix[i]=0.0;
    }

    //create a matrix of crossing angle occurences
    if (SIZE>0) {
        for (int i=0; i<SIZE; i++) {
            for (int j=bound1; j<bound2; j+=binsize)	{
                if (xanglearray[i]==j || xanglearray[i]>j && xanglearray[i]<j+binsize) {
                    xanglematrix[tempR]+=1.0;
                }
                tempR+=1;
            }
            tempR=0;
        }
    }

    //normalize
    for (int i=0; i<R; i++)
    {
        returnxanglematrix[i]=xanglematrix[i]/SIZE;
    }

    //create list to return to the function
    py::list returns(R);

    //fill the return variable
    for (int i=0; i<R; i++) {
                returns[i]=returnxanglematrix[i];
        }

    //return the values
    return_val=returns;

    """

    #variables for weave
    bound1=options["-bounds"].value[0]
    bound2=options["-bounds"].value[1]
    binsize=options["-binsize"].value

    #run weave and return value
    return weave.inline(c_code,['R','xanglearray','bound1','bound2','binsize']) 

'''
contact_matrix

This function iterates through each Calpha atom of one helix and calculates the distance 
between it and all of the other Calpha atoms in the other helix, creating an MxN grid, 
where M and N are the number of residues in each respective helix.

directory: directory the contact matrix is calculated for
'''
def contact_matrix(directory):

    #initialize numpy arrays
    contact_matrix_all=numpy.zeros(shape=(options["-TM1l"].value,options["-TM2l"].value)) 
    contact_matrix_left_handed=numpy.zeros(shape=(options["-TM1l"].value,options["-TM2l"].value)) 
    contact_matrix_right_handed=numpy.zeros(shape=(options["-TM1l"].value,options["-TM2l"].value)) 
    contactmatrixindividualframes=[]

    #initialize xanglearray
    xanglearray=[]

    #variable for normalisation
    numframesconsideredall=0
    numframesconsideredlh=0
    numframesconsideredrh=0

    #define universe for calculations
    universe=MDAnalysis.Universe(str(directory)+"/"+str(options["-g"].value),str(directory)+"/"+str(options["-x"].value))

    #loop over the whole trajectory
    for ts in universe.trajectory[::options["-s"].value]:

        #define temporary contact matrix
        contact_matrix_temp=numpy.zeros(shape=(options["-TM1l"].value,options["-TM2l"].value))

        #get the selections for the CA atoms
        TM1CAs=universe.selectAtoms("bynum %d:%d and (name %s)" % (options["-res"].value[0], options["-res"].value[1], options["-ca"].value))
        TM1CAcoords=TM1CAs.coordinates()
        TM1CA_geometric_center=TM1CAs.centerOfGeometry()
        TM2CAs=universe.selectAtoms("bynum %d:%d and (name %s)" % (options["-res"].value[2],options["-res"].value[3], options["-ca"].value))
        TM2CAcoords=TM2CAs.coordinates()
        TM2CA_geometric_center=TM2CAs.centerOfGeometry()

        #calculate the distance between either the COM, or the minimum distance between the helices
        distance=0
        if (options["-mindist"].value):
            distance=100
            for i in range(len(TM1CAcoords)):
                for j in range(len(TM2CAcoords)):
                    tempdist=math.sqrt((TM1CAcoords[i][0]-TM2CAcoords[j][0])**2+(TM1CAcoords[i][1]-TM2CAcoords[j][1])**2+(TM1CAcoords[i][2]-TM2CAcoords[j][2])**2)
                    if tempdist<distance:
                        distance=tempdist
        else:
            distance = math.sqrt((TM1CA_geometric_center[0]-TM2CA_geometric_center[0])**2+(TM1CA_geometric_center[1]-TM2CA_geometric_center[1])**2+(TM1CA_geometric_center[2]-TM2CA_geometric_center[2])**2)

        #check if distance between helices is less than cutoff
        if (distance<options["-cutoff"].value):

            #loop over all of the residues in the matrix
            for i in range(options["-TM1l"].value):
                for j in range(options["-TM2l"].value):
                    contact_matrix_temp[i][j]+=math.sqrt((TM1CAcoords[i][0]-TM2CAcoords[j][0])**2+(TM1CAcoords[i][1]-TM2CAcoords[j][1])**2+(TM1CAcoords[i][2]-TM2CAcoords[j][2])**2)

            #calculate crossing angle
            xangle=crossing_angle(universe)
            xanglearray.append(float(xangle))
            contactmatrixindividualframes.append(contact_matrix_temp)

            #check if there are any contacts, and sort into left-handed and right-handed matrices
            numframesconsideredall+=1
            contact_matrix_all+=contact_matrix_temp
            if (float(xangle)>0):
                contact_matrix_left_handed+=contact_matrix_temp
                numframesconsideredlh+=1
            else:
                contact_matrix_right_handed+=contact_matrix_temp
                numframesconsideredrh+=1

    #get index for output
    index=directory-1

    #update the progress bar
    if totalforbar.value==options["-numdirs"].value:
        update_progress(1)
        print "\n"
    else:
        modify()
        if totalforbar.value<options["-numdirs"].value:
            update_progress(float(totalforbar.value)/options["-numdirs"].value)

    #print directory

    #return the matrices and the number of frames - will make multiprocessing easier
    return index,contact_matrix_all,contact_matrix_left_handed,contact_matrix_right_handed,numframesconsideredall,numframesconsideredlh,numframesconsideredrh,xanglearray,contactmatrixindividualframes

'''
createreslist

This function creates a residue list from the TM's Calphas for gnuplot, and returns it.

TMCAs: universe containing the TM's CAs.
'''
def createreslist(TMCAs):
    reslist=""
    num=1.5
    offset=options["-offset"].value
    for res in TMCAs:
        reslist=reslist+"\\\""+str(res.resname)+str(offset)+"\\\""+str(num)+", "
        offset+=1
        num+=1
    reslist=reslist[:-2]
    return reslist

'''
A program that finds a particular value in a set of numbers - from the following website:

http://www.phanderson.com/C/find_idx.html

a: the array containing all of the values
num_elements: the length of the array
value: the value you are looking for
'''
def find_index(a,num_elements,value):
    for i in range(num_elements):
        if a[i]==value:
            return(value)
    return -1

'''
checknumcores

This function checks the number of cores available for processing.  If there is only one core, then
this script won't run (unless if not availablecores>1 is changed to if not availablecores>0). This
function returns the number of cores, and takes no arguments.
'''
def checknumcores():
    availablecores=multiprocessing.cpu_count()
    if not availablecores>1:
        sys.exit('Need more than 1 core to exploit multi-core code.')
    return availablecores

'''
process_individual_trajectories

This function does the following for all of the trajectories being considered

	1. Splits each trajectory into 8 pieces for multiprocessing
	2. Calculates all of the crossing angles and distances between contacts for that simulation
	3. Calculates the number of frames that are considered in calculating crossing angles and
	   distance
	4. Puts all of the matrices and arrays into their respective bigger arrays
	5. Removes all of the temporary xtc files
	6. Returns all of the arrays and the CAs for the TMs

This function takes no arguments.

It returns the following:
TM1CAs: the universe object containing TM1's CAs
TM2CAs: the universe object containing TM2's CAs
allcontactmatricesnumdir: array of all of the contact matrices of the number of directories
lhcontactmatricesnumdir: array of all of the contact matrices corresponding to a lh crossing angle
rhcontactmatricesnumdir: array of all of the contact matrices corresponding to a rh crossing angle
arraynumframesconsideredall: array containing the number of frames that was considered in the 
			     crossing angle and distance calculations
arraynumframesconsideredlh: array containing the number of frames that was considered for the left
			    handed contact matrices in the crossing angle and distance 
			    calculations
arraynumframesconsideredrh: array containing the number of frames that was considered for the 
                            right-handed contact matrices in the crossing angle and distance 
                            calculations
xangles: array of xangle arrays
'''
def process_individual_trajectories():

    #initialize arrays 
    allcontactmatricesnumdirtemp=options["-numdirs"].value*[None]
    lhcontactmatricesnumdirtemp=options["-numdirs"].value*[None]
    contactmatricesindividualframestemp=options["-numdirs"].value*[None]
    rhcontactmatricesnumdirtemp=options["-numdirs"].value*[None]
    arraynumframesconsideredall=options["-numdirs"].value*[None]
    arraynumframesconsideredlh=options["-numdirs"].value*[None]
    arraynumframesconsideredrh=options["-numdirs"].value*[None]
    xangles=options["-numdirs"].value*[None]

    listofresults=[]

    def log_results_from_pool(listofresults):
        index=listthing[0]
        allcontactmatricesnumdirtemp[index]=listofresults[1]
        lhcontactmatricesnumdirtemp[index]=listofresults[2]
        rhcontactmatricesnumdirtemp[index]=listofresults[3]
        arraynumframesconsideredall[index]=listofresults[4]
        arraynumframesconsideredlh[index]=listofresults[5]
        arraynumframesconsideredrh[index]=listofresults[6]
        xangles[index]=listofresults[7]
        contactmatricesindividualframestemp[index]=listofresults[8]

    #define pool
    pool=multiprocessing.Pool()

    #define universe
    universe=MDAnalysis.Universe("1/"+str(options["-g"].value),"1/"+str(options["-x"].value))

    #get the length of the helices here
    TM1CAs=universe.selectAtoms("bynum %d:%d and (name %s)" % (options["-res"].value[0],options["-res"].value[1],options["-ca"].value))
    TM2CAs=universe.selectAtoms("bynum %d:%d and (name %s)" % (options["-res"].value[2],options["-res"].value[3],options["-ca"].value))

    #update the progress bar
    update_progress(0)

    #loop over all individual trajectories
    for i in range(1,options["-numdirs"].value+1):

        #do pool thing
        pool.apply_async(contact_matrix,[i],callback=log_results_from_pool)

    #close and join the pool
    pool.close()
    pool.join()

    #set value of these, and make sure they are float32 tupe
    allcontactmatricesnumdir=numpy.array(allcontactmatricesnumdirtemp,dtype="float32")
    lhcontactmatricesnumdir=numpy.array(lhcontactmatricesnumdirtemp,dtype="float32")
    rhcontactmatricesnumdir=numpy.array(rhcontactmatricesnumdirtemp,dtype="float32")

    #return the variables
    return TM1CAs,TM2CAs,allcontactmatricesnumdir,lhcontactmatricesnumdir,rhcontactmatricesnumdir,arraynumframesconsideredall,arraynumframesconsideredlh,arraynumframesconsideredrh,xangles,contactmatricesindividualframestemp

'''
contact_matrix_closest_contacts_thetans

This function averages all of the contact matrices together, and then determines the top X contacts 
for the overall, left-handed and right-handed contact matrices (X is specified by the user).  It 
then writes the contact matrix to a file, and then returns the top 6 contacts.

allcontactmatrices: all of the contact frames for each simulation
lhcontactmatrices: all of the left-handed contact frames for each simulation
rhcontactmatrices: all of the right-handed contact fraomes for each simulation
TM1CAs: Helix 1's CAs
TM2CAs: Helix 2's CAs
'''
def contact_matrix_closest_contacts_thetans(allcontactmatrices,lhcontactmatrices,rhcontactmatrices,TM1CAs,TM2CAs):

    #create an overall contact matrix
    overallCM=numpy.zeros(shape=(options["-TM1l"].value,options["-TM2l"].value))
    lhCM=numpy.zeros(shape=(options["-TM1l"].value,options["-TM2l"].value))
    rhCM=numpy.zeros(shape=(options["-TM1l"].value,options["-TM2l"].value))

    #initialize variables for finding the closest options["-topres"].value contacts
    lowestdistall=options["-topres"].value*[None]
    lowestdistlh=options["-topres"].value*[None]
    lowestdistrh=options["-topres"].value*[None]
    jsoverall=options["-topres"].value*[0]
    ksoverall=options["-topres"].value*[0]
    jslh=options["-topres"].value*[0]
    kslh=options["-topres"].value*[0]
    jsrh=options["-topres"].value*[0]
    ksrh=options["-topres"].value*[0]

    #initialize the distances
    for i in range(options["-topres"].value):
        lowestdistall[i]=9999
        lowestdistlh[i]=9999
        lowestdistrh[i]=9999

    #find the options["-topres"].value closest contacts for each contact matrix and adds it to overall contact matrices
    for j in range(options["-TM1l"].value):
        for k in range(options["-TM2l"].value):
            for l in range(options["-topres"].value):
                if allcontactmatrices[j][k]<lowestdistall[l]:
                    for m in reversed(range(l,options["-topres"].value-1,1)):
                        lowestdistall[m+1]=lowestdistall[m]
                        jsoverall[m+1]=jsoverall[m]
                        ksoverall[m+1]=ksoverall[m]
                    lowestdistall[l]=allcontactmatrices[j][k]
                    jsoverall[l]=j
                    ksoverall[l]=k
                    break
    for j in range(options["-TM1l"].value):
        for k in range(options["-TM2l"].value):
            for l in range(options["-topres"].value):
                if lhcontactmatrices[j][k]<lowestdistlh[l]:
                    for m in reversed(range(l,options["-topres"].value-1,1)):
                        lowestdistlh[m+1]=lowestdistlh[m]
                        jslh[m+1]=jslh[m]
                        kslh[m+1]=kslh[m]
                    lowestdistlh[l]=lhcontactmatrices[j][k]
                    jslh[l]=j
                    kslh[l]=k
                    break
    for j in range(options["-TM1l"].value):
        for k in range(options["-TM2l"].value):
            for l in range(options["-topres"].value):
                if rhcontactmatrices[j][k]<lowestdistrh[l]:
                    for m in reversed(range(l,options["-topres"].value-1,1)):
                        lowestdistrh[m+1]=lowestdistrh[m]
                        jsrh[m+1]=jsrh[m]
                        ksrh[m+1]=ksrh[m]
                    lowestdistrh[l]=rhcontactmatrices[j][k]
                    jsrh[l]=j
                    ksrh[l]=k
                    break
    for j in range(options["-TM1l"].value):
        for k in range(options["-TM2l"].value):
            for l in range(options["-topres"].value):
                if jsoverall[l]==j and ksoverall[l]==k:
                    overallCM[j][k]+=1
                if jslh[l]==j and kslh[l]==k:
                    lhCM[j][k]+=1
                if jsrh[l]==j and ksrh[l]==k:
                    rhCM[j][k]+=1

    #open output files
    allCMfile=open(options["-CMall"].value,"w")
    lhCMfile=open(options["-CMlh"].value,"w")
    rhCMfile=open(options["-CMrh"].value,"w")

    #write out the final contact matrix file
    for i in range(options["-TM1l"].value):
        for j in range(options["-TM2l"].value):
            allCMfile.write(str(overallCM[i][j])+" ")
            lhCMfile.write(str(lhCM[i][j])+" ")	
            rhCMfile.write(str(rhCM[i][j])+" ")	
        allCMfile.write("\n")	
        lhCMfile.write("\n")	
        rhCMfile.write("\n")

    #close output files
    allCMfile.close()
    lhCMfile.close()
    rhCMfile.close()

    #reset the counter
    for i in range(options["-topres"].value):
        lowestdistall[i]=0
        lowestdistlh[i]=0
        lowestdistrh[i]=0
        jsoverall[i]=0
        ksoverall[i]=0
        jslh[i]=0
        kslh[i]=0
        jsrh[i]=0
        ksrh[i]=0

    #find overall closest contacts
    for j in range((options["-TM1l"].value)):
        for k in range(options["-TM2l"].value):
            for l in range(options["-topres"].value):
                if overallCM[j][k]>lowestdistall[l]:
                    for m in reversed(range(l,options["-topres"].value-1,1)):
                        lowestdistall[m+1]=lowestdistall[m]
                        jsoverall[m+1]=jsoverall[m]
                        ksoverall[m+1]=ksoverall[m]
                    lowestdistall[l]=overallCM[j][k]
                    jsoverall[l]=j
                    ksoverall[l]=k
                    break
    for j in range((options["-TM1l"].value)):
        for k in range(options["-TM2l"].value):
            for l in range(options["-topres"].value):
                if lhCM[j][k]>lowestdistlh[l]:
                    for m in reversed(range(l,options["-topres"].value-1,1)):
                        lowestdistlh[m+1]=lowestdistlh[m]
                        jslh[m+1]=jslh[m]
                        kslh[m+1]=kslh[m]
                    lowestdistlh[l]=lhCM[j][k]
                    jslh[l]=j
                    kslh[l]=k
                    break
    for j in range((options["-TM1l"].value)):
        for k in range(options["-TM2l"].value):
            for l in range(options["-topres"].value):
                if rhCM[j][k]>lowestdistrh[l]:
                    for m in reversed(range(l,options["-topres"].value-1,1)):
                        lowestdistrh[m+1]=lowestdistrh[m]
                        jsrh[m+1]=jsrh[m]
                        ksrh[m+1]=ksrh[m]
                    lowestdistrh[l]=rhCM[j][k]
                    jsrh[l]=j
                    ksrh[l]=k
                    break

    #initialize output arrays
    closestcontactsoverall=[]
    closestcontactslh=[]
    closestcontactsrh=[]

    #loop over everything
    for i in range(options["-topres"].value):

        #first the overall closest contacts
        resnum1overall=jsoverall[i]+options["-offset"].value
        resnum2overall=ksoverall[i]+options["-offset"].value
        resid1overall=TM1CAs[jsoverall[i]].resname
        resid2overall=TM2CAs[ksoverall[i]].resname
        res1overall=str(resnum1overall)+str(resid1overall)
        res2overall=str(resnum2overall)+str(resid2overall)
        distoverall=str()
        closestcontactsoverall.append([res1overall,res2overall])
        
        #second, the left-handed closest contacts
        resnum1lh=jslh[i]+options["-offset"].value
        resnum2lh=kslh[i]+options["-offset"].value
        resid1lh=TM1CAs[jslh[i]].resname
        resid2lh=TM2CAs[kslh[i]].resname
        res1lh=str(resnum1lh)+str(resid1lh)
        res2lh=str(resnum2lh)+str(resid2lh)
        closestcontactslh.append([res1lh,res2lh])

        #third, the right-handed closest contacts
        resnum1rh=jsrh[i]+options["-offset"].value
        resnum2rh=ksrh[i]+options["-offset"].value
        resid1rh=TM1CAs[jsrh[i]].resname
        resid2rh=TM2CAs[ksrh[i]].resname
        res1rh=str(resnum1rh)+str(resid1rh)
        res2rh=str(resnum2rh)+str(resid2rh)
        closestcontactsrh.append([res1rh,res2rh])

    #return closest contacts for all, lh and rh
    return closestcontactsoverall,closestcontactslh,closestcontactsrh

'''
contact_matrix_closest_contacts

This function loops over all of the contact matrices (all, lh and rh) and tallies the X closest 
contacts from each of the matrices (X is specified by the user).  It then puts them into an overall 
contact matrix for the overall, lh and rh, and writes it to a file.  It then returns the top X 
contacts for the overall, lh and rh contact matrices.

allcontactmatrices: all of the overall contact matrices
lhcontactmatrices: all of the left-handed contact matrices
rhcontactmatrices: all of the right-handed contact matrices
TM1CAs: residues of TM1
TM2CAs: residues of TM2
'''
def contact_matrix_closest_contacts(allcontactmatrices,lhcontactmatrices,rhcontactmatrices,TM1CAs,TM2CAs):

    #create an overall contact matrix
    overallCM=numpy.zeros(shape=(options["-TM1l"].value,options["-TM2l"].value))
    lhCM=numpy.zeros(shape=(options["-TM1l"].value,options["-TM2l"].value))
    rhCM=numpy.zeros(shape=(options["-TM1l"].value,options["-TM2l"].value))

    #initialize variables for finding the closest options["-topres"].value contacts
    lowestdistall=options["-topres"].value*[None]
    lowestdistlh=options["-topres"].value*[None]
    lowestdistrh=options["-topres"].value*[None]
    jsoverall=options["-topres"].value*[0]
    ksoverall=options["-topres"].value*[0]
    jslh=options["-topres"].value*[0]
    kslh=options["-topres"].value*[0]
    jsrh=options["-topres"].value*[0]
    ksrh=options["-topres"].value*[0]

    #initialize the distances
    for i in range(options["-topres"].value):
        lowestdistall[i]=9999
        lowestdistlh[i]=9999
        lowestdistrh[i]=9999

    #find the X closest contacts for each contact matrix and adds it to the overall contact matrix
    for i in range(len(allcontactmatrices)):
        for j in range(options["-TM1l"].value):
            for k in range(options["-TM2l"].value):
                for l in range(options["-topres"].value):
                    if allcontactmatrices[i][j][k]<lowestdistall[l]:
                        for m in reversed(range(l,options["-topres"].value-1,1)):
                            lowestdistall[m+1]=lowestdistall[m]
                            jsoverall[m+1]=jsoverall[m]
                            ksoverall[m+1]=ksoverall[m]
                        lowestdistall[l]=allcontactmatrices[i][j][k]
                        jsoverall[l]=j
                        ksoverall[l]=k
                        break
        for j in range(options["-TM1l"].value):
            for k in range(options["-TM2l"].value):
                for l in range(options["-topres"].value):
                    if lhcontactmatrices[i][j][k]<lowestdistlh[l]:
                        for m in reversed(range(l,options["-topres"].value-1,1)):
                            lowestdistlh[m+1]=lowestdistlh[m]
                            jslh[m+1]=jslh[m]
                            kslh[m+1]=kslh[m]
                        lowestdistlh[l]=lhcontactmatrices[i][j][k]
                        jslh[l]=j
                        kslh[l]=k
                        break
        for j in range(options["-TM1l"].value):
            for k in range(options["-TM2l"].value):
                for l in range(options["-topres"].value):
                    if rhcontactmatrices[i][j][k]<lowestdistrh[l]:
                        for m in reversed(range(l,options["-topres"].value-1,1)):
                            lowestdistrh[m+1]=lowestdistrh[m]
                            jsrh[m+1]=jsrh[m]
                            ksrh[m+1]=ksrh[m]
                        lowestdistrh[l]=rhcontactmatrices[i][j][k]
                        jsrh[l]=j
                        ksrh[l]=k
                        break
        for j in range(options["-TM1l"].value):
            for k in range(options["-TM2l"].value):
                for l in range(options["-topres"].value):
                    if jsoverall[l]==j and ksoverall[l]==k:
                        overallCM[j][k]+=1
                    if jslh[l]==j and kslh[l]==k:
                        lhCM[j][k]+=1
                    if jsrh[l]==j and ksrh[l]==k:
                        rhCM[j][k]+=1
        for j in range(options["-topres"].value):
            lowestdistall[j]=9999
            lowestdistlh[j]=9999
            lowestdistrh[j]=9999

    #open output files
    allCMfile=open(options["-CMall"].value,"w")
    lhCMfile=open(options["-CMlh"].value,"w")
    rhCMfile=open(options["-CMrh"].value,"w")

    #write out the final contact matrix file
    for i in range(options["-TM1l"].value):
        for j in range(options["-TM2l"].value):
            allCMfile.write(str(overallCM[i][j])+" ")
            lhCMfile.write(str(lhCM[i][j])+" ")
            rhCMfile.write(str(rhCM[i][j])+" ")
        allCMfile.write("\n")
        lhCMfile.write("\n")
        rhCMfile.write("\n")

    #close output files
    allCMfile.close()
    lhCMfile.close()
    rhCMfile.close()

    #reset the counter
    for i in range(options["-topres"].value):
        lowestdistall[i]=0
        lowestdistlh[i]=0
        lowestdistrh[i]=0

    #find overall closest contacts
    for j in range((options["-TM1l"].value)):
        for k in range(options["-TM2l"].value):
            for l in range(options["-topres"].value):
                if overallCM[j][k]>lowestdistall[l]:
                    for m in reversed(range(l,options["-topres"].value-1,1)):
                        lowestdistall[m+1]=lowestdistall[m]
                        jsoverall[m+1]=jsoverall[m]
                        ksoverall[m+1]=ksoverall[m]
                    lowestdistall[l]=overallCM[j][k]
                    jsoverall[l]=j
                    ksoverall[l]=k
                    break
    for j in range((options["-TM1l"].value)):
        for k in range(options["-TM2l"].value):
            for l in range(options["-topres"].value):
                if lhCM[j][k]>lowestdistlh[l]:
                    for m in reversed(range(l,options["-topres"].value-1,1)):
                        lowestdistlh[m+1]=lowestdistlh[m]
                        jslh[m+1]=jslh[m]
                        kslh[m+1]=kslh[m]
                    lowestdistlh[l]=lhCM[j][k]
                    jslh[l]=j
                    kslh[l]=k
                    break
    for j in range((options["-TM1l"].value)):
        for k in range(options["-TM2l"].value):
            for l in range(options["-topres"].value):
                if rhCM[j][k]>lowestdistrh[l]:
                    for m in reversed(range(l,options["-topres"].value-1,1)):
                        lowestdistrh[m+1]=lowestdistrh[m]
                        jsrh[m+1]=jsrh[m]
                        ksrh[m+1]=ksrh[m]
                    lowestdistrh[l]=rhCM[j][k]
                    jsrh[l]=j
                    ksrh[l]=k
                    break

    #initialize output arrays
    closestcontactsoverall=[]
    closestcontactslh=[]
    closestcontactsrh=[]

    #loop over everything
    for i in range(options["-topres"].value):

        #first the overall closest contacts
        resnum1overall=jsoverall[i]+options["-offset"].value
        resnum2overall=ksoverall[i]+options["-offset"].value
        resid1overall=TM1CAs[jsoverall[i]].resname
        resid2overall=TM2CAs[ksoverall[i]].resname
        res1overall=str(resnum1overall)+str(resid1overall)
        res2overall=str(resnum2overall)+str(resid2overall)
        closestcontactsoverall.append([res1overall,res2overall])
        
        #second, the left-handed closest contacts
        resnum1lh=jslh[i]+options["-offset"].value
        resnum2lh=kslh[i]+options["-offset"].value
        resid1lh=TM1CAs[jslh[i]].resname
        resid2lh=TM2CAs[kslh[i]].resname
        res1lh=str(resnum1lh)+str(resid1lh)
        res2lh=str(resnum2lh)+str(resid2lh)
        closestcontactslh.append([res1lh,res2lh])

        #third, the right-handed closest contacts
        resnum1rh=jsrh[i]+options["-offset"].value
        resnum2rh=ksrh[i]+options["-offset"].value
        resid1rh=TM1CAs[jsrh[i]].resname
        resid2rh=TM2CAs[ksrh[i]].resname
        res1rh=str(resnum1rh)+str(resid1rh)
        res2rh=str(resnum2rh)+str(resid2rh)
        closestcontactsrh.append([res1rh,res2rh])

    #return closest contacts for all, lh and rh
    return closestcontactsoverall,closestcontactslh,closestcontactsrh


'''
contact_matrix_all_sims

This function makes the overall contact matrices and overall crossing angle array, to calculate the
Theta(n)'s.  It then returns the allxangle matrix, as well as the overall all, left-handed and 
right-handed contact matrices. 

allcontactmatrices: all of the contact matrices for all of the directories
lhcontactmatrices: all of the left-handed contact matrices for all of the directories
rhcontactmatrices: all of the right-handed contact matrices for all of the directories
arraynumframesconsideredall: array containing the number of frames considered in all matrices for
			     that particular simulation
arraynumframesconsideredlh: array containing the number of left-handed contact marix frames that
			    were considered for each particular simulation
arraynumframesconsideredrh: array containing the number of right-handed contact marix frames that
                            were considered for each particular simulation
xangles: array of crossing angle arrays for each simulation
options["-TM1l"].value: length of TM1
options["-TM2l"].value: length of TM2
'''
def contact_matrix_all_sims(allcontactmatrices,lhcontactmatrices,rhcontactmatrices,arraynumframesconsideredall,arraynumframesconsideredlh,arraynumframesconsideredrh,xangles,contactmatrixindividualframes):

    #initialize array for the whole thing
    allxangle=[]

    #initialize the final matrices
    thetaCMall=numpy.zeros(shape=(options["-TM1l"].value,options["-TM2l"].value))
    thetaCMlh=numpy.zeros(shape=(options["-TM1l"].value,options["-TM2l"].value))
    thetaCMrh=numpy.zeros(shape=(options["-TM1l"].value,options["-TM2l"].value))
    allcontactmatrixindividualframes=[]

    #initialize variables
    numframesconsideredall=0
    numframesconsideredlh=0
    numframesconsideredrh=0

    #create xanglearray for this group
    for i in range(1,options["-numdirs"].value+1):
        for j in xangles[i-1]:
            allxangle.append(float(j))
        for j in contactmatrixindividualframes[i-1]:
            allcontactmatrixindividualframes.append(j)
        if (arraynumframesconsideredall[i-1]>0):
            thetaCMall+=allcontactmatrices[i-1]
            numframesconsideredall+=arraynumframesconsideredall[i-1]
        if (arraynumframesconsideredlh[i-1]>0):
            thetaCMlh+=lhcontactmatrices[i-1]
            numframesconsideredlh+=arraynumframesconsideredlh[i-1]
        if (arraynumframesconsideredrh[i-1]>0):
            thetaCMrh+=rhcontactmatrices[i-1]
            numframesconsideredrh+=arraynumframesconsideredrh[i-1]

    #normalize the contact matrices
    for i in range(len(thetaCMall)):
        for j in range(len(thetaCMall[0])):
            if (numframesconsideredall>0):
                thetaCMall[i][j]=Decimal(thetaCMall[i][j])*(Decimal(1)/Decimal(numframesconsideredall))
            if (numframesconsideredlh>0):
                thetaCMlh[i][j]=Decimal(thetaCMlh[i][j])*(Decimal(1)/Decimal(numframesconsideredlh))
            if (numframesconsideredrh>0):
                thetaCMrh[i][j]=Decimal(thetaCMrh[i][j])*(Decimal(1)/Decimal(numframesconsideredrh))

    #return variables
    return allxangle,thetaCMall,thetaCMlh,thetaCMrh,allcontactmatrixindividualframes,numframesconsideredall,numframesconsideredlh,numframesconsideredrh

'''
getmodestructures

This function looks through each individual frame of the simulations, and determines if it is +/- x 
degrees from one of the modes (x specified by the -degwindow option).  If it is, it is put into the 
group corresponding to that particular mode.  All of the frames are then averaged together, then it 
chooses the top x contacts (x specified by the -topres option).  It write the overall density maps
of each of the modes to a file, and it returns the number of frames that were considered for each 
mode.

thetanlhmode1: mode of left-handed crossing angle with the highest probability
thetanlhmode2: mode of left-handed crossing angle with the second-highest probability
thetanrhmode1: mode of right-handed crossing angle with the highest probability
thetanrhmode2: mode of right-handed crossing angle with the second-highest probability
allxangle: overall crossing angle array of all N frames
contactmatrixindividualframes: all individual frames from all N simulations
'''
def getmodestructures(thetanlhmode1,thetanlhmode2,thetanrhmode1,thetanrhmode2,allxangle,contactmatrixindividualframes):

    #average structure files
    averagethetanlhmode1temp=numpy.zeros(shape=(options["-TM1l"].value,options["-TM2l"].value))
    averagethetanlhmode2temp=numpy.zeros(shape=(options["-TM1l"].value,options["-TM2l"].value))
    averagethetanrhmode1temp=numpy.zeros(shape=(options["-TM1l"].value,options["-TM2l"].value))
    averagethetanrhmode2temp=numpy.zeros(shape=(options["-TM1l"].value,options["-TM2l"].value))

    #counting frames
    lhmode1frames=0
    lhmode2frames=0
    rhmode1frames=0
    rhmode2frames=0

    temp=[]

    #iterate overthe entire xangle array
    for i in range(len(allxangle)):
        if thetanlhmode1!=-99.00:
            if allxangle[i]>=float(thetanlhmode1-options["-degwindow"].value) and allxangle[i]<=float(thetanlhmode1+options["-degwindow"].value):
                averagethetanlhmode1temp+=contactmatrixindividualframes[i]
                lhmode1frames+=1
        if thetanlhmode2!=-99.00:
            if allxangle[i]>=float(thetanlhmode2-options["-degwindow"].value) and allxangle[i]<=float(thetanlhmode2+options["-degwindow"].value):
                averagethetanlhmode2temp+=contactmatrixindividualframes[i]
                lhmode2frames+=1
        if thetanrhmode1!=-99.00:
            if allxangle[i]>=float(thetanrhmode1-options["-degwindow"].value) and allxangle[i]<=float(thetanrhmode1+options["-degwindow"].value):
                averagethetanrhmode1temp+=contactmatrixindividualframes[i]
                rhmode1frames+=1
        if thetanrhmode2!=-99.00:
            if allxangle[i]>=float(thetanrhmode2-options["-degwindow"].value) and allxangle[i]<=float(thetanrhmode2+options["-degwindow"].value):	
                averagethetanrhmode2temp+=contactmatrixindividualframes[i]
                rhmode2frames+=1

    #initialize variables
    thetanlhmode1structure=numpy.zeros(shape=(options["-TM1l"].value,options["-TM2l"].value))
    thetanlhmode2structure=numpy.zeros(shape=(options["-TM1l"].value,options["-TM2l"].value))
    thetanrhmode1structure=numpy.zeros(shape=(options["-TM1l"].value,options["-TM2l"].value))
    thetanrhmode2structure=numpy.zeros(shape=(options["-TM1l"].value,options["-TM2l"].value))

    #loop over and average
    for i in range(options["-TM1l"].value):
        for j in range(options["-TM2l"].value):
            if lhmode1frames!=0:
                thetanlhmode1structure[i][j]=float(averagethetanlhmode1temp[i][j])/float(lhmode1frames)
            if lhmode2frames!=0:
                thetanlhmode2structure[i][j]=float(averagethetanlhmode2temp[i][j])/float(lhmode2frames)
            if rhmode1frames!=0:
                thetanrhmode1structure[i][j]=float(averagethetanrhmode1temp[i][j])/float(rhmode1frames)
            if rhmode2frames!=0:
                thetanrhmode2structure[i][j]=float(averagethetanrhmode2temp[i][j])/float(rhmode2frames)

    #initialize variables
    thetanlhmode1finaltally=numpy.zeros(shape=(options["-TM1l"].value,options["-TM2l"].value))
    thetanlhmode2finaltally=numpy.zeros(shape=(options["-TM1l"].value,options["-TM2l"].value))
    thetanrhmode1finaltally=numpy.zeros(shape=(options["-TM1l"].value,options["-TM2l"].value))
    thetanrhmode2finaltally=numpy.zeros(shape=(options["-TM1l"].value,options["-TM2l"].value))
    lowestdistall=options["-topres"].value*[9999]
    jsoverall=options["-topres"].value*[0]
    ksoverall=options["-topres"].value*[0]

    #find overall closest contacts from average matrix
    for j in range((options["-TM1l"].value)):
        for k in range(options["-TM2l"].value):
            for l in range(options["-topres"].value):
                if thetanlhmode1structure[j][k]<lowestdistall[l]:
                    for m in reversed(range(l,options["-topres"].value-1,1)):
                        lowestdistall[m+1]=lowestdistall[m]
                        jsoverall[m+1]=jsoverall[m]
                        ksoverall[m+1]=ksoverall[m]
                    lowestdistall[l]=thetanlhmode1structure[j][k]
                    jsoverall[l]=j
                    ksoverall[l]=k
                    break
    for j in range(options["-TM1l"].value):
        for k in range(options["-TM2l"].value):
            for l in range(options["-topres"].value):
                if jsoverall[l]==j and ksoverall[l]==k:
                    thetanlhmode1finaltally[j][k]+=1
    for j in range(options["-topres"].value):
        lowestdistall[j]=9999
        jsoverall[j]=0
        ksoverall[j]=0
    for j in range((options["-TM1l"].value)):
        for k in range(options["-TM2l"].value):
            for l in range(options["-topres"].value):
                if thetanlhmode2structure[j][k]<lowestdistall[l]:
                    for m in reversed(range(l,options["-topres"].value-1,1)):
                        lowestdistall[m+1]=lowestdistall[m]
                        jsoverall[m+1]=jsoverall[m]
                        ksoverall[m+1]=ksoverall[m]
                    lowestdistall[l]=thetanlhmode2structure[j][k]
                    jsoverall[l]=j
                    ksoverall[l]=k
                    break
    for j in range(options["-TM1l"].value):
        for k in range(options["-TM2l"].value):
            for l in range(options["-topres"].value):
                if jsoverall[l]==j and ksoverall[l]==k:
                    thetanlhmode2finaltally[j][k]+=1
    for j in range(options["-topres"].value):
        lowestdistall[j]=9999
        jsoverall[j]=0
        ksoverall[j]=0
    for j in range((options["-TM1l"].value)):
        for k in range(options["-TM2l"].value):
            for l in range(options["-topres"].value):
                if thetanrhmode1structure[j][k]<lowestdistall[l]:
                    for m in reversed(range(l,options["-topres"].value-1,1)):
                        lowestdistall[m+1]=lowestdistall[m]
                        jsoverall[m+1]=jsoverall[m]
                        ksoverall[m+1]=ksoverall[m]
                    lowestdistall[l]=thetanrhmode1structure[j][k]
                    jsoverall[l]=j
                    ksoverall[l]=k
                    break
    for j in range(options["-TM1l"].value):
        for k in range(options["-TM2l"].value):
            for l in range(options["-topres"].value):
                if jsoverall[l]==j and ksoverall[l]==k:
                    thetanrhmode1finaltally[j][k]+=1
    for j in range(options["-topres"].value):
        lowestdistall[j]=9999
        jsoverall[j]=0
        ksoverall[j]=0
    for j in range((options["-TM1l"].value)):
        for k in range(options["-TM2l"].value):
            for l in range(options["-topres"].value):
                if thetanrhmode2structure[j][k]<lowestdistall[l]:
                    for m in reversed(range(l,options["-topres"].value-1,1)):
                        lowestdistall[m+1]=lowestdistall[m]
                        jsoverall[m+1]=jsoverall[m]
                        ksoverall[m+1]=ksoverall[m]
                    lowestdistall[l]=thetanrhmode2structure[j][k]
                    jsoverall[l]=j
                    ksoverall[l]=k
                    break
    for j in range(options["-TM1l"].value):
        for k in range(options["-TM2l"].value):
            for l in range(options["-topres"].value):
                if jsoverall[l]==j and ksoverall[l]==k:
                    thetanrhmode2finaltally[j][k]+=1
    for j in range(options["-topres"].value):
        lowestdistall[j]=9999
        jsoverall[j]=0
        ksoverall[j]=0

    #open files
    thetanlhmode1matrix=open("thetanlhmode1matrix.dat","w")
    thetanlhmode2matrix=open("thetanlhmode2matrix.dat","w")
    thetanrhmode1matrix=open("thetanrhmode1matrix.dat","w")
    thetanrhmode2matrix=open("thetanrhmode2matrix.dat","w")

    #write to files
    for i in range(options["-TM1l"].value):
        for j in range(options["-TM2l"].value):
            thetanlhmode1matrix.write(str(thetanlhmode1finaltally[i][j])+" ")
            thetanlhmode2matrix.write(str(thetanlhmode2finaltally[i][j])+" ")
            thetanrhmode1matrix.write(str(thetanrhmode1finaltally[i][j])+" ")
            thetanrhmode2matrix.write(str(thetanrhmode2finaltally[i][j])+" ")
        thetanlhmode1matrix.write("\n")
        thetanlhmode2matrix.write("\n")
        thetanrhmode1matrix.write("\n")
        thetanrhmode2matrix.write("\n")

    #close files
    thetanlhmode1matrix.close()
    thetanlhmode2matrix.close()
    thetanrhmode1matrix.close()
    thetanrhmode2matrix.close()

    #return frames
    return lhmode1frames,lhmode2frames,rhmode1frames,rhmode2frames

'''
calculations

This function calculates probabilities for left and right-handed crossing angles, mean and modes
for right and left-handed crossing angles, and overall, left-handed and right-handed contact matrix.
It writes these to a database file so the RAM doesn't get bogged down.  Function does not return 
anything.

start: starting index of the permutations
finish: last index of the permutations
numberpredicted: the number of total permutations predicted
ssobject: the subsets object that will give the correct permutations
group: number of simulations to be excluded
filename: name of the database file for writing
options["-TM1l"].value: length of TM1
options["-TM2l"].value: length of TM2
allxangle: overall xangle array
allcontactmatricesnumdirin: pointer to array of the overall contact matrices for N simulations (for C code) 
lhcontactmatricesnumdirin: pointer to array of the left-handed contact matrices for N simulations (for C code)
rhcontactmatricesnumdirin: pointer to array of the right-handed contact matrices for N simulations (for C code)
Carraynumframesconsideredall: pointer to array of the number of overall frames considered in each simulation (for C code)
Carraynumframesconsideredlh: pointer to array of the number of left-handed frames considered in each simulation (for C code)
Carraynumframesconsideredrh: pointer to array of the number of right-handed frames considered in each simulation (for C code)
xanglespointer: pointer to array of xangles from each of the N simulations (for C code)
totalforbar: global variable that will count the total permutations done for the bar on the terminal screen
'''
def calculations(start,finish,numberpredicted,ssobject,group,filename,allcontactmatricesnumdirin,lhcontactmatricesnumdirin,rhcontactmatricesnumdirin,Carraynumframesconsideredall,Carraynumframesconsideredlh,Carraynumframesconsideredrh,xanglespointer,totalforbar):

	#open files for PyTables output
    datafile=openFile(filename,mode="w",expectedrows=(finish-start+1))
    grouph5=datafile.createGroup("/",'detector','Detector Information')
    table=datafile.createTable(grouph5,'readout',Data,"Data Readout")
    oneData=table.row
    n=1

	#loop over subset of total permutations
    for i in xrange(start,finish+1,1):

        #generate a random number and check if it is unique
        integer=randint(0,numberpredicted-1)
        boolean=True
        while (boolean):
            if integer in integerarray:
                integer=randint(0,numberpredicted-1)
            else:
                putinarray(integer)
                boolean=False

        #get the subset for the random number
        subset=ssobject.generateSubset(options["-numdirs"].value,group,integer)
        for j in range(len(subset)):
            subset[j]+=1

        #convert subset of simulations from generateSubset to C type
        groupin=(ctypes.c_int * len(subset))(*subset)

        #c code I will run later
        c_code="""

        //initialize the matrices
        float contact_matrix_all[TM1length][TM2length];
        float contact_matrix_left_handed[TM1length][TM2length]; 
        float contact_matrix_right_handed[TM1length][TM2length];
        float contact_matrix_all_temp[TM1length][TM2length];
        float contact_matrix_left_handed_temp[TM1length][TM2length];
        float contact_matrix_right_handed_temp[TM1length][TM2length];
        float allcontactmatricesnumdirtemp[TM1length][TM2length];
        float lhcontactmatricesnumdirtemp[TM1length][TM2length];
        float rhcontactmatricesnumdirtemp[TM1length][TM2length];

        //make sure that they all are initialized at 0
        for (int i=0; i<TM1length; i++) { 
            for (int j=0; j<TM2length; j++) {
                contact_matrix_all[i][j]=0.0;
                contact_matrix_left_handed[i][j]=0.0;
                contact_matrix_right_handed[i][j]=0.0;
                contact_matrix_all_temp[i][j]=0.0;
                contact_matrix_left_handed_temp[i][j]=0.0;
                contact_matrix_right_handed_temp[i][j]=0.0;
                allcontactmatricesnumdirtemp[i][j]=0.0;
                lhcontactmatricesnumdirtemp[i][j]=0.0;
                rhcontactmatricesnumdirtemp[i][j]=0.0; 
            }
        }

        //initialize the frames variables
        int numframesconsideredall=0;
        int numframesconsideredlh=0;
        int numframesconsideredrh=0;

        //initialize total xangles
        int totxangles=0; 

        //loop over all stuff
        py::tuple tempgroup(group);
        for (int i=0; i<group; i++) {
            tempgroup[i]=groupin[i];
        }

        //create the in variables for the find_index function
        py::tuple in(3);
        in[0]=tempgroup;
        in[1]=group;

        //loop over the numdirs
        for (int i=1; i<numdirs+1; i++) {
            in[2]=i;
            int num = (int) find_index.call(in);
            if (num==-1) {
                //first, check if the number of crossing angles for each array is greater than 1
                int temp=Carraynumframesconsideredall[i-1];
                if (Carraynumframesconsideredall[i-1]==1) {
                    totxangles+=1;
                }
                else if (Carraynumframesconsideredall[i-1]>1) {
                    int length=Carraynumframesconsideredall[i-1];
                    totxangles+=length;
                }
                else {
                    totxangles+=0;
                }
                if (Carraynumframesconsideredall[i-1]>0) {
                    for (int j=0; j<TM1length; j++) {
                        for (int k=0; k<TM2length; k++) {
                            float temp=allcontactmatricesnumdirin[i-1][j][k];
                            contact_matrix_all_temp[j][k]+=temp;
                        }
                    }
                    int temp=Carraynumframesconsideredall[i-1];
                    numframesconsideredall+=temp;
                }
                if (Carraynumframesconsideredlh[i-1]>0) {
                    for (int j=0; j<TM1length; j++) {
                        for (int k=0; k<TM2length; k++) {
                            float temp=lhcontactmatricesnumdirin[i-1][j][k];
                            contact_matrix_left_handed_temp[j][k]+=temp;
                        }
                    }
                    int temp=Carraynumframesconsideredlh[i-1];
                    numframesconsideredlh+=temp;
                }
                if (Carraynumframesconsideredrh[i-1]>0) {
                    for (int j=0; j<TM1length; j++) {
                        for (int k=0; k<TM2length; k++) {
                            float temp=rhcontactmatricesnumdirin[i-1][j][k];
                            contact_matrix_right_handed_temp[j][k]+=temp;
                        }
                    }
                    int temp=Carraynumframesconsideredrh[i-1];
                    numframesconsideredrh+=temp;
                }
            }
        }

        //create an array of all crossing angles
        float xanglearray[totxangles];

        //initialize xanglearray
        for (int i=0; i<totxangles; i++) {
            xanglearray[i]=0.0;
        }

        //put the values into one large xanglearray
        int tot=0;
        for (int i=1; i<numdirs+1; i++) {
            in[2]=i;
            int num = (int) find_index.call(in);
            if (num==-1) {
                if (Carraynumframesconsideredall[i-1]==1) {
                    float temp=xanglespointer[i-1][0];
                    xanglearray[tot]+=temp;
                    tot++;
                }
                else if (Carraynumframesconsideredall[i-1]>1) {
                    int tempj=Carraynumframesconsideredall[i-1];
                    for (int j=0; j<tempj; j++) {
                        float temp=xanglespointer[i-1][j];
                        xanglearray[tot]+=temp;
                        tot++;
                    }
                }
            }
        }

        //normalize for the number of total frames considered
        for (int i=0; i<TM1length; i++) {
            for (int j=0; j<TM2length; j++) {
                if (numframesconsideredall>0) {
                    contact_matrix_all[i][j]=contact_matrix_all_temp[i][j]/numframesconsideredall;
                }
                if (numframesconsideredlh>0) {
                    contact_matrix_left_handed[i][j]=contact_matrix_left_handed_temp[i][j]/numframesconsideredlh;
                }
                if (numframesconsideredrh>0) {
                    contact_matrix_right_handed[i][j]=contact_matrix_right_handed_temp[i][j]/numframesconsideredrh;
                }
            }
        }

        //initialize python tuple for returning values
        py::tuple xangles(totxangles);
        for (int i=0; i<totxangles; i++) {
            xangles[i]=xanglearray[i];
        }

        //create python lists for returning matrices
        py::list allCM(TM1length);
        py::list lhCM(TM1length);
        py::list rhCM(TM1length);
        for (int i=0; i<TM1length; i++) {
            py::list allCMtemplist(TM2length);
            py::list lhCMtemplist(TM2length);
            py::list rhCMtemplist(TM2length);
            for (int j=0; j<TM2length; j++) {
                allCMtemplist[j]=contact_matrix_all[i][j];
                lhCMtemplist[j]=contact_matrix_left_handed[i][j];
                rhCMtemplist[j]=contact_matrix_right_handed[i][j];
            }
            allCM[i]=allCMtemplist;
            lhCM[i]=lhCMtemplist;
            rhCM[i]=rhCMtemplist;
        }

        //create some py::tuple objects for frames
        py::tuple allframes(1);
        allframes[0]=numframesconsideredall;
        py::tuple lhframes(1);
        lhframes[0]=numframesconsideredlh;
        py::tuple rhframes(1);
        rhframes[0]=numframesconsideredrh;

        //create the return variable list
        py::list matricesxangleandframes(7);
        matricesxangleandframes[0]=allCM;
        matricesxangleandframes[1]=lhCM;
        matricesxangleandframes[2]=rhCM;
        matricesxangleandframes[3]=xangles;
        matricesxangleandframes[4]=allframes;
        matricesxangleandframes[5]=lhframes;
        matricesxangleandframes[6]=rhframes;

        //return the variables
        return_val=matricesxangleandframes;
        """

        #set values for C code
        TM1length=options["-TM1l"].value
        TM2length=options["-TM2l"].value
        numdirs=options["-numdirs"].value
        bound1=options["-bounds"].value[0]
        bound2=options["-bounds"].value[1]
        binsize=options["-binsize"].value

        #run C code
        contact_matrix_all,contact_matrix_left_handed,contact_matrix_right_handed,xanglearray,numframesconsideredall,numframesconsideredlh,numframesconsideredrh=weave.inline(c_code,['numdirs','group','Carraynumframesconsideredall','Carraynumframesconsideredlh','Carraynumframesconsideredrh','find_index','allcontactmatricesnumdirin','lhcontactmatricesnumdirin','rhcontactmatricesnumdirin','xanglespointer','bound1','bound2','binsize','TM1length','TM2length','groupin'])

        #get the number of frames for all frames, lhframes, and rhframes
        numframesconsideredall=int(numframesconsideredall[0])
        numframesconsideredlh=int(numframesconsideredlh[0])
        numframesconsideredrh=int(numframesconsideredrh[0])

        #test for values
        if len(xanglearray)>0:

            #create xangle matrix for the mode
            R=((math.fabs(options["-bounds"].value[0])+math.fabs(options["-bounds"].value[1]))/options["-binsize"].value)+1
            xanglematrix=put_crossing_angle(xanglearray,R)

            #run the mode calculations - figure out if I should do something with the lh and rh probs
            lhmodes,rhmodes,lhprobs,rhprobs=mode(xanglematrix)

            #calculate all of the probabilities of the crossing angle
            if numframesconsideredlh>0:
                thetalh=float(numframesconsideredlh)/float(numframesconsideredall)
            else:
                thetalh=0
            if numframesconsideredrh>0:
                thetarh=float(numframesconsideredrh)/float(numframesconsideredall)
            else:
                thetarh=0
            thetalhmean,thetarhmean=splitandmean(xanglearray)

        else:

            #set values to -99.00 so functions can exclude them
            thetalh=-99.00
            thetarh=-99.00
            thetalhmean=-99.00
            thetarhmean=-99.00
            lhmodes=[-99.00]
            rhmodes=[-99.00]

        #store data into the PyTables Object
        oneData['lhprob']=thetalh
        oneData['rhprob']=thetarh
        oneData['lhmean']=thetalhmean
        oneData['rhmean']=thetarhmean
        if len(lhmodes)==1:
            oneData['lhmode1']=lhmodes[0]
            oneData['lhmode2']=-99.00
        elif len(lhmodes)==2:
            oneData['lhmode1']=lhmodes[0]
            oneData['lhmode2']=lhmodes[1]
        if len(rhmodes)==1:
            oneData['rhmode1']=rhmodes[0]
            oneData['rhmode2']=-99.00
        elif len(rhmodes)==2:
            oneData['rhmode1']=rhmodes[0]
            oneData['rhmode2']=rhmodes[1]
        oneData['CMall']=contact_matrix_all
        oneData['CMlh']=contact_matrix_left_handed
        oneData['CMrh']=contact_matrix_right_handed
        oneData.append()

        #check if it is time to flush the data into the file
        if ((i+1)/n)==100:
            table.flush()
        elif (i==finish):
            table.flush()

        #check to see about the number of subsets
        value=0
        if (options["-numsubsets"].value==None or options["-numsubsets"].value > numberpredicted):
            value=numberpredicted	
        else:
            value=options["-numsubsets"].value

        #update the progress bar
        if totalforbar.value==value:
            update_progress(1)
            print "\n"
        else:
            modify()
            if totalforbar.value<value:	
                update_progress(float(totalforbar.value)/value)

    #close the datafile when the loop is done
    datafile.close()

'''
calculations_multiprocessing

This function creates C pointers to pass to calculations, and creates X processes to calculate all
possible permutations for this (X specified by user, or default uses all cores on computer).  It 
also calculates the jackknife value from all these permutations, writes these things to a file, and 
creates png files of the contact matrices.  Function does not return anything.

group: number of simulations to be excluded
options["-TM1l"].value: length of TM1
options["-TM2l"].value: length of TM2
TM1reslist: list of residues in TM1
TM2reslist: list of residues in TM2
allxangle: overall xangle array
thetanCMall: overall contact matrix for total N simulations
thetanCMlh: left-handed contact matrix for total N simulations
thetanCMrh: right-handed contact matrix for total N simulations
thetanlh: left-handed probability for total N simulations
thetanrh: right-handed probability for total N simulations
thetanlhmean: left-handed mean of the crossing angle for N simulations
thetanrhmean: right-handed mean of the crossing angle for N simulations
thetanlhmode: left-handed mode of the crossing angle for N simulations
thetanrhmode: right-handed mode of the crossing angle for N simulations
allcontactmatricesnumdir: array of the overall contact matrices for N simulations
lhcontactmatricesnumdir: array of the left-handed contact matrices for N simulations
rhcontactmatricesnumdir: array of the right-handed contact matrices for N simulations
arraynumframesconsideredall: array of the number of overall frames considered in each simulation
arraynumframesconsideredlh: array of the number of left-handed frames considered in each simulation
arraynumframesconsideredrh: array of the number of right-handed frames considered in each simulation
xangles: array of xangles from each of the N simulations
TM1CAs: array of CAs from TM1
TM2CAs: array of CAs from TM2
quartilefileLHprob: file containing data for the boxplots for the LH probability
quartilefileRHprob: file containing data for the boxplots for the RH probability
quartilefileLHmean: file containing data for the boxplots for the LH mean 
quartilefileRHmean: file containing data for the boxplots for the RH mean
quartilefileLHmode: file containing data for the boxplots for the LH mode
quartilefileRHmode: file containing data for the boxplots for the RH mode
'''
def calculations_multiprocessing(group,TM1reslist,TM2reslist,allxangle,thetanCMall,thetanCMlh,thetanCMrh,thetanlh,thetanrh,thetanlhmean,thetanrhmean,thetanlhmode1,thetanrhmode1,thetanlhmode2,thetanrhmode2,allcontactmatricesnumdir,lhcontactmatricesnumdir,rhcontactmatricesnumdir,arraynumframesconsideredall,arraynumframesconsideredlh,arraynumframesconsideredrh,xangles,TM1CAs,TM2CAs,quartilefileLHprob,quartilefileRHprob,quartilefileLHmean,quartilefileRHmean,quartilefileLHmode1,quartilefileRHmode1,quartilefileLHmode2,quartilefileRHmode2):

    #predict number of groups in chosen set
    numberpredicted=factorial(options["-numdirs"].value)/(factorial(options["-numdirs"].value-group)*factorial(group))

    #print information about permutation
    print "groupsize                    = "+str(options["-numdirs"].value-group)
    print "predicted permutations       = "+str(numberpredicted)
    if (options["-numsubsets"].value==None):
        print "number selected permutations = "+str(numberpredicted)
    elif (options["-numsubsets"].value>numberpredicted):
        print "number selected permutations = "+str(options["-numsubsets"].value)
        print 
        print "WARNING: the number of selected permutations you have chosen is more than the number of possible unique permutations - will only use number of predicted permutations\n"
    else:
        print "number selected permutations = "+str(options["-numsubsets"].value)
    print

    #commands to make jackknifing subdirectories, and then move all of the files there for less clutter
    if os.path.isdir("%ichoose%i" % (options["-numdirs"].value,options["-numdirs"].value-group)):
        print "The directory %ichoose%i already exists - will rename all old files with the extension .bak\n" % (options["-numdirs"].value,options["-numdirs"].value-group)
        for path,dirs,filenames in os.walk("%ichoose%i" % (options["-numdirs"].value,options["-numdirs"].value-group)):
            os.chdir("%ichoose%i/" % (options["-numdirs"].value,options["-numdirs"].value-group))
            for filename in filenames:
                if filename[:4]==".nfs":
                    pass
                else:
                    os.system("mv %s %s.bak" % (filename,filename))
    else:
        os.system("mkdir %ichoose%i" % (options["-numdirs"].value,options["-numdirs"].value-group))
        os.chdir("%ichoose%i/" % (options["-numdirs"].value,options["-numdirs"].value-group))

    #initializing arrays to be fed into C
    allcontactmatricesnumdirin=[]
    lhcontactmatricesnumdirin=[]
    rhcontactmatricesnumdirin=[]

    #change all variables to C type variables	
    for i in range(len(allcontactmatricesnumdir)):
        matrix=[array.ctypes.data_as(ctypes.POINTER(ctypes.c_float)) for array in allcontactmatricesnumdir[i]]
        allcontactmatricesnumdirin.append(matrix)
    for i in range(len(lhcontactmatricesnumdir)):
        matrix=[array.ctypes.data_as(ctypes.POINTER(ctypes.c_float)) for array in lhcontactmatricesnumdir[i]]
        lhcontactmatricesnumdirin.append(matrix)
    for i in range(len(rhcontactmatricesnumdir)):
        matrix=[array.ctypes.data_as(ctypes.POINTER(ctypes.c_float)) for array in rhcontactmatricesnumdir[i]]
        rhcontactmatricesnumdirin.append(matrix)

    #create the C pointers for the number of frames arrays
    Carraynumframesconsideredall=(ctypes.c_int * len(arraynumframesconsideredall))(*arraynumframesconsideredall)
    Carraynumframesconsideredlh=(ctypes.c_int * len(arraynumframesconsideredlh))(*arraynumframesconsideredlh)
    Carraynumframesconsideredrh=(ctypes.c_int * len(arraynumframesconsideredrh))(*arraynumframesconsideredrh)

    #create xanglespointer
    xanglespointer=[(ctypes.c_float * len(array)) (*array) for array in xangles]

    #create the subsets object to get all the permutations
    ssobject=subsets.EnumeratedSubsets()

    #make sure the division of labour is correct
    filenames=[]
    starts=[]
    ends=[]
    if (options["-numsubsets"].value==None or options["-numsubsets"].value>numberpredicted):
        dnum=numberpredicted/options["-numcores"].value
        remainder=numberpredicted%options["-numcores"].value
    else:
        dnum=options["-numsubsets"].value/options["-numcores"].value
        remainder=options["-numsubsets"].value%options["-numcores"].value
    start=0
    end=dnum

    #create and start all of the processes
    for i in range(options["-numcores"].value):
        filenames.append("datafile%i.h5" % (i+1))
        starts.append(start)
        ends.append(end)
        start=end+1
        end+=dnum
        if i==options["-numcores"].value-2:
            end=end+remainder-1

    #import a manager for the processes
    manager=multiprocessing.Manager()

    #run calculations
    jobs=[]
    for i in range(options["-numcores"].value):
        process=multiprocessing.Process(target=calculations,args=(starts[i],ends[i],numberpredicted,ssobject,group,filenames[i],allcontactmatricesnumdirin,lhcontactmatricesnumdirin,rhcontactmatricesnumdirin,Carraynumframesconsideredall,Carraynumframesconsideredlh,Carraynumframesconsideredrh,xanglespointer,totalforbar))
        jobs.append(process)
        process.start()

    #don't want control flow to continue in the parent until all of the children are finished
    for job in jobs:
        job.join()

    #restart the count for the next group
    restart()

    #reset the integerarray for the next group
    reset()

    #create a link to the file containing all of the data
    tables=[]
    datafiles=[]
    for i in range(options["-numcores"].value):
        datafile=openFile(filenames[i],mode="r")
        datafiles.append(datafile)
        tables.append(datafile.root.detector.readout)

    #retrieve all of the data
    if thetanlh!=-99.00:
        lhprobthetas=[]
        for table in tables:
            for x in table.iterrows():
                lhprobthetas.append(x['lhprob'])
    if thetanrh!=-99.00:
        rhprobthetas=[]
        for table in tables:
            for x in table.iterrows():
                rhprobthetas.append(x['rhprob'])
    if thetanlhmean!=-99.00:
        lhmeanthetas=[]
        for table in tables:
            for x in table.iterrows():
                lhmeanthetas.append(x['lhmean'])
    if thetanrhmean!=-99.00:
        rhmeanthetas=[]
        for table in tables:
            for x in table.iterrows():
                rhmeanthetas.append(x['rhmean'])
    if thetanlhmode1!=-99.00:
        lhmode1thetas=[]
        for table in tables:
            for x in table.iterrows():
                lhmode1thetas.append(x['lhmode1'])
    if thetanrhmode1!=-99.00:
        rhmode1thetas=[]
        for table in tables:
            for x in table.iterrows():
                rhmode1thetas.append(x['rhmode1'])
    if thetanlhmode2!=-99.00:
        lhmode2thetas=[]
        for table in tables:
            for x in table.iterrows():
                lhmode2thetas.append(x['lhmode2'])
    if thetanrhmode2!=-99.00:
        rhmode2thetas=[]
        for table in tables:
            for x in table.iterrows():
                rhmode2thetas.append(x['rhmode2'])

	#get the quartiles and write to their respective file
    minLHprob,firstquartileLHprob,medianLHprob,thirdquartileLHprob,maxLHprob=boxplotstats(lhprobthetas)
    quartilefileLHprob.write(str(options["-numdirs"].value-group)+" "+str(minLHprob)+" "+str(firstquartileLHprob)+" "+str(medianLHprob)+" "+str(thirdquartileLHprob)+" "+str(maxLHprob)+"\n")
    minRHprob,firstquartileRHprob,medianRHprob,thirdquartileRHprob,maxRHprob=boxplotstats(rhprobthetas)
    quartilefileRHprob.write(str(options["-numdirs"].value-group)+" "+str(minRHprob)+" "+str(firstquartileRHprob)+" "+str(medianRHprob)+" "+str(thirdquartileRHprob)+" "+str(maxRHprob)+"\n")
    minLHmean,firstquartileLHmean,medianLHmean,thirdquartileLHmean,maxLHmean=boxplotstats(lhmeanthetas)
    quartilefileLHmean.write(str(options["-numdirs"].value-group)+" "+str(minLHmean)+" "+str(firstquartileLHmean)+" "+str(medianLHmean)+" "+str(thirdquartileLHmean)+" "+str(maxLHmean)+"\n")
    minRHmean,firstquartileRHmean,medianRHmean,thirdquartileRHmean,maxRHmean=boxplotstats(rhmeanthetas)
    quartilefileRHmean.write(str(options["-numdirs"].value-group)+" "+str(minRHmean)+" "+str(firstquartileRHmean)+" "+str(medianRHmean)+" "+str(thirdquartileRHmean)+" "+str(maxRHmean)+"\n")
    if thetanlhmode1!=-99.00:
        minLHmode1,firstquartileLHmode1,medianLHmode1,thirdquartileLHmode1,maxLHmode1=boxplotstats(lhmode1thetas)
        quartilefileLHmode1.write(str(options["-numdirs"].value-group)+" "+str(minLHmode1)+" "+str(firstquartileLHmode1)+" "+str(medianLHmode1)+" "+str(thirdquartileLHmode1)+" "+str(maxLHmode1)+"\n")
    if thetanrhmode1!=-99.00:
        minRHmode1,firstquartileRHmode1,medianRHmode1,thirdquartileRHmode1,maxRHmode1=boxplotstats(rhmode1thetas)
        quartilefileRHmode1.write(str(options["-numdirs"].value-group)+" "+str(minRHmode1)+" "+str(firstquartileRHmode1)+" "+str(medianRHmode1)+" "+str(thirdquartileRHmode1)+" "+str(maxRHmode1)+"\n")
    if thetanlhmode2!=-99.00:
        minLHmode2,firstquartileLHmode2,medianLHmode2,thirdquartileLHmode2,maxLHmode2=boxplotstats(lhmode2thetas)
        quartilefileLHmode2.write(str(options["-numdirs"].value-group)+" "+str(minLHmode2)+" "+str(firstquartileLHmode2)+" "+str(medianLHmode2)+" "+str(thirdquartileLHmode2)+" "+str(maxLHmode2)+"\n")
    if thetanrhmode2!=-99.00:
        minRHmode2,firstquartileRHmode2,medianRHmode2,thirdquartileRHmode2,maxRHmode2=boxplotstats(rhmode2thetas)
        quartilefileRHmode2.write(str(options["-numdirs"].value-group)+" "+str(minRHmode2)+" "+str(firstquartileRHmode2)+" "+str(medianRHmode2)+" "+str(thirdquartileRHmode2)+" "+str(maxRHmode2)+"\n")

    #retrieve all of the contact matrices
    allcontactmatrices=[]
    for table in tables:
        for x in table.iterrows():
            allcontactmatrices.append(x['CMall'])
    lhcontactmatrices=[]
    for table in tables:
        for x in table.iterrows():
            lhcontactmatrices.append(x['CMlh'])
    rhcontactmatrices=[]
    for table in tables:
        for x in table.iterrows():
            rhcontactmatrices.append(x['CMrh'])

	#write the contact matrices and get closest contacts
    closestcontactsoverall,closestcontactslh,closestcontactsrh=contact_matrix_closest_contacts(allcontactmatrices,lhcontactmatrices,rhcontactmatrices,TM1CAs,TM2CAs)

	#print contacts to a file
    contactsfile=open(options["-outtopres"].value,"w")
    contactsfile.write("Top contacts for %i simulations\n\n" % (options["-numdirs"].value))
    contactsfile.write("Top overall contacts:\n\n")
    for i in range(options["-topres"].value):
        contactsfile.write("\t%s-%s\n" % (closestcontactsoverall[i][0],closestcontactsoverall[i][1]))
    contactsfile.write("\nTop left-handed contacts:\n\n")
    for i in range(options["-topres"].value):
        contactsfile.write("\t%s-%s\n" % (closestcontactslh[i][0],closestcontactslh[i][1]))
    contactsfile.write("\nTop right-handed contacts:\n\n")
    for i in range(options["-topres"].value):
        contactsfile.write("\t%s-%s\n" % (closestcontactsrh[i][0],closestcontactsrh[i][1]))
    contactsfile.close()

    #set the pallete for the contact matrix
    if (options["-numsubsets"].value==None):
        var1=10*(numberpredicted/100)
        var2=20*(numberpredicted/100)
        var3=30*(numberpredicted/100)
    else:
        var1=10*(options["-numsubsets"].value/100)
        var2=20*(options["-numsubsets"].value/100)
        var3=30*(options["-numsubsets"].value/100)
    palette="0 \\\"yellow\\\", "+str(var1)+" \\\"purple\\\", "+str(var2)+" \\\"red\\\", "+str(var3)+" \\\"black\\\""
    if (options["-numsubsets"].value==None):
        tempcbrange=50*(float(numberpredicted)/float(100))
    else:
        tempcbrange=50*(float(options["-numsubsets"].value)/float(100))
    cbrange="0:"+str(tempcbrange)

    #make the gnuplot script, then run gnuplot on all of the contact matrix files to make all of the plots
    os.system("cp /sansom/sb8/bioc1030/scripts/make_contact_matrix.gnu tmp.gnu")
    os.system("sed -e \"s|XTICS|"+str(TM1reslist)+"|g;s|YTICS|"+str(TM2reslist)+"|g;s|PALETTE|"+str(palette)+"|g;s|CBRANGE|"+str(cbrange)+"|g\" tmp.gnu > make_contact_matrix.gnu")
    os.system("rm tmp.gnu")
    subsetsize=options["-numdirs"].value-group
    if (options["-numsubsets"].value==None):
        os.system("/usr/bin/gnuplot -e \"f1='"+str(options["-CMall"].value)+"';f2='"+str(options["-CMlh"].value)+"';f3='"+str(options["-CMrh"].value)+"';numdirs='"+str(options["-numdirs"].value)+"';subset='"+str(subsetsize)+"';subsets='All';max='"+str(options["-TM1l"].value+1)+"'\" make_contact_matrix.gnu")
    else:
        os.system("/usr/bin/gnuplot -e \"f1='"+str(options["-CMall"].value)+"';f2='"+str(options["-CMlh"].value)+"';f3='"+str(options["-CMrh"].value)+"';numdirs='"+str(options["-numdirs"].value)+"';subset='"+str(subsetsize)+"';subsets='"+str(options["-numsubsets"].value)+"';max='"+str(options["-TM1l"].value+1)+"'\" make_contact_matrix.gnu")

    #close the data file, remove it, and change into another directory, in case there is another group
    for datafile in datafiles:
        datafile.close()	
    os.system("rm *.h5")
    os.chdir("../")

'''
Option class
 
borrowed from martinize.py - this gives each Option object a type that the function is, e.g. 
boolean or string), a number (the number of arguments expected for the option), the value 
of the option (can be a default, or specified by the user), and the description of the 
option (what the option is for, how to use it, when is it best used, etc.)
'''
class Option:
    def __init__(self,func=str,num=1,default=None,description=""):
        self.func        = func
        self.num         = num
        self.value       = default
        self.description = description
    def __nonzero__(self):
        if self.func == bool:
            return self.value != False
        return bool(self.value)
    def __str__(self):
        return self.value and str(self.value) or ""
    def setvalue(self,v):
        if len(v) == 1:
            self.value = self.func(v[0])
        else:
            self.value = [ self.func(i) for i in v ]

'''
options

Structure borrowed from martinize.py - list of all the options names, as well as their respective 
option objects (with type, number of arguments, default values, and explanations).
'''
options=[
    ("-g",          Option(str,  1,None,                 "Input file (gro)")),
    ("-x",          Option(str,  1,None,                 "Input file (xtc)")),
    ("-s",          Option(int,  1,1,                    "Number of frames to skip (default=1)")),
    ("-ts",         Option(int,  1,200,                  "Timestep used in the simulation (default=200)")),
    ("-res",        Option(int,  4,None,                 "Selections for helices - takes 4 arguments:\n\t\t\t\t       arg1: first residue of first helix\n\t\t\t\t             arg2: last residue of first helix\n\t\t\t\t             arg3: first residue last helix\n\t\t\t\t             arg4: last residue of second helix")),
    ("-bounds",     Option(float,2,None,                 "Lower and upper bound for crossing angle \n\t\t\t\t       takes 2 arguments:\n\t\t\t\t             arg1: lower bound xangle\n\t\t\t\t             arg2: upper bound xangle")),
    ("-binsize",    Option(float,1,1,                    "Number of bins for the crossing angle \n\t\t\t\t       probability (default=1)")),
    ("-numdirs",    Option(int,  1,5,                    "Number of directories you want to loop over")),
    ("-ca",         Option(str,  1,"CA",                 "Selections for backbone residues \n\t\t\t\t             BOND: use CA (default)\n\t\t\t\t             MARTINI2.1: use B*\n\t\t\t\t             MARTINI2.2: use BB")),
    ("-totaltime",  Option(int,  1,1000000,              "Length of simulation in ps (default=1000000)")),
    ("-offset",     Option(int,  1,1,                    "Residue offset you want to use (default=1)")),
    ("-eigen",      Option(bool, 0,False,                "Use the eigenvector option when calculating the \n\t\t\t\t       vector through helices")),
    ("-mindist",    Option(bool, 0,False,                "Option for calculating the distance between \n\t\t\t\t       helices using the minimum distance, not \n\t\t\t\t       the distance between the COM of the helices")),
    ("-TM1l",       Option(int,  1,33,		     "Length of TM1 (default=33)")),
    ("-TM2l",       Option(int,  1,33,                   "Length of TM2 (default=33)")),
    ("-dstart",     Option(int,  1,None,                 "d that you want to start with (i.e. if you have \n\t\t\t\t       50 simulations, and you wanted to only \n\t\t\t\t       choose 1 simulation to star with, then your \n\t\t\t\t       first d would be 49)")),
    ("-dend",       Option(int,  1,None,                 "d that you want to end with (i.e. if you have \n\t\t\t\t       50 simulations, and you want to only \n\t\t\t\t       choose (n-1) simulations, your last d would be 1)")),
    ("-inc",        Option(int,  1,1,                    "increment for the groups (default=1)")),
    ("-degwindow",  Option(float,1,1.0,                  "degree window for mode analysis (default=1.0)")),
    ("-numsubsets", Option(int,  1,1000,                 "Number of random subsets to use in calculations \n\t\t\t\t       (default=1000)")),
    ("-numcores",   Option(int,  1,None,                 "Number of cores you want to use (default is \n\t\t\t\t       using all cores on your machine)")),
    ("-topres",     Option(int,  1,10,                   "Number of residues to be written out to the \n\t\t\t\t       output file (default=6)")),
    ("-cutoff",     Option(int,  1,10,                   "Cutoff for determining the crossing angle of the helix")),
    ("-outtopres",  Option(str,  1,"topcontacts.dat",    "File that contains the top X residues of the \n\t\t\t\t       jackknife calculation (X set by the options \"-topres\")")),
    ("-thetans",    Option(str,  1,"thetans.dat",        "File that has the overall Theta(n)'s for the  \n\t\t\t\t       simulations")),
    ("-CMall",      Option(str,  1,"CMoverall.dat",      "Name of the file with the overall contact matrix \n\t\t\t\t       data for each jackknife run")),
    ("-CMlh",       Option(str,  1,"CMlh.dat",           "Name of the file with the left-handed contact \n\t\t\t\t       matrix data for each jackknife run")),
    ("-CMrh",       Option(str,  1,"CMrh.dat",           "Name of the file with the right-handed contact \n\t\t\t\t       matrix data for each jackknife run")),
    ("-h",          Option(bool, 0,False,                "Display this help screen and exit"))
]

'''
help

prints the help screen for this script and exits.  Will list the option tag, the type it is/takes,
the number of arguments it takes, and the description of what it is/what it does.  
'''
def help():
    print
    print "\t\t\t      ==> jackknife.py <=="
    print
    print "\t\t\t     Written by Amanda Buyan "
    print
    print "\tDESCRIPTION"
    print "\t--------------------------------------------------------------------"
    print
    print "\t jackknife.py calculations the crossing angles and the mean contact "
    print "\t   matrix for each simulation, and for each d provided, calculates "
    print "\t        the following for each randomly selected subset:"
    print
    print "\t 	    1. Probability of LH crossing angle "
    print "\t 	    2. Probability of RH crossing angle "
    print "\t 	    3. Mean LH crossing angle "
    print "\t 	    4. Mean RH crossing angle "
    print "\t 	    5. LH crossing angle mode "
    print "\t 	    6. RH crossing angle mode "
    print "\t 	    7. Overall Contact Matrix "
    print "\t       8. Contacts associated with LH crossing angle "
    print "\t       9. Contacts associated with RH crossing angle "
    print
    print "\t  It does this for the list of subsets specified by the user, and "
    print "\t creates the contact matrix for each one.  After looping over all "
    print "\t  specified subsets, it creates png files of the contact matrix. "
    print "\t After doing this for all specified d's, it then creates separate "
    print "\t  boxplots of 1-6, and writes out the overall statistics of the N "
    print "\t              simulations in a separate directory."
    print
    print "\t    The directory will be called \"Nsimswithnsubsets\", with the "
    print "\t     subdirectories named \"Nchoose(n-d)\" for every d specified, "
    print "\t    \"thetans\" for the overall statistics of N simulations, and "
    print "\t     \"overall_data_and_plots\" where all the boxplots will be."
    print
    print "\t     Option   Type  Num Args   Description"
    print "\t--------------------------------------------------------------------"
    print
    for item in options:
            if type(item) == str:
                    print item
    for item in options:
            if type(item) != str:
                    print "\t%11s  %5s     %s       %s"%(item[0],item[1].func.__name__,item[1].num,item[1].description)
    print
    sys.exit()

'''
option_parser

creates a library of options, which the function then returns for use by the program.

args: arguments from the options
options: options library that is returned
'''
def option_parser(args,options):

    # Check whether there is a request for help
    if '-h' in args or '--help' in args:
        help()

    # Convert the option list to a dictionary, discarding all comments
    options = dict([i for i in options if not type(i) == str])

    # get the arguments in the args
    while args:
        ar = args.pop(0)
        if ar=="-eigen" or ar=="-pick" or ar=="-mindist":
            options[ar].setvalue([True])
        else:
            options[ar].setvalue([args.pop(0) for i in range(options[ar].num)])

	#return options
    return options

'''
update_progress

This function prints a bar tracking your progress through the permutations, so it doesn't look like
it is hanging.  

progress: a number denoting the percentage of progress that has happened (between 0 and 1)
'''
def update_progress(progress):
        barLength=50 
        status=""
        if isinstance(progress,int):
            progress=float(progress)
        if not isinstance(progress,float):
            progress=0
            status="error: progress var must be float\r\n"
        if progress<0:
            progress=0
            status="Halt...\r\n"
        if progress >=1:
            progress=1
            status="Done...\r\n"
        block=int(round(barLength*progress))
        text="\rPercent: [{0}] {1:.1f}% {2}".format("#"*(block)+" "*((barLength-block)),progress*100,status) #{1}%
        sys.stdout.write(text)
        sys.stdout.flush()

'''
declare global variables and check variables
'''
args=sys.argv[1:]
options=option_parser(args,options)
getcontext().prec = 25
cores=checknumcores()
if (options["-numcores"].value==None):
	options["-numcores"].value=cores
totalforbar=Value('i',1)
integerarray=Array('i',range(options["-numsubsets"].value))
for i in range(len(integerarray)):
	integerarray[i]=-1

'''
modify()

Function from 

http://stackoverflow.com/questions/10588317/python-function-global-variables

Modifies the variable totalforbar by increasing it by 1
'''
def modify():
	global totalforbar
	totalforbar.value+=1

'''
restart()

Function that resets the totalforbar value to 1 for the next group
'''
def restart():
	global totalforbar
	totalforbar.value=1

'''
putinarray

Function that puts a value into the global integerarray

value: value to go into the integerarray
'''
def putinarray(value):
	global integerarray
	for i in range(len(integerarray)):
		if integerarray[i]==-1:
			integerarray[i]=value
			break

'''
reset()

Function that resets the integerarray to []
'''
def reset():
	global integerarray
	for i in range(len(integerarray)):
		integerarray[i]=-1

'''
Data

A class object holding all of the data that will be generated from the crossing angle and contact
matrix calculations for each permutation.  This is from the PyTables extension, and how to use the
PyTables extension is 
'''
class Data(IsDescription):
	lhprob=Float64Col()
	rhprob=Float64Col()
	lhmean=Float64Col()
	rhmean=Float64Col()
	lhmode1=Float64Col()
	rhmode1=Float64Col()
	lhmode2=Float64Col()
	rhmode2=Float64Col()
	CMall=Float64Col(shape=(options["-TM1l"].value,options["-TM2l"].value))
	CMlh=Float64Col(shape=(options["-TM1l"].value,options["-TM2l"].value))
	CMrh=Float64Col(shape=(options["-TM1l"].value,options["-TM2l"].value))

'''
lhmode=Float64Col()
rhmode=Float64Col()
'''

'''
main

This function does the following:

	1. Reads all trajectories, and analyses them for crossing angles and contacts
	2. Calculates overall statistics and writes them to its own separate directory
	3. Calculates the number of unique permutations of N samples and group size d
	4. Chooseis a subset of all the unique permutations, and calculate overall crossing angle
	   LH and RH probability, mean and mode, as well as the overall contacts for each chosen
	   permutation
	5. Iterate over all d's specified by the user
	6. Generate graphs for overall statistics

This function takes no arguments, and returns none. 
'''
def main():

    #beginning message
    print "\nReading in all the individual trajectories...\n"

    #process each individual simulation first to get all of the data
    TM1CAs,TM2CAs,allcontactmatricesnumdir,lhcontactmatricesnumdir,rhcontactmatricesnumdir,arraynumframesconsideredall,arraynumframesconsideredlh,arraynumframesconsideredrh,xangles,contactmatricesindividualframes=process_individual_trajectories()

    #restart value for totalforbar for the permutations
    restart()

    #get reslist for both TMs
    TM1reslist=createreslist(TM1CAs)
    TM2reslist=createreslist(TM2CAs)

    #calculate overall contact matrices and crossing angle array 	
    allxangle,thetanCMall,thetanCMlh,thetanCMrh,allcontactmatricesindividualframes,numframesconsideredall,numframesconsideredlh,numframesconsideredrh=contact_matrix_all_sims(allcontactmatricesnumdir,lhcontactmatricesnumdir,rhcontactmatricesnumdir,arraynumframesconsideredall,arraynumframesconsideredlh,arraynumframesconsideredrh,xangles,contactmatricesindividualframes)

    #create overall xanglematrix
    R=((math.fabs(options["-bounds"].value[0])+math.fabs(options["-bounds"].value[1]))/options["-binsize"].value)+1
    xanglematrix=put_crossing_angle(allxangle,R)

    #calculate thetan for probability, mean and mode
    thetanlh=float(numframesconsideredlh)/float(numframesconsideredall)
    thetanrh=float(numframesconsideredrh)/float(numframesconsideredall)
    thetanlhmean,thetanrhmean=splitandmean(allxangle)
    thetanlhmodes,thetanrhmodes,thetanlhprobs,thetanrhprobs=mode(xanglematrix)
    if len(thetanlhmodes)==1:
        thetanlhmode1=thetanlhmodes[0]
        thetanlhprob1=thetanlhprobs[0]
        thetanlhmode2=-99.00
        thetanlhprob2=0
    elif len(thetanlhmodes)==2 or len(thetanlhmodes)>2:
        thetanlhmode1=thetanlhmodes[0]
        thetanlhprob1=thetanlhprobs[0]
        thetanlhmode2=thetanlhmodes[1]
        thetanlhprob2=thetanlhprobs[1]
    if len(thetanrhmodes)==1:
        thetanrhmode1=thetanrhmodes[0]
        thetanrhprob1=thetanrhprobs[0]
        thetanrhmode2=-99.00
        thetanrhprob2=0
    elif len(thetanrhmodes)==2 or len(thetanrhmodes)>2:
        thetanrhmode1=thetanrhmodes[0]
        thetanrhprob1=thetanrhprobs[0]
        thetanrhmode2=thetanrhmodes[1]
        thetanrhprob2=thetanrhprobs[1]

    #check to see if the number of subsets is specified before making names of files
    if (options["-numsubsets"].value==None):
        if os.path.isdir("%isimswithAllsubsets" % (options["-numdirs"].value)):
            print "The directory %isimswithAllsubsets exists - will rename all old files with the extension .bak\n"
            for path,dirs,filenames in os.walk("%isimswithAllsubsets" % (options["-numdirs"].value)): 
                os.chdir("%isimswithAllsubsets" % (options["-numdirs"].valuei))
                for filename in filenames:
                    if filename[:4]==".nfs":
                        pass
                    else:
                        os.system("mv %s %s.bak" % (filename,filename))
        else:
            os.system("mkdir %isimswithAllsubsets" % (options["-numdirs"].value))
            os.chdir("%isimswithAllsubsets" % (options["-numdirs"].valuei))
        name1="boxplotdatafor"+str(options["-numdirs"].value)+"_Allsubsets_LHprob.dat"
        name2="boxplotdatafor"+str(options["-numdirs"].value)+"_Allsubsets_RHprob.dat"
        name3="boxplotdatafor"+str(options["-numdirs"].value)+"_Allsubsets_LHmean.dat"
        name4="boxplotdatafor"+str(options["-numdirs"].value)+"_Allsubsets_RHmean.dat"
        name5="boxplotdatafor"+str(options["-numdirs"].value)+"_Allsubsets_LHmode1.dat"
        name6="boxplotdatafor"+str(options["-numdirs"].value)+"_Allsubsets_RHmode1.dat"
        name7="boxplotdatafor"+str(options["-numdirs"].value)+"_Allsubsets_LHmode2.dat"
        name8="boxplotdatafor"+str(options["-numdirs"].value)+"_Allsubsets_RHmode2.dat"
    else:
        if os.path.isdir("%isimswith%isubsets" % (options["-numdirs"].value,options["-numsubsets"].value)):
            print "The directory %isimswith%isubsets exists - will rename all old files with the extension .bak\n" % (options["-numdirs"].value,options["-numsubsets"].value)
            for path,dirs,filenames in os.walk("%isimswith%isubsets" % (options["-numdirs"].value,options["-numsubsets"].value)):
                os.chdir("%isimswith%isubsets" % (options["-numdirs"].value,options["-numsubsets"].value))
                for filename in filenames: 
                    if filename[:4]==".nfs":
                        pass
                    else:
                        os.system("mv %s %s.bak" % (filename,filename))
        else:
            os.system("mkdir %isimswith%isubsets" % (options["-numdirs"].value,options["-numsubsets"].value))
            os.chdir("%isimswith%isubsets" % (options["-numdirs"].value,options["-numsubsets"].value))
        name1="boxplotdatafor"+str(options["-numdirs"].value)+"_subsets"+str(options["-numsubsets"].value)+"_LHprob.dat"
        name2="boxplotdatafor"+str(options["-numdirs"].value)+"_subsets"+str(options["-numsubsets"].value)+"_RHprob.dat"
        name3="boxplotdatafor"+str(options["-numdirs"].value)+"_subsets"+str(options["-numsubsets"].value)+"_LHmean.dat"
        name4="boxplotdatafor"+str(options["-numdirs"].value)+"_subsets"+str(options["-numsubsets"].value)+"_RHmean.dat"
        name5="boxplotdatafor"+str(options["-numdirs"].value)+"_subsets"+str(options["-numsubsets"].value)+"_LHmode1.dat"
        name6="boxplotdatafor"+str(options["-numdirs"].value)+"_subsets"+str(options["-numsubsets"].value)+"_RHmode1.dat"
        name7="boxplotdatafor"+str(options["-numdirs"].value)+"_subsets"+str(options["-numsubsets"].value)+"_LHmode2.dat"
        name8="boxplotdatafor"+str(options["-numdirs"].value)+"_subsets"+str(options["-numsubsets"].value)+"_RHmode2.dat"

    #open writing files
    quartilefileLHprob=open(name1,"w")
    quartilefileRHprob=open(name2,"w")
    quartilefileLHmean=open(name3,"w")
    quartilefileRHmean=open(name4,"w")
    quartilefileLHmode1=open(name5,"w")
    quartilefileRHmode1=open(name6,"w")
    quartilefileLHmode2=open(name7,"w")
    quartilefileRHmode2=open(name8,"w")

    #make the thetans directory
    if os.path.isdir("thetans"):
        print "The directory thetans exists - will rename all old files with the extension .bak\n"
        for path,dirs,filenames in os.walk("thetans"):
            os.chdir("thetans")
            for filename in filenames:
                if filename[:4]==".nfs":
                    pass
                else:
                    os.system("mv %s %s.bak" % (filename,filename))	
    else:
        os.system("mkdir thetans")
        os.chdir("thetans")

	#first, write out the data file for the overall contact matrix and get overall contacts
    closestcontactsoverall,closestcontactslh,closestcontactsrh=contact_matrix_closest_contacts_thetans(thetanCMall,thetanCMlh,thetanCMrh,TM1CAs,TM2CAs)

	#get the contact map or heat map for the modes (+/-5 degrees)
    framesLHmode1,framesLHmode2,framesRHmode1,framesRHmode2=getmodestructures(thetanlhmode1,thetanlhmode2,thetanrhmode1,thetanrhmode2,allxangle,allcontactmatricesindividualframes)

	#second, write out all the thetan's in a big data file - LHprob,RHprob,LHmean,RHmean
    thetanoutfile=open(options["-thetans"].value,"w")
    if (options["-numsubsets"].value==None):
        thetanoutfile.write("Theta(n)'s of the following run of "+str(options["-numdirs"].value)+" with All subsets:\n\n")
    else:
        thetanoutfile.write("Theta(n)'s of the following run of "+str(options["-numdirs"].value)+" with "+str(options["-numsubsets"].value)+" subsets:\n\n")
    thetanoutfile.write("Probability:\n\n")
    if thetanlh!=-99.00:
        thetanoutfile.write("\tLeft-handed:  %0.1f \n" % (thetanlh))
    else:
        thetanoutfile.write("\tLeft-handed:  No left-handed crossing angles detected\n")
    if thetanrh!=-99.00:
        thetanoutfile.write("\tRight-handed: %0.1f \n\n" % (thetanrh))
    else:
        thetanoutfile.write("\tRight-handed: No right-handed crossing angles detected\n\n")
	thetanoutfile.write("Mean:\n\n")
    if thetanlhmean!=-99.00:
        thetanoutfile.write("\tLeft-handed:  %i \n" % (thetanlhmean))
    else:
        thetanoutfile.write("\tLeft-handed:  No left-handed crossing angles detected\n")
    if thetanrhmean!=-99.00:
        thetanoutfile.write("\tRight-handed: %i \n\n" % (thetanrhmean))
    else:
        thetanoutfile.write("\tRight-handed: No right-handed crossing angles detected\n\n")
	thetanoutfile.write("Mode(s):\n\n")
    if thetanlhmode1!=-99.00 and thetanlhmode2==-99.00:
        thetanoutfile.write("\tLeft-handed:  %i\t%0.3f\n" % (thetanlhmode1,thetanlhprob1))
    elif thetanlhmode1!=-99.00 and thetanlhmode2!=-99.00:
        thetanoutfile.write("\tFirst Left-handed Mode :  %i\t%0.3f \n" % (thetanlhmode1,thetanlhprob1))
        thetanoutfile.write("\tSecond Left-handed Mode:  %i\t%0.3f \n\n" % (thetanlhmode2,thetanlhprob2))
    else:
        thetanoutfile.write("\tLeft-handed:  No left-handed crossing angles detected\n")
    if thetanrhmode1!=-99.00 and thetanrhmode2==-99.00:
        thetanoutfile.write("\tRight-handed: %i\t%0.3f \n\n" % (thetanrhmode1,thetanrhprob1))
    elif thetanrhmode1!=-99.00 and thetanrhmode2!=-99.00:
        thetanoutfile.write("\tFirst Right-handed Mode :  %i\t%0.3f \n" % (thetanrhmode1,thetanrhprob1))
        thetanoutfile.write("\tSecond Right-handed Mode:  %i\t%0.3f \n" % (thetanrhmode2,thetanrhprob2))
    else:
        thetanoutfile.write("\tRight-handed: No right-handed crossing angles detected\n\n")
    thetanoutfile.write("\n\nThe "+str(options["-topres"].value)+" closest contacts for all frames are: \n\n")
    for i in range(options["-topres"].value):
        thetanoutfile.write("\t%i: %s-%s\n" % (i+1,closestcontactsoverall[i][0],closestcontactsoverall[i][1]))
    if thetanlh!=-99.00:
        thetanoutfile.write("\nThe "+str(options["-topres"].value)+" closest contacts for left-handed frames are: \n\n")
        for i in range(options["-topres"].value):
            thetanoutfile.write("\t%i: %s-%s\n" % (i+1,closestcontactslh[i][0],closestcontactslh[i][1]))
    else:
        thetanoutfile.write("\nThere are no left-handed frames in the simulations. \n\n")
    if thetanrh!=-99.00:
        thetanoutfile.write("\nThe "+str(options["-topres"].value)+" closest contacts for right-handed frames are: \n\n")
        for i in range(options["-topres"].value):
            thetanoutfile.write("\t%i: %s-%s\n" % (i+1,closestcontactsrh[i][0],closestcontactsrh[i][1]))
    else:
        thetanoutfile.write("\nThere are no right-handed frames in the simulations. \n\n")
    thetanoutfile.close()

	#set the pallete for the contact matrix
    palette="0 \\\"yellow\\\", 0.3 \\\"purple\\\", 0.6 \\\"red\\\", 1.0 \\\"black\\\""
    cbrange="0:1"

    #make the gnuplot script, then run gnuplot on all of the contact matrix files to make all of the plots
    os.system("cp /sansom/sb8/bioc1030/scripts/make_contact_matrix_thetans.gnu tmp.gnu")
    os.system("sed -e \"s|XTICS|"+str(TM1reslist)+"|g;s|YTICS|"+str(TM2reslist)+"|g;s|PALETTE|"+str(palette)+"|g;s|CBRANGE|"+str(cbrange)+"|g\" tmp.gnu > make_contact_matrix_thetans.gnu")
    os.system("rm tmp.gnu")
    os.system("/usr/bin/gnuplot -e \"f1='"+str(options["-CMall"].value)+"';f2='"+str(options["-CMlh"].value)+"';f3='"+str(options["-CMrh"].value)+"';numdirs='"+str(options["-numdirs"].value)+"';subsets='"+str(options["-numsubsets"].value)+"';maxX='"+str(options["-TM1l"].value+1)+"';maxY='"+str(options["-TM2l"].value)+"'\" make_contact_matrix_thetans.gnu")

	#make the mode plots here
    os.system("cp /sansom/sb8/bioc1030/scripts/make_contact_matrix_modes.gnu tmp.gnu")
    os.system("sed -e \"s|XTICS|"+str(TM1reslist)+"|g;s|YTICS|"+str(TM2reslist)+"|g;s|PALETTE|"+str(palette)+"|g;s|CBRANGE|"+str(cbrange)+"|g;s|LHMODE1|"+str(framesLHmode1)+"|g;s|LHMODE2|"+str(framesLHmode2)+"|g;s|RHMODE1|"+str(framesRHmode1)+"|g;s|RHMODE2|"+str(framesRHmode2)+"|g;\" tmp.gnu > make_contact_matrix_modes.gnu")
    os.system("rm tmp.gnu")
    os.system("/usr/bin/gnuplot -e \"f1='thetanlhmode1matrix.dat';f2='thetanlhmode2matrix.dat';f3='thetanrhmode1matrix.dat';f4='thetanrhmode2matrix.dat';lhmode1='"+str(thetanlhmode1)+"';lhmode2='"+str(thetanlhmode2)+"';rhmode1='"+str(thetanrhmode1)+"';rhmode2='"+str(thetanrhmode2)+"';maxX='"+str(options["-TM1l"].value+1)+"';maxY='"+str(options["-TM2l"].value+1)+"'\" make_contact_matrix_modes.gnu")

    #change into the above directory
    os.chdir("../")

    #print the warning that this will start the permutations
    print
    print "Starting permutations..."
    print 

    #loop over all of the groups specified, and run the calculations
    for group in reversed(range(options["-dend"].value,options["-dstart"].value+1,options["-inc"].value)):
        calculations_multiprocessing(group,TM1reslist,TM2reslist,allxangle,thetanCMall,thetanCMlh,thetanCMrh,thetanlh,thetanrh,thetanlhmean,thetanrhmean,thetanlhmode1,thetanrhmode1,thetanlhmode2,thetanrhmode2,allcontactmatricesnumdir,lhcontactmatricesnumdir,rhcontactmatricesnumdir,arraynumframesconsideredall,arraynumframesconsideredlh,arraynumframesconsideredrh,xangles,TM1CAs,TM2CAs,quartilefileLHprob,quartilefileRHprob,quartilefileLHmean,quartilefileRHmean,quartilefileLHmode1,quartilefileRHmode1,quartilefileLHmode2,quartilefileRHmode2)

	#close the writing files
    quartilefileLHprob.close()
    quartilefileRHprob.close()
    quartilefileLHmean.close()
    quartilefileRHmean.close()
    quartilefileLHmode1.close()
    quartilefileRHmode1.close()
    quartilefileLHmode2.close()
    quartilefileRHmode2.close()

	#make a new directory, move files, and change into that directory
    if os.path.isdir("overall_data_and_plots"):
        print "The directory overall_data_and_plots exists - will rename all old files with the extension .bak\n"
        for path,dirs,filenames in os.walk("overall_data_and_plots"):
            os.chdir("overall_data_and_plots")
            for filename in filenames:
                if filename[:4]==".nfs":
                    pass
                else:
                    os.system("mv %s %s.bak" % (filename,filename))
    else:
        os.system("mkdir overall_data_and_plots")
        os.chdir("overall_data_and_plots")

    #move the dat files into the correct folder
    os.system("mv ../*.dat .")

    #copy all of the scripts for the boxplot into the correct directory
    os.system("cp /sansom/sb8/bioc1030/scripts/make_boxplot_*.gnu .")

    #set tics
    xtics=""
    ticsLHmeanandmode=""
    ticsRHmeanandmode=""
    for i in frange(0,options["-numdirs"].value+1,float(options["-numdirs"].value)/float(10)):
        xtics+=str(i)+","
    for i in frange(options["-bounds"].value[0],5,options["-binsize"].value*5):
        ticsRHmeanandmode+=str(i)+","
    for i in frange(0,options["-bounds"].value[1]+5,options["-binsize"].value*5):
        ticsLHmeanandmode+=str(i)+","
    xtics=xtics[:-1]
    ticsLHmeanandmode=ticsLHmeanandmode[:-1]
    ticsRHmeanandmode=ticsRHmeanandmode[:-1]
    yrangeLH="0:"+str(options["-bounds"].value[1])
    yrangeRH=str(options["-bounds"].value[0])+":0"

    #edit gnuplot scripts
    os.system("sed -e \"s|XTICS|"+str(xtics)+"|g\" make_boxplot_prob.gnu > tmp.gnu")
    os.system("mv tmp.gnu make_boxplot_prob.gnu")
    os.system("sed -e \"s|XTICS|"+str(xtics)+"|g\" make_boxplot_mean.gnu > tmp.gnu")
    os.system("sed -e \"s|YTICSLH|"+str(ticsLHmeanandmode)+"|g\" tmp.gnu > tmp2.gnu")
    os.system("sed -e \"s|YTICSRH|"+str(ticsRHmeanandmode)+"|g\" tmp2.gnu > tmp3.gnu")
    os.system("sed -e \"s|YRANGELH|"+str(yrangeLH)+"|g\" tmp3.gnu > tmp4.gnu")
    os.system("sed -e \"s|YRANGERH|"+str(yrangeRH)+"|g\" tmp4.gnu > make_boxplot_mean.gnu")
    os.system("rm tmp*.gnu")
    os.system("sed -e \"s|XTICS|"+str(xtics)+"|g\" make_boxplot_mode.gnu > tmp.gnu")
    os.system("sed -e \"s|YTICSLH|"+str(ticsLHmeanandmode)+"|g\" tmp.gnu > tmp2.gnu")
    os.system("sed -e \"s|YTICSRH|"+str(ticsRHmeanandmode)+"|g\" tmp2.gnu > tmp3.gnu")
    os.system("sed -e \"s|YRANGELH|"+str(yrangeLH)+"|g\" tmp3.gnu > tmp4.gnu")
    os.system("sed -e \"s|YRANGERH|"+str(yrangeRH)+"|g\" tmp4.gnu > make_boxplot_mode.gnu")
    os.system("rm tmp*.gnu")
    os.system("sed -e \"s|XTICS|"+str(xtics)+"|g\" make_boxplot_twomodes.gnu > tmp.gnu")
    os.system("sed -e \"s|YTICSLH|"+str(ticsLHmeanandmode)+"|g\" tmp.gnu > tmp2.gnu")
    os.system("sed -e \"s|YTICSRH|"+str(ticsRHmeanandmode)+"|g\" tmp2.gnu > tmp3.gnu")
    os.system("sed -e \"s|YRANGELH|"+str(yrangeLH)+"|g\" tmp3.gnu > tmp4.gnu")
    os.system("sed -e \"s|YRANGERH|"+str(yrangeRH)+"|g\" tmp4.gnu > make_boxplot_twomodes.gnu")
    os.system("rm tmp*.gnu")

    #make the box plots
    if (options["-numsubsets"].value==None):
        os.system("/usr/bin/gnuplot -e \"f1='"+str(name1)+"';f2='"+str(name2)+"';numdirs='"+str(options["-numdirs"].value)+"';subsets='All'\" make_boxplot_prob.gnu")
        os.system("/usr/bin/gnuplot -e \"f1='"+str(name3)+"';f2='"+str(name4)+"';numdirs='"+str(options["-numdirs"].value)+"';subsets='All'\" make_boxplot_mean.gnu")
        if thetanlhmode1!=-99.00 and thetanlhmode2!=-99.00 and thetanrhmode1!=-99.00 and theatnrhmode2!=-99.00:
            os.system("/usr/bin/gnuplot -e \"f1='"+str(name6)+"';f2='"+str(name6)+"';f3='"+str(name7)+"';f4='"+str(name8)+"';numdirs='"+str(options["-numdirs"].value)+"';subsets='All'\" make_boxplot_twomodes.gnu")
        else:
            os.system("/usr/bin/gnuplot -e \"f1='"+str(name6)+"';f2='"+str(name6)+"';numdirs='"+str(options["-numdirs"].value)+"';subsets='All'\" make_boxplot_mode.gnu")
    else:
        os.system("/usr/bin/gnuplot -e \"f1='"+str(name1)+"';f2='"+str(name2)+"';numdirs='"+str(options["-numdirs"].value)+"';subsets='"+str(options["-numsubsets"].value)+"'\" make_boxplot_prob.gnu")
        os.system("/usr/bin/gnuplot -e \"f1='"+str(name3)+"';f2='"+str(name4)+"';numdirs='"+str(options["-numdirs"].value)+"';subsets='"+str(options["-numsubsets"].value)+"'\" make_boxplot_mean.gnu")
        if thetanlhmode1!=-99.00 and thetanlhmode2!=-99.00 and thetanrhmode1!=-99.00 and thetanrhmode2!=-99.00:
            os.system("/usr/bin/gnuplot -e \"f1='"+str(name5)+"';f2='"+str(name6)+"';f3='"+str(name7)+"';f4='"+str(name8)+"';numdirs='"+str(options["-numdirs"].value)+"';subsets='"+str(options["-numsubsets"].value)+"'\" make_boxplot_twomodes.gnu")
        else:
            os.system("/usr/bin/gnuplot -e \"f1='"+str(name5)+"';f2='"+str(name6)+"';numdirs='"+str(options["-numdirs"].value)+"';subsets='"+str(options["-numsubsets"].value)+"'\" make_boxplot_mode.gnu")

    #change into the above directory
    os.chdir("../")

'''
execute the program
'''
if __name__=="__main__":
	main()


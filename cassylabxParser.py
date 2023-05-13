'''cassylabxParser.py

   read files in xml-format produced with the Leybold Cassy system
   (very simple example for sensor box 524032)

   contains:

     * print header information and extract lists of data

     * peak search to extract envelope 

.. moduleauthor:: Guenter Quast <g.quast@kit.edu>

'''

def labxParser(file, prlevel=3):
  '''   
  read files in xml-format produced with the Leybold Cassy sytem
   
  Args:
     * file:  input data in .labx format

  Returns:
     * 2d list:         list of measurement vectors read from file 
     * list of strings: names of measurement
  '''
# --------------------------------------------------------------------
# dependencies: xml.etree.ElementTree
#
#  30-Oct-16  initial version
# changes :
# --------------------------------------------------------------------
  import xml.etree.ElementTree as ET
  import numpy as np, matplotlib.pyplot as plt
  import sys

  root = ET.parse(file).getroot()
  if root.tag != 'cassylab':
    print " !!! only cassylab supported - exiting (1) !"
    sys.exit(1)    
  else:
    if(prlevel): print "\n\n*==* labxParser: name of XML root object:\n",\
    ' ', root.tag, root.attrib
#
# some sanity checks wether we got the right kind of input
  if not root.findall('cassys'):
    print " !!! no tag 'casssys' found  - exiting (2) !"
    sys.exit(2)    
  if not root.findall('ios'):
    print "! !!! tag 'ios' not found exiting (3) !"
    sys.exit(3)    
#
# print header info if requested 
  if (prlevel):
    childvals=[]
    childtags=[]
    for child in root:
      childtags.append(child.tag)
      childvals.append(child.attrib)
    print "    %i children found, " %(len(childtags)), 
    print "tags and attributes: \n",
    for i in range(len(childtags)):
      print '   ', childtags[i],' : ', childvals[i]

  if(prlevel>1):
    print '\n *==*  Details:' 
    print " ** found tag 'ios', configuration:"
    for ios in root.findall('ios'):
      print '   ', ios.tag, ios.attrib
    print "   measurement settings:"
    i=0
    for io in ios.findall('io'): 
      i+=1      
      print "  --> io %i:"%i, io.attrib
      for qs in io.iter('quantities'): print '   ', qs.tag, qs.attrib
      for q in qs.iter('quantity'): print '   ', q.tag, q.attrib

  if(prlevel>2):
    if root.findall('calcs'):
      print "\n ** found tag 'calcs', calc settings:"
      for calcs in root.findall('calcs'):
        i=0
        for calc in calcs.findall('calc'): 
          i+=1      
          print "  --> calc %i:"%i, calc.attrib

  print '    - - - - - - - - - - - - - - - - - - - - - - - - -\n ' 

  # cassylab stores data under the tag "channels:channel:values", 
  #    search for and extract data from xml structure
  varray=[]
  vnames=[]
  iv=0
  ic=0
  for clist in root.iter('channels'):
    for c in clist:
      ic+=1
      quantity=c.find('quantity').text
      vnames.append(quantity)
      if(prlevel>1): 
        print "   --> new channel found %i, quantity %s"\
          %(ic, quantity)
        if(prlevel>2): print "     ", c.attrib
      values=c.find('values')
      if(prlevel>2): print "     number of values: ", values.attrib
      varray.append([])
      for v in values:
        varray[iv].append(np.float32(v.text))
      iv+=1

  if (prlevel): 
    print "*==* labxParser:  %i value lists found"%iv
    print "Quantities: ", vnames
    print "\n\n"

  return varray, vnames


def Fourier_fft(t, a):
  ''' 
  Fourier transform of the amplitude spectrum a(t) using numpy.fft

    Args:
      * t: np-array of time values
      * a: np-array amplidude a(t)
 
    Returns:
      * arrays f, a_f: frequencies and amplites
  '''
# -----------------------------------------------
  from numpy.fft import fft, fftfreq

  n = len(t)
  dt = (t[-1]-t[0])/(n-1.) # time step
  freq = fftfreq(n, dt)
  amp = abs(fft(a))

  return freq, amp


def FourierSpectrum(t, a):
  ''' 
  Fourier transform of the amplitude spectrum a(t)
   (an example of a calculation "by hand")

    Args:
      * t: np-array of time values
      * a: np-array amplidude a(t)
 
    Returns:
      * arrays freq, amp: frequencies and amplites
  '''
# -----------------------------------------------
  import numpy as np

  n = len(t)
  dt = (t[-1]-t[0])/(n-1.) # time step
  fmax = 0.5/dt
  df = fmax/n
  freq = np.linspace(0., fmax, n/2)
  amp = np.zeros(len(freq))

  i=0
  for f in freq:
    omega = 2. * np.pi * f
    s=sum(a * np.sin(omega * t))/n
    c=sum(a * np.cos(omega * t))/n
    amp[i] = np.sqrt(s**2 + c**2)
    i+=1

  return freq, amp


def peakfinder(x,y, th=0.):
  ''' 
  find positions of all Peaks and Dips in data
    x-coordinates are determined from weighted average over 3 data points

    Args:
      * x: np-array of values
      * y: np-array of values at postion x
      * th: float, threshold for peaks
 
    Returns:
      * np-arrays xpeak, ypeak, xdip, ydip
  '''
# -----------------------------------------------
  xpeak=[]
  ypeak=[]
  xdip=[]
  ydip=[]
  if y[0]-y[1]>th and y[0]-y[2]>th:
    xpeak.append(x[0])    
    ypeak.append(y[0])    
  if y[0]-y[1]<th and y[0]-y[2]<th:
    xdip.append(x[0])    
    ydip.append(y[0])    
  for i in range(1,len(x)-1):
    if y[i]-y[i-1]>=0. and y[i]-y[i+1]>th: 
      xpeak.append(sum(x[i-1:i+1]*y[i-1:i+1])/sum(y[i-1:i+1]))    
      ypeak.append(y[i])    
    if y[i]-y[i-1]<=0. and y[i]-y[i+1]<-th: 
      xdip.append(sum(x[i-1:i+1]*y[i-1:i+1])/sum(y[i-1:i+1]))    
      ydip.append(y[i])    
  if y[-1]-y[-2]>th and y[-1]-y[-3]>th:
    xpeak.append(x[-1])    
    ypeak.append(y[-1])    
  if y[-1]-y[-2]<th and y[-1]-y[-3]<th:
    xdip.append(x[-1])    
    ydip.append(y[-1])    

  return xpeak, ypeak, xdip, ydip


# -----example Code illustrating usage --------------------
if __name__ == "__main__":
  import numpy as np, matplotlib.pyplot as plt
  from PhyPraKit import kFit

  filename="CassyExample.labx"
  values, names = labxParser(filename, prlevel=1)

  # interpret returned vectors:
  i_n=0
  i_t=1
  i_dummy=2
  i_phic=3
  i_phi=4
  i_omega=5
  nl=len(names)

  if nl>=4:
  # convert to numpy arrays
    n=np.array(values[i_n][:])
    t=np.array(values[i_t][:])
    phic=np.array(values[i_phic][:])
    phi=np.array(values[i_phi][:])
    print "time=", len(phic), t
    print "phi_count=", len(phic), phic
    print "phi=", len(phi), phi

# ----  some examples of analysis -----

# calculate fourier spectrum 
  freq, amp = FourierSpectrum(t, phi)

# run peak finder
  tp, phip, td, phid = peakfinder(t, phi, th=0.01)

# define functional form of envelope curves ...
  def envelope_p(x, a=1., b=1., c=1.):
    return a*x**2+b*x+c

  def envelope_d(x, a=-1., b=1., c=1.):
    return a*x**2+b*x+c

# ... and fit parameters with kafe (kFit from PhyPraPit)  
  parp, parpe, corp, chi2p =\
    kFit(envelope_p, tp[:-10], phip[:-10], 0., 0.03, 
     plot=False, quiet=True)
  print "fit of positive envelope, chi2/df= %.2g"%(chi2p/(len(tp[:-10])-3.))
  np.set_printoptions(precision=3)
  print " -> parameters:   ", parp
  np.set_printoptions(precision=2)
  print " -> uncertainties:", parpe 
  print " -> correlations:\n", corp

  pard, parde, cord, chi2d =\
    kFit(envelope_d, td[:-10], phid[:-10], 0., 0.03, 
     plot=False, quiet=True)
  print "fit of negative envelope, chi2/df= %.2g"%(chi2d/(len(td[:-10])-3.))
  np.set_printoptions(precision=3)
  print " -> parameters:   ", pard
  np.set_printoptions(precision=2)
  print " -> uncertainties:", parde 
  print " -> correlations:\n", cord

# plot data and analysis results
  fig=plt.figure(figsize=(15., 10.))

  ax1=fig.add_subplot(2,1,1)
  ax1.plot(t,phi)
  ax1.plot(tp, phip, 'rx', label='peak')
  ax1.plot(td, phid, 'gx', label='dip')
  x=np.linspace(0., max(tp),100)
  ax1.plot(x, envelope_p(x, parp[0], parp[1], parp[2]), 
    'r-', label='positive envelope fit')
  ax1.plot(x, envelope_d(x, pard[0], pard[1], pard[2]), 
    'g-', label='negative envelope fit')
  ax1.set_xlabel('$time$ $t$ (s)', size='x-large')
  ax1.set_ylabel('$amplitude$  $\phi$', size='x-large')
  ax1.legend(loc='best',numpoints=1)
  ax1.grid()

  ax2=fig.add_subplot(2,1,2)
  ax2.plot(freq, amp, 'b.')
  ax2.set_xlabel('$frequency$ $f$ (Hz)', size='x-large')
  ax2.set_ylabel('$amplitude$', size='x-large')
  ax2.set_yscale('log')
  ax2.grid()

  plt.show()

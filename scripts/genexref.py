
# Abstract the process of crossreferencing gene names, e.g. unigene x
# refgene or ucsc x refgene.  Usually this is just a dictionary lookup
# based on data read from a file, but it could also be based on data
# from a database.
#

import csv


class GeneXref:
  """
  callback for matplotlib to display an annotation when points are clicked on.  The
  point which is closest to the click and within xtol and ytol is identified.
    
  Register this function like this:
    
  scatter(xdata, ydata)
  af = AnnoteFinder(xdata, ydata, annotes)
  connect('button_press_event', af)
  """

  def __init__(self, xrefFile):
    self._refmap = dict()
    f = open(xrefFile, 'rt')
    try:
        refmap_reader = csv.DictReader(f, delimiter='\t')
        for row in refmap_reader:
            self._refmap[row['locusLinkId']] = row['mrnaAcc']
        
    finally:
        f.close()

  # Read the refseq_reflink.txt
  # name		product			mrnaAcc		protAcc		geneName	prodName	locusLinkId	omimId
  # MIR375		NR_035777		390601		0		100316383	0
  # Os08g0544800	hypothetical protein	NM_001068944	NP_001062409	235443		97		4346214	0
  # Os04g0476700	hypothetical protein	NM_001059617	NP_001053082	265619		97		4336157	0
  # Os01g0170300	hypothetical protein	NM_001048674	NP_001042139	233119		97		4327633	0
  #

  def unigene2refgene(self, unigene):
      return self._refmap[unigene]
    

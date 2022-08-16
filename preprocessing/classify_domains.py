#!/usr/bin/env python

import re
import sys
import xml.etree.ElementTree as ET


inputFileName = '../../../data/interpro77.0/interpro.xml'


""" Remove citations
"""
print('Remove citations in square brackets... ', file = sys.stderr)
with open(inputFileName , 'r') as inputFile, \
     open('../output/classify_domains/interpro_nocit.xml', 'w') as outputFile:
    inputXML = inputFile.read()
    outputXML = re.sub(r"\[[^\[\]]*\]", \
                       "", \
                       inputXML)
    outputFile.write(outputXML)
print('Done', file = sys.stderr)


""" Remove the beginning part and 
    correct formatting errors
"""
print('Correct formatting errors... ', file = sys.stderr)
with open('../output/classify_domains/interpro_nocit.xml', 'r') as inputFile, \
     open('../output/classify_domains/interpro_nocit_corrected.xml', 'w') as outputFile:
    for lineNumber, line in enumerate(inputFile, start = 1):
        if (lineNumber == 2) or \
           ((lineNumber >= 4) and (lineNumber <= 26)):
            continue
        elif lineNumber == 41700:
            replacementLine = '<p>The cation dependent mannose-6-phosphate (man-6-P) \
receptor is one of two transmembrane proteins involved \
in the transport of lysosomal enzymes from the Golgi \
complex and the cell surface to lysosomes . Lysosomal \
enzymes bearing phosphomannosyl residues bind specifically \
to man-6-P receptors in the Golgi apparatus and the resulting \
receptor-ligand complex is transported to an acidic prelyosomal \
compartment, where the low pH mediates dissociation of the \
complex. Binding is optimal in the presence of divalent cations.</p> \
<p>The amino acid sequence is a single polypeptide chain that \
contains a putative signal sequence and a transmembrane domain . \
The cation-dependent mannose 6-phosphate (M6P)'
            outputFile.write(replacementLine + '\n')
        elif lineNumber == 41718:
            continue
        elif lineNumber == 1000964:
            replacementLine = '<p>This group of metallopeptidases \
belong to the MEROPS peptidase family M11 \
(gametolysin family, clan MA(M)). The protein fold \
of the peptidase domain for members of this family \
resembles that of thermolysin, the type example for \
clan MA and the predicted active site residues for \
members of this family and thermolysin occur in \
the motif HEXXH . The type example is \
gametolysin from the  unicellular \
biflagellated alga, Chlamydomonas \
reinhardtii Gametolysin is a zinc-containing metallo-protease, \
which is responsible for the degradation of the cell wall. \
Homologues of gametolysin have also been reported in \
the simple multicellular organism, Volvox .</p>'
            outputFile.write(replacementLine + '\n')
        else:
            outputFile.write(line)
print('Done', file = sys.stderr)


""" Parse InterPro entries 
"""
print('Parse InterPro entries... ', file = sys.stderr)
tree = ET.parse('../output/classify_domains/interpro_nocit_corrected.xml')
root = tree.getroot()
print('Done', file = sys.stderr)


""" Classify domains into ' ', which stands for 'putative DBD',
    and 'Other'
"""
print('Classify domains... ', file = sys.stderr)
iprClassification = dict()
for ipr in root:
    if ipr.tag == 'interpro':
        iprAccession = ipr.get('id')
        iprName = ipr.find('name').text
        if ipr.find('abstract'):
            iprAbstract = '\n'.join(p.text for p in ipr.find('abstract').iter(tag = 'p'))
        else:
            iprAbstract = ''
        if ipr.find('class_list'):
            iprGOTerms = ipr.find('class_list')
            goDescriptions = '; '. join(goTermDescription.text for goTermDescription in iprGOTerms.iter('description'))
        else:
            goDescriptions = ''
        if (not re.search('DNA', iprName)) and \
           (not re.search('DNA', iprAbstract)) and \
           (not re.search('DNA', goDescriptions)):
            iprClassification[iprAccession] = 'Other'
        else:
            iprClassification[iprAccession] = ' '
print('Done', file = sys.stderr)


""" Output the classification table
"""
print('Output classification table... ', file = sys.stderr)
outputTSV = '\n'.join(iprAccession + '\t' + iprClassification[iprAccession] for iprAccession in iprClassification) + '\n'
outputTSV = 'ipr_accession' + '\t' + 'category' + '\n' + outputTSV 
with open('../output/classify_domains/interpro_domain_binary_classification.tsv', 'w') as outputFile:
    outputFile.write(outputTSV)
print('Done', file = sys.stderr)

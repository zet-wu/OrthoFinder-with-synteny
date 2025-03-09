from scipy import sparse
import numpy as np

class IdFromAccession():
    def __init__(self, idsFilename):
        self.nameToIdDicts_in_species = []
        with open(idsFilename, 'r') as idsFile:
            for line in idsFile:
                id, rest = line.split(": ", 1)
                accession = rest.split(None, 1)[0]
                iSp = int(id.split("_")[0])
                if (cur_len := len(self.nameToIdDicts_in_species)) < iSp + 1:
                    self.nameToIdDicts_in_species.extend([{}] * (iSp + 1 - cur_len))
                
                if accession in self.nameToIdDicts_in_species[iSp]:
                    raise RuntimeError("ERROR: A duplicate id was found in the fasta files: % s" % id)
                self.nameToIdDicts_in_species[iSp][accession] = id   

def syntenyMatrix(seqsInfo, syntenyFilename:str, idFromName:IdFromAccession):
    if(syntenyFilename.endswith('.anchors')):
        format = 'anchors'
    elif(syntenyFilename.endswith('.collinearity')):
        format = 'collinearity'
    with open(syntenyFilename, 'r') as syntenyFile: 
        iSp = -1
        jSp = -1
        syntenyiijj = []
        for line in syntenyFile:
            if(line.startswith('#')):
                continue
            values = line.split('\t')
            if(format == 'anchors'):
                if(len(values) != 3):
                    continue
                iName, jName, score = values
            elif(format == 'collinearity'):
                if(len(values) != 4):
                    continue
                values = list(map(str.strip, values))
                _, iName, jName, score = values
            if(iSp < 0):
                for i, dict in enumerate(idFromName.nameToIdDicts_in_species):
                    if(iName in dict):
                        iSp = i
                        break
                if(iSp < 0):
                    continue
            if(jSp < 0):
                for j, dict in enumerate(idFromName.nameToIdDicts_in_species):
                    if(jName in dict):
                        jSp = j
                        break
                if(jSp < 0):
                    continue
            try:
                iNum = int(idFromName.nameToIdDicts_in_species[iSp][iName].split('_')[1])
                jNum = int(idFromName.nameToIdDicts_in_species[jSp][jName].split('_')[1])
            except KeyError:
                continue
            syntenyiijj.append((iNum, jNum))
        if(syntenyiijj):
            ii, jj = zip(*syntenyiijj)
        else:
            ii = ()
            jj = ()
        onesArray = np.ones(len(syntenyiijj))
        mat = sparse.csr_matrix((onesArray,  (ii, jj)),
                                shape=(seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[iSp]], seqsInfo.nSeqsPerSpecies[seqsInfo.speciesToUse[jSp]]))
        return mat, (iSp, jSp)

    
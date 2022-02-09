import sys
import pandas as pd

input_file = sys.argv[1]
out_file = sys.argv[2]

with open(input_file, 'r') as fopen:
    data = fopen.readlines()[44:]
saved_data = []


showEnergy = ["Backbone", "Rama", "Contact", "Fragment", "Qdiff", "Beta", "Pap", "DNABond", "DNAAngle", "DNAStacking", "DNADihe", "DNABP",\
              "DNACS", "ExclusionD", "ElectroD", "ExclusionPD", "ElectroPD", "ATPBond", "ATPBias", "ExclusionPL", "ElectroPL", "Total"]
line = " ".join(["{0:<8s}".format(i) for i in ["Steps"] + showEnergy])
#print(line)
step_count = 0
for each_line in data:
    data = each_line.split()
    if data[0] == 'TotalEnergy':
        saved_data.append({"Steps":0, "Total":float(data[-2])})
        step_count += 1
    if data[0] == 'Ligandbone':
        saved_data[-1].update({"ATPBond":float(data[-2])})
    if data[0] == 'Bond':
        saved_data[-1].update({"DNABond":float(data[-2])})
    if data[0] == 'Angle':
        saved_data[-1].update({"DNAAngle":float(data[-2])})
    if data[0] == 'Stacking':
        saved_data[-1].update({"DNAStacking":float(data[-2])})
    if data[0] == 'BasePair':
        saved_data[-1].update({"DNADP":float(data[-2])})
    if data[0] == 'Dihedral':
        saved_data[-1].update({"DNADihe":float(data[-2])})
    if data[0] == 'CrossStacking':
        saved_data[-1].update({"DNACS":float(data[-2])})
    if data[0] == 'Exclusion':
        saved_data[-1].update({"ExclusionD":float(data[-2])})
    if data[0] == 'Electrostatics':
        saved_data[-1].update({"ElectroD":float(data[-2])})
    if data[0] == 'Connectivity':
        saved_data[-1].update({"Backbone":float(data[-2])})
    if data[0] == 'rama':
        saved_data[-1].update({"Rama":float(data[-2])})
    if data[0] == 'contact':
        saved_data[-1].update({"Contact":float(data[-2])})
    if data[0] == 'beta1':
        saved_data[-1].update({"Beta":float(data[-2])})
    if data[0] == 'pap1':
        saved_data[-1].update({"Pap":float(data[-2])})
    if data[0] == 'fragment':
        saved_data[-1].update({"Fragment":float(data[-2])})
    if data[0] == 'qdiff':
        saved_data[-1].update({"Qdiff":float(data[-2])})
    if data[0] == 'LigandDistanceRestraint1':
        saved_data[-1].update({"ATPBias":float(data[-2])})
    if data[0] == 'ExclusionProteinLigand':
        saved_data[-1].update({"ExclusionPL":float(data[-2])})
    if data[0] == 'ElectrostaticsProteinLigand':
        saved_data[-1].update({"ElectroPL":float(data[-2])})
    if data[0] == 'ExclusionProteinDNA':
        saved_data[-1].update({"ExclusionPD":float(data[-2])})
    if data[0] == 'ElectrostaticsProteinDNA':
        saved_data[-1].update({"ElectroPD":float(data[-2])})
    if len(data) == 1:
        saved_data.append({"Steps":step_count, "Total":float(data[0])})
        step_count += 1
            
    #print(saved_data)
df = pd.DataFrame.from_dict(saved_data)
#print(df)
df.to_csv(sys.argv[2], index=None, sep='\t', float_format='%8f', mode='w')

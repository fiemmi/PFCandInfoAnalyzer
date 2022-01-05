import argparse
import os
import sys

parser=argparse.ArgumentParser()
parser.add_argument("--nlists", type=int, default=16, help="Number of output lists")
parser.add_argument("--base_dir", type=str, default="/afs/cern.ch/work/f/fiemmi/private/CMSSW_10_6_20/src/PFCandInfo/PFCandInfoAnalyzer/sorting_framework", help="Base directory")
parser.add_argument("--list1", type=str, default="/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/events_EpsilonPU_EXT80k_v9-v1.txt", help="Path to 1st .txt file to open")
parser.add_argument("--list2", type=str, default="/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/events_PU_EXT80k_v9-v1.txt", help="Path to 2nd .txt file to open")

FLAGS=parser.parse_args()
nlists=FLAGS.nlists
base_dir=FLAGS.base_dir
list1=FLAGS.list1
list2=FLAGS.list2

sys.path.append(base_dir)

file1 = open(list1, "r")
file1_l = []
for line in file1 :
    file1_l.append(line)

file2 = open(list2, "r")
file2_l = []
lines_read = 0
for line in file2 :
    file2_l.append(line)


#print(file1_l)
#print(file2_l)

indices = []
indices_read = 0
print("Finding indices...")
for element in file1_l :
    indices_read += 1
    if indices_read % 10000 == 0 :
            print("ievent = {}".format(indices_read))
    indices.append(file2_l.index(element))

#print(indices)

range_ = len(file1_l)/nlists
lists = [indices[i:i+range_] for i in range(0,len(file1_l), range_)]
#print(lists)


count_events = 0
if not os.path.exists("lists_of_indices") :
    os.mkdir("lists_of_indices")

for i in range(nlists) :
    file_ = open("lists_of_indices/list_"+str(i)+".txt", "w")
    #file_ = open("list_"+str(i)+".txt", "w")
    for j in range(len(lists[i])) :
        count_events += 1
        if count_events % 10000 == 0 :
            print("i-event = {}".format(count_events))
        file_.write("%s\n" % lists[i][j])

print("Done. Processed {} events.".format(count_events))

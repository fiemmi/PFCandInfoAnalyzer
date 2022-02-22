import argparse
import os
import sys

parser=argparse.ArgumentParser()
parser.add_argument("--evts_per_list", type=int, default=16, help="Number of events per output list")
parser.add_argument("--base_dir", type=str, default="/afs/cern.ch/work/f/fiemmi/private/CMSSW_10_6_20/src/PFCandInfo/PFCandInfoAnalyzer/sorting_framework", help="Base directory")
parser.add_argument("--list1", type=str, default="/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/events_EpsilonPU_EXT80k_v9-v1.txt", help="Path to 1st .txt file to open")
parser.add_argument("--list2", type=str, default="/eos/user/f/fiemmi/JetMET/ntuplize/CMSSW_10_6_16/src/events_PU_EXT80k_v9-v1.txt", help="Path to 2nd .txt file to open")

FLAGS=parser.parse_args()
evts_per_list=FLAGS.evts_per_list
base_dir=FLAGS.base_dir
list1=FLAGS.list1
list2=FLAGS.list2

sys.path.append(base_dir)

with open(list1, "r") as file1:
    file1_l = []
    for line in file1 :
        file1_l.append(line)

with open(list2, "r") as file2:
    file2_l = []
    lines_read = 0
    for line in file2 :
        file2_l.append(line)

indices = []
indices_read = 0
print("Finding indices...")
for element in file1_l :
    indices_read += 1
    if indices_read % 10000 == 0 :
        print("ievent = {}".format(indices_read))
    try:
        indices.append(file2_l.index(element))
    except: #if event not found in second list, skip it
        continue

range_ = evts_per_list
lists = [indices[i:i+range_] for i in range(0,len(file2_l), range_)]

count_events = 0
if not os.path.exists("lists_of_indices") :
    os.mkdir("lists_of_indices")

for i in range(len(lists)) :
    with open("lists_of_indices/list_"+str(i)+".txt", "w") as file_:
        for j in range(len(lists[i])) :
            count_events += 1
            if count_events % 10000 == 0 :
                print("i-event = {}".format(count_events))
            file_.write("%s\n" % lists[i][j])

print("Done. Processed {} events.".format(count_events))

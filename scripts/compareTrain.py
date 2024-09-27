import sys

filename_1 = sys.argv[1]

with open(filename_1, 'r') as file1:
    content1_ = file1.readlines()
    content1  = content1_[1:]
    header    = content1_[0]
    with open(filename_1, 'r') as file2:
        content2 = file2.readlines()[1:]
    count_duplicates = 0
    new_content = list(header)
    idx = []
    for i in range(len(content1)):
        for j in range(i+1,len(content2)):
            if content1[i]==content2[j]:
                count_duplicates+=1
                idx.append(j)
                
for i in range(len(content1)):
    if i not in idx:
        new_content.append(content1[i])
print(len(idx))
print(count_duplicates)
outfile = filename_1[:-4] + "_new.csv"

with open(outfile,"w") as file_:
     file_.writelines(new_content)


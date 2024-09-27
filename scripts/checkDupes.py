import sys

filename_read = sys.argv[1]
if len(sys.argv) > 2:
    option = sys.argv[2]

file_read = open(filename_read, 'r+')
data_read = file_read.readlines()

counter = 0
for i in range(len(data_read)):
    line_check_1 = data_read[i]
    # print("Checking line " + str(i) + "\n")
    for j in range(len(data_read)):
        # print("Against line " + str(j) + "\n")
        line_check_2 = data_read[j]
        if i == j:
            continue
        if line_check_1 == line_check_2:
            if len(sys.argv) > 2 and sys.argv[2] == "v":
                print("MATCH LINE " + str(i) + " AND " + str(j))
                print(line_check_1)
                print(line_check_2)
            counter += 1

new_count=int(counter/2)
print(str(new_count))



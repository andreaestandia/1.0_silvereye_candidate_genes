import sys

# Load dictionary from map file
dictionary = {}
with open(sys.argv[2]) as map_file:
    for line in map_file:
        fields = line.strip().split()
        dictionary[fields[0]] = fields[1]

# Update population labels in BA3 file
with open(sys.argv[1]) as ba3_file:
    for line in ba3_file:
        fields = line.strip().split()
        if fields[0] in dictionary:
            fields[1] = dictionary[fields[0]]
        print('\t'.join(fields))

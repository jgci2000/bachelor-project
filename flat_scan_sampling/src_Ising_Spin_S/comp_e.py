
NN_table = list()

with open("./neighbour_tables/neighbour_table_2D_SS_4NN_L4.txt") as f:
    for line in f:
        NN_table.append(line.strip().split(" "))
    
spins = [-2, -2, 2, -2, -2, -2, -2, -2, -2, -2, -2, -2, 0, -2, -2, -2]
E = 0

for i in range(len(spins)):
    E += - spins[i] * (spins[int(NN_table[i][0])]+ spins[int(NN_table[i][1])] + spins[int(NN_table[i][2])] + spins[int(NN_table[i][3])])

print(E/2)

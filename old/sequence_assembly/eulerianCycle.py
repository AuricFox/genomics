# https://towardsdatascience.com/russian-bridges-eulerian-circuits-and-genome-assembly-9fd1d84832e9

import copy
import random
import sys

def EulerianCycle(strings, format=True):
  ##Ignore this formatting block, it's only for a desired input
    if format:
        graph = [i.split(' -> ') for i in strings]
        graph = dict(graph)
        for (key, val) in graph.items():
            val = val.split(',')
            graph[key] = val

    copy_graph = copy.deepcopy(dict(graph))

    #We need to calculate the length of our graph so that we know when to stop the algorithm. This means we
    #need to calculate the number of edges in our graph.
    l = 1
    values = []
    for val in copy_graph.values():
        for i in val:
            l += 1
            values.append(i)

    #validate is set as default by true, you don't need this
    validate = True
    degrees = {key: [] for key in copy_graph.keys()}

    #Here, we're calculating the degree of the graph to determine whether each node has even degrees.
    for key in degrees.keys():
        degrees[key].append(values.count(key))
        degrees[key].append(len(copy_graph[key]))
        degrees[key].append(degrees[key][1] - degrees[key][0])

    #The final sequence is what will be our answer
    final_sequence = []

    #The cycle list will be the current cycle we are on
    cycle = []

    #Initialize a Eulerian cycle with a random node
    cn = random.choice([key for key in copy_graph.keys()])

    while len(final_sequence) != l:

        #The way we're going to pass through our graph is by removing nodes as we pass over all of their respective edges.
        #We keep on adding to the cycle as a result and the node we removed becomes our next node.

        if copy_graph[cn] != []:
            cycle.append(cn) #Add the start node to our cycle
            next_possibles = copy_graph[cn]
            new_cn = random.choice(next_possibles) #Using our dictionary, we find what nodes are connected to our current
                                                   #starting node and pick a random node to walk through next.
            #Remove the edge we've already walked through
            copy_graph[cn].remove(new_cn)
            cn = new_cn #our new node becomes our next starting node

        #A dead end of the graph is found if the node we are on has no more edges to walk through. Thus our cycle is complete.
        elif copy_graph[cn] == []:
            #We add this end node at the end of our final sequence.
            final_sequence.insert(0, cn)
            if len(cycle) == 0: #We're checking to see if we've passed through all possible cycles here
                # final_sequence.append(final_sequence[0])
                break

            #If there are still more cycles, we backtrack to the previous node and continue our walk from there
            else:
                cn = cycle[-1] #Our new starting node is the previous node -> this new node is later checked to make sure
                               #that it has unvisited edges, otherwise this is also definitely part of our final sequence
                               #and then we backtrack to the next previous node.
                #To prevent reusing this node, we remove it from the current cycle
                cycle.pop()

    return final_sequence


def main():
    kmer = sys.argv[1]
    f = open('output/temp/spike_protein_directed_graph_10.txt', 'r')

    text = []
    for line in f:
        text.append(line.replace('\n', ''))
    f.close()

    ec = EulerianCycle(text)
    print('->'.join(ec))

    filename = 'output/eulerian/eulerianCycle_' + kmer + '.txt'
    with open(filename, 'w') as f:
        f.write(''.join(ec))
    f.close()

if __name__ == '__main__':
    main()

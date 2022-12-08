# https://towardsdatascience.com/russian-bridges-eulerian-circuits-and-genome-assembly-9fd1d84832e9

import random
import copy

def EulerianPath(strings, format=True):
    #Similar formatting for turning txt file to DeBruijn graph in python dict form
    if format:
        graph = [i.split(' -> ') for i in strings]
        graph = dict(graph)
        for (key, val) in graph.items():
            val = val.split(',')
            graph[key] = val
        copy_graph = copy.deepcopy(dict(graph))
    else:
        graph = strings
        copy_graph = copy.deepcopy(graph)

    #Calculate the length the same way we did for the Eulerian cycle
    l = 1
    values = []
    for val in copy_graph.values():
        for i in val:
            l += 1
            values.append(i)

    validate = True

    #For a Eulerian path to be true, the end node must have a degree of 1, and it most cases the end node does
    #not connect to any other node with an edge, hence the graph is nearly balanced.

    end = [val for val in values if val not in copy_graph.keys()] #check if end node is connected to any other nodes

    #In some cases, this is not always true so we can check whether there is an end node or not
    if end == []:
        pass

    #If the condition is true, then we need to format our final dictionary such that we add a key in our dictionary
    #where the end node is assigned an empty list. If you remember the cycle loop in the Eulerian cycle code, this
    #empty list here will help us identify whether we have reached the end of the cycle or not.

    else:
        copy_graph[end[0]] = []
    degrees = {key: [] for key in copy_graph.keys()}

    #Regardless, we can verify our start and end nodes by checking their in-out degrees
    for key in degrees.keys():
        degrees[key].append(values.count(key))
        degrees[key].append(len(copy_graph[key]))
        degrees[key].append(degrees[key][1] - degrees[key][0])

    final_sequence = []
    cycle = []

    #cn or our current start node is determined if its degree is odd (or 1 in this case) TRY
    cn = -1
    try:
        cn = [key for key in degrees.keys() if degrees[key][2] == abs(1)][0]
    except Exception:
        pass
    finally:
        # test on other indexs
        cn = list(degrees.keys())[2:3]
        cn = str(cn[0])

    #We implement the exact same loop here as in the Eulerian cycle only with a difference where
    #we already know the starting node
    while len(final_sequence) != l:
        if copy_graph[cn] != []:
            cycle.append(cn)
            next_possibles = copy_graph[cn]
            new_cn = random.choice(next_possibles)
            copy_graph[cn].remove(new_cn)
            cn = new_cn
        elif copy_graph[cn] == []:
            final_sequence.insert(0, cn)

            #The length of our current cycle will be 0 or cycle=[] if there are no more nodes to backtrack to, therefore
            #the current node that we are on right now is the start node since it has no more outgoing edges nor incoming
            #ones that can allow us to backtrack from.

            if len(cycle) == 0:
                break
            else:
                cn = cycle[-1]
                cycle.pop()

    return final_sequence

def main():
    f = open('output/spike_protein_directed_graph.txt', 'r')

    text = []
    for line in f:
        text.append(line.replace('\n', ''))

    f.close()

    print('->'.join(EulerianPath(text)))

if __name__ == '__main__':
    main()

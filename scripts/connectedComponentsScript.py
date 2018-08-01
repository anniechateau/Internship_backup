import networkx as nx
import matplotlib.pyplot as plt
import pydot
import os


#create subdirectory structure if it doesn't already exist:



from networkx.drawing.nx_pydot import write_dot
# Dependencies to install: NetworkX + pydot + matplotlib


def seqInfoExtractor(seq):
    # function to extract team name, sequence length and
    # sequence name from the csv line
    # and fuse it into a valid NetworkX node
    rest,seqLen = seq.split("_length=")
    # joins with "-" the first two elements gained by splitting on "-"
    # little trickery needed since "-" is present in the team name XS
    seqTeam = "-".join(rest.split("-",2)[0:2])
    seqName = rest.split("-",2)[2]


    return seqName,int(seqLen),seqTeam


def statWriter(filename,statlist):
    # Function for writing a dictionary of statistic and information to csv.
    # Made to write statistics like # of teams, avg. sequence length,
    # max and min sequence lengths, ...
    # as well as binary qualifiers like if the component is a clique,
    # or contains all the teams.

    # organization:
    ## components  ,stat1,stat2,stat3
    ## component001, ... , ... , ...

    f = open(filename, "a+")

    f.write(";".join(statlist))
    f.write("\n")


    f.close()



def upSetWriter(filename,componentID,cur_teamList,teamList):
    # Function to write team list data in the UpSet compatible binary format
    # AKA., 0 for team not in said component, 1 for team in said component.
    # 1 line per component.
    f = open(filename, "a")
    f.write(str(componentID) + ";")
    data = ""
    for team in sorted(teamList):
        if team in cur_teamList:
            data += "1;"
        else:
            data += "0;"
    f.write(data.strip(";"))
    f.write("\n")
    f.close()

######## COMPLETES ###########
# make dirs if they don't exist
dir = "./connectedComponents/complete/"

if not os.path.exists(dir):
    os.makedirs(dir)


G0 = nx.Graph()


complete = "./graph/complete.csv"
inclusions = "./graph/inclusions.csv"
TotalSeqCountComplete = 0
TotalSeqCountInclusions = 0


with open(complete,"rb") as blastFile:
    for comparison in blastFile:
        # get the two node names, keeping the common stats separated
        # to be reused as edge information later on.
        seqAInfo, seqBInfo ,edge_info = comparison.split("\t",2)
        # get ind. sequence information from sequence identifiers:
        seqAName, seqALen,seqATeam = seqInfoExtractor(seqAInfo)
        seqBName, seqBLen,seqBTeam = seqInfoExtractor(seqBInfo)




        # add sequences as nodes by their names (W/O team and length)
        # nodes being a set, there's no need to check for repeat elements
        G0.add_node(seqAName)

        G0.add_node(seqBName)

        # Add associated data to each node
        G0.node[seqAName]["Team"]   = seqATeam
        G0.node[seqAName]["Length"] = seqALen

        G0.node[seqBName]["Team"]   = seqBTeam
        G0.node[seqBName]["Length"] = seqBLen

        # add edge information

        # extract edge info we want to keep (% ident, coverage, e-val, bitscore)
        edgeMetadata = edge_info.strip("\n").split("\t")
        id,cov = edgeMetadata[0:2]
        eval,bitsc =edgeMetadata[8:11]

        G0.add_edge(seqAName,seqBName,identity=float(id),coverage = int(cov),e_value = float(eval),bitscore= float(bitsc))




# Create subgraphs for each connected component, representing a set of sequences
# identified as equal.

graphs = list(nx.connected_component_subgraphs(G0, copy=True))

# number of cliques (every node is connected to every other)
cliqueCount = 0

# number of connected components (including non-cliques)
numberOfConnectedComponents_complete = len(graphs)

print("Writing Dot files for each Connected Components -- Complete\n")

# get all available teams:
teams = sorted(set(nx.get_node_attributes(G0,"Team").values()))

#create first line of UpSet data file for Complete graphs
upsetFile = "connectedComponents/upsetDataComplete.csv"
f = open(upsetFile,"w+")
teamNames = ";".join(teams)
f.write("Row;" + teamNames + "\n")
f.close()


# write header of Statistics file:

upsetstats = "connectedComponents/statsComplete.csv"
f = open(upsetstats,"w+")
statsHeader = "component;num.Nodes;num.Edges;avg.length;avg.Score;minimal.Score;maximal.Score;Clique\n"
f.write(statsHeader)
f.close()


# Write each connected component in separate DOT files.
#for every graph in the connected components list:

for i in range(0,len(graphs)):
    print("processing connected complete component #" + str(i) + "\n")

    write_dot(graphs[i], "connectedComponents/complete/ConnectedComponent%04d.dot" % i )

    # list of sequence lengths for subgraph i
    lenGraph = nx.get_node_attributes(graphs[i],"Length").values()
    # list of unique teams present in the subgraph
    cur_Teams = set(nx.get_node_attributes(graphs[i],"Team").values())
    id = "connectedComponent%04d" % i

    upSetWriter(upsetFile,id,cur_Teams,teams)


    # Determining if this specific subgraph is a clique:
    # AKA., is (1/2) * (NumNodes * (NumNodes - 1)) == NumEdges ?
    # which is the formulae that every clique fulfills.
    n = graphs[i].number_of_nodes()
    cliqueBol = graphs[i].number_of_edges() == ((n * (n - 1 ) ) / 2)
    #print("Number of edges:"+ str(graphs[i].number_of_edges()) + " versus expected for a clique: " + str( (n * (n-1)) / 2  ))
    #print("Is a clique ? :" + str(cliqueBol))

    bitscores = nx.get_edge_attributes(graphs[i],"bitscore").values()
    coverage  = nx.get_edge_attributes(graphs[i],"coverage").values()

    avgBitscore = sum(bitscores) / len(bitscores)
    maxBitscore = max(bitscores)
    minBitscore = min(bitscores)

    #statsHeader = "component;num.Nodes;num.Edges;avg.length;avg.Score;minimal.Score;maximal.Score;Clique\n"
    stats = [str(id), str(len(graphs[i].nodes)), str(len(graphs[i].edges)), str(sum(lenGraph)/len(lenGraph)), str(avgBitscore), str(minBitscore), str(maxBitscore), str(cliqueBol)]
    statWriter(upsetstats,stats)


    ## keep a sub-folder for copies of connected components that are not cliques
    ## plus another folder for those that have more than 5 sequences in them.

    if cliqueBol != True:
        dir = "./connectedComponents/complete/notclique"

        if not os.path.exists(dir):
            os.makedirs(dir)
        write_dot(graphs[i], "connectedComponents/complete/notclique/ConnectedComponent%04d.dot" % i )
    if len(lenGraph) > 5:
        dir = "./connectedComponents/complete/highSeqCount"

        if not os.path.exists(dir):
            os.makedirs(dir)

        write_dot(graphs[i], "connectedComponents/complete/highSeqCount/ConnectedComponent%04d.dot" % i )




    if len(lenGraph) < 15:


        #print("\nSet positions for graph#" + str(i) + "\n")
        ## !!! I need to study which layout is best !!!
        pos = nx.spring_layout(graphs[i])

        #print("Draw graph #"+str(i)+"\n")
        drawn = nx.draw(graphs[i], pos)
        #print("setting edge and node labels for graph #" + str(i) + "\n")
        labels_node = nx.get_node_attributes(graphs[i],'Team')
        nx.draw_networkx_labels(graphs[i], pos, labels = labels_node)
        labels_edge = nx.get_edge_attributes(graphs[i],'Identity')
        nx.draw_networkx_edge_labels(graphs[i], pos, labels_edge)

        plt.savefig("connectedComponents/complete/ConnectedComponent%04d.png" % i)
        plt.close()
    else:
        print("\nskipping visualization, too many nodes\n")




    print("\nAverage length of graph #" + str(i)  + " : " + str(sum(lenGraph)/len(lenGraph)))
    print("\nNumber of teams represented in this graph: " + str(len(cur_Teams)))
    print("\nTeams represented in this graph: " + str(cur_Teams))



print("# of connected components " + str(len(graphs)) + "\n")
#print("# of cliques " + str(cliques) + "\n")





#######  INCLUSIONS  #######

G0.clear()
G0 = nx.DiGraph()

# make dirs if they don't exist
dir = "./connectedComponents/inclusions/"

if not os.path.exists(dir):
    os.makedirs(dir)



with open(inclusions,"rb") as blastFile:
    for comparison in blastFile:
        # get the two node names, keeping the common stats separated
        # to be reused as edge information later on.
        seqAInfo, seqBInfo ,edge_info = comparison.split("\t",2)
        # get ind. sequence information from sequence identifiers:
        seqAName, seqALen,seqATeam = seqInfoExtractor(seqAInfo)
        seqBName, seqBLen,seqBTeam = seqInfoExtractor(seqBInfo)




        # add sequences as nodes by their names (W/O team and length)
        # nodes being a set, there's no need to check for repeat elements
        G0.add_node(seqAName)

        G0.add_node(seqBName)

        # Add associated data to each node
        G0.node[seqAName]["Team"]   = seqATeam
        G0.node[seqAName]["Length"] = seqALen

        G0.node[seqBName]["Team"]   = seqBTeam
        G0.node[seqBName]["Length"] = seqBLen

        if seqBLen > seqALen:
            G0.node[seqAName]["Type"] = "inclusion"
            incl  = seqAName
            G0.node[seqBName]["Type"] = "main"
            main = seqBName
        else:
            G0.node[seqBName]["Type"] = "inclusion"
            incl  = seqBName
            G0.node[seqAName]["Type"] = "main"
            main = seqAName



        # add edge information

        # extract edge info we want to keep (% ident, coverage, e-val, bitscore)
        edgeMetadata = edge_info.strip("\n").split("\t")
        id,cov = edgeMetadata[0:2]
        eval,bitsc =edgeMetadata[8:11]

        G0.add_edge(incl,main,Identity=float(id),Coverage = int(cov),e_value = float(eval),bitscore= float(bitsc))



# Create subgraphs for each connected component, representing a set of sequences
# identified as equal.

graphs = list(nx.weakly_connected_component_subgraphs(G0, copy=True))

numberOfWeaklyConnectedComponents_inclusions = len(graphs)


print("Writing Dot files for each Connected Components -- Inclusions\n")

# Write each connected component in separate DOT files.

#create first line of UpSet data file
teams = set(nx.get_node_attributes(G0,"Team").values())

upsetFile = "connectedComponents/upsetDataInclusions.csv"
f = open(upsetFile,"w+")
teamNames = ";".join(teams)
f.write("Row;" + teamNames + "\n")
f.close()



#for every graph in the connected components list:

for i in range(0,len(graphs)):
    print("processing connected inclusion component #" + str(i) + "\n")
    write_dot(graphs[i], "connectedComponents/inclusions/ConnectedComponent%04d.dot" % i )

    lenGraph = nx.get_node_attributes(graphs[i],"Length").values()
    uniqueTeams = set(nx.get_node_attributes(graphs[i],"Team").values())

    n = graphs[i].number_of_nodes()
    cliqueBol = graphs[i].number_of_edges() == ((n * (n - 1 ) ) / 2)

    print("Number of edges:"+ str(graphs[i].number_of_edges()) + " versus expected for a clique: " + str( (n * (n - 1)) /2  ))
    print("Is a clique ? :" + str(cliqueBol))

    # list of unique teams present in the subgraph
    cur_Teams = set(nx.get_node_attributes(graphs[i],"Team").values())
    id = "connectedComponent%04d" % i

    upSetWriter(upsetFile,id,cur_Teams,teams)



    bitscores = nx.get_edge_attributes(graphs[i],"bitscore").values()
    coverage  = nx.get_edge_attributes(graphs[i],"coverage").values()

    avgBitscore = sum(bitscores) / len(bitscores)
    maxBitscore = max(bitscores)
    minBitscore = min(bitscores)

    #statsHeader = "component;num.Nodes;num.Edges;avg.length;avg.Score;minimal.Score;maximal.Score;Clique\n"
    stats = [str(id), str(len(graphs[i].nodes)), str(len(graphs[i].edges)), str(sum(lenGraph)/len(lenGraph)), str(avgBitscore), str(minBitscore), str(maxBitscore), str(cliqueBol)]
    statWriter(upsetstats,stats)



    if len(lenGraph) < 30:


        #print("\nSet positions for graph#" + str(i) + "\n")
        ## !!! I need to study which layout is best !!!
        pos = nx.spring_layout(graphs[i])

        #print("Draw graph #"+str(i)+"\n")
        drawn = nx.draw(graphs[i], pos)
        #print("setting edge and node labels for graph #" + str(i) + "\n")
        labels_node = nx.get_node_attributes(graphs[i],'Team')
        nx.draw_networkx_labels(graphs[i], pos, labels = labels_node)
        labels_edge = nx.get_edge_attributes(graphs[i],'Identity')
        nx.draw_networkx_edge_labels(graphs[i], pos, labels_edge)

        plt.savefig("connectedComponents/inclusions/ConnectedComponent%04d.png" % i)
        plt.close()
    else:
        print("skipping visualization, too many nodes\n")



    ## keep a sub-folder for copies of connected components that are not cliques
    ## plus another folder for those that have more than 5 sequences in them.

    if cliqueBol != True:

        # make dirs if they don't exist
        dir = "./connectedComponents/inclusions/notClique/"

        if not os.path.exists(dir):
            os.makedirs(dir)

        write_dot(graphs[i], "connectedComponents/inclusions/notClique/ConnectedComponent%04d.dot" % i )


    if len(lenGraph) > 5:

        # make dirs if they don't exist
        dir = "./connectedComponents/inclusions/highSeqCount"

        if not os.path.exists(dir):
            os.makedirs(dir)

        write_dot(graphs[i], "connectedComponents/inclusions/highSeqCount/ConnectedComponent%04d.dot" % i )




    print("Average length of graph #" + str(i)  + " : " + str(sum(lenGraph)/len(lenGraph)))
    print("\nNumber of teams represented in this graph: " + str(len(uniqueTeams)))
    print("\nTeams represented in this graph: " + str(uniqueTeams))


print("# of connected components " + str(len(graphs)) + "\n")
#print("# of cliques " + str(cliques) + "\n")

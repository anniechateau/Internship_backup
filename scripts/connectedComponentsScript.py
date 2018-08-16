import networkx as nx  # necesitates networkX install
import matplotlib.pyplot as plt
import pydot
import os
from itertools import combinations  # easy combinatorics, default library.

# necessitates biopython install (using pip)
from Bio.Blast.Applications import NcbiblastpCommandline
from StringIO import StringIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

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


def completeGraph(completeCSV):
    # function takes the CSV file for the long matches ("complete")
    # and fills the graph using the data.
    # then returns said graph

    # initialize graph
    G0 = nx.Graph()

    # open the CSV file and process it into a Graph object of the correct type.
    with open(completeCSV,"rb") as blastFile:
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

    return G0

def inclusionsGraph(inclusionsCSV):
        # function takes the CSV file for the long matches ("complete")
        # and fills the graph using the data.
        # then returns said graph

    # initialize graph
    G0 = nx.DiGraph()



    with open(inclusionsCSV,"rb") as blastFile:
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

    return G0

def drawGraph(subgraph,outputDir,id):
    #This function calculates the optimal positions for each node,
    # then adds extra information like team names and identity % to the graph.
    # and saves it to a file

    pos = nx.spring_layout(subgraph)

    #print("Draw graph #"+str(i)+"\n")
    drawn = nx.draw(subgraph, pos)
    #print("setting edge and node labels for graph #" + str(i) + "\n")
    labels_node = nx.get_node_attributes(subgraph,'Team')
    nx.draw_networkx_labels(subgraph, pos, labels = labels_node)
    labels_edge = nx.get_edge_attributes(subgraph,'Identity')

    nx.draw_networkx_edge_labels(subgraph, pos, labels_edge)

    plt.savefig(outputDir + "/"  + id + ".png")
    plt.close()

def nodes_connected(subgraph,u, v):
    # checks if nodes u and v are neighbours.
    # AKA., do they share an edge ?
    return u in subgraph.neighbors(v)

def findFastaSeqforID(fastaFile, id,seqName):
    # parses a Fasta file, looking for a partial match in the ID,
    # returning the corresponding sequence and length, otherwise returning False
    # if no match is found.

    for seq_record in SeqIO.parse(fastaFile, "fasta"):
        match = False
        curID = seq_record.id
        if id in curID:
            match = True
            SeqIO.write(seq_record, (seqName + ".fasta"), "fasta")
            break


def listFastaFiles(dir):
    # returns a list containing the tagged fasta files we want.
    fileList = []
    for root, dirs, files in os.walk("."):
        for filename in files:
            if "_tagged.fasta" in filename:
                if not "nhr" in filename and not "nin" in filename and not "nsq" in filename:
                    fileList.append("./tagged_fasta/" + filename)
    return fileList



def blast(u,v,fastaDir):
    # takes the names of the two nodes,
    # finds their corresponding sequences in the initial fasta files
    # and BLASTs one sequence against the other,
    # (blasting the smallest on the biggest sequence).

    # returns the resulting blast results in a biopython object
    fastaFiles = listFastaFiles(fastaDir)
    seqName1 = "sequence1"
    seqName2 = "sequence2"
    for fastaFile in fastaFiles:
        findFastaSeqforID(fastaFile, u,seqName1)
        findFastaSeqforID(fastaFile, v,seqName2)

    output = NcbiblastpCommandline(query= (seqName1 + ".fasta"), subject= (seqName2 + ".fasta"), outfmt=5, out="comparison.xml")()[0]
    blast_result_record = NCBIXML.parse(open('comparison.xml', 'r'))
    return blast_result_record

def addEdgesToSubgraph(subgraph,fastaDir):
    # find every unique 2-tuple combination of nodes in the subgraph.

    # adds an edge between all combinations U,V that do not already posses one,
    # using Blast to get the informtion necessary to establish the clique.
    # Essentially completing the subgraph and making it a Clique.
    # (if a valid alignment is available between the two entities, of course)
    combinationsList = [[u,v] for u,v in combinations(subgraph.nodes, 2)]

    for u,v in combinationsList:
        if not nodes_connected(subgraph,u,v):
            blast_result_record = blast(u,v,fastaDir)
            for query in blast_result_record:
                for alignment in query.alignments:
                    for hsp in alignment.hsps:
                        id = float(hsp.identities) / float(hsp.align_length)
                        cov = float(hsp.align_length) / float(query.query_length)
                        eval = hsp.expect
                        bitsc = hsp.score
                        subgraph.add_edge(u,v,Identity=float(id),Coverage = int(cov),e_value = float(eval),bitscore= float(bitsc))






def processInclusionSubgraps(graphs,teams, inclusionDir,fastaDir):
    for i in range(0,len(graphs)):
        print("processing connected inclusion component #" + str(i) + "\n")
        id = "connectedComponent%04d" % i
        write_dot(graphs[i], (inclusionDir + "/" + id + ".dot"))

        lenGraph = nx.get_node_attributes(graphs[i],"Length").values()
        uniqueTeams = set(nx.get_node_attributes(graphs[i],"Team").values())

        n = graphs[i].number_of_nodes()
        cliqueBol = graphs[i].number_of_edges() == ((n * (n - 1 ) ) / 2)

        print("Number of edges:"+ str(graphs[i].number_of_edges()) + " versus expected for a clique: " + str( (n * (n - 1)) /2  ))
        print("\nIs a clique ? :" + str(cliqueBol) + "\n")

        # list of unique teams present in the subgraph
        cur_Teams = set(nx.get_node_attributes(graphs[i],"Team").values())

        upSetWriter(upsetFile,id,cur_Teams,teams)



        bitscores = nx.get_edge_attributes(graphs[i],"bitscore").values()
        coverage  = nx.get_edge_attributes(graphs[i],"coverage").values()

        avgBitscore = sum(bitscores) / len(bitscores)
        maxBitscore = max(bitscores)
        minBitscore = min(bitscores)

        #statsHeader = "component;num.Nodes;num.Edges;avg.length;avg.Score;minimal.Score;maximal.Score;Clique\n"
        stats = [str(id), str(len(graphs[i].nodes)), str(len(graphs[i].edges)), str(sum(lenGraph)/len(lenGraph)), str(avgBitscore), str(minBitscore), str(maxBitscore), str(cliqueBol)]
        statWriter(upsetstats,stats)



        if len(lenGraph) < 15:

            drawGraph(graphs[i],inclusionDir,id)

        else:

            print("skipping visualization, too many nodes\n")



        ## keep a sub-folder for copies of connected components that are not cliques
        ## plus another folder for those that have more than 5 sequences in them.

        if cliqueBol != True:

            # make dirs if they don't exist
            createDir(inclusionDir + "/notClique/")

            write_dot(graphs[i], (inclusionDir +"/notClique/" + id + "_before.dot") )


            addEdgesToSubgraph(graphs[i],fastaDir)

            write_dot(graphs[i], (inclusionDir +"/notClique/" + id + "_after.dot") )


        if len(lenGraph) > 5:

            # make dirs if they don't exist
            createDir(inclusionDir + "/highSeqCount/")

            write_dot(graphs[i], (inclusionDir +"/highSeqCount/" + id + ".dot") )




        print("Average length of graph #" + str(i)  + " : " + str(sum(lenGraph)/len(lenGraph)))
        print("\nNumber of teams represented in this graph: " + str(len(uniqueTeams)))
        print("\nTeams represented in this graph: " + str(uniqueTeams))


    print("# of connected components " + str(len(graphs)) + "\n")
    #print("# of cliques " + str(cliques) + "\n")


def processCompleteSubgraphs(graphs,teams,outputDir,fastaDir):

    for i in range(0,len(graphs)):
        print("processing connected complete component #" + str(i) + "\n")
        id = "connectedComponent%04d" % i
        write_dot(graphs[i], (outputDir + "/" + id + ".dot") )

        # list of sequence lengths for subgraph i
        lenGraph = nx.get_node_attributes(graphs[i],"Length").values()
        # list of unique teams present in the subgraph
        cur_Teams = set(nx.get_node_attributes(graphs[i],"Team").values())

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




        if len(lenGraph) < 15:

            drawGraph(graphs[i],outputDir,id)

        else:

            print("\nskipping visualization, too many nodes\n")


        ## keep a sub-folder for copies of connected components that are not cliques
        ## plus another folder for those that have more than 5 sequences in them.
        if cliqueBol != True:
            createDir(outputDir + "/notclique")
            write_dot(graphs[i], outputDir + "/notclique/" + id + ".dot")
            # do a BLAST of the nodes that aren't connected to each other,
            # and save the result in a new DOT file, as well as create a new image of it.
            drawGraph(graphs[i],outputDir + "/notclique/",id + "_before")
            addEdgesToSubgraph(graphs[i],fastaDir)
            drawGraph(graphs[i],outputDir + "/notclique/",id + "_after")


        if len(lenGraph) > 5:
            createDir(outputDir + "/highSeqCount")
            write_dot(graphs[i], (outputDir + "/highSeqCount/" + id + ".dot" ))





        print("\nAverage length of graph #" + str(i)  + " : " + str(sum(lenGraph)/len(lenGraph)))
        print("\nNumber of teams represented in this graph: " + str(len(cur_Teams)))
        print("\nTeams represented in this graph: " + str(cur_Teams))



    print("# of connected components " + str(len(graphs)) + "\n")
    #print("# of cliques " + str(cliques) + "\n")




def writeHeader(csvFile,data):
    f = open(csvFile,"w+")
    teamNames = data
    f.write(teamNames)
    f.close()

def createDir(dir):
    # create the directory and sub-directories if they don't exist.
    if not os.path.exists(dir):
        os.makedirs(dir)

#def connectedComponentsWriter()

######## COMPLETES ###########
# make dirs if they don't exist
outputDir = "./connectedComponents/complete/"
fastaDir = "./tagged_fasta"
# if path doesn't exist, create it
createDir(outputDir)


completeCSV = "./graph/complete.csv"

G0 = completeGraph(completeCSV)


# Create subgraphs for each connected component, representing a set of sequences
# identified as equal.
graphs = list(nx.connected_component_subgraphs(G0, copy=True))


# number of cliques (every node is connected to every other)
cliqueCount = 0

# get all available teams:
teams = sorted(set(nx.get_node_attributes(G0,"Team").values()))

#create first line of UpSet data file for Complete graphs
upsetFile = "connectedComponents/upsetDataComplete.csv"
header = "Row;" + ";".join(teams) + "\n"
writeHeader(upsetFile,header)


# write header of Statistics file:

upsetstats = "connectedComponents/statsComplete.csv"
statsHeader = "component;num.Nodes;num.Edges;avg.length;avg.Score;minimal.Score;maximal.Score;Clique\n"
writeHeader(upsetstats,statsHeader)



# number of connected components (including non-cliques)
numberOfConnectedComponents_complete = len(graphs)

print("Writing Dot files for each Connected Components -- Complete\n")

# Write each connected component in separate DOT files.
#for every graph in the connected components list:

processCompleteSubgraphs(graphs,teams,outputDir,fastaDir)





#######  INCLUSIONS  #######


inclusionsCSV = "./graph/inclusions.csv"


# make dirs if they don't exist
inclusionDir = "./connectedComponents/inclusions/"

createDir(inclusionDir)

# purge the Complete Graph data.
G0.clear()

# input the data for the Inclusions graph.
G0 = inclusionsGraph(inclusionsCSV)






# Create subgraphs for each connected component, representing a set of sequences
# identified as equal.

graphs = list(nx.weakly_connected_component_subgraphs(G0, copy=True))

numberOfWeaklyConnectedComponents_inclusions = len(graphs)


print("Writing Dot files for each Connected Components -- Inclusions\n")

# Write each connected component in separate DOT files.

#create first line of UpSet data file
teams = set(nx.get_node_attributes(G0,"Team").values())

upsetFile = "connectedComponents/upsetDataInclusions.csv"
header = "Row;" + ";".join(teams) + "\n"
writeHeader(upsetFile,header)


# write stats statsHeader
statFile = "connectedComponents/inclusionsStats.csv"
header = "component;num.Nodes;num.Edges;avg.length;avg.Score;minimal.Score;maximal.Score;Clique\n"
# adjust stats for inclusions specific charaters
writeHeader(statFile,header)


#for every graph in the connected components list:

processInclusionSubgraps(graphs,teams,inclusionDir,fastaDir)

import glob, os

# get fasta files with extension for snakemake
SAMPLES = [os.path.basename(f) for f in glob.glob('initialFastaFiles/*')]

# extract the names themselves from the generated list
names = []
for name in SAMPLES:
    names.append(name.split(".")[0])

# use said generated list to generate the names of our BLast result files,
#so that snakemake knows what should be expected in the results folder.
blastNames = []

for i in range(0,len(names)):
    name1 = names[i]
    for j in range(i,len(names)):
        name2 = names[j]
        if (name1 != name2):
            csvName = name1 + "_vs_" + name2 + ".csv"
            blastNames.append(csvName)




# Rule all needs to have as input the final output files, so that
# Snakemake knows how to organize the execution of the rules to produce
# the wanted output files.
rule all:
    input:
        expand("preprocessed_fasta/{ID}", ID=SAMPLES),
        expand("tagged_fasta/{ID}", ID=SAMPLES),
        expand("tagged_fasta/{TAGGED}_tagged.fasta",TAGGED=names),
        expand("tagged_fasta/{ID}_tagged.fasta.nhr", ID=names),
        expand("tagged_fasta/{ID}_tagged.fasta.nin", ID=names),
        expand("tagged_fasta/{ID}_tagged.fasta.nsq", ID=names),
        expand("blastResults/{BLASTNAMES}",BLASTNAMES=blastNames),
        #"graph/filteredComplete.csv"
        "graph/complete.csv",
        "graph/inclusions.csv",
        "graph_step.txt"
    message:
        "Checking that all final output files exist.\n"


rule preprocess_fasta:
    input:
        "initialFastaFiles/{names}.fasta"
    output:
        "preprocessed_fasta/{names}.fasta"
    message:
        "Preprocessing {names} file"
    shell:
        "sed 's/\t/_/g' {input} > {output}"




rule fasta:
    input:
        expand("preprocessed_fasta/{FASTA}.fasta", FASTA=names)

rule copyFasta:
    input:
        "preprocessed_fasta/{FASTA}.fasta"
    output:
        "tagged_fasta/{FASTA}.fasta"
    shell:
        "cp {input} {output}"



rule tagged:
    input:
        expand("tagged_fasta/{TAGGED}.fasta", TAGGED=names)

rule prefix_ID_to_Sequences:
    input:
        "tagged_fasta/{TAGGED}.fasta"
    output:
        "tagged_fasta/{TAGGED}_pretagged.fasta"
    shell:
        # wildcards.tagged is the syntax to access one of the wildcards inside a shell command.
        "sed 's/>/>{wildcards.TAGGED}-/g' {input} > {output}"


rule addSequenceLengths:
    input:
        "tagged_fasta/{TAGGED}_pretagged.fasta"
    output:
        "tagged_fasta/{TAGGED}_tagged.fasta"
    shell:
        "cat {input} | seqkit fx2tab --length | awk -F \"\t\" '{{print $1\"-sequence_length=\"$4\"\t\"$2}}' | seqkit tab2fx > {output}"



rule formatdb:
    input:
        "tagged_fasta/{names}_tagged.fasta"
    output:
        "tagged_fasta/{names}_tagged.fasta.nhr",
        "tagged_fasta/{names}_tagged.fasta.nin",
        "tagged_fasta/{names}_tagged.fasta.nsq"
    message:
        "Creating database files for {input}\n"

    shell:
        "formatdb -p F -i {input}"

# Given we want to run the Megablast command for every possible combination of
# Database and Subject in the preprocessed_fasta directory,
# We create an anonymous snakemake rule that uses said combinations to run the megablast.
for i in range(0,len(names)):
    A = names[i]
    for j in range(i,len(names)):
        B = names[j]
        if (A != B):


            # use of an anonymous rule to avoid problems when rerunning the rule
            # for each blast of A vs B.
            rule :
                # input needs all DB-related files for A.fasta as well as the subject B.fasta file.
                input:
                    expand("tagged_fasta/{db}_tagged.fasta", db=A),
                    expand("tagged_fasta/{db}_tagged.fasta.nhr", db=A),
                    expand("tagged_fasta/{db}_tagged.fasta.nin", db=A),
                    expand("tagged_fasta/{db}_tagged.fasta.nsq", db=A),
                    expand("tagged_fasta/{subject}_tagged.fasta", subject=B)
                output:
                    expand("blastResults/{db}_vs_{subject}.csv", db=A, subject=B)
                message:
                    "megablasting {input[0]} vs {input[4]}\n"
                shell:
                    # A double free / Memory coruption error occurs when number of processors used > 1
                    "megablast -d {input[0]} -i {input[4]}  -o {output} -q 5 -r -4 -D 3 -a 1 -W 32 -m 8"


rule Filter_Blast_Results:
    input:
        expand("blastResults/{BLASTNAMES}", BLASTNAMES=blastNames)
    output:
        "graph/complete.csv",
        "graph/inclusions.csv"
    message:
        "Filtering blast results to keep those with coverage > 75%, identity > 95% and subj/query ~= 1 "
    shell:
        # The conditions are:
        # AAA || BBB   <==> if subj.len / query.len close to 1 (AKA., it's not an inclusion)
        # AND Coverage % (coverage.len/Query.len) > 0.75
        # AND identity % > 95%
        # AND E.Val < 1*10^-15
        #"awk '/^[^#]/ {{if ( (substr($1,index($1,\"sequence_length=\")+16) / (substr($2,index($2,\"sequence_length=\")+16))) < 1.5 || (substr($1,index($1,\"sequence_length=\")+16) / (substr($2,index($2,\"sequence_length=\")+16))) > 0.5 && $4/(substr($1,index($1,\"sequence_length=\")+16)) > 0.75 && $3 > 95 && $11 < 0.000000000000001) print $0}}' {input} > {output}"
        "./scripts/blastFilter.awk {input}"



# new rule to execute the Graph generation and analysis process
rule graphProcess:
    input:
        "graph/complete.csv",
        "graph/inclusions.csv"
    output:
        "graph_step.txt"
    shell:
        "python ./scripts/connectedComponentsScript.py;touch graph_step.txt"

#download WikiPathways human subset from https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#C2
#download WikiPathways mouse subset from https://www.gsea-msigdb.org/gsea/msigdb/mouse/collections.jsp?targetSpeciesDB=Mouse#M2
# run script as python map_incytr_2_wikipathway.py human allpath.csv

import sys
import os

def get_WikiPathways(inputFile):
    """
    Reads a file containing MSigDB pathways and returns a dictionary
    mapping pathway names to sets of genes.
    
    Returns:
        dict: A dictionary where keys are pathway names and values are sets of genes.
    """
    #get MSigDB pathway
    pathways={}
    try:
         # Open the file and read lines
         with open(inputFile, 'r') as msigdb:
             for line in msigdb:
                 elements=line.strip().split("\t")
                 pathway_name=elements[0] # First column is the pathway name
                 genelist=elements[2:] # Genes start from the third column
                 pathways[pathway_name]=set(genelist)  # Store as a set for fast lookup
    except FileNotFoundError:
        print(f"Error: File '{inputFile}' not found.")
        raise
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        raise
    return pathways
    
def map_2_pathway(pathways, Incytrpath):
    path=Incytrpath.split("*") 
    myset=set()   
    for gene in path:
        myset.add(gene) 
    mappedPath=""   
    for pathway, genes in pathways.items():
        if myset.issubset(genes):
            if not mappedPath:
                mappedPath=pathway
                    #print(mappedPath)
            else:
                mappedPath = mappedPath + "; " + pathway  
    return mappedPath; 
    
def get_column_index(columnname, title_line):
    columnList=title_line.replace("\"", "").split("\t")
    if columnname in columnList:
        indx=columnList.index(columnname) 
        return indx 
    else:
        raise ValueError(f"Column '{columnname}' not found in the title line.")     
       
def main():
    inputSpecies=sys.argv[1].strip()
    incytrFile=sys.argv[2].strip()
    #pathwayFile
    if inputSpecies.lower() == 'human' or inputSpecies.lower() == 'h':
        pathwayFile='c2.cp.wikipathways.v2024.1.Hs.symbols.gmt'    
    elif inputSpecies.lower() == 'mouse' or inputSpecies.lower() == 'm':
        pathwayFile='m2.cp.wikipathways.v2024.1.Mm.symbols.gmt'    
    #read in WikiPathways
    pathways=get_WikiPathways(pathwayFile)
    try:
        #read in Incytr pathways
        with open(incytrFile, 'r') as Incytr:
            # Split the file path into base name and extension
            base_name, extension = os.path.splitext(incytrFile)
            outfile_name=base_name + "_WikiPathways.tsv"
            # Open the output file
            with open(outfile_name, "w") as outfile:
                # Read the title line
                title_line=next(Incytr).strip().replace(",", "\t")
                outfile.write(title_line + "\tWikiPathways\n")
                # Get the index of the "Path" column
                path_index=get_column_index("Path", title_line)
                # Process the rest of the lines
                for line in Incytr: 
                    line=line.strip().replace(",", "\t")
                    elements=line.replace("\"", "").split("\t")
                    path=elements[path_index].strip()
                    mappedPath=map_2_pathway(pathways, path)
                    outfile.write(line + "\t" + mappedPath + "\n") 
    except FileNotFoundError:
        print(f"Error: File '{incytrFile}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")           
               
    Incytr.close() 
    outfile.close()
    
if __name__ == "__main__":
    main()    
  
   

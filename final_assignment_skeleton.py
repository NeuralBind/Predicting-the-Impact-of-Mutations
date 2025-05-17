#!/usr/bin/python
# final_assignment_skeleton.py alignments/rno-mmu.sif GO/rno.go GO/mmu.go mapping/rno.map mapping/mmu.map as an example on how to run in the command line
#here we iport the sys arguments lib
import sys


def get_mapping(map_file):
    # Open the file.
    f = open(map_file, "r")

    # Result is a list of dictionaries.
    header = []
    words=[]
    # Skip the header on the first line.
    header = f.readline()
    mapping_list = []
    
    header= header.split("\t")
    
    # Checking how many columns there are and then append the same number of dictionaries -1
    for x in range(len(header)-1):
        
        mapping_list.append({})
    # Read the 1st line of our file
    for line in f:
        # Take out new line characters
        words = line.strip() 
        # Split tab values
        words = words.split("\t") 
        
        # Skip headers
        
        if (line == header):
            continue
        else:
                
                # Add the ENS values to the keys in the mapping list
                for i in range (1,len(words)):
                    
                    mapping_list[i-1][words[i]]=words[0]            
            
    # Use .pop to remove the ENS values without keys 
    for k in range (len(mapping_list)):
        
        mapping_list[k].pop("",None)

    # Remember to close the file after we're done.
    f.close()

    return mapping_list

def get_go_terms(mapping_list, go_file):
    # Open the file.
    f = open(go_file, "r")
    # This will be the dictionary that this function returns.
    # Entries will have as a key an EnS ID and the value will
    # be a set of GO terms.
    go_dict = dict()
    info = "!" # For the comments in the GO file
    go_set = {} # Initialise and then make it a set 
    go_set = set()
    words = ""
    go_list=[] # Save the tab delimited values in this list
    id_checker= [] #I made this to check each ID of each line of the GO file
    j=0 # Works as a counter
    ENSid = "" # Save ENS IDs
    flag="green" 
    
    for line in f:
        # Check if the line starts with ! which is a comment on the GO file
        if (line[0] == info):
            continue
            
        else:
            
            # Strip the new line character and split the text info by sorted by tab inside a list
            words = line.strip()
            go_list = (words.split("\t"))                
            
            # I append the IDs into a list
            id_checker.append(go_list[1])

            # We check the current line with the previous if they share the same ID, if not process with the sets
            if (id_checker[j] != id_checker[j-1]):
            # A way to check if we found matches and our sets are not empty & also skip the first time running
                    if (go_set != set()): 
                        
                        if (flag != "red"):
                            # Initialise the dict's key and value
                            go_dict[ENSid]=(go_set)
                            # Make and empty set and ENSid and move to make the next set
                            go_set= set()
                            ENSid = ""
                        elif (flag == "red"):
                            # For the IDS in the file that are not the common IDS or for same ids that are not in order and try to overwrite my set
                               go_dict[ENSid].update(go_set)
                               flag = "green" 
                               go_set= set()
                               ENSid = ""
                               
            # Search through the mapping least if we can find the ID from the GO file, if we do then add go term to the set
            for i in range (len(mapping_list)):
                
                if go_list[1] in mapping_list[i]:
                    # Here we save the ENS ID related to matching ID of the GO and the Mapping file
                    ENSid = mapping_list[i][go_list[1]]
                    # Add the new match in the set (the number after the GO:)
                    go_set.add(go_list[4][3:])
            # If there is already a value in this key, we make flag= "red" in order for the script to understand that we have to update the old set with the new match and not overwrite it       
            if go_dict.get(ENSid):
                flag = "red"     
          
            # Add j+1 to the counter and go to the next line 
            j+=1  
    # I sorted the items with the keys because i checked thats how the codegrade wants it, but i guess it is not needed      
    # go_dict = {k: v for k, v in sorted(list(go_dict.items()))}   
     
    # Remember to close the file after we're done.
    f.close()
    # Return the value
    return go_dict


def compute_score(alignment_file, go_one_dict, go_two_dict):
    # Open the file.
    f = open(alignment_file, "r")

    # Keep track of the number of proteins we can't map to GO terms
    # and the score.
    unmappable_one = 0
    unmappable_two = 0
    score = 0.0
    ID1 = ''
    ID2 = ''
    words = ''
    Not_match = []
    ID_list = []
    
    for line in f:
        # Seperate the id of the columns and the strip the new line character
        words = line.strip()
        ID_list = (words.split("\t"))
        # From the 1st map/go
        ID1 = ID_list[0] 
        # 2nd map/go
        ID2 = ID_list[1] 
        
        # We have to check if the IDs of the alignment file exist in our GO dictionaries
        if (ID1 in go_one_dict) and (ID2 in go_two_dict):
            # Calculate Union of the dicts' sets
            Not_match = go_one_dict[ID1] | go_two_dict[ID2] 
            # Give the matching set values (Intersection)
            match = go_one_dict[ID1] & go_two_dict[ID2] 
            # Calculate the score (Jaccard index)
            score = (len(match)/(len(Not_match))) + score 
        # Here we check  and trace the unmappables  
        else:     
            if ID1 not in go_one_dict:
                unmappable_one+=1
            if ID2 not in go_two_dict:
                unmappable_two+=1
        
    # Remember to close the file after we're done.
    f.close()

    # Return the unmappables and the score back so the main code
    return unmappable_one, unmappable_two, score


def main():
     
    #Here use the arguments of the script call in the command line, each of them is a path for alignent file, maping and GO files
    align_path = sys.argv[1]
    go_path1 = sys.argv[2] 
    go_path2 = sys.argv[3]
    map_path1 = sys.argv[4] 
    map_path2 = sys.argv[5]
    # I do a small check to see if our arguments are correct
    if ('.map' not in (map_path1 and map_path2)) or ('.go' not in (go_path1 and go_path2)) or ('.sif' not in align_path) or (len(sys.argv)!=5) :
        
        sys.exit(r'ERROR: filename in the input file path arguments should contain argv[1]:Path .sif for Alignment, argv[2],[3] Path .map for the Mapping and argv[4][5] Path .go  for the GO file extensions!!!                                   !!!The first Map file should match the 1st GO file!!!')
    
    else:
        # Not much here, get the returned values and use them for the next step/function
        mapping_list = get_mapping(map_path1)
        go_one_dict = get_go_terms(mapping_list, (go_path1))
        mapping_list = get_mapping(map_path2)
        go_two_dict = get_go_terms(mapping_list, (go_path2))
        unmap1, unmap2, score_rno_mmu = compute_score(align_path, go_one_dict, go_two_dict)
        # Print the alignment the unmapped and the Score
        print(align_path,"\nUnmappable 1: ", unmap1, "\nUnmappable 2: ", unmap2, "\nScore: ",score_rno_mmu)
        
if __name__ == '__main__':
    main()



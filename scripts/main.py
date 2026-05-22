import random
#פונקציות#
#------------------------------------------------
def seq_to_string(seq):#הופך את הליסט לסטרינג
    seq_string = "".join(seq)
    return seq_string
#------------------------------------------------

def is_conserved(amino_acid_list):
    #מיפוי תכונות של חומצות אמינו
    #amino_groups = [['R','H','K'],['D','E'],['S', 'T', 'N', 'Q'],
    #                ['P', 'C', 'G'],['A', 'V', 'I','L', 'M', 'F','Y', 'W']]
    
    #counter_list = [0, 0, 0, 0, 0]
    
    """
    #1
    Electric_charged_sidechains_P = ['R','H','K']
    #2
    Electric_charged_sidechains_N = ['D','E']
    #3
    Electric_polar_uncharged_sidechains = ['S', 'T', 'N', 'Q']
    #4
    special_cases = ['P', 'C', 'G']
    #5
    Hydrophobic_sidechains = ['A', 'V', 'I','L', 'M', 'F','Y', 'W']
    """
    
    cnt = 0
    for i in range(amino_acid_list):
        if amino_acid_list[i] == "-":
            cnt += 1

    if (len(set(amino_acid_list)) <= 2 or cnt <=3):
        return 1
    else:
        return 0
        
#------------------------------------------------
def max_seq (zero_or_one, list):
    max_len = 0
    start = -1
    stop = -1

    i = 0
    while i < len(list):
        if list[i] == zero_or_one:
            temp_start = i
            cnt  = 0
            
            while i < len(list) and list[i] and list[i] == zero_or_one:
                cnt += 1
                i += 1

            if cnt > max_len:
                max_len = cnt
                start = temp_start
                stop = i - 1
        else:
            i += 1


    return max_len, start, stop

#------------------------------------------------

#------------------------------------------------
def compare(zero_and_one_list, animle_list):
    len_conserved, start1, stop1 = max_seq(zero_and_one_list, 1)
    conserved_seq_original = animle_list[start1: stop1]
    num_mutation1 = int(len_conserved * 0.2)
    conserved_seq_mutated = Mutate_protein(conserved_seq_original, num_mutation1)


    len_non_conserved, start0, stop0 = max_seq(zero_and_one_list, 0)
    non_conserved_seq_original = animle_list[start0: stop0]
    num_mutation0 = int(len_non_conserved * 0.2)
    non_conserved_seq_mutated = Mutate_protein(non_conserved_seq_original, num_mutation0)
#------------------------------------------------
def Mutate_protein(seq, num_mutation):#מקבלת רצף חומצות אמינו של חיה ומחזירה את אותו רצף רק עם מוטצות נוקדתיות של החלפה
  '''
  '''
  for i in range(num_mutation):
    amini_acids_list = ['F', 'S', 'Y', 'C', 'L', '*', 'W',
                        'P', 'H', 'R', 'Q',
                        'I', 'T', 'N', 'K', 'M',
                        'V', 'A', 'D', 'G', 'E']

    rand_acid = random.choice(amini_acids_list)
    rand_num = random.randrange(0,len(seq))
    
    if seq[rand_num] != rand_acid:
        change_protein = seq[0:rand_num]+ rand_acid + seq[(rand_num+1):]
    
    else:
        amini_acids_list.remove(rand_acid)
        rand_acid = random.choice(amini_acids_list)
        change_protein = seq[0:rand_num]+ rand_acid + seq[(rand_num+1):]
  return change_protein
#------------------------------------------------
def file_to_list(file):
    curr_seq = ""
    seq_list = []
    
    for line in file:
        line = line.rstrip('\r\n')

        # רצף החומצות אמינו מופיע בשורות שאינן מתחילות בסימן "<" לכן "נדלג" על שורה זו
        if line == "" or line[0] == ">":
            if curr_seq != "":
                seq_list.append(curr_seq)
            curr_seq = ""

        else:
            curr_seq += line

    if curr_seq != "":
        seq_list.append(curr_seq)
    
    return seq_list
#------------------------------------------------
def conserved_list(seq_str):#
    one_zero_list=[]
    seq_length=len(seq_str)
    for i in range (seq_length):
        x=is_conserved(seq_str)
        one_zero_list.append(x)
    return one_zero_list

#------------------------------------------------
def mutate_ustage(organizim_list):# מעביר את הרצף 70 מוטציות ושומר את זה במשתנה חדש
    mutated_list=[]
    for seq in organizim_list:
        mutated_seq=seq
        for i in range (70):
            mutated_seq=Mutate_protein(mutated_seq)
        mutated_list.append(mutated_seq)
    return mutated_list 

#------------------------------------------------
def conserved_list(organizm_list):#שימוש של ה פונקציה של conserved
    place_number=[]
    is_conserved_all_organizams=[]
    counter=0
    seq_length=len(organizm_list[0])
    for h in range(seq_length):
        for i in range(8):
            x=organizm_list[i][h]
            place_number.append(x)
        is_conserved_all_organizams.append(is_conserved(place_number))
    return is_conserved_all_organizams
#------------------------------------------------
#תוכנית ראשית#
# פתיחת הקבצים
GAPDH_file = open('data/GAPDH_MSA.fasta', 'r')
RBP1_file = open('data/RBP1_MSA.fasta', 'r')

#הגדרת שתי רשימות המכילות את רצפי החלבונים שבקצים, כל רשימה מייצגת קובץ אחר.
GAPDH_list = []
RBP1_list = []

#ליסט של שתי הקבצים במקום קובץ
GAPDH_list = file_to_list(GAPDH_file)
RBP1_list = file_to_list(RBP1_file)

#העברת כל אחד מהרצפים בקבצים 70 מוטציות
mutated_GAPDH_list =mutate_ustage(GAPDH_list)
mutated_RBP1_list=mutate_ustage(RBP1_list)


"""
def ustage_conserved_list(organizam_list):
    for seq in organizam_list:
        list_organizem_concerved=[]
        str_seq=seq_to_string(seq)
        yes_or_not_conserved=conserved_list(str_seq)
        list_organizem_concerved.append(yes_or_not_conserved)
    return list_organizem_concerved

print(ustage_conserved_list(mutated_GAPDH_list))
"""



GAPDH_list_conserved=conserved_list(GAPDH_list)
mutated_GAPDH_list_conserved=conserved_list(mutated_GAPDH_list)
RBP1_list_conserved=conserved_list(RBP1_list)
mutated_RBP1_list_conserved=conserved_list(mutated_RBP1_list)


#שימוש בפונקציה compare שמשווה בין שתי הגרסאות של הרצפים
comapred_GAPDH_list=compare(GAPDH_list_conserved, mutated_GAPDH_list_conserved)
comapred_RBP1_list=compare(RBP1_list_conserved, mutated_RBP1_list_conserved)



x = [0,0,0,1,0,1,1,1,0,1,1,1,1,1,1]
max, sta, sto = max_seq(1, x)
print(max)
print(sta)
print(sto)
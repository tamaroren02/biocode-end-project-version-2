import random
#פונקציות#
#------------------------------------------------
def seq_to_string(seq):#הופך את הליסט לסטרינג
    seq_string = "".join(seq)
    return seq_string
#------------------------------------------------

def is_conserved(amino_list):
    #מיפוי תכונות של חומצות אמינו
    amino_groups = [['R','H','K'],['D','E'],['S', 'T', 'N', 'Q'],
                    ['P', 'C', 'G'],['A', 'V', 'I','L', 'M', 'F','Y', 'W']]
    
    counter_list = [0, 0, 0, 0, 0]
    
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
    #מחזירה 0/1 (שמור או לא שמור)
    for i in range(len(amino_list)):
        curr_acid = amino_list[i]

        for h in range(len(amino_groups)):
            for j in range(len(amino_groups[h])):
                if curr_acid == amino_groups[h][j]:
                    counter_list[h] += 1

    for i in range(len(counter_list)):
        if counter_list[i] >= 4:
            return 1
    return 0   

#------------------------------------------------

#משווה בין שתי הרצפים (הישן והחדש אחרי מוטציות) ומוצאת איפה שמור ואיפה לא
#לא נכון עד הסוף צריך לסיים אותה!!!
def compre(original_list, mutation_list):
    is_conserved_list = []
    if original_list % 10 != 0:
        left_region = original_list - (original_list % 10)

        for i in range(0, (original_list - left_region), 10):
            for h in range(10):
                original_sum = 0
                original_sum = original_sum + original_list[h]

                mutation_sum = 0
                mutation_sum = mutation_sum + mutation_list[h]

            if max(original_sum, mutation_sum) - min(original_sum, mutation_sum) <= 5:
                is_conserved_list.append("saved")
            else:
                is_conserved_list.append("not saved")

    if original_list % 10 != 0:
        for i in range((original_list - left_region), original_list):
            for h in range(10):
                original_sum = 0
                original_sum = original_sum + original_list[h]

                mutation_sum = 0
                mutation_sum = mutation_sum + mutation_list[h]

        if max(original_sum, mutation_sum) - min(original_sum, mutation_sum) <= 5:
            is_conserved_list.append("saved")
        else:   
            is_conserved_list.append("not saved")

    return is_conserved_list    
#------------------------------------------------
def Mutate_protein(seq):#מקבלת רצף חומצות אמינו של חיה ומחזירה את אותו רצף רק עם מוטצות נוקדתיות של החלפה
  '''
  '''
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
def conserved_list(seq_str):
    one_zero_list=[]
    seq_length=len(seq_str)
    for i in range (seq_length):
        x=is_conserved(seq_str)
        one_zero_list.append(x)
    return one_zero_list

#תוכנית ראשית#
# פתיחת הקבצים
GAPDH_file = open('data/GAPDH_MSA.fasta', 'r')
RBP1_file = open('data/RBP1_MSA.fasta', 'r')

#הגדרת שתי רשימות המכילות את רצפי החלבונים שבקצים, כל רשימה מייצגת קובץ אחר.
GAPDH_list = []
RBP1_list = []

GAPDH_list = file_to_list(GAPDH_file)
RBP1_list = file_to_list(RBP1_file)


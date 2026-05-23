import random
#פונקציות#

def is_conserved(amino_acid_list):
    '''
    הפונקציה בודקת האם מיקום בחלבון שמור/ לא שמור/ לא זה ולא זה.
    מקבלת: amino_acid_list
    מחזירה: 1/ 0 / '-'.
    '''
    cnt = 0
    for i in range(len(amino_acid_list)):
        if amino_acid_list[i] == "-":
            cnt += 1

    if cnt >= 3:
        return '-'

    if len(set(amino_acid_list)) <= 2:
        return 1
    else:
        return 0
        
#------------------------------------------------
def max_seq (zero_or_one, list):
    '''
    הפונקציה מקבלת רשימה ומספר ומחזירה את אורכו הרצף הכי ארוך ברשימה, את המיקום ברשימה בו הוא מתחיל ומסתיים.
    מקבלת: zero_or_one, list
    מחזירה: max_len, start, stop
    '''
    max_len = 0
    start = -1
    stop = -1

    i = 0
    while i < len(list):
        if list[i] == zero_or_one:
            temp_start = i
            cnt  = 0
            
            while i < len(list) and list[i] == zero_or_one:
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
def amino_acid_groups(original_seq, mutated_seq):
    '''
    הפונקציה מקבלת רצף בצורתו המקורית ואחרי מוטציות, היא בודקת כמה שינויים יש בין שני הרצפים לפי המיפוי של החומצות אמינו לפי התכונות שלהן.
    '''
    #מיפוי תכונות של חומצות אמינו
    groups = [['R','H','K'],['D','E'],['S', 'T', 'N', 'Q'],
                    ['P', 'C', 'G'],['A', 'V', 'I','L', 'M', 'F','Y', 'W']]
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

    for i in range(len(original_seq)):
        if original_seq[i] != mutated_seq[i]:
            for j in range(len(groups)):
                if mutated_seq[i] in groups[j] and original_seq[i] in groups[j]:
                    cnt += 1
    return cnt

#------------------------------------------------

def Mutate_protein(seq, num_mutation):
  '''
  הפונקציה מקבלת רצף חומצות אמינו ומספר המסמל את כמות המוטציות שהרצף שהתקבל יעבור. הפונקציה מבצעת רק מוטציות נקודתיות של החלפה.
  מקבלת: seq, num_mutation
  מחזירה: seq
  '''
  amino_acids_list = ['F', 'S', 'Y', 'C', 'L', '*', 'W',
                        'P', 'H', 'R', 'Q',
                        'I', 'T', 'N', 'K', 'M',
                        'V', 'A', 'D', 'G', 'E']

  for i in range(num_mutation):
    
    
    rand_acid = random.choice(amino_acids_list)
    
    rand_num = random.randrange(0,len(seq))
    
    while seq[rand_num] == rand_acid:
        rand_acid = random.choice(amino_acids_list)
    
    seq = seq[0:rand_num]+ rand_acid + seq[(rand_num+1):]
    
  return seq

#------------------------------------------------
def compare(zero_and_one_list, animle_list):
    '''
    הפונקציה מעבירה את האזור השמור והלא שמור הארוך ביותר בחלבון מסוים מוטציות ומחזירה באחוזים כמה כל אחד מהם השתנה.
    מקבלת: zero_and_one_list, animle_list
    מחזירה: conserved_percentage, non_conserved_percentage
    '''
    len_conserved, start1, stop1 = max_seq(1, zero_and_one_list)
    conserved_seq_original = animle_list[start1: stop1 + 1]


    num_mutation1 = int(len_conserved * 0.2)
    conserved_mutated_seq = Mutate_protein(conserved_seq_original, num_mutation1)
    
    conserved_cnt = amino_acid_groups(conserved_seq_original, conserved_mutated_seq)

    if len_conserved != 0 :
        # חישוב בכמה אחוזים הרצף ישתנה
        conserved_percentage = 100 * ( conserved_cnt / len_conserved)
    
    else:
        conserved_percentage = 0


    len_non_conserved, start0, stop0 = max_seq(0, zero_and_one_list)
    non_conserved_seq_original = animle_list[start0: stop0 + 1]

    num_mutation0 = int(len_non_conserved * 0.2)
    non_conserved_seq_mutated = Mutate_protein(non_conserved_seq_original, num_mutation0)

    non_conserved_cnt = amino_acid_groups(non_conserved_seq_original, non_conserved_seq_mutated)

    if len_non_conserved != 0:
        # חישוב בכמה אחוזים הרצף ישתנה
        non_conserved_percentage = 100 * ( non_conserved_cnt / len_non_conserved)
        
    else:
        non_conserved_percentage = 0


    return conserved_percentage, non_conserved_percentage

#------------------------------------------------
def file_to_list(file):
    '''
    הפונקציה מקבלת קובץ שבו רצפים של חומצות אמינו, מכניסה את רצפים אלו לרשימה ומחזירה את הרשימה.
    מקבלת: file
    מחזירה: seq_list
    '''
    curr_seq = ""
    seq_list = []
    
    for line in file:
        line = line.rstrip('\r\n')

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
def position(protein_list):
    '''
    הפונקציה יוצרת רשימה של חומצות אמינו ואו '-' אשר נמצאים באותו מיקום בכל הרצפים בקובץ מסוים.
    מקבלת: protein_list
    מחזירה: zero_one_list
    '''
    zero_one_list = []
    for i in range(len(protein_list[0])):
        position_list = []
        for h in range(len(protein_list)):
            position_list.append((protein_list[h])[i])

        zero_one_list.append(is_conserved(position_list))

    return zero_one_list

#------------------------------------------------
#תוכנית ראשית#

# פתיחת הקבצים
GAPDH_file = open('data/GAPDH_MSA.fasta', 'r')
RBP1_file = open('data/RBP1_MSA.fasta', 'r')

#הגדרת שתי רשימות ריקות שיכילו את רצפי החלבונים שבקבצים, כל רשימה מייצגת קובץ אחר.
GAPDH_list = []
RBP1_list = []

# הגדרת משתנים
GAPDH_conserved = 0
GAPDH_non_conserved = 0
RBP1_conserved = 0
RBP1_non_conserved = 0


GAPDH_list = file_to_list(GAPDH_file)
zero_one_GAPDH_list = position(GAPDH_list)
#print (zero_one_GAPDH_list)

GAPDH_conserved, GAPDH_non_conserved = compare(zero_one_GAPDH_list, GAPDH_list[0])

print("The protein: GAPDH")
print(f"{GAPDH_conserved:.2f}% percent of the longest conserved region in the protein changed")
print(f"{GAPDH_non_conserved:.2f}% percent of the longest non-conserved region in the protein changed\n")


RBP1_list = file_to_list(RBP1_file)
zero_one_RBP1_list = position(RBP1_list)
#print (zero_one_RBP1_list)

RBP1_conserved, RBP1_non_conserved = compare(zero_one_RBP1_list, RBP1_list[0])

print("The protein: RBP1")
print(f"{RBP1_conserved:.2f}% percent of the longest conserved region in the protein changed")
print(f"{RBP1_non_conserved:.2f}% percent of the longest non-conserved region in the protein changed")









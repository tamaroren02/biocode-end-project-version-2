import random
#פונקציות#
def is_conserved(amino_acid_list):
    '''
    הפונקציה בודקת האם מיקום בחלבון שמור/ לא שמור/ לא זה ולא זה.
    מקבלת: amino_acid_list
    מחזירה: 1/ 0 / '-'.
    '''
    
    cnt_dict ={}
    for i in range(len(amino_acid_list)):
        if amino_acid_list[i] in cnt_dict:
            cnt_dict[amino_acid_list[i]] += 1
        else:
            cnt_dict[amino_acid_list[i]] = 1

    if cnt_dict.get("-", 0) >= 3:
        return '-'

    if len(cnt_dict) <= 2:
        max_value = max(cnt_dict.values())
        num_amino_acides= 8 - cnt_dict.get('-', 0)
        max_amino_acid = max_value / num_amino_acides

        if max_amino_acid >= 0.75:
            return 1
        else:
            return 0
    
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
    cnt_svd_grp = 0
    nt_cnt_svd_grp=0

    for i in range(len(original_seq)):
        if original_seq[i] != mutated_seq[i]:
            same_group = False

            for group in groups:
                if original_seq[i] in group and mutated_seq[i] in group:
                    same_group = True
                
            if same_group:
                cnt_svd_grp += 1
            else:
                nt_cnt_svd_grp += 1

    return cnt_svd_grp, nt_cnt_svd_grp

#------------------------------------------------

def Mutate_protein(seq, num_mutation):
  '''
  הפונקציה מקבלת רצף חומצות אמינו ומספר המסמל את כמות המוטציות שהרצף שהתקבל יעבור. הפונקציה מבצעת רק מוטציות נקודתיות של החלפה.
  מקבלת: seq, num_mutation
  מחזירה: seq
  '''
  amino_acids_list = ['F', 'S', 'Y', 'C', 'L', 'W',
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
    #print("zero_and_one_list",zero_and_one_list)
    len_conserved, start1, stop1 = max_seq(1, zero_and_one_list)
    conserved_seq_original = animle_list[start1: stop1 + 1]

    num_mutation1 = int(len_conserved * 0.2)
    conserved_mutated_seq = Mutate_protein(conserved_seq_original, num_mutation1)
    
    c_saved ,c_not_saved = amino_acid_groups(conserved_seq_original, conserved_mutated_seq)
    print("len_conserved",len_conserved)
    if len_conserved != 0 :
        # חישוב בכמה אחוזים הרצף ישתנה
        conserved_percentage = 100 * ( c_not_saved / int(len_conserved * 0.2))
    
    else:
        conserved_percentage = 0


    len_non_conserved, start0, stop0 = max_seq(0, zero_and_one_list)
    non_conserved_seq_original = animle_list[start0: stop0 + 1]

    num_mutation0 = int(len_non_conserved * 0.2)
    print("len_non_conserved",len_non_conserved)

    non_conserved_seq_mutated = Mutate_protein(non_conserved_seq_original, num_mutation0)

    nc_saved,nc_not_saved  = amino_acid_groups(non_conserved_seq_original, non_conserved_seq_mutated)

    if len_non_conserved != 0:
        # חישוב בכמה אחוזים הרצף ישתנה
        non_conserved_percentage = 100 * ( nc_not_saved / int(len_non_conserved * 0.2))
        
    else:
        non_conserved_percentage = 0

    #להלן פירוט על חלק משמות המשתנים:
    #nc_not_saved = non_conserved_cnt_not_saved  ,c_not_saved = conserved_cnt_not_saved
    #nc_saved = non_conserved_cnt_saved  , c_saved = conserved_cnt_saved
    return conserved_percentage, non_conserved_percentage, len_conserved, nc_not_saved, c_not_saved, nc_saved, c_saved

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

def seq_lengths(seq_list,organism_file):
    #מקבלת את הליסט של האורגניזמים, הקובץ,ואת הקובץ שרושמים אליו את התשובות הסופיות
    #הפונקציה מוצאת את אורך כל אחד מהאורגניזמים ומכניסה לקובץ התוצאות הסופיות את השם של האורגניזם ומה האורך שלו
    organism_file.seek(0)
    seq_lengths_list = []
    organism_names_list = []

    for seq in seq_list:
        seq_lengths_list.append(len(seq))
    
    for line in organism_file:
        name=line
        if line.startswith(">"):
            organism_name = line[1:line.find("_")]
            organism_names_list.append(organism_name)
    return organism_names_list,seq_lengths_list
#------------------------------------------------
    
def columns_count(zero_one_list,file_for_results,protein_name):
    #מקבל את הרשימה של 1(שמור), 0(לא שמור) ,-(פסול) 
    #את הקובץ לתוצאות ואת שם החלבון
    #מה היא עושה:מדפיסה אל תוך קובץ התוצאות אחוז של כמה כל אחד מהם הופיע ברשימה
    columns_dict={}
    zero_one_list_amount=len(zero_one_list)
    for i in range(zero_one_list_amount):
        if zero_one_list[i] in columns_dict:
            columns_dict[zero_one_list[i]] += 1
        else:
            columns_dict[zero_one_list[i]] = 1
    
    for y in columns_dict:
        columns_precentage= (columns_dict[y]/zero_one_list_amount)*100
        if y == 1:
            file_for_results.write(f"precent of concerved positions in {protein_name} protein= {columns_precentage:.2f}%\n")
        elif y == 0:
            file_for_results.write(f"precent of non concerved positions in {protein_name} protein= {columns_precentage:.2f}%\n")
        elif y == "-":
            file_for_results.write(f"precent of disqualified positions in {protein_name} protein= {columns_precentage:.2f}%\n")

#------------------------------------------------
def writing_to_file_per_org(seq_list,organism_file,file_for_results):
    #כותב לתוך הקובץ את השם של כל חיה : מה האורך שלה וכמה עמודות פסלנו בה
    organism_names_list,seq_lengths_list=seq_lengths(seq_list,organism_file)
    for i in range(8):
        file_for_results.write(f"{organism_names_list[i]},")
    file_for_results.write(f"\ntheir mutual length={seq_lengths_list[i]}\n")
#------------------------------------------------

def graph_changes (zero_and_one_list, organism_list,graph_file):
    #קבצים שנכניס לאקסל ונעשה גרפים#
    #אחוזים של שמור ולא שמור#
    num_reps=100
    total_nc_not_saved=0
    total_c_not_saved=0
    total_nc_saved=0
    total_c_saved=0
    
    graph_file.write(f"conserved,non_conserved\n")

    for i in range(num_reps):  
        conserved, non_conserved, conserved_amount, nc_not_saved, c_not_saved, nc_saved, c_saved= compare(zero_and_one_list,organism_list[0])
        
        graph_file.write(f"{conserved:.2f}%,{non_conserved:.2f}%\n")

        total_nc_not_saved += nc_not_saved
        total_c_not_saved += c_not_saved
        total_nc_saved += nc_saved
        total_c_saved += c_saved

    return total_nc_not_saved, total_c_not_saved, total_nc_saved, total_c_saved

#-----------------------------------------------
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
n_c_s=0
c_s=0
nc_not_saved=0
c_not_saved=0


GAPDH_list = file_to_list(GAPDH_file)
zero_one_GAPDH_list = position(GAPDH_list)

GAPDH_conserved, GAPDH_non_conserved, GAPDH_len_conserved, nc_not_saved, c_not_saved, nc_saved, c_saved = compare(zero_one_GAPDH_list, GAPDH_list[0])

RBP1_list = file_to_list(RBP1_file)
zero_one_RBP1_list = position(RBP1_list)
print(zero_one_RBP1_list)
print(len(zero_one_RBP1_list))


RBP1_conserved, RBP1_non_conserved, RBP1_len_conserved, nc_not_saved, c_not_saved, nc_saved, c_saved = compare(zero_one_RBP1_list, RBP1_list[0])


###🥳🥳תוצאות סופיות🥳🥳###
results_file=open('results/results_file', 'w')

results_file.write("Results:\n")
results_file.write("\n")

#GAPDH:
#השם של כל חיה ואז האורך של הרצף שלה ומה האחוז עמודות שפסלנו בה
results_file.write("GAPDH:\n")
GAPDH_data = writing_to_file_per_org(GAPDH_list,GAPDH_file,results_file)
results_file.write("\n")

#אחוזי השינוי בחלק הכי ארוך השמור והלא שמור#
results_file.write(f"{GAPDH_conserved:.2f}% percent of the longest conserved region in the protein changed\n")
results_file.write(f"{GAPDH_non_conserved:.2f}% percent of the longest non-conserved region in the protein changed\n")
results_file.write("\n")

#הרצף השמור הארוך ביותר#
results_file.write(f" longest conserved sequence length={GAPDH_len_conserved}\n")
results_file.write("\n")

#אחוזים של כמה פעמים מופיע 1(שמור), 0(לא שמור) ו- (פסול) ומעביר לאחוזים
GAPDH="GAPDH"
GAPDH_columns_count=columns_count(zero_one_GAPDH_list,results_file,GAPDH)

###RBP1###
results_file.write("\n")
results_file.write("\n")
results_file.write("RBP1:\n")

#השם של כל חיה ואז האורך של הרצף שלה ומה האחוז עמודות שפסלנו בה
RBP1_data=writing_to_file_per_org(RBP1_list,RBP1_file,results_file)
results_file.write("\n")

#אחוזי השינוי בחלק הכי ארוך השמור והלא שמור#
results_file.write(f"{RBP1_conserved:.2f}% percent of the longest conserved region in the protein changed\n")
results_file.write(f"{RBP1_non_conserved:.2f}% percent of the longest non-conserved region in the protein changed\n")
results_file.write("\n")

#הרצף השמור הארוך ביותר#
results_file.write(f" longest conserved sequence length={RBP1_len_conserved}\n")
results_file.write("\n")

#אחוזים של כמה פעמים מופיע 1(שמור), 0(לא שמור) ו- (פסול) ומעביר לאחוזים
RBP1="RBP1"
RBP1_columns_count=columns_count(zero_one_RBP1_list,results_file,RBP1)
results_file.write("\n")

#------------------------------------------------------------------------------

#קבצים שנכניס לאקסל ונעשה גרפים#
#אחוזים של שמור ולא שמור#
GAPDH_graph_file = open('results/GAPDH_graph', 'w')
RBP1_graph_file = open('results/RBP1_graph', 'w')


GAPDH_nc_not_saved, GAPDH_c_not_saved, GAPDH_nc_saved, GAPDH_c_saved = graph_changes (zero_one_GAPDH_list, GAPDH_list, GAPDH_graph_file)

RBP1_nc_not_saved, RBP1_c_not_saved, RBP1_nc_saved, RBP1_c_saved = graph_changes (zero_one_RBP1_list, RBP1_list, RBP1_graph_file)

GAPDH_file.close()
RBP1_file.close()
results_file.close()
GAPDH_graph_file.close()
RBP1_graph_file.close()

x=1.6
print(int(x))
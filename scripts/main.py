def is_conserved(amino_list):
    #מיפוי תכונות של חומצות אמינו
    amino_groups = [['R','H','K'],['D','E'],['S', 'T', 'N', 'Q'],['P', 'C', 'G'],['A', 'V', 'I','L', 'M', 'F','Y', 'W']]
    
    counter_list = [0, 0, 0, 0, 0, 0]
    

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



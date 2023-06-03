class AminoAcidProfile:
    def __init__(self, sequences):
        self.sequences = sequences
        self.length = len(sequences[0])
        self.num_sequences = len(sequences)
        
    def get_amino_acid(self, position):
        amino_acids = [seq[position] for seq in self.sequences]
        return amino_acids
    
    def get_profile(self, position):
        profile = {}
        amino_acids = self.get_amino_acid(position)
        total = len(amino_acids)
        
        for aa in amino_acids:
            if aa in profile:
                profile[aa] += 1
            else:
                profile[aa] = 1
        
        for aa in profile:
            profile[aa] /= total
        
        return profile

blosum62 = {
    "A": {"A": 4, "R": -1, "N": -2, "D": -2, "C": 0, "Q": -1, "E": -1, "G": 0, "H": -2, "I": -1, "L": -1, "K": -1, "M": -1, "F": -2, "P": -1, "S": 1, "T": 0, "W": -3, "Y": -2, "V": 0},
    "R": {"A": -1, "R": 5, "N": 0, "D": -2, "C": -3, "Q": 1, "E": 0, "G": -2, "H": 0, "I": -3, "L": -2, "K": 2, "M": -1, "F": -3, "P": -2, "S": -1, "T": -1, "W": -3, "Y": -2, "V": -3},
    "N": {"A": -2, "R": 0, "N": 6, "D": 1, "C": -3, "Q": 0, "E": 0, "G": 0, "H": 1, "I": -3, "L": -3, "K": 0, "M": -2, "F": -3, "P": -2, "S": 1, "T": 0, "W": -4, "Y": -2, "V": -3},
    "D": {"A": -2, "R": -2, "N": 1, "D": 6, "C": -3, "Q": 0, "E": 2, "G": -1, "H": -1, "I": -3, "L": -4, "K": -1, "M": -3, "F": -3, "P": -1, "S": 0, "T": -1, "W": -4, "Y": -3, "V": -3},
    "C": {"A": 0, "R": -3, "N": -3, "D": -3, "C": 9, "Q": -3, "E": -4, "G": -3, "H": -3, "I": -1, "L": -1, "K": -3, "M": -1, "F": -2, "P": -3, "S": -1, "T": -1, "W": -2, "Y": -2, "V": -1},
    "Q": {"A": -1, "R": 1, "N": 0, "D": 0, "C": -3, "Q": 5, "E": 2, "G": -2, "H": 0, "I": -3, "L": -2, "K": 1, "M": 0, "F": -3, "P": -1, "S": 0, "T": -1, "W": -2, "Y": -1, "V": -2},
    "E": {"A": -1, "R": 0, "N": 0, "D": 2, "C": -4, "Q": 2, "E": 5, "G": -2, "H": 0, "I": -3, "L": -3, "K": 1, "M": -2, "F": -3, "P": -1, "S": 0, "T": -1, "W": -3, "Y": -2, "V": -2},
    "G": {"A": 0, "R": -2, "N": 0, "D": -1, "C": -3, "Q": -2, "E": -2, "G": 6, "H": -2, "I": -4, "L": -4, "K": -2, "M": -3, "F": -3, "P": -2, "S": 0, "T": -2, "W": -2, "Y": -3, "V": -3},
    "H": {"A": -2, "R": 0, "N": 1, "D": -1, "C": -3, "Q": 0, "E": 0, "G": -2, "H": 8, "I": -3, "L": -3, "K": -1, "M": -2, "F": -1, "P": -2, "S": -1, "T": -2, "W": -2, "Y": 2, "V": -3},
    "I": {"A": -1, "R": -3, "N": -3, "D": -3, "C": -1, "Q": -3, "E": -3, "G": -4, "H": -3, "I": 4, "L": 2, "K": -3, "M": 1, "F": 0, "P": -3, "S": -2, "T": -1, "W": -3, "Y": -1, "V": 3},
    "L": {"A": -1, "R": -2, "N": -3, "D": -4, "C": -1, "Q": -2, "E": -3, "G": -4, "H": -3, "I": 2, "L": 4, "K": -2, "M": 2, "F": 0, "P": -3, "S": -2, "T": -1, "W": -2, "Y": -1, "V": 1},
    "K": {"A": -1, "R": 2, "N": 0, "D": -1, "C": -3, "Q": 1, "E": 1, "G": -2, "H": -1, "I": -3, "L": -2, "K": 5, "M": -1, "F": -3, "P": -1, "S": 0, "T": -1, "W": -3, "Y": -2, "V": -2},
    "M": {"A": -1, "R": -1, "N": -2, "D": -3, "C": -1, "Q": 0, "E": -2, "G": -3, "H": -2, "I": 1, "L": 2, "K": -1, "M": 5, "F": 0, "P": -2, "S": -1, "T": -1, "W": -1, "Y": -1, "V": 1},
    "F": {"A": -2, "R": -3, "N": -3, "D": -3, "C": -2, "Q": -3, "E": -3, "G": -3, "H": -1, "I": 0, "L": 0, "K": -3, "M": 0, "F": 6, "P": -4, "S": -2, "T": -2, "W": 1, "Y": 3, "V": -1},
    "P": {"A": -1, "R": -2, "N": -2, "D": -1, "C": -3, "Q": -1, "E": -1, "G": -2, "H": -2, "I": -3, "L": -3, "K": -1, "M": -2, "F": -4, "P": 7, "S": -1, "T": -1, "W": -4, "Y": -3, "V": -2},
    "S": {"A": 1, "R": -1, "N": 1, "D": 0, "C": -1, "Q": 0, "E": 0, "G": 0, "H": -1, "I": -2, "L": -2, "K": 0, "M": -1, "F": -2, "P": -1, "S": 4, "T": 1, "W": -3, "Y": -2, "V": -2},
    "T": {"A": 0, "R": -1, "N": 0, "D": -1, "C": -1, "Q": -1, "E": -1, "G": -2, "H": -2, "I": -1, "L": -1, "K": -1, "M": -1, "F": -2, "P": -1, "S": 1, "T": 5, "W": -2, "Y": -2, "V": 0},
    "W": {"A": -3, "R": -3, "N": -4, "D": -4, "C": -2, "Q": -2, "E": -3, "G": -2, "H": -2, "I": -3, "L": -2, "K": -3, "M": -1, "F": 1, "P": -4, "S": -3, "T": -2, "W": 11, "Y": 2, "V": -3},
    "Y": {"A": -2, "R": -2, "N": -2, "D": -3, "C": -2, "Q": -1, "E": -2, "G": -3, "H": 2, "I": -1, "L": -1, "K": -2, "M": -1, "F": 3, "P": -3, "S": -2, "T": -2, "W": 2, "Y": 7, "V": -1},
    "V": {"A": 0, "R": -3, "N": -3, "D": -3, "C": -1, "Q": -2, "E": -2, "G": -3, "H": -3, "I": 3, "L": 1, "K": -2, "M": 1, "F": -1, "P": -2, "S": -2, "T": 0, "W": -3, "Y": -1, "V": 4}
}

def blosum62_amino_score(amino1, amino2, gap_penalty):
    if amino1 == "-" or amino2 == "-":
        return gap_penalty
    else:
        return blosum62[amino1][amino2]

def align_sequences_with_profile(sequence, profile, gap_penalty):
    #プロファイル一致率計算での浮動小数点誤差対策
    epsilon = 1e-8

    # 配列の長さを取得
    m = len(sequence)
    n = profile.length
    
    # スコア行列を初期化
    score_matrix = [[0] * (n + 1) for _ in range(m + 1)]
    
    # 初期化ステップ：境界条件の設定
    for i in range(m + 1):
        score_matrix[i][0] = 0
        
    for j in range(n + 1):
        score_matrix[0][j] = 0
    
    # 動的計画法の実行
    for i in range(1, m + 1):
        for j in range(1, n + 1):

            # プロファイルとシーケンスのアミノ酸一致率を計算
            #Blosum62ように改良
            score = [blosum62_amino_score(sequence[i - 1], amino_acid, gap_penalty) for amino_acid in profile.get_amino_acid(j - 1)]

            # 斜め方向のスコア
            match = score_matrix[i - 1][j - 1] + sum(score)
            
            # 上方向のスコア
            if j == n:
                delete = score_matrix[i - 1][j]             
            else:
                delete = score_matrix[i - 1][j] + gap_penalty * profile.num_sequences
            
            # 左方向のスコア
            if i == m:
                insert = score_matrix[i][j - 1]
            else:
                insert = score_matrix[i][j - 1] + gap_penalty * profile.num_sequences
            
            # 最大スコアを選択
            score_matrix[i][j] = max(match, delete, insert)

            #デバッグ用
            #print('\n'.join([' '.join([str(element) for element in row]) for row in score_matrix]))    
    # アラインメントの復元
    aligned_sequence = ""
    aligned_profile = [""] * profile.num_sequences
    i = m
    j = n
    
    while i > 0 and j > 0:        
        if j == n:
            mismatch_penaltyI = 0
        else:
            mismatch_penaltyI = gap_penalty * profile.num_sequences
        if i == m:
            mismatch_penaltyJ = 0
        else:
            mismatch_penaltyJ = gap_penalty * profile.num_sequences        

        #上方向に一致
        if abs(score_matrix[i][j] - (score_matrix[i - 1][j] + mismatch_penaltyI)) < epsilon:
            aligned_sequence = sequence[i - 1] + aligned_sequence
            for k in range(profile.num_sequences):
                aligned_profile[k] = "-" + aligned_profile[k]
            i -= 1
        #横方向に一致
        elif abs(score_matrix[i][j] - (score_matrix[i][j - 1] + mismatch_penaltyJ)) < epsilon:
            aligned_sequence = "-" + aligned_sequence
            for k in range(profile.num_sequences):
                aligned_profile[k] = profile.sequences[k][j - 1] + aligned_profile[k]
            j -= 1
        #斜め方向に一致
        else:
            aligned_sequence = sequence[i - 1] + aligned_sequence
            for k in range(profile.num_sequences):
                aligned_profile[k] = profile.sequences[k][j - 1] + aligned_profile[k]
            i -= 1
            j -= 1
#        print(aligned_sequence)
#        for p in aligned_profile:
#            print(p)

    
    while i > 0:
        aligned_sequence = sequence[i - 1] + aligned_sequence
        for k in range(profile.num_sequences):
            aligned_profile[k] = "-" + aligned_profile[k]
        i -= 1

    
    while j > 0:
        aligned_sequence = "-" + aligned_sequence
        for k in range(profile.num_sequences):
            aligned_profile[k] = profile.sequences[k][j - 1] + aligned_profile[k]
        j -= 1
    
    alignment_score = score_matrix[m][n]
    return aligned_sequence, aligned_profile, alignment_score


# アミノ酸配列
sequence = "ACRCT"

# プロファイル
profile_sequences = ["A--TACGT",
           "ACGTAC-G",
           "ACGTA-GT"]

profile = AminoAcidProfile(profile_sequences)

gap_penalty = -8

# アラインメントの実行
aligned_sequence, aligned_profile, alignment_score = align_sequences_with_profile(sequence, profile, gap_penalty)

# 結果の表示
print(aligned_sequence)
for p in aligned_profile:
    print(p)
print(alignment_score)

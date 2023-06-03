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





def align_sequences_with_profile(sequence, profile, match_score, mismatch_penalty):
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
            """
            #Blosum62ように改良
            score = 0
            for amino_acid in profile.get_amino_acid(j - 1):
                score += blosum62_amino_score(sequence[i - 1], amino_acid, gap_penalty)
            """
            profile_profile = profile.get_profile(j - 1)
            try:
                percentage_of_match = profile_profile[sequence[i - 1]]
            except KeyError:
                percentage_of_match = 0

            # 斜め方向のスコア
            point = match_score * percentage_of_match + mismatch_penalty * (1 - percentage_of_match)
            match = score_matrix[i - 1][j - 1] + point
            
            # 上方向のスコア
            delete = score_matrix[i - 1][j] + mismatch_penalty
            
            # 左方向のスコア
            insert = score_matrix[i][j - 1] + mismatch_penalty
            
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
        profile_profile = profile.get_profile(j - 1)
        try:
            percentage_of_match = profile_profile[sequence[i - 1]]
        except KeyError:
            percentage_of_match = 0
        point = match_score * percentage_of_match + mismatch_penalty * (1 - percentage_of_match)
        #if score_matrix[i][j] == score_matrix[i - 1][j - 1] + point:
        #斜め方向に一致
        if sequence[i - 1] in profile.get_amino_acid(j - 1) and (abs(score_matrix[i][j] - (score_matrix[i - 1][j - 1] + point)) < epsilon):
            aligned_sequence = sequence[i - 1] + aligned_sequence
            for k in range(profile.num_sequences):
                aligned_profile[k] = profile.sequences[k][j - 1] + aligned_profile[k]
            i -= 1
            j -= 1
        #上方向に一致
        #elif score_matrix[i][j] == score_matrix[i - 1][j] + mismatch_penalty:
        elif abs(score_matrix[i][j] - (score_matrix[i - 1][j] + mismatch_penalty)) < epsilon:
            aligned_sequence = sequence[i - 1] + aligned_sequence
            for k in range(profile.num_sequences):
                aligned_profile[k] = "-" + aligned_profile[k]
            i -= 1
        #横方向に一致
        else:
            aligned_sequence = "-" + aligned_sequence
            for k in range(profile.num_sequences):
                aligned_profile[k] = profile.sequences[k][j - 1] + aligned_profile[k]
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
    
    return aligned_sequence, aligned_profile


# アミノ酸配列
sequence = "ACTACT"

# プロファイル
profile_sequences = ["A--TACGT",
           "ACGTAC-G",
           "ACGTA-GT"]

profile = AminoAcidProfile(profile_sequences)

# アラインメントの実行
aligned_sequence, aligned_profile = align_sequences_with_profile(sequence, profile, 2, -1)

# 結果の表示
print(aligned_sequence)
for p in aligned_profile:
    print(p)


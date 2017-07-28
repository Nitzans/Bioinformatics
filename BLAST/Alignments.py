#!/usr/bin/python


def fasta_parser(fasta_file):
    with open(fasta_file, mode='r') as fasta_file:
        fasta_dict = {}
        seq_name = ''
        for fasta_line in fasta_file:
            if fasta_line.startswith('>'):
                seq_name = fasta_line[1:-1]
                fasta_dict[seq_name] = ''
            else:
                fasta_dict[seq_name] += fasta_line.replace('\n', '')
        fasta_dict = sorted(fasta_dict.items())  # return list of tuples
    return fasta_dict


def matrix_parser(matrix_file):
    with open(matrix_file, mode='r') as score_file:
        score_file = score_file.readlines()
        matrix = score_file[:][7:]
        matrix = [f.split() for f in matrix]
        matrix[0].insert(0, ' ')  # align the column
        score_matrix = {}
        for x in range(1, len(matrix[0])):
            for y in range(1, len(matrix[0])):
                score_matrix[(matrix[0][x], matrix[y][0])] = matrix[y][x]  # set score to tuple of two characters
    return score_matrix


def local_alignment(seq1, seq2, scoring_matrix):
    """ Perform local alignment & backtracking """
    # Init matrices.
    s = [[0] * (len(seq2) + 1) for i in range(0, len(seq1) + 1)]  # fill the matrix s with zeroes
    backtrack = [[0] * (len(seq2) + 1) for i in range(0, len(seq1) + 1)]  # fill the reconstruction matrix in zeroes

    best = [[0] * 3] * 2
    # Fill out s & backtrack matrices.
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            scores = [
                s[i - 1][j] + int(scoring_matrix[(seq1[i - 1]), '*']),
                s[i][j - 1] + int(scoring_matrix[('*', seq2[j - 1])]),
                s[i - 1][j - 1] + int(scoring_matrix[(seq1[i - 1]), seq2[j - 1]]),
                0
            ]
            # Get & save max score.
            s[i][j] = max(scores)

            # Store location & value of 2 best scoring cell so far.
            if best[0][0] < s[i][j]:
                best[1] = best[0][:]  # move previous best score to be sub-optimal score
                best[0][0], best[0][1], best[0][2] = s[i][j], i, j  # and now update the new best score

            backtrack[i][j] = scores.index(s[i][j])

    # Get the position of the two highest scoring cell and its score.
    i, j = best[0][1], best[0][2]
    score = s[i][j]

    i2, j2 = best[1][1], best[1][2]
    score2 = s[i2][j2]

    # We only want seq1 and seq2 up to seq1[i] and seq2[j].
    seq1_aligned_best, seq2_aligned_best = seq1[:i], seq2[:j]

    # Backtrack local alignment starting at the highest scoring cell.
    while backtrack[i][j] != 3 and i * j != 0:  # backtrack[i][j] == 3 means 0 was the maximum score
        if backtrack[i][j] == 0:  # up
            i -= 1
            seq2_aligned_best = seq2_aligned_best[:j] + '_' + seq2_aligned_best[j:]
        elif backtrack[i][j] == 1:  # left
            j -= 1
            seq1_aligned_best = seq1_aligned_best[:i] + '_' + seq1_aligned_best[i:]
        elif backtrack[i][j] == 2:  # diagonal
            i -= 1
            j -= 1

    # Discard anything before seq1[i] and seq2[j]
    # NOTICE: i & j were updated during the backtrack.
    seq1_aligned_best = seq1_aligned_best[i:]
    seq2_aligned_best = seq2_aligned_best[j:]

    seq1_aligned_2nd, seq2_aligned_2nd = seq1[:i2], seq2[:j2]
    # Backtrack local alignment starting at the second best scoring cell.
    while backtrack[i2][j2] != 3 and i2 * j2 != 0:
        if backtrack[i2][j2] == 0:
            i2 -= 1
            seq2_aligned_2nd = seq2_aligned_2nd[:j2] + '_' + seq2_aligned_2nd[j2:]
        elif backtrack[i2][j2] == 1:
            j2 -= 1
            seq1_aligned_2nd = seq1_aligned_2nd[:i2] + '_' + seq1_aligned_2nd[i2:]
        elif backtrack[i2][j2] == 2:
            i2 -= 1
            j2 -= 1

    # Discard anything before seq1[i] and seq2[j]
    # NOTICE: i & j were updated during the backtrack.
    seq1_aligned_2nd = seq1_aligned_2nd[i2:]
    seq2_aligned_2nd = seq2_aligned_2nd[j2:]

    return score, seq1_aligned_best, seq2_aligned_best, score2, seq1_aligned_2nd, seq2_aligned_2nd


def best_locals():
    texts = fasta_parser('text[1].fasta')  # Create fastas list
    queries = fasta_parser('queries[1].fasta')  # Create fastas list
    matrix = matrix_parser('Score.matrix')  # Create score matrix

    for q, Qseq in queries:
        print "Running {}...".format(q)
        best = ("", 0)
        semi = ("", 0)
        for t, Tseq in texts:
            locals_dict = local_alignment(Qseq, Tseq, matrix)
            max_score = locals_dict[0]
            sub_score = locals_dict[3]
            if best[1] < max_score:
                semi = best
                best = (t, max_score, sub_score)
            elif semi[1] < max_score:
                semi = (t, max_score, sub_score)
        print ('Query {}\n\tBest match: {}, Score: {},{}\n\tSemi-best match: {}, Score: {},{}\n\n'
               .format(q, best[0], best[1], best[2], semi[0], semi[1], semi[2]))

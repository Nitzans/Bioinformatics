#!/usr/bin/python
import sys #for args

def fasta_parser(file):
    with open(file, mode='r') as fasta_file:
        fasta_dict = {}
        current_seq = ''
        for fasta_line in fasta_file:
            if fasta_line.startswith('>'):
                current_seq = fasta_line[1:-1]
                fasta_dict[current_seq] = ''
            else:
                fasta_dict[current_seq] += fasta_line.replace('\n', '')
        fasta_dict = sorted(fasta_dict.items()) #return list of tuples
    return fasta_dict

def matrix_parser(file):
    with open(file, mode='r') as score_file:
        score_file = score_file.readlines()
        matrix = score_file[:][7:]
        matrix = [f.split() for f in matrix]
        matrix[0].insert(0, ' ') #aligne the column
        score_matrix = {}
        for x in range(1, len(matrix[0])):
            for y in range(1, len(matrix[0])):
                score_matrix[(matrix[0][x], matrix[y][0])] = matrix[y][x] #create value of tuple and put in it his score
    return score_matrix

def global_alignment(seq1, seq2, scoring_matrix):
    """ Perform global alignment & backtracking """
    # Init matrices.
    s = [[0] * (len(seq2) + 1) for i in range(0, len(seq1) + 1)] #fill the matrix s with zeroes
    backtrack = [[0] * (len(seq2) + 1) for i in range(0, len(seq1) + 1)] #fill the reconstruction matrix in zeroes

    # Save location & value of best scoring cell.
    best = [0] * 3
    # Fill in the Score and Backtrack matrices.
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            scores = [
                s[i - 1][j] + int(scoring_matrix[(seq1[i - 1]), '*']),
                s[i][j - 1] + int(scoring_matrix[('*', seq2[j - 1])]),
                s[i - 1][j - 1] + int(scoring_matrix[(seq1[i - 1]), seq2[j - 1]])
            ]
            s[i][j] = max(scores)

            # Store location & value of best scoring cell so far.
            if i == len(seq1) and j == len(seq2):
                if best[0] < s[i][j]:
                    best[0], best[1], best[2] = s[i][j], i, j

            backtrack[i][j] = scores.index(s[i][j]) #get the direction 0, 1 or 2

    # Make copies of seq1 & seq2.
    seq1_aligned, seq2_aligned = seq1, seq2

    # Get the position of the highest scoring cell and its score.
    i, j = best[1], best[2]
    max_score = s[i][j]

    # Backtrack to the edge of the matrix starting at the highest scoring cell.
    while i * j != 0:
        if backtrack[i][j] == 0 and j != 0: #up
            i -= 1
            seq2_aligned = seq2_aligned[:j] + '_' + seq2_aligned[j:]
        elif backtrack[i][j] == 1 and i != 0: #left
            j -= 1
            seq1_aligned = seq1_aligned[:i] + '_' + seq1_aligned[i:]
        else:
            i -= 1
            j -= 1

    # insert indels till the begining to align the sequences
    for repeat in range(i):
        seq2_aligned = seq2_aligned[:0] + '_' + seq2_aligned[0:]
    for repeat in range(j):
        seq1_aligned = seq1_aligned[:0] + '_' + seq1_aligned[0:]

    return max_score, seq1_aligned, seq2_aligned

def local_alignment(seq1, seq2, scoring_matrix):
    """ Perform local alignment & backtracking """
    # Init matrices.
    s = [[0] * (len(seq2) + 1) for i in range(0, len(seq1) + 1)] #fill the matrix s with zeroes
    backtrack = [[0] * (len(seq2) + 1) for i in range(0, len(seq1) + 1)] #fill the reconstruction matrix in zeroes

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
                best[1] = best[0][:] #move previous best score to be sub-optimal score
                best[0][0], best[0][1], best[0][2] = s[i][j], i, j #and now update the new best score

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
        if backtrack[i][j] == 0: #up
            i -= 1
            seq2_aligned_best = seq2_aligned_best[:j] + '_' + seq2_aligned_best[j:]
        elif backtrack[i][j] == 1: #left
            j -= 1
            seq1_aligned_best = seq1_aligned_best[:i] + '_' + seq1_aligned_best[i:]
        elif backtrack[i][j] == 2: #diagonal
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

def free_ends_alignment(seq1, seq2, scoring_matrix):
    """ Perform free-ends alignment & backtracking """
    # Init matrices.
    s = [[0] * (len(seq2) + 1) for i in range(0, len(seq1) + 1)] #fill the matrix with zeroes
    backtrack = [[0] * (len(seq2) + 1) for i in range(0, len(seq1) + 1)] #fill the reconstruction matrix in zeroes

    # Save location & value of best scoring cell.
    best = [0] * 3
    # Fill in the Score and Backtrack matrices.
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            scores = [
                s[i - 1][j] + int(scoring_matrix[(seq1[i - 1]), '*']),
                s[i][j - 1] + int(scoring_matrix[('*', seq2[j - 1])]),
                s[i - 1][j - 1] + int(scoring_matrix[(seq1[i - 1]), seq2[j - 1]])
            ]
            s[i][j] = max(scores)

            # Store location & value of best scoring cell so far.
            if i == len(seq1) or j == len(seq2): #last row or last column
                if best[0] < s[i][j]:
                    best[0], best[1], best[2] = s[i][j], i, j

            backtrack[i][j] = scores.index(s[i][j])

    # Make copies of seq1 & seq2.
    seq1_aligned, seq2_aligned = seq1, seq2

    # Get the position of the highest scoring cell and its score.
    i, j = best[1], best[2]
    max_score = s[i][j]

    # Backtrack to the edge of the matrix starting at the highest scoring cell.
    while i * j != 0:
        if backtrack[i][j] == 0:
            i -= 1
            seq2_aligned = seq2_aligned[:j] + '_' + seq2_aligned[j:]
        elif backtrack[i][j] == 1:
            j -= 1
            seq1_aligned = seq1_aligned[:i] + '_' + seq1_aligned[i:]
        else:
            i -= 1
            j -= 1

    # insert indels till the begining to align the sequences
    for repeat in range(i):
        seq2_aligned = seq2_aligned[:0] + '_' + seq2_aligned[0:]
    for repeat in range(j):
        seq1_aligned = seq1_aligned[:0] + '_' + seq1_aligned[0:]

    return max_score, seq1_aligned, seq2_aligned

"""fastas = fasta_parser(sys.argv[3])
for i in fastas:
    print i[1]+"\n"

matrix = matrix_parser(sys.argv[2])
print matrix
"""

fastas1 = fasta_parser(sys.argv[3]) #Create fastas list
fastas2 = fasta_parser(sys.argv[4]) #Create fastas list
matrix = matrix_parser(sys.argv[2]) #Create score matrix

output_file = open("Align.fasta", "w")
if sys.argv[1]=='-g':

    print '|  /^_| _ |_  _ |       /\ |. _  _  _ _  _  _ _|_  |'
    print '|  \_/|(_)|_)(_||      /~~\||(_|| || | |(/_| | |   |'
    print '|____________________________ _|___________________|'

    for k1, v1 in fastas1:
        for k2, v2 in fastas2:
            if k1 != k2 and k1[6] < k2[6]:  # ignore duplicates compares
                output_file.write('\nOptimal Score: {}\n{key1}:\n{}\n{key2}:\n{}\n'
                                  .format(*global_alignment(v1, v2, matrix), key1=k1, key2=k2))
                print('\nOptimal Score: {}\n{key1}:\n{}\n{key2}:\n{}\n'.format(*global_alignment(v1, v2, matrix), key1=k1, key2=k2))

elif sys.argv[1]=='-l':
    print '|  |  _  _ _ |       /\ |. _  _  _ _  _  _ _|_  |'
    print '|  |_(_)(_(_||      /~~\||(_|| || | |(/_| | |   |'
    print '|_________________________ _|___________________|'

    for k1, v1 in fastas1:
        for k2, v2 in fastas2:
            if k1!=k2 and k1[6]<k2[6]: #ignore duplicates compares
                output_file.write('\nOptimal Score: {}\n{key1}:\n{}\n{key2}:\n{}\n\nSub-Optimal Score: {}\n{key1}:\n{}\n{key2}:\n{}\n'
                    .format(*local_alignment(v1, v2, matrix), key1=k1, key2=k2))
                print('\nOptimal Score: {}\n{key1}:\n{}\n{key2}:\n{}\n\nSub-Optimal Score: {}\n{key1}:\n{}\n{key2}:\n{}\n'.format(*local_alignment(v1, v2, matrix), key1=k1, key2=k2))

elif sys.argv[1]=='-o':
    print '|  /^\   _  _  | _  _        /\ |. _  _  _ _  _  _ _|_  |'
    print '|  \_/\/(/_|   |(_||_)      /~~\||(_|| || | |(/_| | |   |'
    print '|__________________|______________ _|___________________|'
    for k1, v1 in fastas1:
        for k2, v2 in fastas2:
            if k1!=k2 and k1[6]<k2[6]: #ignore duplicates compares
                output_file.write('\nOptimal Score: {}\n{key1}:\n{}\n{key2}:\n{}\n'
                        .format(*free_ends_alignment(v1, v2, matrix), key1=k1, key2=k2))
                print('\nOptimal Score: {}\n{key1}:\n{}\n{key2}:\n{}\n'.format(*free_ends_alignment(v1, v2, matrix), key1=k1, key2=k2))
else:
    print('You can use only flags -g, -l, -o right after the program name\n')

print 'Fasta file was created\n'

"""
scores_matrix = matrix_parser(sys.argv[2])
print global_alignment("AAAC","AGC",scores_matrix)
print local_alignment("AAAC", "AGC", scores_matrix)
print free_ends_alignment("AAAC", "AGC", scores_matrix)
"""
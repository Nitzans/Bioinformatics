import time
from Alignments import global_alignment, local_alignment, free_ends_alignment, matrix_parser, fasta_parser

fastas1 = fasta_parser('sequences.fasta')
fastas2 = fasta_parser('sequences.fasta')
score_matrix = matrix_parser('Score_matrix.txt')

print 'Starting benchmark...\n'

print 'Running Global Alignments..\n'
times = list(range(5))
for t in range(5):
    start = time.clock()
    for k1, v1 in fastas1:
        for k2, v2 in fastas2:
            global_alignment(v1, v2, score_matrix)

    end = time.clock()
    times[t] = end - start
    print 'Global alignment #{} ran for: {} sec'.format(t, times[t])

print '\n'
print 'Average time: {} sec'.format(sum(times)/len(times))
print 'Best time: {} sec'.format(str(min(times)))
print 'Worst time: {} sec'.format(str(max(times)))

print('----------------------------------\n')

print('Running Local Alignments...\n')
times = list(range(5))
for t in range(5):
    start = time.clock()
    for k1, v1 in fastas1:
        for k2, v2 in fastas2:
            local_alignment(v1, v2, score_matrix)

    end = time.clock()
    times[t] = end - start
    print 'Local alignment #{} ran for: {} s'.format(t, times[t])

print '\n'
print 'Average time: {} sec'.format(sum(times)/len(times))
print 'Best time: {} sec'.format(str(min(times)))
print 'Worst time: {} sec'.format(str(max(times)))

print('----------------------------------\n')

print('Running Free-ends Alignments..\n')
times = list(range(5))
for t in range(5):
    start = time.clock()
    for k1, v1 in fastas1:
        for k2, v2 in fastas2:
            free_ends_alignment(v1, v2, score_matrix)

    end = time.clock()
    times[t] = end - start
    print('Free-ends alignment #{} ran for: {} s'.format(t, times[t]))

print '\n'
print 'Average time: {} sec'.format(sum(times)/len(times))
print 'Best time: {} sec'.format(str(min(times)))
print 'Worst time: {} sec'.format(str(max(times)))

print '\n'
print 'Done!'
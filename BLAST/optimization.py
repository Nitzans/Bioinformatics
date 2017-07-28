#!/usr/bin/python
from blast import *
import time

"""
This part running differents combination for the three parameters
w, T, X to help us find the optimized trio
"""

WMER = 6
WORD_THRESHOLD = 30
EXCEEDED_THRESHOLD = 4

texts = fasta_parser('text[1].fasta')
queries = fasta_parser('queries[1].fasta')
matrix = matrix_parser('Score.matrix')

times = []
texts_hash = text_pre_process(texts, WMER)

for WMER in range(6, 11 + 1):
    for WORD_THRESHOLD in range(WMER * 5 - 18, WMER * 5 + 1, 9):  # only two mismatches allowed
        for EXCEEDED_THRESHOLD in range(4, 12 + 1, 4):
            print('\n============================================\nWMER: {}\tWORD_THRESHOLD: {}\tEXCEEDED_THRESHOLD: {}'
                  .format(WMER, WORD_THRESHOLD, EXCEEDED_THRESHOLD))
            results = {}
            for qid, query in texts.items():
                start = time.clock()
                for tid, text in queries.items():
                    blast = Blast(text,
                                  query,
                                  texts_hash[tid],
                                  matrix,
                                  WMER,
                                  WORD_THRESHOLD,
                                  EXCEEDED_THRESHOLD)

                    blast.extend_hits()
                    query_coverage, text_coverage = blast.compute_coverage(2)
                    results[(query_coverage, text_coverage)] = ('{} vs {}\tquery_coverage:\t{}%\ttext_coverage:\t{}%'
                                                                .format(qid, tid, query_coverage, text_coverage))

            for qid, tid in sorted(results):
                print(results[(qid, tid)])

            end = time.clock()
            times.append((qid, end - start))

            print('\n---------------------------------------------------------')
            total_run = 0
            for t in times:
                total_run += t[1]
                print "Query {} against all texts ran for: {} sec".format(t[0], t[1])
            print "~~~~~ For parameters: W = {}, T = {}, X = {} ~~~~~".format(WMER, WORD_THRESHOLD, EXCEEDED_THRESHOLD)
            print "\nTotal running time of blast algorithm: {} minutes".format(round(total_run / 60, 1))

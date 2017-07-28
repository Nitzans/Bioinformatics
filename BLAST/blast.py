#!/usr/bin/python
from itertools import islice, product
import time
from Alignments import *

""" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BLAST Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ """


class Blast:
    def __init__(self, text, query, text_hash, scores, wmer, word_threshold, exceeded_threshold):
        self.text = text
        self.query = query
        self.scores = scores
        self.wmer = wmer  # w
        self.word_threshold = int(word_threshold)  # T
        self.exceeded_threshold = int(exceeded_threshold)  # X
        self.query_ranges = []
        self.text_ranges = []
        self.text_hash = text_hash
        self.query_hash = create_hash(create_sliding_window(wmer, query))
        self.query_neighbors_table = create_query_neighbors(self.query_hash,
                                                            self.wmer,
                                                            self.word_threshold,
                                                            self.scores,
                                                            self.text)
        self.hits_table = find_hits(self.text_hash, self.query_neighbors_table)

    # extend hits wrapper
    def extend_hits(self):
        for hit, values in self.hits_table.items():  # iterate for hits and their value (index,score)
            query_indices = self.query_hash[hit]
            for value in values:
                score = value[1]
                for text_index in value[0]:
                    for query_index in query_indices:
                        text_left, text_right, query_left, query_right, score = self.extend_hit(text_index, query_index,
                                                                                                score, score)
                        # after extending update the ranges
                        self.query_ranges.append((int(query_left), int(query_right), score))
                        self.text_ranges.append((int(text_left), int(text_right), score))

    # Receives text indices, query indices and score and try to extend left and right as long as it is possible
    def extend_hit(self, text_index, query_index, score, best_score):
        # update indices for current segment
        text_left, text_right = text_index, text_index + self.wmer - 1
        query_left, query_right = query_index, query_index + self.wmer - 1

        # We continue extending until: 1) out of bounds 2) exceeded X points below current score
        while True:
            left_score, best_l = self.should_extend(text_left, query_left, score, best_score, 'left')
            right_score, best_r = self.should_extend(text_right, query_right, score, best_score, 'right')
            best_score = max(best_l, best_r)
            # print('left score is {} and right score is {}'.format(left_score, right_score))
            if left_score is None:
                if right_score is None:
                    # cannot extend more
                    return text_left, text_right, query_left, query_right, score

                # We can extend only right
                else:
                    score = right_score
                    query_right += 1
                    text_right += 1

            # We can extend only left
            elif right_score is None:
                score = left_score
                query_left -= 1
                text_left -= 1

            # Can extend in both directions, we choose the higher among them
            else:
                if left_score >= right_score:
                    score = left_score
                    query_left -= 1
                    text_left -= 1
                else:
                    score = right_score
                    query_right += 1
                    text_right += 1

    # Check that we don't exceed string bounds
    def should_extend(self, text_index, query_index, current_score, best_score, direction):
        if direction is 'left':
            text_index -= 1
            query_index -= 1
            if text_index < 0 or query_index < 0:
                return None, None
        elif direction is 'right':
            text_index += 1
            query_index += 1
            if text_index >= len(self.text) or query_index >= len(self.query):
                return None, None

        new_score = int(self.scores[(self.text[text_index], self.query[query_index])]) + int(current_score)
        if best_score < new_score:  # if the score improved - save it
            best_score = new_score
        if new_score >= best_score - self.exceeded_threshold:  # check if new score doesn't drop X below max score
            return new_score, best_score

        return None, None

    # coverage wrapper (for the final stage)
    def compute_coverage(self, precision):
        final_query = merge_overlap_seg(self.query_ranges)  # query after merge overlapped
        final_text = merge_overlap_seg(self.text_ranges)  # text after merge overlapped
        query_coverage, query_normalized = compute_coverage(final_query, self.query)
        text_coverage, text_normalized = compute_coverage(final_text, self.text)
        return round(float(query_coverage), precision), query_normalized, round(float(text_coverage),
                                                                                precision), text_normalized

""" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Helper Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ """

""" Part 2.1 - Text Pre-process"""


# This function gets a list of sliding windows and mapped each to its index in the text
def create_hash(windows):
    result = {}
    for i, slidewind in enumerate(windows):  # enumerate attach index i for every window in the list (tuple(i, window))
        result.setdefault(slidewind, []).append(i)  # if slidwind doesn't exist in dict then setdefault add empty list,
                                                    # afterwards append the index of the sliding window to this list.
    return result


# Iterates the text and saves in a list any w-length word
def create_sliding_window(wmer, text):
    window = (islice(text, i, None) for i in range(wmer))  # islice split string to chars in range (i,i+wmer) from text
    return (''.join(i) for i in list(zip(*window)))  # unzip the tuple's list and create list of ('X','Y','Z') windows


# Pre-processing all the texts using the two previous functions
def text_pre_process(texts, wmer):
    return {key: create_hash(create_sliding_window(wmer, txt)) for key, txt in texts.items()}


""" Part 2.2 - Query pre-process"""


# Check which neighbors have score above threshold
def create_query_neighbors(query_hash, wmer, threshold, score_matrix, text):
    if text is not None:
        permutations = list(create_sliding_window(wmer, text))
    else:
        permutations = [''.join(p) for p in product(['A', 'C', 'G', 'T'], repeat=wmer)]
    result = {}
    for query_window in query_hash.keys():
        for neighbor in permutations:
            score = similarity_score(query_window, neighbor, score_matrix)
            if score >= threshold:
                result.setdefault(query_window, []).append([neighbor, score])
    return result


# This function gets query window and some neighbor and calculate the differences score base on the matrix score
def similarity_score(q_window, neighbor, scores):
    result = 0
    for queryChar, neighborChar in zip(q_window, neighbor):  # for each char sum the match score
        result += int(scores[queryChar, neighborChar])
    return result


""" Part 2.3 - Finding hits """


# For every query finds its neighbors in the text and insert them as ([indices in the text], scores)
def find_hits(text_hash, query_hash):
    result = {}
    for query_key, similar_words in query_hash.items():
        result[query_key] = [(text_hash[word[0]], word[1]) for word in similar_words if word[0] in text_hash]
        #  if neighbor is in text - add it to the list
    return result


""" Extra Part -  """


# After extending we merge overlapped segments in order to filtered unnecessary segments and create monotonically ranges
def merge_overlap_seg(ranges):
    ranges = iter(sorted(ranges))
    current_start, current_stop, current_score = next(ranges)  # init first iteration
    for start, stop, score in ranges:
        if start > current_stop:  # Gap between segments: add the current segment and start a new one.
            yield current_start, current_stop, (current_score * (current_stop - current_start))  # adding the new range
            current_start, current_stop, current_score = start, stop, score  # jump to next segment
        else:  # Segments adjacent or overlapping - merged them
            if current_stop < stop:  # the range start after and end after the cu
                current_stop = stop
                current_score += score * (stop - current_stop)
    yield current_start, current_stop, current_score


def compute_coverage(ranges, text):
    count = 0
    total_score = 0
    for r1, r2, score in ranges:
        count += float(r2 - r1)
        total_score += score
    length = float(len(text))
    cover = float(count / length) * 100
    normalized = round(total_score * float(count / length), 2)
    return cover, normalized


""" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Parsing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ """


def fasta_parser(fasta_file):
    with open(fasta_file, mode='r') as fasta_file:
        fasta_dict = {}
        current_seq = ''
        for fasta_line in fasta_file:
            if fasta_line.startswith('>'):
                current_seq = fasta_line[1:-1]
                fasta_dict[current_seq] = ''
            else:
                fasta_dict[current_seq] += fasta_line.replace('\n', '')
                # fasta_dict = sorted(fasta_dict.items()) #return list of tuples
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
                score_matrix[(matrix[0][x], matrix[y][0])] = matrix[y][x]  # match for each substitution its score
    return score_matrix


""" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Set optimized parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ """

WMER = 5
WORD_THRESHOLD = 0
EXCEEDED_THRESHOLD = 20

texts = fasta_parser('text[1].fasta')
queries = fasta_parser('queries[1].fasta')
matrix = matrix_parser('Score.matrix')

""" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Running Part 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ """
# Print the two best texts for each query in local alignment
"""" Remove comment to activate """
# best_locals()

""" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Running Part 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ """
output_file = open("output.txt", "w")
times = []
texts_hash = text_pre_process(texts, WMER)  # sliding window for the text

print('\n=====================================================\n'
      'WMER: {}\tWORD_THRESHOLD: {}\tEXCEEDED_THRESHOLD: {}\n'
      .format(WMER, WORD_THRESHOLD, EXCEEDED_THRESHOLD))
output_file.write('\n=====================================================\n'
                  'WMER: {}\tWORD_THRESHOLD: {}\tEXCEEDED_THRESHOLD: {}\n'
                  .format(WMER, WORD_THRESHOLD, EXCEEDED_THRESHOLD))
results = {}

for qid, Qseq in sorted(queries.items()):
    start = time.clock()
    best_cover = ("", 0, 0)
    semi_best_cover = ("", 0, 0)
    for tid, Tseq in sorted(texts.items()):
        # print('Calculating {} vs {}...'.format(qid, tid))
        blast = Blast(Tseq,
                      Qseq,
                      texts_hash[tid],
                      matrix,
                      WMER,
                      WORD_THRESHOLD,
                      EXCEEDED_THRESHOLD)

        blast.extend_hits()
        query_coverage, query_normalized, text_coverage, text_normalized = blast.compute_coverage(2)
        results[qid, tid] = ('{} vs {}\tquery coverage:\t{}%\ttext coverage:\t{}%\n\t\t\tnormalized score with text: {}'
                             .format(qid, tid, query_coverage, text_coverage, text_normalized))
        if best_cover[1] < text_coverage:
            semi_best_cover = best_cover
            best_cover = (tid, text_coverage, text_normalized)
        elif semi_best_cover[1] < text_coverage:
            semi_best_cover = (tid, text_coverage, text_normalized)

    print ("Query {} two best matches:\n\t{} with {}% coverage and normalized score of {}\n\t"
           "{} with {}% coverage and normalized score of: {}\n"
           .format(qid, best_cover[0], best_cover[1], best_cover[2], semi_best_cover[0], semi_best_cover[1],
                   semi_best_cover[2]))
    output_file.write("Query {} two best matches:\n\t{} with {}% coverage and normalized score of {}\n\t"
                      "{} with {}% coverage and normalized score of: {}\n"
                      .format(qid, best_cover[0], best_cover[1], best_cover[2], semi_best_cover[0], semi_best_cover[1],
                              semi_best_cover[2]))

    end = time.clock()
    times.append((qid, end - start))

print('\n---------------------------------------------------------')
output_file.write('\n---------------------------------------------------------\n')
for qid, tid in sorted(results):
    print(results[(qid, tid)])

total_run = 0
for t in times:
    total_run += t[1]
    print "Query {} against all texts ran for: {} sec".format(t[0], t[1])
    output_file.write("Query {} against all texts ran for: {} sec\n".format(t[0], t[1]))

print "\nTotal running time of blast algorithm: {} minutes".format(round(total_run / 60, 1))
output_file.write("\nTotal running time of blast algorithm: {} minutes\n".format(round(total_run / 60, 1)))

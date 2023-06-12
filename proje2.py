import Bio
from Bio.Align import PairwiseAligner
aligner = PairwiseAligner()

#aligner.mode = 'local' 
#'local' for local alignment
#aligner.match_score = 2
#aligner.mismatch_score = -1
#aligner.open_gap_score = -0.5
#aligner.extend_gap_score = -0.1
"""
seq1 = 'ACCGT'
seq2 = 'ACG'
alignment = aligner.align(seq1, seq2)

#1. örnek aligment score 3
for aln in alignment:
    print('Alignment score:', aln.score)
    print('Sequence 1:', aln.aligned[0])
    print('Sequence 2:', aln.aligned[1])


#2.örnek global aligment
aligner.mode = 'global'
seq1 = 'ACCGT'
seq2 = 'ACG'
alignment = aligner.align(seq1, seq2)

for aln in alignment:
    print('Alignment score:', aln.score)
    print('Sequence 1:', aln.aligned[0])
    print('Sequence 2:', aln.aligned[1])


#import Bio
#print(Bio._version_)

#3. örnek local aligment 
aligner.mode='local'
seq1 = 'ACCGT'
seq2 = 'ACG'
alignment = aligner.align(seq1, seq2)

for aln in alignment:
    print('Alignment score:', aln.score)
    print('Sequence 1:', aln.aligned[0])
    print('Sequence 2:', aln.aligned[1])


#4. örnek global aligment 
aligner.mode='global'
aligner.match_score = 2 # her match için 2 puan
aligner.mismatch_score = -1 # her mismatch için -1 puan
seq1 = 'ACCGT'
seq2 = 'ACG'
alignment = aligner.align(seq1, seq2)

for aln in alignment:
    print('Alignment score:', aln.score)
    print('Sequence 1:', aln.aligned[0])
    print('Sequence 2:', aln.aligned[1])


#5. örnek 
aligner.mode='global'
aligner.match_score = 2 # her match için 2 puan
aligner.mismatch_score = -1 # her mismatch için -1 puan
aligner.open_gap_score = -0.5 #açık boşluk cezası yeni bir boşluğun oluşturulması için bir cezadır
aligner.extend_gap_score = -0.1 #uzatma boşluğu cezası ise mevcut bir boşluğu uzatmak için cezadır

seq1 = 'ACCGT'
seq2 = 'ACG'
alignment = aligner.align(seq1, seq2)

for aln in alignment:
    print('Alignment score:', aln.score)
    print('Sequence 1:', aln.aligned[0])
    print('Sequence 2:', aln.aligned[1])


#6. örnek
aligner.mode='global'
aligner.match_score = 5 # her match için 5 puan
aligner.mismatch_score = -4 # her mismatch için -4 puan
aligner.open_gap_score = -1 # açık boşluk için -1 puan
aligner.extend_gap_score = -0.1 # boşluk uzatma için -0.1 puan

seq1 = 'A'
seq2 = 'T'
alignment = aligner.align(seq1, seq2)

for aln in alignment:
    print('Alignment score:', aln.score)
    print('Sequence 1:', aln.aligned[0])
    print('Sequence 2:', aln.aligned[1])

#7. örnek
aligner.mode='global'
aligner.match_score = 5 # her match için 5 puan
aligner.mismatch_score = -4 # her mismatch için -4 puan
aligner.open_gap_score = -3 # açık boşluk için -3 puan
aligner.extend_gap_score = -0.1 # boşluk uzatma için -0.1 puan

seq1 = 'A'
seq2 = 'T'
alignment = aligner.align(seq1, seq2)

for aln in alignment:
    print('Alignment score:', aln.score)
    print('Sequence 1:', aln.aligned[0])
    print('Sequence 2:', aln.aligned[1])

# 8. örnek
from Bio.Align import substitution_matrices
matrix = substitution_matrices.load("BLOSUM62")
aligner.mode = 'global'
aligner.substitution_matrix=matrix
seq1 = 'KEVLA'
seq2 = 'EVL'

alignment = aligner.align(seq1, seq2)

for aln in alignment:
    print('Alignment score:', aln.score)
"""
"""
#9.örnek
from math import log
def gap_function(x, y):  # x is gap position in seq, y is gap length
    if y == 0:  # No gap
        return 0
    elif y == 1:  # Gap open penalty
        return -2
    return - (2 + y/4.0 + log(y)/2.0)
aligner.match_score = 5
aligner.mismatch_score = -4
aligner.open_gap_score = gap_function()
aligner.extend_gap_score = gap_function()
#globalmc fonksiyonu Align modülü içinde bulunmuyor
alignment = aligner.align.globalmc("ACCCCCGT", "ACG")

for aln in alignment:
    print('Alignment score:', aln.score)
"""
from math import log
from pairwise2 import align, format_alignment

# Define custom gap function
def gap_function(x, y):  # x is gap position in seq, y is gap length
    if y == 0:  # No gap
        return 0
    elif y == 1:  # Gap open penalty
        return -2
    return - (2 + y/4.0 + log(y)/2.0)

# Set match and mismatch scores
match_score = 5
mismatch_score = -4

# Align sequences using custom gap function and match/mismatch scores
alignment = align.globalmc("ACCCCCGT", "ACG", match_score, mismatch_score, gap_function, gap_function)

# Print alignment
for aln in alignment:
    print(format_alignment(*aln))

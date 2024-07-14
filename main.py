from generation import *
from alignment import *
from hmm import *


alphabet = ["A", "C", "G", "T"]
patterns = ["AATTGA", "CGCTTAT", "GGACTCAT", "TTATTCGTA"]
similarity = 1
dissimilarity = -0.5
penalty = -1

print("---= Part 1: Dataset Generation =---")
# Generate the dataset
dataset = generate_dataset(alphabet, patterns, 50)

# Out of the dataset, pick the first 15 samples for subsetA and the rest for  subsetB
subsetA = dataset[:15]
subsetB = dataset[15:]

print("The two subsets are:")
print("Subset A:")
for index, sequence in enumerate(subsetA):
    print(f"Sequence {index + 1}/{len(subsetA)}:{' ' if index + 1 < 10 else ''} {sequence}")
print("Subset B:")
for index, sequence in enumerate(subsetB):
    print(f"Sequence {index + 1}/{len(subsetB)}:{' ' if index + 1 < 10 else ''} {sequence}")

print("---= Part 2: Multiple Sequence Alignment =---")
# Perform multiple sequence alignment on subsetA
aligned_A = align(subsetA, [similarity, dissimilarity, penalty])

print("The aligned sequences of subset A are:")
for index, sequence in enumerate(aligned_A):
    print(f"Sequence {index + 1}/{len(aligned_A)}:{' ' if index + 1 < 10 else ''} {sequence}")

print("---= Part 3: Hidden Markov Model Profiling =---")
# Instantiate the HMM model object and find and print the HMM sequence as well as the alignment scores and the alignment paths of subsetB based on the HMM
HMM = HMM(alphabet + ['_'], aligned_A)
hmm_sequence = HMM.calculate_hmm_sequence()
# Parallel lists for each sequence
scores, paths = HMM.align_dataset(subsetB)

print(f"The Hidden Markov Model Profile is: {hmm_sequence}")
for index, sequence in enumerate(subsetB):
    print(f"""Sequence {index + 1}/{len(subsetB)}: {sequence}, 
        Score: {scores[index]}, 
        Path: {paths[index]}""")

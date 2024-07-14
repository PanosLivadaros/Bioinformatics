def align(dataset, scores):
    aligned_sequences = []

    base_sequence = ""
    # Find and set the biggest sequence as the base for the multiple sequence alignment process
    for sequence in dataset:
        if len(base_sequence) < len(sequence):
            base_sequence = sequence
    aligned_sequences.append(base_sequence)

    # Perform multiple sequence alignment for the rest of the sequences in the dataset
    for sequence in dataset:
        if sequence not in aligned_sequences:
            aligned_sequences.append(calculate_grid(aligned_sequences[0], sequence, scores))

    return aligned_sequences


# Calculate the alignment grid between 2 sequences
def calculate_grid(sequenceA, sequenceB, scores):
    # Initialize the 2D Matrix
    grid = [[0] * (len(sequenceA) + 1) for _ in range(len(sequenceB) + 1)]

    # Set the initial gap penalties on the first horizontal line and on the first vertical line
    for i in range(1, len(sequenceB) + 1):
        grid[i][0] = grid[i - 1][0] + scores[2]
    for j in range(1, len(sequenceA) + 1):
        grid[0][j] = grid[0][j - 1] + scores[2]

    # Calculate cell scores
    for i in range(1, len(sequenceB) + 1):
        for j in range(1, len(sequenceA) + 1):
            top_left = grid[i - 1][j - 1] + (scores[0] if sequenceB[i - 1] == sequenceA[j - 1] else scores[1])
            left = grid[i - 1][j] + scores[2]
            top = grid[i][j - 1] + scores[2]
            grid[i][j] = max(top_left, left, top)

    # After filling in the grid determine the best path for the sequence alignment
    path = find_path(grid, sequenceA, sequenceB, scores)
    aligned_sequence = ["_"] * len(sequenceA)

    # If moving diagonally add indexed letter from sequenceB to the aligned sequence, else if moving vertically or horizontally add "_" to the aligned sequence
    for p in path:
        if p[1] > 0:
            aligned_sequence[p[1] - 1] = sequenceB[p[0] - 1] if p[0] > 0 else "_"

    return "".join(aligned_sequence)


# Finds the best path for sequence alignment and returns a list of the path's cells' coordinates
def find_path(grid, sequenceA, sequenceB, scores):
    path = []
    i, j = len(sequenceB), len(sequenceA)

    while i > 0 or j > 0:
        path.append((i, j))
        if grid[i][j] == grid[i - 1][j - 1] + (scores[0] if sequenceA[j - 1] == sequenceB[i - 1] else scores[1]):
            i -= 1
            j -= 1
        elif grid[i][j] == grid[i - 1][j] + scores[2]:
            i -= 1
        else:
            j -= 1

    path.append((0, 0))
    path.reverse()
    return path

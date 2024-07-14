import numpy


class HMM:
    alphabet = []
    states = ['M', 'D']
    aligned_dataset = []
    transition_probs = []
    emission_probs = []

    # Constructor
    def __init__(self, alphabet, aligned_sequences):
        self.alphabet = alphabet
        self.aligned_dataset = aligned_sequences
        self.transition_probs, self.emission_probs = self.calculate_probabilities(self.aligned_dataset)

    # Calculate the emission probabilities and the transition probabilities for the aligned dataset
    def calculate_probabilities(self, aligned_dataset):
        # Initialize transition and emission probabilities
        transition_counts = numpy.zeros((len(self.states), len(self.states)))
        emission_counts = numpy.zeros((len(self.states), len(self.alphabet)))

        for seq in aligned_dataset:
            prev_state = None
            for symbol in seq:
                if symbol == "_":
                    state = "D"
                else:
                    state = "M"

                state_idx = self.states.index(state)
                symbol_idx = self.alphabet.index(symbol)

                if prev_state is not None:
                    prev_state_idx = self.states.index(prev_state)
                    transition_counts[prev_state_idx][state_idx] += 1

                emission_counts[state_idx][symbol_idx] += 1
                prev_state = state

        # Convert counts to probabilities
        transition_probs = transition_counts / transition_counts.sum(axis=1, keepdims=True)
        emission_probs = emission_counts / emission_counts.sum(axis=1, keepdims=True)
        return transition_probs, emission_probs

    # Find, print and return the HMM profile sequence of states
    def calculate_hmm_sequence(self):
        sequence = []
        for i in range(len(self.aligned_dataset[0])):
            cols = []
            for j in range(len(self.aligned_dataset)):
                cols.append(self.aligned_dataset[j][i])
            if "_" not in cols:
                sequence.append([f"M{i + 1}"])
            elif cols.count("_") == len(self.aligned_dataset):
                sequence.append([f"D{i + 1}"])
            else:
                sequence.append([f"M{i + 1}", f"D{i + 1}"])

        return sequence

    # Calculate the alignment scores and the alignment paths
    def align_dataset(self, dataset):
        alignment_scores = self.forward_algorithm(dataset)
        alignment_paths = self.viterbi_algorithm(dataset)
        return alignment_scores, alignment_paths

    # Calculate the alignment scores for the given dataset, based on the emission probabilities and the transition probabilities
    def forward_algorithm(self, dataset):
        alignment_scores = []
        for sequence in dataset:
            sequence = list(sequence)
            state_count = len(self.states)
            seq_len = len(sequence)

            # Initialize forward probabilities matrix
            fwd = numpy.zeros((state_count, seq_len + 1))
            fwd[:, 0] = 1.0 / state_count  # Starting probability is uniformly distributed

            for i in range(1, seq_len + 1):
                symbol = sequence[i - 1]
                symbol_idx = self.alphabet.index(symbol)

                for state_idx in range(state_count):
                    for prev_state_idx in range(state_count):
                        fwd[state_idx, i] += fwd[prev_state_idx, i - 1] * self.transition_probs[
                            prev_state_idx, state_idx] * self.emission_probs[state_idx, symbol_idx]

            alignment_scores.append(numpy.sum(fwd[:, seq_len]))
        return alignment_scores

    # Calculate the alignment paths for the given dataset, based on the emission probabilities and the transition probabilities
    def viterbi_algorithm(self, dataset):
        alignment_path = []
        for sequence in dataset:
            sequence = list(sequence)
            state_count = len(self.states)
            seq_len = len(sequence)

            # Initialize viterbi matrix and backpointer matrix
            viterbi = numpy.zeros((state_count, seq_len + 1))
            backpointer = numpy.zeros((state_count, seq_len + 1), dtype=int)

            viterbi[:, 0] = 1.0 / state_count  # Starting probability is uniformly distributed

            for i in range(1, seq_len + 1):
                symbol = sequence[i - 1]
                symbol_idx = self.alphabet.index(symbol)

                for state_idx in range(state_count):
                    max_prob = -1
                    max_state = 0

                    for prev_state_idx in range(state_count):
                        prob = viterbi[prev_state_idx, i - 1] * self.transition_probs[prev_state_idx, state_idx] * \
                               self.emission_probs[state_idx, symbol_idx]

                        if prob > max_prob:
                            max_prob = prob
                            max_state = prev_state_idx

                    viterbi[state_idx, i] = max_prob
                    backpointer[state_idx, i] = max_state

            # Backtrack to find the optimal path
            best_path = []
            best_last_state = numpy.argmax(viterbi[:, seq_len])
            best_path.append(best_last_state)

            for i in range(seq_len, 0, -1):
                best_last_state = backpointer[best_last_state, i]
                best_path.append(best_last_state)

            best_path.reverse()
            # Append the states to the path and enumerate them
            alignment_path.append([f"{self.states[state]}{index + 1}" for index, state in enumerate(best_path[1:])])
        return alignment_path

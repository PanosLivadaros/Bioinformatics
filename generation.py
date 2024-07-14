import random


# Generates a dataset of strings based on given alphabet and patterns
def generate_dataset(alphabet, patterns, dataset_size):
    dataset = []

    for _ in range(dataset_size):

        # Pick 1-3 random letters from the alphabet and add them to the start of the sequence
        new_seq = ""
        for _ in range(random.randint(1, 3)):
            new_seq += random.choice(alphabet)

        # Pick each patterns once modify it max 2 times and then add them to the new sequence string
        for pattern in patterns:
            new_seq += modify(pattern, alphabet + ["_"])

        # Pick 1-2 random letters from the alphabet and add them to the end of the sequence
        for _ in range(random.randint(1, 2)):
            new_seq += random.choice(alphabet)

        dataset.append(new_seq)

    return dataset


# Changes the pattern in 1-2 places and returns it
def modify(pattern, alphabet):
    pattern = list(pattern)
    restricted_position = -1

    for i in range(random.randint(1, 2)):

        # Initialize the restricted position to avoid replacing the same symbol twice
        if i == 1:
            restricted_position = position

        position = random.randint(0, len(pattern) - 1)
        # Check to avoid replacing the same symbol twice
        while position == restricted_position:
            position = random.randint(0, len(pattern) - 1)

        element = random.choice(alphabet)
        # Check to avoid replacing a symbol with itself
        while pattern[position] == element:
            element = random.choice(alphabet)

        pattern[position] = element

    return "".join(pattern)

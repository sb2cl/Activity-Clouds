import itertools
import statistics

POSSIBLE_CHARACTERS = ['A', 'C', 'G', 'T', '_']
RBS_CORE_LENGTH = 6


class RBSSequence:
    def __init__(self, sequence, mean=None):
        self.sequence = sequence
        self.mean = mean
        self.std = 0
        self.cv = 0
        self.children = []

    def generate_children(self):
        """
        Generate all possible children sequences of this sequence.
        :return:
        """

        children = []

        if '_' in self.sequence:

            for i in range(RBS_CORE_LENGTH):

                if self.sequence[i] == '_':
                    for letter in "ACGT":
                        new_sequence = self.sequence[:i] + letter + self.sequence[i + 1:]
                        children.append(new_sequence)

        return children

    def calculate_statistics(self):
        """
        Update all the statistics of this sequence.
        :return:
        """
        if self.children:
            means = [child.mean for child in self.children]
            self.mean = statistics.mean(means)
            if len(means) > 1:
                self.std = statistics.stdev(means)
            else:
                self.std = 0
            if self.mean != 0:
                self.cv = self.std / self.mean
            else:
                self.cv = 0
        else:
            self.std = 0
            self.cv = 0

    def is_specific(self):
        """
        Returns True if this core RBS is specific (i.e. no wildcard is present).
        :return: boolean
        """
        return '_' not in self.sequence

    def __repr__(self):
        return (f"RBSSequence(sequence='{self.sequence}', mean={self.mean}, "
                f"std={self.std}, cv={self.cv})")


def collect_node_data(node):
    """
    Collect the data of all the nodes and produce a dict of {sequence: cv}.
    :param node:
    :return:
    """
    # Initialize an empty dictionary to store the results
    result = {}

    # Recursive helper function to traverse the tree
    def traverse(current_node):
        # Add the current node's sequence and cv to the dictionary
        result[current_node.sequence] = current_node.cv

        # Recursively traverse each child node
        for child in current_node.children:
            traverse(child)

    # Start traversing from the root node
    traverse(node)

    return result


def get_all_paths(root):
    paths = []

    def dfs(node, path):
        if node is None:
            return

        # Add the current node to the path
        path.append(node.sequence)

        # If the current node is a leaf, save the path
        if not node.children:
            paths.append(path[:])
        else:
            # Continue the DFS for each child
            for child in node.children:
                dfs(child, path)

        # Backtrack: remove the current node from the path
        path.pop()

    # Initialize the DFS from the root with an empty path
    dfs(root, [])

    return paths


def generate_all_sequences(length):
    """
    Generate all possible sequences of a given length using POSSIBLE_CHARACTERS.

    Args:
        length (int): The length of the sequences to generate.

    Yields:
        str: The next sequence of the specified length.
    """
    for combination in itertools.product(POSSIBLE_CHARACTERS, repeat=length):
        yield ''.join(combination)
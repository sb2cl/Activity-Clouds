import itertools
import statistics
import subprocess
import re

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

# Anti-Shine-Dalgarno sequence (16S rRNA E. coli)
ANTI_SD_SEQ = "ACCUCCUUA"

def extract_number(line):
    """
    Extracts a number from a string provided by RNAFold in parenthesis.

    Args:
        line (str): The string to extract the number from.

    Returns:
        float: The extracted number, or None if not found.
    """
    match = re.search(r'\(\s*([-+]?[0-9]*\.?[0-9]+)\s*\)\s*$', line)
    return float(match.group(1)) if match else None

def compute_dG_SD_interaction(sd_candidate):
    """
    Compute the Gibbs free energy (ΔG) of SD-aSD interaction using RNAcofold.

    Returns:
        dG (float): The Gibbs free energy (ΔG) of SD-aSD interaction.
    """
    input_data = f"{sd_candidate}&{ANTI_SD_SEQ}"

    # Run RNAcofold
    result = subprocess.run(["RNAcofold", "−−noPS"], input=input_data, text=True, capture_output=True)
    output = result.stdout.strip().split("\n")

    # Get the numerical value from the parenthesis from the second line
    if len(output) > 1:
        dG = extract_number(output[1])
        return dG

    return None

def find_best_SD(mRNA_seq, search_range=(-14, -5), start_codon="ATG"):
    """
    Go through all possible 6-nt sequences in the 5' UTR and select the one with the strongest binding to 16S rRNA.

    Returns:
        best_sd (str): The best SD sequence.
        best_spacer (int): The distance from the end of SD to ATG.
        best_dG (float): The Gibbs free energy (ΔG) of SD-aSD interaction.
        aug_pos (int): The position of the AUG codon.
    """
    # Find the position of the start codon
    aug_pos = mRNA_seq.rfind(start_codon)

    if aug_pos == -1:
        return None, None, None  # Start codon not found

    # Trim the sequence to the 5' UTR region (typically -14 to -4 nt from ATG)
    search_start = max(0, aug_pos + search_range[0])
    search_end = aug_pos + search_range[1]
    search_region = mRNA_seq[search_start:search_end]

    best_sd = None
    best_dG = float("inf")
    best_spacer = None

    # Go through all possible 6-nt sequences in the given range from right to left
    for i in range(len(search_region) - 6, -1, -1):
        candidate_sd = search_region[i: i + 6]
        dG = compute_dG_SD_interaction(candidate_sd)

        if dG is not None and dG < best_dG:  # Looking for the strongest binding (lowest ΔG)
            best_sd = candidate_sd
            best_dG = dG
            best_spacer = (aug_pos - (search_start + i + 6))  # Distance from the end of SD to ATG

    return best_sd, best_spacer, best_dG, aug_pos

def run_RNAfold(mRNA_seq):
    """
    Run RNAfold and get an mRNA structure (dot bracket notation) and Gibbs free energy.

    Returns:
        structure (str): The mRNA structure.
        dG (float): The Gibbs free energy.
    """
    result = subprocess.run(["RNAfold", "--MEA", "−−noPS"], input=mRNA_seq, text=True, capture_output=True)
    output = result.stdout.strip().split("\n")

    structure = output[1].split(" ")[0]  # Dot-bracket structure
    dG = float(output[1].split(" ")[-1].strip("()"))  # Gibbs energy
    return structure, dG
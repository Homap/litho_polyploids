import sys
import re

def remove_branch_lengths(input_tree, output_tree):
    with open(input_tree, 'r') as infile:
        tree_content = infile.read()
    cleaned_tree = re.sub(r":\d*\.?\d*", "", tree_content)
    with open(output_tree, 'w') as outfile:
        outfile.write(cleaned_tree)

input_tree = sys.argv[1]
output_tree = sys.argv[2]

remove_branch_lengths(sys.argv[1], sys.argv[2])
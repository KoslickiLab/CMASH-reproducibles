import sys

# test import
def test_run():
    print("Happy coding : )")

# read file and build TST tree
def build_ref(input_file):
    f = open(input_file, 'r')
    query_list = [x.strip() for x in list(f)]
    f.close()

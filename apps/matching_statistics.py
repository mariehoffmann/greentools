import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compute matching statistics of two strings.')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--file', dest='file_name', action='store_const', help='give file name of file storing both strings in fasta format')
    group.add_argument('--strings', dest='strings', type=string, nargs=2, help='give directly two input strings')

    args = parser.parse_args()
    print args

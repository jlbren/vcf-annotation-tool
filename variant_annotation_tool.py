from vcf_annotate import VCFAnnotate
from argparse import ArgumentParser
import os

def _parse_arguments():
    arg_parser = ArgumentParser(description='VCF Annotation Tool',
                                          prog='VCFAnnotate')
    arg_parser.add_argument('input_file', nargs=1,
                             help='Input VCF file to be annotated.')
    arg_parser.add_argument('--local-only', action='store_true', default=False,
                             help='Optional flag to skip calls to ExAC API.\n'
                                  'Missing values will be marked as NA.')
    arg_parser.add_argument('-o', '--output', default=os.getcwd(),
                             help='Optional flag to specify output directory.\n'
                                  'Defaults to current working directory.')

    return arg_parser.parse_args()

def main():
    args = _parse_arguments()
    # Create a new VCFAnnotate object.
    va = VCFAnnotate(args.input_file[0], args.local_only, args.output)
    # Parse VCF file. 
    va.parse()
    va.annotate()
main()


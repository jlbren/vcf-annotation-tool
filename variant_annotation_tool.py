from vcf_parser import VCFParser
import argparser 

  def _parse_arguments(self):
        arg_parser = argparse.ArgumentParser(description='VCF Annotation Tool',
                                             prog='VCFAnnotate')
        arg_parser.add_argument('input_file', nargs=1, required=True,
                                help='Input VCF file to be annotated.')
        arg_parser.add_argument('--local-only', action='store_true', default=False,
                                help='Optional flag to skip calls to ExAC API.\n'
                                     'Missing values will be marked as NA.')
        arg_parser.add_argument('-o', '--output', default=os.getcwd()
                                help='Optional flag to specify output directory.\n'
                                     'Defaults to current working directory.')

        return arg_parser.parse_args()

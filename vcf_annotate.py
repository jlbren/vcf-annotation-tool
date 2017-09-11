from vcf_parser import VCFParser
import argparse
import os

class VCFAnnotate:
    def __init__(self, vcf_file, local_only):
        self.file_name = vcf_file
        self.skip_api = local_only


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
    def parse(self):
        # Create a VCFParser object from the input file. 
        self.vcf_data = VCFParser(self.file_name)
    
    def annotate(self):
            

    def split_info_field(self, line):
        """Method for splitting a VCF info field into a dict of its keys and values.
        Arguments:
            line (str): Single info entry to be split into a dictionary.
        Returns:
            info (dict): Dictionary containing each info field as a key 
                         with its corresponding value.
        """
        #Get key-value pairs. 
        pairs = line.split(';')
        info = {}
        for pair in pairs:
            field = pair.split('=')
            # Add key-value to dictionary 
            info[field[0]] = field[1]

        return info


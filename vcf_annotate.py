from vcf_parser import VCFParser
import os
import pandas
import requests

class VCFAnnotate:
    def __init__(self, vcf_file, local, output):
        self.file_name = vcf_file
        self.local_only = local
        self.out_path = output

    def parse(self):
        """Wrapper method for creating a VCFParser object"""
        # Create a VCFParser object from the input file.
        self.vcf_data = VCFParser(self.file_name)

    def annotate(self):
        out_name = 'annotated_' + self.vcf_data.meta_data.file_name
        self.out_file = os.path.join(self.out_path, out_name)
        variant_type = []
        depth_per_base = []
        alt_obs = []
        percent_variant = []

        for row in self.vcf_data.vcf_df['INFO']:
            fields = self.split_info_field(row)
            variant_type.append(fields['TYPE'])
            depth_per_base.append(fields['DPB'])
            var_alleles = fields['AO']
            ref_alleles = int(fields['RO'])
            # Sum observations for multiple variants.
            if var_alleles.find(',') > -1:
                # Split list of observations and map to int.
                var_list = map(int, var_alleles.split(','))
                # Sum observations.
                var_alleles = sum(var_list)
            else:
                var_alleles = int(var_alleles)
            alt_obs.append(var_alleles)
            # Catch divide by zero errors where there are no reference observations
            try:
                var_vs_ref = float(var_alleles / ref_alleles)
            except ZeroDivisionError:
                # Percent variant is 100% where no ref obs where found.
                var_vs_ref = 1.00

            # Convert back to string and format to two decimals.
            percent_variant.append( '%.2f' % var_vs_ref)

        fields_df = pandas.DataFrame({'TYPE': variant_type,
                                      'DBP': depth_per_base,
                                      'AO': alt_obs,
                                      'AO/RO': percent_variant})
        print(fields_df)

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

    def get_exac_api_request(self, variant):
        request_url = 'http://exac.hms.harvard.edu/rest/variant/variant/'
                      + variant 
        response = requests.get(url=request_url).json() # TODO check other args 





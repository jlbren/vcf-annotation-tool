from vcf_parser import VCFParser
from collections import OrderedDict
import os
import pandas
import requests

class VCFAnnotate:
    """Class for initiating and managing the VCF annotation pipeline.
       Makes use of VCFParser to read an input VCF file and produces an
       output table with annotations.
    Attributes:
        input_file (str): Full path to input VCF file.
        local_only (bool): Flag to disable sending requests to the ExAC API.
        out_path (str): Path to output directory.
        vcf_data (VCFParser): VCFParser object generated from input file.
        out_file (str): Full path to output file in the form of annotated_<input_file_name>.
        out_table (pandas.DataFrame): Pandas dataframe containing annotation output data.
    """
    def __init__(self, vcf_file, local, output):
        self.input_file = vcf_file
        self.local_only = local
        self.out_path = output

    def parse(self):
        """Wrapper method for creating a VCFParser object"""
        # Create a VCFParser object from the input file.
        self.vcf_data = VCFParser(self.input_file)

    def annotate(self):
        """Public method for starting the annotation process."""
        # Create output file.
        out_name = 'annotated_' +os.path.basename( self.input_file)
        self.out_file = os.path.join(self.out_path, out_name)
        # Extract relevant data from VCFParser object
        self.out_table = self.parse_info_col()
        # Get info from ExAC if local_only was not passed at runtime.
        if not self.local_only:
            self.get_api_data()

    def write_out(self):
        """Write output dataframe to file. Represent missing vals as NA.
           Do not include quotes or row numbers.
        Outputs:
            annotated_input_file_name.vcf: Annotated output file.
        """
        self.out_table.to_csv(path_or_buf=self.out_file, na_rep='NA',
                              index=False)

    def parse_info_col(self):
        """Method for parsing relavent info from the INFO col of the VCF dataframe.
        Returns:
            fields_df (pandas.DataFrame): Returns dataframe consisting of mutation type,
                                          Depth of coverage per base, number of alternate
                                          observations, and alternate obs/reference obs
        """
        variant_type = []
        depth_per_base = []
        alt_obs = []
        percent_variant = []
        # Walk through each row of the INFO col.
        for row in self.vcf_data.vcf_df['INFO']:
            # Split fields string into a dictionary.
            fields = self.split_info_field(row)
            # Split TYPE field into list to check for multiple variants.
            var_type = fields['TYPE'].split(',')
            if len(var_type) > 1:
                # Get most harmful variant from list.
                var_type = self.rank_mutations(var_type)
            else:
                # Convert singleton list back to string.
                var_type = var_type[0]
            # Get TYPE annotation and append to output df.
            type_annotation = self.get_type_annotation(var_type)
            variant_type.append(type_annotation)
            depth_per_base.append(fields['DPB'])
            var_alleles = fields['AO'].split(',')
            ref_alleles = int(fields['RO'])
            # Sum observations for multiple variants.
            if len(var_alleles) > 1:
                # Convert all values to int.
                var_list = map(int, var_alleles)
                # Sum observations.
                var_alleles = sum(var_list)
            else:
                var_alleles = int(var_alleles[0])
            alt_obs.append(var_alleles)
            # Catch divide by zero errors where there are no reference observations
            try:
                var_vs_ref = float(var_alleles / ref_alleles)
            except ZeroDivisionError:
                # Percent variant is 100% where no ref obs where found.
                var_vs_ref = 1.00
            # Convert back to string and format to two decimals.
            percent_variant.append( '%.2f' % var_vs_ref)
        # Create dataframe using OrderedDict to preserve col order. 
        fields_df = pandas.DataFrame(OrderedDict({'CHROM': self.vcf_data.vcf_df['CHROM'],
                                                  'POS': self.vcf_data.vcf_df['POS'],
                                                  'VARIANT': self.vcf_data.vcf_df['ALT'],
                                                  'REF': self.vcf_data.vcf_df['REF'],
                                                  'TYPE': variant_type,
                                                  'DBP': depth_per_base,
                                                  'AO': alt_obs,
                                                  'AO/RO': percent_variant}))
        return fields_df

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

    def rank_mutations(self, mutations):
        """Method for determining the most deleterious mutation where multiple variants
           are present.
        Arguments:
            mutations (list(str)): List of all mutation types from multiple variants.
        Returns:
            worst (string): Single most harmful mutation from list of variants.
        """
        worst = ''
        for mutation in mutations:
            # Complex mutations are the most harmful due to compounding effects.
            if mutation == 'complex':
                return mutation
            # Indels are more harmful than substituions due to effects on
            # downstream reading frames.
            elif mutation == 'del' or mutation == 'ins':
                worst = mutation
            # Substitutions (SNPs, MNPs) are least harmful.
            elif mutation == '':
                worst = mutation
        return worst

    def get_type_annotation(self, variant_type):
        """Method that accepts a variant type field and returns an expanded annotation string.
        Arguments:python past values of arrays
            variant_type (str): Single variant type field.
        Returns:
            (str): Expanded annotation for variant type.
        """

        if variant_type == 'snp' or variant_type == 'mnp':
            return 'Substitution'
        elif variant_type == 'del':
            return 'Deletion'
        elif variant_type == 'ins':
            return 'Insertion'
        elif variant_type =='complex':
            return 'Complex'
        # Default case where type does not match known keys.
        else:
            return variant_type

    def get_api_data(self):
        """Method for passing request keys and inserting corresponding
           parsed responses into the output dataframe.
        """
        # Get list of list from selected df cols.
        request_keys = self.get_request_keys()
        allele_freq = []
        genes = []
        consequences = []
        for key in request_keys:
            # Send request to ExAC API for each key.
            exac_response = self.get_exac_api_request(key)
            # Check if response is populated.
            if len(exac_response) != 0:
                parsed_response = self.parse_exac_response(exac_response)
                allele_freq.append(parsed_response['FREQ'])
                genes.append(parsed_response['GENES'])
                consequences.append(parsed_response['CONSEQUENCE'])
            else:
                # Append missing values.
                allele_freq.append('NA')
                genes.append('NA')
                consequences.append('NA')

        exac_df = pandas.DataFrame({'ALLELE_FREQ': allele_freq,
                                    'GENE': genes,
                                    'CONSEQUENCE': consequences})
        # Add cols to output dataframe.
        self.out_table = pandas.concat([self.out_table, exac_df], axis=1)

    def get_request_keys(self):
        """Method that pastes relevant dataframe cols together to form exac api request keys.
        Returns:
            request_keys (list(str)): List of variant request keys for the ExAC API.
        """
        # Slice rows from dataframe needed to form request key.
        key_vals = self.vcf_data.vcf_df[['CHROM', 'POS', 'REF', 'ALT']].values
        request_keys = []
        for row in key_vals:
            # Convert CHROM and POS to strings.
            row = [str(x) for x in row]
            # Join list of strings and append.
            request_keys.append('-'.join(row))
        return request_keys

    def get_exac_api_request(self, key):
        """Method which takes in a request key, querries the ExAC API, and
           returns the corresponding response.
        Arguments:
            key (str): ExAC variant request key in the form CHROM-POS-VARIANT-REF.
        Returns:
           response (dict): ExAC response as JSON serialized into a dictionary.
                            Empty if response was not received.
        """
        # Form request from base URL and key.
        request_url = ('http://exac.hms.harvard.edu/rest/variant/' + key)
        try:
            response = requests.get(url=request_url).json()
        # Catch broken responses and return an empty dict.
        except requests.exceptions.RequestException as e:
            print('ExAC API connection timeout for key: %s.\n'
                  'Missing values will be marked as NA.'% key)
            #print(e)
            return {}
        return response

    def parse_exac_response(self, exac_response):
        """Method for taking in a complete ExAC response and returning a reduced
           dictionary consisting of relevant information.
        Arguemnts:
            exac_response (dict): Complete response as a serialized JSON dictionary.
        Returns:
            parsed_values (dict): Dictionary consisting only of the response fields of interest.
        """
        parsed_values = {}
        # Try-catch blocks for selecting each desired response field.
        # Except KeyErrors to return values as NA
        try:
            variant = exac_response['variant']
        except KeyError:
            parsed_values['FREQ'] = 'NA'
            parsed_values['GENES'] = 'NA'
            return parsed_values
        try:
            parsed_values['FREQ'] = variant['allele_freq']
        except KeyError:
            parsed_values['FREQ'] = 'NA'
        try:
            # Join all genes into single string.
            parsed_values['GENES'] = ';'.join(variant['genes'])
        except KeyError:
            parsed_values['GENES'] = 'NA'
        try:
            consequence = exac_response['consequence']
            # Make sure the consequence field is not empty.
            if consequence is not None:
                # Join all consequences into a single string.
                parsed_values['CONSEQUENCE'] = ';'.join(consequence.keys())
            else:
                parsed_values['CONSEQUENCE'] = 'NA'
        except KeyError:
            parsed_values['CONSEQUENCE'] = 'NA'

        return parsed_values


